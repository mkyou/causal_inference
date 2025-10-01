library(data.table)
library(ggplot2)
library(pheatmap)

check_causal_structure_design = function(scenario_config) {
  is_zero_func = function(f) {
    is.null(f) || identical(body(f), 0) || identical(body(f), 0L)
  }
  feedback_presente = !all(sapply(scenario_config$gamma_matrix, is_zero_func))
  confounding_markovian = any(scenario_config$theta_matrix[, -1] != 0)
  confounding_historic = !is.null(scenario_config$zeta_funcs_matrix)
  
  summary_dt = data.table(
    Verificacao = c(
      "Feedback (A -> X) Desenhado",
      "Confundimento Markoviano (X(t-1) -> A(t)) Desenhado",
      "Confundimento Histórico (∫X(u)du -> A(t)) Desenhado"
    ),
    Status = c(
      ifelse(feedback_presente, "Presente", "Ausente"),
      ifelse(confounding_markovian, "Presente", "Ausente"),
      ifelse(confounding_historic, "Presente", "Ausente")
    )
  )
  return(summary_dt)
}

plot_feedback_counterfactual = function(processed_data) {
  cov_names = grep(
    "^X[0-9]+$", names(processed_data$functional_vars), value = TRUE
  )
  plot_list = list()
  for (cov_name in cov_names) {
    a_cov_name = paste0("a_", cov_name)
    b_cov_name = paste0("b_", cov_name)
    if (!all(c(a_cov_name, b_cov_name) %in%
             names(processed_data$functional_vars))) next
    mean_a_traj = colMeans(processed_data$functional_vars[[a_cov_name]])
    mean_b_traj = colMeans(processed_data$functional_vars[[b_cov_name]])
    plot_df = data.table(
      t = processed_data$time_grid,
      Regime_A = mean_a_traj,
      Regime_B = mean_b_traj
    )
    p = ggplot(melt(plot_df, id.vars = "t"), aes(x=t, y=value, color=variable)) +
      geom_line(linewidth = 1.2) +
      labs(
        title = paste("Efeito Causal do Regime em", cov_name),
        subtitle = "Diferença E[Xk | Regime A] - E[Xk | Regime B] (Feedback)",
        x = "Tempo", y = paste("Valor Médio de", cov_name), color = "Regime"
      ) +
      theme_bw()
    plot_list[[cov_name]] = p
  }
  return(plot_list)
}

plot_confounding_factual = function(processed_data) {
  treatment_times = as.numeric(
    gsub("A_s", "", names(processed_data$scalar_vars)[-1])
  )
  time_grid = processed_data$time_grid
  K = length(grep("^X[0-9]+$", names(processed_data$functional_vars), value = TRUE))
  plot_list = list()
  for (s_time in treatment_times) {
    dt_step = time_grid[2] - time_grid[1]
    prev_t_idx = which.min(abs(time_grid - (s_time - dt_step)))
    plot_df = data.frame(A = processed_data$scalar_vars[[paste0("A_s", s_time)]])
    for (k in 1:K) {
      plot_df[[paste0("X", k)]] =
        processed_data$functional_vars[[paste0("X", k)]][, prev_t_idx]
    }
    plot_df_melt = melt(as.data.table(plot_df), id.vars = "A")
    p = ggplot(plot_df_melt, aes(x = value, y = A)) +
      geom_jitter(aes(color = variable), height = 0.05, alpha = 0.3) +
      geom_smooth(method = "glm", method.args = list(family = "binomial")) +
      facet_wrap(~variable, scales = "free_x") +
      labs(
        title = paste("Confundimento: Relação de X(t-1) com A(", s_time, ")"),
        x = "Valor da Covariável Suavizada em t-1",
        y = "Prob. de Tratamento",
        color = "Covariável"
      ) +
      theme_bw()
    plot_list[[as.character(s_time)]] = p
  }
  return(plot_list)
}

check_positivity_smooth = function(processed_data) {
  treatment_props = colMeans(processed_data$scalar_vars[, -1])
  
  # MUDANÇA CRÍTICA: Extrai o número do nome da coluna antes de converter.
  treatment_times = as.numeric(gsub("A_s", "", names(treatment_props)))
  
  time_grid = processed_data$time_grid
  cov_names = grep("^X[0-9]+$", names(processed_data$functional_vars), value = TRUE)
  K = length(cov_names)
  plot_list = list()
  
  # Adiciona uma verificação para o caso de treatment_times ser NA
  if (any(is.na(treatment_times))) {
    warning("Não foi possível extrair tempos de tratamento válidos.")
    return(list(summary_proportions = data.table(), plots = list()))
  }
  
  for (s_time in treatment_times) {
    dt_step = time_grid[2] - time_grid[1]
    prev_t_idx = which.min(abs(time_grid - (s_time - dt_step)))
    model_df = data.frame(A = processed_data$scalar_vars[[paste0("A_s", s_time)]])
    for (k in 1:K) {
      cov_name = cov_names[k]
      model_df[[paste0(cov_name, "_lag1")]] =
        processed_data$functional_vars[[cov_name]][, prev_t_idx]
    }
    fit = glm(A ~ ., data = model_df, family = "binomial")
    prop_scores = predict(fit, type = "response")
    p = ggplot(data.frame(prop_scores), aes(x = prop_scores)) +
      geom_histogram(
        aes(y = after_stat(density)), bins = 30,
        fill = "dodgerblue", color = "black", alpha = 0.7
      ) +
      geom_density(color = "blue", linewidth = 1.1) +
      labs(
        title = paste("Distribuição dos Escores de Propensão (t=", s_time, ")"),
        subtitle = "Valores próximos de 0 ou 1 indicam positividade fraca",
        x = "Probabilidade Ingênua de Tratamento (Escore de Propensão)",
        y = "Densidade"
      ) +
      coord_cartesian(xlim = c(0, 1)) +
      theme_bw()
    plot_list[[as.character(s_time)]] = p
  }
  return(
    list(
      summary_proportions = data.table(
        treatment_time = treatment_times,
        proportion = as.vector(treatment_props)
      ),
      plots = plot_list
    )
  )
}

check_collinearity_smooth = function(processed_data) {
  cov_names = grep("^X[0-9]+$", names(processed_data$functional_vars), value = TRUE)
  mean_levels = sapply(
    cov_names,
    function(name) rowMeans(processed_data$functional_vars[[name]])
  )
  colnames(mean_levels) = cov_names
  cor_matrix = cor(mean_levels)
  return(list(correlation_matrix = cor_matrix))
}

analyze_scenario_processed = function(scenario_object, scenario_name = "Cenário") {
  scenario_data = scenario_object$data
  scenario_config = scenario_object$config
  
  if (!is.list(scenario_data) ||
      !all(c("functional_vars", "scalar_vars", "diagnostics") %in%
           names(scenario_data))) {
    stop("O elemento '$data' não parece estar no formato processado esperado.")
  }
  
  results = list(
    design_summary = check_causal_structure_design(scenario_config),
    smoothing_diagnostics = scenario_data$diagnostics,
    positivity_summary = check_positivity_smooth(scenario_data),
    feedback_plots = plot_feedback_counterfactual(scenario_data),
    confounding_plots = plot_confounding_factual(scenario_data),
    collinearity_results = check_collinearity_smooth(scenario_data)
  )
  
  return(results)
}