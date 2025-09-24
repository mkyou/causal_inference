# ===================================================================
# Suíte de Análise para Cenários em Formato Processado
#
# Funções para analisar e validar a estrutura causal e o pré-processamento
# diretamente nos dados suavizados e processados.
# ===================================================================

library(data.table)
library(ggplot2)
library(pheatmap)

# --- 1. Exibir Diagnósticos da Suavização ---
print_smoothing_diagnostics = function(processed_data) {
  if (is.null(processed_data$diagnostics)) {
    cat("Diagnósticos de suavização não encontrados.\n")
    return(invisible(NULL))
  }
  cat("--- (Check 1/6) Diagnósticos da Etapa de Suavização ---\n")
  print(processed_data$diagnostics)
  cat("Pseudo-R2: % da variância dos dados mantida na curva suavizada.\n\n")
}

# --- 2. Checar Desenho Causal (via Config) ---
check_causal_structure_design = function(scenario_config) {
  cat("--- (Check 2/6) Verificando Estrutura Causal do Desenho ---\n")
  
  feedback_presente = !all(sapply(scenario_config$gamma_matrix, is.null))
  
  # alpha_funcs_matrix (antiga zeta) em vez de theta_integral
  confounding_markovian = any(scenario_config$theta_matrix[, -1] != 0)
  confounding_historic = !is.null(scenario_config$alpha_funcs_matrix)
  
  summary_dt = data.table(
    Verificacao = c("Feedback (A -> X) Desenhado", 
                    "Confundimento Markoviano (X(t-1) -> A(t)) Desenhado",
                    "Confundimento Histórico (∫X(u)du -> A(t)) Desenhado"),
    Status = c(ifelse(feedback_presente, "Presente", "Ausente"),
               ifelse(confounding_markovian, "Presente", "Ausente"),
               ifelse(confounding_historic, "Presente", "Ausente"))
  )
  print(summary_dt)
  cat("\n")
  return(summary_dt)
}

# --- 3. Plotar Feedback (A -> X) via Contrafactuais ---
plot_feedback_counterfactual = function(processed_data) {
  cat("--- (Check 3/6) Verificando Feedback (A -> X) via Contrafactuais ---\n")
  
  cov_names = grep(
    "^X[0-9]+$", 
    names(processed_data$functional_vars), value = TRUE
  )
  plot_list = list()
  
  for (cov_name in cov_names) {
    a_cov_name = paste0("a_", cov_name)
    b_cov_name = paste0("b_", cov_name)
    
    if (
      !all(c(a_cov_name, b_cov_name) %in% 
           names(processed_data$functional_vars))
    ) next
    
    mean_a_traj = colMeans(processed_data$functional_vars[[a_cov_name]])
    mean_b_traj = colMeans(processed_data$functional_vars[[b_cov_name]])
    
    plot_df = data.table(
      t = processed_data$time_grid, 
      Regime_A = mean_a_traj, 
      Regime_B = mean_b_traj
    )
    
    p = ggplot(melt(plot_df, id.vars="t"), aes(x=t, y=value, color=variable)) +
      geom_line(linewidth=1.2) +
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

# --- 4. Plotar Confundimento (X -> A) via Dados Factuais ---
plot_confounding_factual = function(processed_data) {
  cat("--- (Check 4/6) Verificando Confundimento (X -> A) via Dados ---\n")
  
  treatment_times = as.numeric(
    gsub("A_s", "", names(processed_data$scalar_vars)[-1])
  )
  time_grid = processed_data$time_grid
  K = length(
    grep("^X[0-9]+$", names(processed_data$functional_vars), value = TRUE)
  )
  plot_list = list()
  
  for (s_time in treatment_times) {
    prev_t_idx = which.min(
      abs(time_grid - (s_time - (time_grid[2]-time_grid[1])))
    )
    
    plot_df = data.frame(
      A = processed_data$scalar_vars[[paste0("A_s", s_time)]]
    )
    for(k in 1:K) {
      plot_df[[paste0("X", k)]] = processed_data$functional_vars[[
        paste0("X", k)
        ]][, prev_t_idx]
    }
    plot_df_melt = melt(as.data.table(plot_df), id.vars="A")
    
    p = ggplot(plot_df_melt, aes(x = value, y = A)) +
      geom_jitter(aes(color = variable), height = 0.05, alpha = 0.3) +
      geom_smooth(method = "glm", method.args = list(family = "binomial")) +
      facet_wrap(~variable, scales = "free_x") +
      labs(
        title = paste("Confundimento: Relação de X(t-1) com A(", s_time, ")"),
        x = "Valor da Covariável Suavizada em t-1", 
        y = "Prob. de Tratamento", 
        color="Covariável"
      ) +
      theme_bw()
    plot_list[[as.character(s_time)]] = p
  }
  return(plot_list)
}

# --- 5. Checar Positividade ---
check_positivity_smooth = function(processed_data) {
  cat("--- (Check 5/6) Verificando Positividade ---\n")
  
  # --- a. Proporções marginais ---
  # Calcula a proporção de tratados em cada tempo de tratamento
  treatment_props = colMeans(processed_data$scalar_vars[, -1])
  cat("Proporção de Indivíduos Tratados em cada tempo:\n")
  print(treatment_props)
  
  # --- b. Gráfico da distribuição dos escores de propensão ---
  # Extrai metadados necessários do objeto de dados processados
  treatment_times = as.numeric(
    gsub("A_s", "", names(processed_data$scalar_vars)[-1])
  )
  time_grid = processed_data$time_grid
  cov_names = grep(
    "^X[0-9]+$", 
    names(processed_data$functional_vars), value = TRUE
  )
  K = length(cov_names)
  
  plot_list = list()
  # Itera sobre cada tempo de tratamento para criar um gráfico de diagnóstico
  for (s_time in treatment_times) {
    # Encontra o índice de tempo anterior ao tratamento (s-1)
    # Necessário para pegar o valor das covariáveis que afetam o tratamento s
    dt = time_grid[2] - time_grid[1]
    prev_t_idx = which.min(abs(time_grid - (s_time - dt)))
    
    # Monta um data.frame para ajustar um modelo de propensão simples
    model_df = data.frame(
      A = processed_data$scalar_vars[[paste0("A_s", s_time)]]
    )
    for(k in 1:K) {
      cov_name = cov_names[k]
      # Adiciona o valor suavizado da covariável em s-1 como preditor
      model_df[[paste0(cov_name, "_lag1")]] = processed_data$functional_vars[[
        cov_name]][, prev_t_idx]
    }
    
    # Ajusta um modelo logístico simples para estimar os escores de propensão
    fit = glm(A ~ ., data = model_df, family = "binomial")
    prop_scores = predict(fit, type = "response")
    
    # Gera o gráfico de distribuição (histograma + densidade)
    p = ggplot(
      data.frame(prop_scores), aes(x = prop_scores)
      ) +
      geom_histogram(
        aes(y = ..density..), 
        bins = 30, 
        fill = "dodgerblue", 
        color = "black", 
        alpha = 0.7
      ) +
      geom_density(color = "blue", linewidth = 1.1) +
      labs(
        title = paste(
          "Distribuição dos Escores de Propensão Estimados em t =", s_time
        ),
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
        treatment = names(treatment_props), 
        proportion = treatment_props
      ),
      plots = plot_list
    )
  )
}

# --- 6. Checar Colinearidade ---
check_collinearity_smooth = function(processed_data) {
  cat("--- (Check 6/6) Verificando Colinearidade Funcional ---\n")
  
  cov_names = grep(
    "^X[0-9]+$", 
    names(processed_data$functional_vars), value = TRUE
  )
  mean_levels = sapply(cov_names, function(name) 
    rowMeans(processed_data$functional_vars[[name]]))
  colnames(mean_levels) = cov_names
  cor_matrix = cor(mean_levels)
  
  p = pheatmap(
    cor_matrix, 
    display_numbers = TRUE, 
    main = "Correlação Média entre Covariáveis"
  )
  
  return(list(correlation_matrix = cor_matrix, plot = p))
}

# --- Ponto de Entrada Principal ---
analyze_scenario_processed = function(
    scenario_data,
    scenario_config,
    scenario_name = "Cenário"
) {
  cat(
    paste(
      "\n--- Analisando Cenário (Dados Processados):", scenario_name, "---\n"
    )
  )
  
  if (
    !is.list(scenario_data) || 
    !all(c("functional_vars", "scalar_vars", "diagnostics") %in% 
         names(scenario_data)
    )
  ) {
    stop("Os dados não parecem estar no formato processado esperado.")
  }
  
  # 1. Diagnósticos da Suavização
  print_smoothing_diagnostics(scenario_data)
  
  # 2. Estrutura Causal do Desenho
  design_summary = check_causal_structure_design(scenario_config)
  
  # 3. Feedback (via Contrafactuais)
  feedback_plots = plot_feedback_counterfactual(scenario_data)
  
  # 4. Confundimento (via Factuais)
  confounding_plots = plot_confounding_factual(scenario_data)
  
  # 5. Positividade
  positivity_summary = check_positivity_smooth(scenario_data)
  
  # 6. Colinearidade
  collinearity_results = check_collinearity_smooth(scenario_data)
  
  results = list(
    design_summary = design_summary,
    smoothing_diagnostics = scenario_data$diagnostics,
    positivity_summary = positivity_summary,
    feedback_plots = feedback_plots,
    confounding_plots = confounding_plots,
    collinearity_results = collinearity_results
  )
  
  cat(paste("--- Análise para", scenario_name, "concluída. ---\n"))
  return(results)
}