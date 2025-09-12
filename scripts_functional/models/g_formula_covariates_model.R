# ===================================================================
# G-FORMULA FUNCIONAL: MODELAGEM E VALIDAÇÃO DAS COVARIÁVEIS (VERSÃO 6.0)
#
# Abordagem com Ajuste para Confundidores Variantes no Tempo:
# - Adiciona os valores das covariáveis pré-tratamento como preditores
#   no modelo GAM para controlar o viés de seleção.
# - A função de simulação foi reescrita para ser sequencial, refletindo
#   a dependência temporal do novo modelo.
# ===================================================================

library(mgcv)
library(ggplot2)
library(data.table)

#' @title Prepara Dados para Modelagem da G-Formula (VERSÃO FINAL)
#' @description Prepara os dados adicionando colunas para o tratamento (A_s)
#'   e para o valor das covariáveis que antecedem cada tratamento (Xk_pre_s).
prepare_gformula_data = function(dt_long, treatment_times) {
  dt_prepared = copy(dt_long)
  dt_prepared[, id := as.factor(id)]
  
  all_cov_names = grep("^X[0-9]+$", names(dt_prepared), value = TRUE)
  time_grid = sort(unique(dt_prepared$t))
  
  # Adiciona os valores das covariáveis pré-tratamento (confundidores)
  for (s in treatment_times) {
    s_idx = which.min(abs(time_grid - s))
    if (s_idx > 1) {
      confounder_time = time_grid[s_idx - 1]
      confounders_at_s_minus_1 = dt_prepared[
        t == confounder_time, .SD, .SDcols = c("id", all_cov_names)
      ]
      new_conf_names = paste0(all_cov_names, "_pre_s", s)
      setnames(confounders_at_s_minus_1, all_cov_names, new_conf_names)
      dt_prepared = merge(dt_prepared, confounders_at_s_minus_1, by = "id", all.x = TRUE)
    }
  }
  
  # Adiciona as colunas de tratamento A_s
  for (s in treatment_times) {
    time_idx = which.min(abs(unique(dt_prepared$t) - s))
    actual_time = unique(dt_prepared$t)[time_idx]
    treatment_at_s = dt_prepared[t == actual_time, .SD, .SDcols = c("id", "A")]
    setnames(treatment_at_s, "A", paste0("A_s", s))
    dt_prepared = merge(dt_prepared, treatment_at_s, by = "id", all.x = TRUE)
  }
  
  return(dt_prepared)
}


#' @title Ajusta Modelos GAM para Covariáveis (VERSÃO FINAL COM INTERAÇÃO)
#' @description Usa interações com `(t >= s)` para garantir que os confundidores
#'   só afetem a trajetória APÓS o momento em que são medidos.
fit_covariate_models = function(
    prepared_data, 
    treatment_times, 
    cov_to_model
) {
  fitted_models = list()
  time_grid = sort(unique(prepared_data$t))
  
  cat(paste0("Ajustando modelos GAM (com ajuste de confundidor específico) para: ", 
             paste(cov_to_model, collapse=", "), "\n"))
  
  for (cov_name in cov_to_model) {
    
    # --- INÍCIO DA MUDANÇA PRINCIPAL ---
    
    # 1. Termo de base para a tendência de tempo
    base_formula = 's(t, bs = "cr", k = 20)'
    
    # 2. Termos de efeito do tratamento (feedback)
    treatment_effects_formula = paste(
      sapply(treatment_times, function(s) {
        paste0('s(t, bs = "tp", k = 10, by = I(A_s', s, ' * (t >= ', s, ')))')
      }), collapse = " + ")
    
    # 3. Termos de ajuste de confundimento (AGORA CONDICIONAIS NO TEMPO)
    confounder_terms_list = sapply(treatment_times, function(s) {
      s_idx = which.min(abs(time_grid - s))
      if (s_idx > 1) {
        # A fórmula agora só aplica o efeito do confundidor para t >= s
        # Usamos o nome da própria covariável para o ajuste
        conf_var = paste0(cov_name, "_pre_s", s)
        paste0(conf_var, " * I(t >= ", s, ")")
      } else { NULL }
    })
    confounder_terms = paste(unlist(confounder_terms_list), collapse = " + ")
    
    # --- FIM DA MUDANÇA PRINCIPAL ---
    
    model_formula_str = paste(cov_name, "~", base_formula, "+", treatment_effects_formula)
    if (nchar(confounder_terms) > 0) {
      model_formula_str = paste(cov_name, "~", confounder_terms, "+", base_formula, "+", treatment_effects_formula)
    }
    
    fitted_models[[cov_name]] = bam(
      as.formula(model_formula_str), data = prepared_data, discrete = TRUE, select = TRUE, method = "fREML"
    )
  }
  cat("Ajuste de modelos GAM concluído.\n")
  return(fitted_models)
}


#' @title Simula Trajetórias Contrafactuais (VERSÃO NÃO-SEQUENCIAL FINAL)
#' @description Versão final que corrige a extração de nomes de colunas
#'   para a criação da tabela de predição.
simulate_counterfactual_trajectories = function(
    fitted_models,
    original_data,
    treatment_times,
    regime,
    simple_covs_data = NULL
) {
  
  # 1. Cria a grade de predição completa
  pred_data = CJ(id = unique(original_data$id), t = sort(unique(original_data$t)))
  pred_data[, id := as.factor(id)]
  
  # 2. Adiciona as colunas do regime de tratamento
  for (i in seq_along(treatment_times)) {
    pred_data[, (paste0("A_s", treatment_times[i])) := regime[i]]
  }
  
  # 3. Adiciona as colunas de confundidores
  all_cov_names = names(fitted_models)
  time_grid = sort(unique(pred_data$t))
  
  for (s in treatment_times) {
    s_idx = which.min(abs(time_grid - s))
    if (s_idx > 1) {
      confounder_time = time_grid[s_idx - 1]
      
      # --- INÍCIO DA CORREÇÃO ---
      # Linha antiga (com bug):
      # confounder_cols_to_get = gsub("_pre_s[0-9\\.]+$", "", names(fitted_models[[1]]$coefficients)[grep("_pre_s", names(fitted_models[[1]]$coefficients))]) |> unique()
      
      # Linha nova (corrigida e simples):
      # Pega o nome de todas as covariáveis que estão sendo modeladas.
      confounder_cols_to_get = names(fitted_models)
      # --- FIM DA CORREÇÃO ---
      
      confounder_data = original_data[
        t == confounder_time, .SD, .SDcols = c("id", confounder_cols_to_get)
      ]
      
      confounder_data[, id := as.factor(id)]
      
      new_conf_names = paste0(confounder_cols_to_get, "_pre_s", s)
      setnames(confounder_data, confounder_cols_to_get, new_conf_names)
      
      pred_data = merge(pred_data, confounder_data, by = "id", all.x = TRUE)
    }
  }
  
  # 4. Faz a predição em um único passo
  simulated_dt = pred_data[, .(id, t)]
  for (cov_name in names(fitted_models)) {
    predictions = predict(fitted_models[[cov_name]], newdata = pred_data)
    simulated_dt[, (paste0(cov_name, "_sim")) := predictions]
  }
  
  simulated_dt[, id := as.numeric(as.character(id))]
  return(simulated_dt)
}


#' @title Executa a Suíte de Validação da G-Formula (VERSÃO FINAL E ROBUSTA)
run_gformula_validation_suite = function(
    scenario_data, 
    scenario_config, 
    static_model_type = "mean_trajectory"
) {
  
  dt_conf = scenario_data$confounders
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(dt_conf), value = TRUE)
  
  # --- Lógica de decisão híbrida (sem alterações) ---
  model_strategies = list()
  for (i in seq_along(cov_names)) {
    cov_name = cov_names[i]
    is_static = isTRUE(all.equal(
      scenario_config$mean_baseline_covariates[[i]](1),
      scenario_config$mean_baseline_covariates[[i]](10)
    ))
    has_feedback = !all(sapply(
      scenario_config$effects_covariates_from_A[[i]],
      function(f) f(5, 1) == 0
    ))
    if (!is_static || has_feedback) {
      model_strategies[[cov_name]] = "gam"
    } else {
      model_strategies[[cov_name]] = static_model_type
    }
  }
  
  # --- Ajuste dos modelos (sem alterações) ---
  covs_requiring_gam = names(model_strategies)[model_strategies == "gam"]
  if (length(covs_requiring_gam) > 0) {
    prepared_data = prepare_gformula_data(dt_conf, treatment_times)
    fitted_gam_models = fit_covariate_models(
      prepared_data, treatment_times, covs_requiring_gam
    )
  } else {
    fitted_gam_models = list()
  }
  
  # --- Pré-cálculo para as covariáveis simples (sem alterações) ---
  covs_requiring_mean = names(model_strategies)[model_strategies != "gam"]
  mean_trajectories_data = NULL
  if (length(covs_requiring_mean) > 0) {
    mean_trajectories_data = dt_conf[
      ,
      lapply(.SD, mean), by = id, .SDcols = covs_requiring_mean
    ]
  }
  
  # --- SIMULAÇÃO (sem alterações) ---
  cat("\nSimulando para o regime 'a_seq'...\n")
  simulated_a = simulate_counterfactual_trajectories(
    fitted_models = fitted_gam_models,
    original_data = dt_conf,
    treatment_times = treatment_times,
    regime = scenario_config$a_seq,
    simple_covs_data = mean_trajectories_data
  )
  
  cat("\nSimulando para o regime 'another_a_seq'...\n")
  simulated_another_a = simulate_counterfactual_trajectories(
    fitted_models = fitted_gam_models,
    original_data = dt_conf,
    treatment_times = treatment_times,
    regime = scenario_config$another_a_seq,
    simple_covs_data = mean_trajectories_data
  )
  
  # --- INÍCIO DA CORREÇÃO NA SEÇÃO DE VALIDAÇÃO ---
  
  # 1. Renomeia as colunas de predição ANTES de qualquer merge
  # Para simulated_a: X1_sim -> a_X1_pred
  sim_names_a = names(simulated_a)[grep("_sim$", names(simulated_a))]
  pred_names_a = paste0("a_", gsub("_sim$", "", sim_names_a), "_pred")
  setnames(simulated_a, sim_names_a, pred_names_a)
  
  # Para simulated_another_a: X1_sim -> another_a_X1_pred
  sim_names_another_a = names(simulated_another_a)[grep("_sim$", names(simulated_another_a))]
  pred_names_another_a = paste0("another_a_", gsub("_sim$", "", sim_names_another_a), "_pred")
  setnames(simulated_another_a, sim_names_another_a, pred_names_another_a)
  
  # 2. Faz o merge das duas tabelas de predição
  sim_data = merge(simulated_a, simulated_another_a, by = c("id", "t"))
  
  # 3. Faz o merge final com os dados verdadeiros
  validation_data = merge(dt_conf, sim_data, by = c("id", "t"))
  
  # --- FIM DA CORREÇÃO ---
  
  plot_data = validation_data[, lapply(.SD, mean, na.rm=TRUE), by = t]
  
  all_metrics = list()
  all_plots = list()
  
  for (cov_name in cov_names) {
    true_a_col = paste0("a_", cov_name)
    pred_a_col = paste0("a_", cov_name, "_pred")
    true_another_a_col = paste0("another_a_", cov_name)
    pred_another_a_col = paste0("another_a_", cov_name, "_pred")
    
    # Adicionando uma verificação para garantir que as colunas existem antes de calcular
    if (!all(c(true_a_col, pred_a_col, true_another_a_col, pred_another_a_col) %in% names(validation_data))) {
      warning(paste("Colunas de validação ausentes para a covariável:", cov_name))
      next # Pula para a próxima covariável
    }
    
    model_used_label = model_strategies[[cov_name]]
    if(model_used_label == "gam") {
      model_used_label = "GAM-simplified"
    }
    
    metrics = data.table(
      covariate = cov_name, model_used = model_used_label,
      RMSE_a = sqrt(mean( (validation_data[[true_a_col]] - validation_data[[pred_a_col]])^2, na.rm = TRUE )),
      Bias_a = mean( validation_data[[pred_a_col]] - validation_data[[true_a_col]], na.rm = TRUE ),
      Cor_a = tryCatch(cor( plot_data[[true_a_col]], plot_data[[pred_a_col]], use = "complete.obs" ), error = function(e) NA_real_),
      RMSE_another_a = sqrt(mean( (validation_data[[true_another_a_col]] - validation_data[[pred_another_a_col]])^2, na.rm = TRUE )),
      Bias_another_a = mean( validation_data[[pred_another_a_col]] - validation_data[[true_another_a_col]], na.rm = TRUE ),
      Cor_another_a = tryCatch(cor( plot_data[[true_another_a_col]], plot_data[[pred_another_a_col]], use = "complete.obs" ), error = function(e) NA_real_)
    )
    all_metrics[[cov_name]] = metrics
    
    plot_title = paste("Validação da G-Formula para", cov_name, "| Modelo:", model_used_label)
    subtitle_text = paste0("RMSE: ", round(metrics$RMSE_a, 3), " (Sempre) / ", round(metrics$RMSE_another_a, 3), " (Nunca) | ", "Corr: ", round(metrics$Cor_a, 3), " (Sempre) / ", round(metrics$Cor_another_a, 3), " (Nunca)")
    
    plot = ggplot(plot_data, aes(x = t)) +
      geom_line(aes(y = .data[[true_a_col]], color = "Verdadeiro (Sempre Tratar)"), linetype = "dashed", size = 1) +
      geom_line(aes(y = .data[[pred_a_col]], color = "Simulado (Sempre Tratar)"), size = 1) +
      geom_line(aes(y = .data[[true_another_a_col]], color = "Verdadeiro (Nunca Tratar)"), linetype = "dashed", size = 1) +
      geom_line(aes(y = .data[[pred_another_a_col]], color = "Simulado (Nunca Tratar)"), size = 1) +
      scale_color_manual(values = c("Verdadeiro (Sempre Tratar)"="blue", "Simulado (Sempre Tratar)"="lightblue", "Verdadeiro (Nunca Tratar)"="red", "Simulado (Nunca Tratar)"="pink")) +
      labs(title = plot_title, subtitle = subtitle_text, x = "Tempo", y = "Valor da Covariável", color = "Trajetória") +
      theme_minimal()
    
    all_plots[[cov_name]] = plot
  }
  
  return(
    list(
      metrics = rbindlist(all_metrics), 
      plots = all_plots,
      simulations = list(
        a_seq = sim_data[, .SD, .SDcols = c("id", "t", pred_names_a)],
        another_a_seq = sim_data[, .SD, .SDcols = c("id", "t", pred_names_another_a)]
      )
    )
  )
}


