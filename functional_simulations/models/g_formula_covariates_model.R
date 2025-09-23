library(mgcv)
library(data.table)

#' @title Ajusta os modelos para cada covariável
fit_covariate_models_mgcv = function(prepared_data, treatment_times) {
  fitted_models = list()
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  fit_data = na.omit(prepared_data)
  fit_data[, id := as.factor(id)]
  
  for (k_name in cov_names) {
    base_formula = 's(t, k = 20)'
    treatment_formula = paste(
      sapply(
        treatment_times, 
        function(s) paste0('s(t, by = A_s', s, ', k = 12)')
      ), 
      collapse = " + "
    )
    other_lagged = setdiff(lagged_names, paste0(k_name, "_lag1"))
    history_formula = paste(
      sapply(
        other_lagged, 
        function(lag) paste0('te(t, ', lag, ', k = c(12, 12))')
      ), collapse = " + "
    )
    model_formula = as.formula(
      paste(
        k_name, "~", base_formula, "+", treatment_formula, "+", 
        history_formula
      )
    )
    
    fitted_models[[k_name]] = bam(
      model_formula, data = fit_data, discrete = TRUE, select = TRUE
    )
  }
  return(fitted_models)
}

run_gformula_covariates_for_simulation = function(
    prepared_data,
    scenario_config,
    regime_a,
    regime_b
) {
  
  cat("--- Iniciando G-Formula (Etapa de Covariáveis) ---\n")
  
  # --- 1. AJUSTE DOS MODELOS DAS COVARIÁVEIS ---
  cat("  - Ajustando modelos para as covariáveis...\n")
  fitted_models = fit_covariate_models_mgcv(
    prepared_data, scenario_config$treatment_times
  )
  
  # --- 2. SIMULAÇÃO SEQUENCIAL ---
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  ids = unique(prepared_data$id)
  time_grid = sort(unique(prepared_data$t))
  
  # Esta função interna executa a simulação iterativa para um regime
  simulate_one_regime = function(regime_to_sim, regime_name) {
    cat(
      paste0(
        "  -- Iniciando simulação sequencial para o Regime: ", 
        regime_name, 
        " --\n"
      )
    )
    
    ids = unique(fit_data$id)
    n_subjects = length(ids)
    time_grid = sort(unique(fit_data$t))
    
    # Prepara o data.table que armazenará os resultados da simulação
    simulated_dt = CJ(id = ids, t = time_grid)
    for (k in cov_names) {
      simulated_dt[, paste0(k, "_sim") := NA_real_]
    }
    
    # Inicializa em t=0 com os valores factuais
    first_time_vals = prepared_data[
      t == time_grid[1], 
      .SD, 
      .SDcols = c("id", cov_names)
    ]
    
    for(k in cov_names) {
      simulated_dt[
        t == time_grid[1], 
        paste0(k, "_sim") := first_time_vals[
          match(id, first_time_vals$id), get(k)
        ]
      ]
    }
    
    # O CORAÇÃO DA G-FORMULA: O LOOP SEQUENCIAL
    for (i in 2:length(time_grid)) {
      current_t = time_grid[i]
      previous_t = time_grid[i-1]
      
      # Pega o histórico SIMULADO do passo anterior
      history = simulated_dt[
        t == previous_t, 
        .SD, 
        .SDcols = c("id", paste0(cov_names, "_sim"))
      ]
      setnames(history, paste0(cov_names, "_sim"), lagged_names)
      
      # Prepara o 'newdata' para a predição em 'current_t'
      pred_data_t = data.table(id = ids, t = current_t)
      for (j in 1:length(treatment_times)) {
        s_time = treatment_times[j]
        pred_data_t[, paste0("A_s", s_time) := regime_to_sim[j]]
      }
      pred_data_t = merge(pred_data_t, history, by = "id")
      
      # Prediz o valor de cada covariável e armazena
      for (k in cov_names) {
        predictions = predict(fitted_models[[k]], newdata = pred_data_t)
        simulated_dt[t == current_t, (paste0(k, "_sim")) := predictions]
      }
    }
    cat(
      paste0("  -- Simulação para o Regime: ", regime_name, " concluída --\n")
    )
    return(simulated_dt)
  }
  
  sim_a = simulate_one_regime(regime_a, "Regime A")
  sim_b = simulate_one_regime(regime_b, "Regime B")
  
  # Renomeia colunas para o formato final
  setnames(sim_a, paste0(cov_names, "_sim"), paste0("sim_a_", cov_names))
  setnames(sim_b, paste0(cov_names, "_sim"), paste0("sim_b_", cov_names))
  
  simulation_results = merge(sim_a, sim_b, by = c("id", "t"))
  
  cat("--- Estágio de Covariáveis Concluído ---\n")
  return(list(
    models = fitted_models,
    simulated_data = simulation_results
  ))
}