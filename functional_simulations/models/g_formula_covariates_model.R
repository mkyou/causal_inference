library(mgcv)
library(data.table)

source("functional_simulations/utils/revert_to_long_format.R")

prepare_covariate_data = function(factual_data, treatment_times) {
  dt = copy(factual_data)
  setkey(dt, id, t)
  
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  
  dt[, (lagged_names) := shift(.SD, n = 1, type = "lag"), 
     by = id, .SDcols = cov_names]
  
  for (s_time in treatment_times) {
    A_s_name = paste0("A_s", s_time)
    treatment_status_at_s = factual_data[t == s_time, .(id, A)]
    setnames(treatment_status_at_s, "A", A_s_name)
    
    dt = merge(dt, treatment_status_at_s, by = "id", all.x = TRUE)
  }
  return(dt)
}


fit_covariate_models_mgcv = function(prepared_data, treatment_times) {
  
  fitted_models = list()
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  
  fit_data = na.omit(prepared_data)
  fit_data[, id := as.factor(id)]
  
  for (k_name in cov_names) {
    
    # Termo de base para o intercepto funcional f(t)
    base_formula = '1 + s(t, k = 15)'
    
    # Efeito do tratamento (coeficiente funcional para A_s)
    treatment_formula = paste(
      sapply(treatment_times, function(s) {
        paste0('s(t, by = A_s', s, ', k = 8)')
      }), collapse = " + "
    )
    
    other_lagged = setdiff(lagged_names, paste0(k_name, "_lag1"))
    history_formula = paste(
      sapply(other_lagged, function(lag) {
        paste0('te(t, ', lag, ', k = c(8, 8))')
      }),
      collapse = " + "
    )

    model_formula = as.formula(
      paste(k_name, "~", base_formula, "+", 
            treatment_formula, "+", history_formula)
    )
    
    cat("\n--- Fórmula para a covariável:", k_name, "---\n")
    print(model_formula)
    
    fitted_models[[k_name]] = bam(
      model_formula, 
      data = fit_data,
      discrete = TRUE,
      select = TRUE
    )
  }
  
  return(fitted_models)
}


run_gformula_covariates_for_assessment = function(
    scenario_data, 
    scenario_config, 
    regime_a, 
    regime_b
) {
  
  cat("--- Iniciando G-Formula (mgcv) para Covariáveis ---\n")
  
  # Converte os dados do formato wide (entrada) para o formato long,
  # que é o formato esperado pelo resto do código.
  scenario_data_long = revert_to_long_format(scenario_data)
  
  treatment_times = scenario_config$treatment_times
  
  # 1. Preparar os dados factuais para o ajuste
  # Usamos o novo objeto 'scenario_data_long' que está no formato correto.
  factual_data = scenario_data_long[, .SD, .SDcols = c(
    "id", "t", "A", grep("^X[0-9]+$", names(scenario_data_long), value = TRUE)
  )]
  prepared_data = prepare_covariate_data(factual_data, treatment_times)
  
  # 2. Ajustar os modelos de covariáveis
  covariate_models = fit_covariate_models_mgcv(prepared_data, treatment_times)
  
  # 3. Simular para cada regime
  # Coloque esta função dentro de run_gformula_covariates_mgcv,
  # substituindo a original.
  
  simulate_one_regime_debug = function(regime_to_sim, regime_name) {
    
    cat(paste0("  -- Iniciando simulação para o Regime: ", regime_name, " --\n"))
    
    cov_names = names(covariate_models)
    ids = unique(factual_data$id)
    n_subjects = length(ids)
    id_levels = levels(as.factor(ids))
    time_grid = sort(unique(factual_data$t))
    
    sim_names = paste0(cov_names, "_sim")
    
    # A correção para o colapso da variância.
    residual_sds = list()
    for (k in cov_names) {
      residual_sds[[k]] = sd(residuals(covariate_models[[k]]), na.rm = TRUE)
      cat(
        paste0(
          "   - [Setup] Desvio Padrão dos Resíduos para ", 
          k, ": ", 
          round(residual_sds[[k]], 4), "\n"
        )
      )
    }
    
    cov_ranges = list()
    fit_data_for_ranges = na.omit(prepared_data)
    for (k in cov_names) {
      cov_ranges[[k]] = range(fit_data_for_ranges[[k]], na.rm = TRUE)
    }
    
    simulated_dt = CJ(id = ids, t = time_grid)
    for (k in cov_names) {
      simulated_dt[, paste0(k, "_sim") := NA_real_]
    }
    
    first_time = factual_data[
      t == time_grid[1], .SD, .SDcols = c("id", cov_names)
    ]
    setnames(first_time, cov_names, sim_names)
    first_time[, t := time_grid[1]]
    simulated_dt[
      first_time, on=.(id,t), 
      (sim_names) := mget(paste0("i.", sim_names))
    ]
    
    for (i in 2:length(time_grid)) {
      current_t = time_grid[i]
      previous_t = time_grid[i-1]
      
      lagged_names = paste0(cov_names, "_lag1")
      history = simulated_dt[t == previous_t, .SD, .SDcols = c("id", sim_names)]
      setnames(history, sim_names, lagged_names)
      
      # LOG de depuração: Verificar a variância ANTES da predição
      if(i %% 20 == 0) {
        cat(
          paste0(
            "   - [", regime_name, "] t = ", 
            round(current_t, 2), ". Verificando variância de entrada:\n"
          )
        )
        print(summary(history[, .SD, .SDcols = lagged_names]))
      }
      
      if(any(!is.finite(as.matrix(history[, .SD, .SDcols = lagged_names])))) {
        cat("Dados não finitos encontrados na 'history' antes da predição!\n")
      }
      
      pred_data_t = data.table(id = ids, t = current_t)
      for (j in 1:length(treatment_times)) {
        s_time = treatment_times[j]
        pred_data_t[, paste0("A_s", s_time) := regime_to_sim[j]]
      }
      pred_data_t = merge(pred_data_t, history, by = "id")
      pred_data_t[, id := factor(id, levels = id_levels)]
      
      for (k in cov_names) {
        predictions_mean = predict(
          covariate_models[[k]], 
          newdata = as.data.frame(pred_data_t)
        )
        
        if (any(!is.finite(predictions_mean))) {
          cat(
            paste0(
              "ALERTA: Previsão MÉDIA não finita para ", k, 
              " em t=", current_t, "!\n"
            )
          )
          browser()
        }
        
        # CORREÇÃO: Adicionando ruído para prevenir o colapso da variância
        noise = rnorm(n_subjects, mean = 0, sd = residual_sds[[k]])
        predictions_stochastic = predictions_mean + noise
        
        min_val = cov_ranges[[k]][1]
        max_val = cov_ranges[[k]][2]
        predictions_final = pmax(min_val, pmin(max_val, predictions_stochastic))
        
        simulated_dt[t == current_t, (paste0(k, "_sim")) := predictions_final]
      }
    }
    
    cat(
      paste0("  -- Simulação para o Regime: ", regime_name, " concluída --\n")
    )
    return(simulated_dt[, .SD, .SDcols = c("id", "t", sim_names)])
  }
  
  cat("  - Simulando trajetórias para o Regime A...\n")
  sim_a = simulate_one_regime_debug(regime_a, "Regime A (Treatment)")
  
  cat("  - Simulando trajetórias para o Regime B...\n")
  sim_b = simulate_one_regime_debug(regime_b, "Regime B (No Treatment)")
  
  # 4. Combina os resultados
  sim_cols = names(sim_a)[grepl("_sim$", names(sim_a))]
  
  # Cria os novos nomes de forma limpa: sim_a_X1, sim_a_X2, etc.
  new_names_a = paste0("sim_a_", gsub("_sim", "", sim_cols))
  new_names_b = paste0("sim_b_", gsub("_sim", "", sim_cols))
  
  setnames(sim_a, old = sim_cols, new = new_names_a)
  setnames(sim_b, old = sim_cols, new = new_names_b)
  
  # Usamos 'scenario_data_long' para o merge final
  results_dt = merge(scenario_data_long, sim_a, by = c("id", "t"))
  results_dt = merge(results_dt, sim_b, by = c("id", "t"))
  
  cat("--- Estágio de Covariáveis Concluído ---\n")
  return(results_dt)
}