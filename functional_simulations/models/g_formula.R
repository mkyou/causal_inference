library(mgcv)
library(data.table)

# Carrega as funções que já criamos e depuramos para a etapa das covariáveis
source("functional_simulations/models/g_formula_covariates_model.R")
source("functional_simulations/utils/revert_to_long_format.R")

run_gformula_covariates_for_simulation = function(
    scenario_data, 
    scenario_config, 
    regime_a, 
    regime_b
) {
  
  cat("--- Iniciando G-Formula (Covariates for SIMULATION) ---\n")
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
  simulate_one_regime = function(regime_to_sim, regime_name) {
    
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
  
  sim_a = simulate_one_regime(regime_a, "Regime A (Treatment)")
  sim_b = simulate_one_regime(regime_b, "Regime B (No Treatment)")
  
  sim_cols = names(sim_a)[grepl("_sim$", names(sim_a))]
  new_names_a = paste0("sim_a_", gsub("_sim", "", sim_cols))
  new_names_b = paste0("sim_b_", gsub("_sim", "", sim_cols))
  setnames(sim_a, old = sim_cols, new = new_names_a)
  setnames(sim_b, old = sim_cols, new = new_names_b)
  
  # Junta apenas os resultados das duas simulações
  simulation_results = merge(sim_a, sim_b, by = c("id", "t"))
  
  cat("--- Estágio de Covariáveis (Simulation) Concluído ---\n")
  
  # Retorna APENAS o id, tempo e as trajetórias simuladas
  return(simulation_results)
}

run_gformula_full = function(
    scenario_data,
    scenario_config,
    regime_a,
    regime_b
) {
  
  cat("--- Iniciando G-Formula Completa ---\n")
  
  # --- ETAPA 1: SIMULAÇÃO DAS COVARIÁVEIS ---
  cat("  - Etapa 1: Simulando trajetórias das covariáveis...\n")
  covariate_sims = run_gformula_covariates_for_simulation(
    scenario_data, scenario_config, regime_a, regime_b
  )
  
  # --- ETAPA 2: AJUSTE DO MODELO DE DESFECHO ---
  cat("  - Etapa 2: Preparando e ajustando modelo de y...\n")
  scenario_data_long = revert_to_long_format(scenario_data)
  
  # CORREÇÃO: Usar um nome de variável consistente
  full_data_for_fitting = merge(
    scenario_data_long, covariate_sims, by = c("id", "t")
  )
  
  treatment_times = scenario_config$treatment_times
  cat("    - Criando colunas de tratamento (A_s)...\n")
  for (s_time in treatment_times) {
    A_s_name = paste0("A_s", s_time)
    treatment_status_at_s = full_data_for_fitting[t == s_time, .(id, A)]
    setnames(treatment_status_at_s, "A", A_s_name)
    full_data_for_fitting = merge(full_data_for_fitting,
                                  treatment_status_at_s, by = "id", all.x = TRUE)
  }
  a_s_cols = paste0("A_s", treatment_times)
  full_data_for_fitting[, (a_s_cols) := lapply(
    .SD,
    function(x) ifelse(is.na(x), 0, x)
    ), .SDcols = a_s_cols
  ]
  
  # CORREÇÃO: Usar a variável com o nome correto
  fit_data_y = na.omit(full_data_for_fitting)
  fit_data_y[, id := as.factor(id)]
  
  cov_names = grep("^X[0-9]+$", names(fit_data_y), value = TRUE)
  base_formula_y = "s(t, k = 15)"
  cov_form_y = paste(
    sapply(
      cov_names,
      function(k) paste0("s(", k, ") + ti(t, ", k, ")")
    ), collapse = " + "
  )
  
  treat_form_y = paste(
    sapply(
      treatment_times,
      function(s) paste0("s(t, by = A_s", s, ", k=8)")
    ), collapse = " + "
  )
  
  model_formula_y = as.formula(
    paste("y ~", base_formula_y, "+", cov_form_y, "+", treat_form_y)
  )
  
  cat("    - Fórmula do modelo de desfecho:\n")
  print(model_formula_y)
  
  outcome_model = bam(model_formula_y, data = fit_data_y,
                      discrete = TRUE)
  
  # --- ETAPA 3: PREDIÇÃO ---
  cat("  - Etapa 3: Prevendo desfechos contrafactuais...\n")
  newdata_a = copy(full_data_for_fitting)
  sim_a_cols = names(newdata_a)[grepl("^sim_a_X", names(newdata_a))]
  factual_cov_cols = gsub("sim_a_", "", sim_a_cols)
  setnames(newdata_a, sim_a_cols, factual_cov_cols)
  for (i in 1:length(treatment_times)) {
    s_time = treatment_times[i]
    newdata_a[, paste0("A_s", s_time) := regime_a[i]]
  }
  
  newdata_b = copy(full_data_for_fitting)
  sim_b_cols = names(newdata_b)[grepl("^sim_b_X", names(newdata_b))]
  setnames(newdata_b, sim_b_cols, factual_cov_cols)
  for (i in 1:length(treatment_times)) {
    s_time = treatment_times[i]
    newdata_b[, paste0("A_s", s_time) := regime_b[i]]
  }
  
  Y_a_sim = predict(outcome_model, newdata = as.data.frame(newdata_a))
  Y_b_sim = predict(outcome_model, newdata = as.data.frame(newdata_b))
  
  # --- ETAPA 4: EMPACOTANDO OS RESULTADOS ---
  cat("  - Etapa 4: Empacotando resultados...\n")
  
  cov_results = copy(full_data_for_fitting)
  
  outcome_results = full_data_for_fitting[, .(id, t, a_y, b_y)]
  outcome_results[, y_a_sim := Y_a_sim]
  outcome_results[, y_b_sim := Y_b_sim]
  
  causal_effect_results = full_data_for_fitting[, .(id, t)]
  causal_effect_results[, effect_estimate := Y_a_sim - Y_b_sim]
  
  cat("--- G-Formula Completa Concluída ---\n")
  return(list(
    covariates = cov_results,
    outcomes = outcome_results,
    causal_effect = causal_effect_results
  ))
}