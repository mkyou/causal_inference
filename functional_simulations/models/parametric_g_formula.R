# ===================================================================
# MODELO 2: G-Formula Funcional
#
# Contém a implementação completa do pipeline da G-formula:
# - Preparação de dados (criação de lags, etc.)
# - Ajuste dos modelos de covariáveis
# - Ajuste do modelo de desfecho
# - Simulação sequencial de trajetórias contrafactuais
# - Predição dos desfechos contrafactuais
# ===================================================================

library(mgcv)
library(data.table)

# --- Função de Preparação de Dados ---
# Usada tanto na modelagem quanto na simulação, por isso definida primeiro.
prepare_gformula_data = function(dt_long, treatment_times) {
  dt_prepared = copy(dt_long)
  
  setkey(dt_prepared, id, t)
  
  all_cov_names = grep("^X[0-9]+$", names(dt_prepared), value = TRUE)
  lagged_cov_names = paste0(all_cov_names, "_lag1")
  
  dt_prepared[, (lagged_cov_names) := shift(.SD, n=1, type="lag"), 
              by = id, .SDcols = all_cov_names]
  
  # Adiciona colunas A_s para cada tempo de tratamento
  for (s_time in treatment_times) {
    A_s_name = paste0("A_s", s_time)
    
    # Pega o valor de A no tempo de tratamento exato
    treatment_data = dt_long[t == s_time, .(id, A)]
    setnames(treatment_data, "A", A_s_name)
    
    # Garante consistência de tipos para o merge
    dt_prepared[, id := as.character(id)]
    treatment_data[, id := as.character(id)]
    
    dt_prepared = merge(dt_prepared, treatment_data, by = "id", all.x = TRUE)
    
    # Preenche NAs para indivíduos que não tinham obs naquele tempo
    if (any(is.na(dt_prepared[[A_s_name]]))) {
      dt_prepared[is.na(get(A_s_name)), (A_s_name) := 0] 
    }
  }
  
  return(dt_prepared)
}

#' @title Simula Trajetórias Contrafactuais (Sequencialmente)
#' @description Implementa a simulação sequencial para a arquitetura de
#'   histórico lagado. Em cada passo t, usa os valores simulados em t-1
#'   como preditores para o estado em t.
simulate_counterfactual_trajectories = function(
    fitted_models,
    original_data,
    treatment_times,
    regime
) {
  
  dt = copy(original_data)
  setkey(dt, id, t)
  
  dt[, id := as.factor(id)]
  ids = unique(dt$id) 
  
  covs_gam = names(fitted_models)
  all_covs = covs_gam
  time_grid = sort(unique(dt$t))
  
  original_data[, id := as.factor(id)]
  id_levels = levels(original_data$id)
  
  simulated_dt = CJ(id = ids, t = time_grid)
  for (cov_name in all_covs) {
    simulated_dt[, (paste0(cov_name, "_sim")) := NA_real_]
  }
  
  sim_names_gam = paste0(covs_gam, "_sim")
  
  first_time_data = dt[t == time_grid[1], 
                       .SD, .SDcols = c("id", covs_gam)]
  first_time_data[, t := time_grid[1]]
  setnames(first_time_data, covs_gam, sim_names_gam)
  
  simulated_dt[first_time_data, on = .(id, t),
               (sim_names_gam) := mget(paste0("i.", sim_names_gam))]
  
  # Loop de simulação sequencial (t = t_2, ..., t_T)
  for (i in 2:length(time_grid)) {
    current_t = time_grid[i]
    previous_t = time_grid[i-1]
    
    history_sim_names = paste0(all_covs, "_sim")
    lagged_names = paste0(all_covs, "_lag1")
    history_data = simulated_dt[t == previous_t, 
                                .SD, .SDcols = c("id", history_sim_names)]
    setnames(history_data, history_sim_names, lagged_names)
    
    pred_data_t = data.table(id = ids, t = current_t)
    for (j in seq_along(treatment_times)) {
      pred_data_t[, paste0("A_s", treatment_times[j]) := regime[j]]
    }
    
    pred_data_t = merge(pred_data_t, history_data, by = "id")
    
    for (cov_name in covs_gam) {
      predictions = predict(
        fitted_models[[cov_name]], 
        newdata = as.data.frame(pred_data_t)
      )
      simulated_dt[t == current_t, 
                   (paste0(cov_name, "_sim")) := predictions]
    }
  }
  
  simulated_dt[, id := as.numeric(as.character(id))]
  return(simulated_dt)
}



#' @title Ajusta modelos de covariáveis para G-formula
fit_covariate_models_gformula = function(
    prepared_data, # Recebe dados já preparados
    treatment_times, 
    covs_to_model
) {
  fitted_models = list()
  
  prepared_data[, id := as.factor(id)]
  model_fitting_data = na.omit(prepared_data)
  
  # O grep agora funciona corretamente pois os lags já existem
  all_cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  lagged_cov_names = paste0(all_cov_names, "_lag1")
  
  for (cov_name in covs_to_model) {
    base_formula = 's(t, k = 15)'
    treatment_effects_formula = paste(
      sapply(treatment_times, function(s) {
        paste0('s(t, bs="tp", k=8, by=I(A_s', s, ' * (t >= ', s, ')))')
      }), collapse = " + "
    )
    history_terms = paste(
      sapply(lagged_cov_names, function(lag_name) {
        paste0('s(', lag_name, ',  k = 8)')
      }), collapse = " + "
    )
    model_formula_str = paste(
      cov_name, "~", base_formula, "+", 
      treatment_effects_formula, "+", history_terms
    )
    
    # Usa os dados já limpos
    fitted_models[[cov_name]] = bam(
      as.formula(model_formula_str), data = model_fitting_data, 
      discrete = TRUE, select = TRUE
    )
  }
  return(fitted_models)
}


#' @title Ajusta modelo de desfecho para G-formula
fit_outcome_model_gformula = function(
    prepared_data, # Recebe dados já preparados
    treatment_times
) {
  prepared_data[, id := as.factor(id)]
  model_fitting_data = na.omit(prepared_data)
  
  cov_names = grep("^X[0-9]+$", names(model_fitting_data), value = TRUE)
  
  treatment_effects_formula = paste(
    sapply(treatment_times, function(s) {
      paste0('s(t, bs="tp", k=8, by=I(A_s', s, ' * (t >= ', s, ')))')
    }), collapse = " + "
  )
  covariate_effects = paste(
    sapply(cov_names, function(cov) paste0("s(t, by = ", cov, ", k=8)")),
    collapse = " + "
  )
  
  model_formula_str = paste(
    "y ~ s(t, k = 15) +", treatment_effects_formula, "+", covariate_effects
  )
  
  # Usa os dados já limpos
  outcome_model = bam(
    as.formula(model_formula_str), data = model_fitting_data, 
    discrete = TRUE, select = TRUE
  )
  return(outcome_model)
}


predict_potential_outcomes_gformula = function(
    outcome_model, simulated_covariates, treatment_times, regime
) {
  pred_data = copy(simulated_covariates)
  
  id_levels = levels(outcome_model$model$id)
  pred_data[, id := factor(id, levels = id_levels)]
  
  sim_names = names(pred_data)[grep("_sim$", names(pred_data))]
  new_names = gsub("_sim$", "", sim_names)
  setnames(pred_data, sim_names, new_names)
  
  for (i in seq_along(treatment_times)) {
    pred_data[, paste0("A_s", treatment_times[i]) := regime[i]]
  }
  
  predictions = predict(outcome_model, newdata = pred_data)
  return(data.table(id = pred_data$id, t = pred_data$t, y_pred = predictions))
}


# --- Função de "Entry Point" para o Pipeline de Avaliação ---

fit_g_formula = function(
    scenario_data, scenario_config, regime_a, regime_b
) {
  cat("--- Ajustando modelo: G-Formula Funcional ---\n")
  
  treatment_times = scenario_config$treatment_times
  
  # 1. Preparação dos dados (cria lags e colunas A_s)
  prepared_data = prepare_gformula_data(scenario_data, treatment_times)
  
  # 2. Ajuste dos modelos (Estágio I)
  cov_names = grep("^X[0-9]+$", names(scenario_data), value = TRUE)
  
  cov_models = fit_covariate_models_gformula(
    prepared_data, treatment_times, cov_names
  )
  outcome_model = fit_outcome_model_gformula(prepared_data, treatment_times)
  
  # 3. Simulação e Predição (Estágio II)
  
  # Simula trajetórias X para o regime A
  sim_covs_a = simulate_counterfactual_trajectories(
    fitted_models = cov_models,
    original_data = scenario_data,
    treatment_times = treatment_times,
    regime = regime_a
  )
  # Prediz o desfecho Y para o regime A
  pred_y_a = predict_potential_outcomes_gformula(
    outcome_model, sim_covs_a, treatment_times, regime_a
  )
  
  # Repete para o regime B
  sim_covs_b = simulate_counterfactual_trajectories(
    fitted_models = cov_models,
    original_data = scenario_data,
    treatment_times = treatment_times,
    regime = regime_b
  )
  pred_y_b = predict_potential_outcomes_gformula(
    outcome_model, sim_covs_b, treatment_times, regime_b
  )
  
  # 4. Calcula o efeito e retorna
  effect_dt = merge(
    pred_y_a, pred_y_b, by = c("id", "t"), suffixes = c("_a", "_b")
  )
  effect_dt[, effect_estimate := y_pred_a - y_pred_b]
  
  cat("--- Estimação com G-Formula concluída ---\n")
  return(effect_dt[, .(id, t, effect_estimate)])
}