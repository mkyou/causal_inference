# ===================================================================
# G-FORMULA: ESTÁGIO DE MODELAGEM DAS COVARIÁVEIS
#
# Este script contém as funções para:
# 1. Preparar os dados (lags e colunas de tratamento).
# 2. Ajustar os modelos GAM para cada covariável.
# 3. Simular as trajetórias contrafactuais das covariáveis.
# ===================================================================

library(mgcv)
library(data.table)

#' @title Prepara Dados para Modelos de Covariáveis
prepare_covariate_data = function(factual_data, treatment_times) {
  dt = copy(factual_data)
  setkey(dt, id, t)
  
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  
  dt[, (lagged_names) := shift(.SD, n = 1, type = "lag"), 
     by = id, .SDcols = cov_names]
  
  # Itera sobre os tempos de tratamento para criar as colunas A_s
  for (s_time in treatment_times) {
    A_s_name = paste0("A_s", s_time)
    
    # Pega o status do tratamento no tempo exato
    treatment_status_at_s = factual_data[t == s_time, .(id, A)]
    
    # --- LÓGICA DE MERGE CORRIGIDA E ROBUSTA ---
    # Renomeia a coluna 'A' ANTES do merge para evitar conflitos de nome
    setnames(treatment_status_at_s, "A", A_s_name)
    
    # Junta a nova coluna A_s ao data.table principal.
    # Como não há mais conflito de nomes, o merge é limpo.
    dt = merge(dt, treatment_status_at_s, by = "id", all.x = TRUE)
  }
  
  return(dt)
}

#' @title Ajusta Modelos GAM para as Covariáveis
#' @description Ajusta um modelo GAM para cada covariável, usando uma
#'   especificação de fórmula funcional e numericamente estável.
#'
#' @param prepared_data O data.table da função `prepare_covariate_data`,
#'   já contendo as colunas de lag e de tratamento (`A_s`).
#' @param treatment_times Vetor numérico com os tempos de tratamento.
#'
#' @return Uma lista nomeada com os modelos GAM ajustados para cada covariável.

fit_covariate_models = function(prepared_data, treatment_times) {
  fitted_models = list()
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  
  fit_data = na.omit(prepared_data)
  fit_data[, id := as.factor(id)]
  
  for (k_name in cov_names) {
    base_formula = 's(t, k = 15)'
    treatment_formula = paste(sapply(treatment_times, function(s) {
      paste0('s(t, bs="tp", k=8, by=I(A_s', s, ' * (t >= ', s, ')))')
    }), collapse = " + ")
    
    # Lógica Adaptativa para a História
    history_terms = c()
    for (lag in lagged_names) {
      if (var(fit_data[[lag]], na.rm = TRUE) < 1e-6) {
        history_terms = c(history_terms, lag)
      } else {
        history_terms = c(
          history_terms, paste0('s(t, by=', lag, ', k=8)')
        )
      }
    }
    history_formula = paste(history_terms, collapse = " + ")
    
    model_formula = as.formula(
      paste(k_name, "~", base_formula, "+", 
            treatment_formula, "+", history_formula)
    )
    
    fitted_models[[k_name]] = bam(
      model_formula, 
      data = fit_data,
      discrete = TRUE # Reativado conforme sua preferência
    )
  }
  return(fitted_models)
}


#' @title (Refatorada) Simula Trajetórias Contrafactuais de Covariáveis
#' @description Usando modelos ajustados, simula sequencialmente as
#'   trajetórias das covariáveis sob um regime de tratamento fixo.
#'
#' @param fitted_models Uma lista nomeada de modelos GAM ajustados.
#' @param factual_data O data.table contendo os dados factuais, usado
#'   para obter a lista de indivíduos e os valores iniciais (t=0).
#' @param treatment_times Vetor numérico com os tempos de tratamento.
#' @param regime Vetor numérico definindo o regime de tratamento fixo.
#'
#' @return Um data.table com `id`, `t`, e as trajetórias simuladas
#'   (ex: `X1_sim`, `X2_sim`).
simulate_covariate_trajectories = function(
    fitted_models,
    factual_data,
    treatment_times,
    regime
) {
  cov_names = names(fitted_models)
  setkey(factual_data, id, t)
  ids = unique(factual_data$id)
  id_levels = levels(as.factor(ids))
  time_grid = sort(unique(factual_data$t))
  
  simulated_dt = CJ(id = ids, t = time_grid)
  for (k in cov_names) {
    simulated_dt[, paste0(k, "_sim") := NA_real_]
  }
  
  sim_names = paste0(cov_names, "_sim")
  first_time_data = factual_data[t == time_grid[1], 
                                 .SD, .SDcols = c("id", cov_names)]
  setnames(first_time_data, cov_names, sim_names)
  first_time_data[, t := time_grid[1]]
  
  simulated_dt[, id := as.character(id)]
  first_time_data[, id := as.character(id)]
  simulated_dt[first_time_data, on = .(id, t), (sim_names) := mget(
    paste0("i.", sim_names)
  )]
  
  for (i in 2:length(time_grid)) {
    current_t = time_grid[i]
    previous_t = time_grid[i-1]
    
    pred_data_t = data.table(id = ids, t = current_t)
    
    for (j in seq_along(treatment_times)) {
      s_time = treatment_times[j]
      pred_data_t[, paste0("A_s", s_time) := regime[j]]
    }
    
    lagged_names = paste0(cov_names, "_lag1")
    history_data = simulated_dt[t == previous_t, 
                                .SD, .SDcols = c("id", sim_names)]
    setnames(history_data, sim_names, lagged_names)
    
    pred_data_t[, id := as.character(id)]
    history_data[, id := as.character(id)]
    
    pred_data_t = merge(pred_data_t, history_data, by = "id")
    pred_data_t[, id := factor(id, levels = id_levels)]
    
    for (k in cov_names) {
      predictions = predict(
        fitted_models[[k]], 
        newdata = as.data.frame(pred_data_t)
      )
      simulated_dt[t == current_t, (paste0(k, "_sim")) := predictions]
    }
  }
  
  final_cols = c("id", "t", paste0(cov_names, "_sim"))
  simulated_dt[, id := as.numeric(as.character(id))]
  
  return(simulated_dt[, ..final_cols])
}


#' @title Ponto de Entrada para o Estágio das Covariáveis da G-Formula
#' @description Versão com a lógica de renomeação de colunas corrigida.
run_gformula_covariates = function(
    scenario_data, scenario_config, regime_a, regime_b
) {
  cat("--- Iniciando Estágio de Covariáveis da G-Formula ---\n")
  
  factual_data = scenario_data[, .SD, .SDcols = c(
    "id", "t", "A", grep("^X[0-9]+$", names(scenario_data), value = TRUE)
  )]
  
  prepared_data = prepare_covariate_data(
    factual_data, scenario_config$treatment_times
  )
  
  covariate_models = fit_covariate_models(
    prepared_data, scenario_config$treatment_times
  )
  
  cat("  - Simulando trajetórias para o Regime A...\n")
  sim_a = simulate_covariate_trajectories(
    covariate_models, factual_data, scenario_config$treatment_times, regime_a
  )
  cat("  - Simulando trajetórias para o Regime B...\n")
  sim_b = simulate_covariate_trajectories(
    covariate_models, factual_data, scenario_config$treatment_times, regime_b
  )
  
  sim_a_old_names = grep("_sim$", names(sim_a), value = TRUE)
  sim_a_new_names = paste0("sim_a_", gsub("_sim$", "", sim_a_old_names))
  setnames(sim_a, sim_a_old_names, sim_a_new_names)
  
  sim_b_old_names = grep("_sim$", names(sim_b), value = TRUE)
  sim_b_new_names = paste0("sim_b_", gsub("_sim$", "", sim_b_old_names))
  setnames(sim_b, sim_b_old_names, sim_b_new_names)
  
  results_dt = merge(scenario_data, sim_a, by = c("id", "t"))
  results_dt = merge(results_dt, sim_b, by = c("id", "t"))
  
  cat("--- Estágio de Covariáveis da G-Formula Concluído ---\n")
  return(results_dt)
}