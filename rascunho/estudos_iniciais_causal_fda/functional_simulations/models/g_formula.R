library(mgcv)
library(data.table)

source("functional_simulations/models/g_formula_covariates_model.R")

run_gformula_full = function(
    prepared_data,
    scenario_config,
    regime_a,
    regime_b
) {
  
  cat("--- Iniciando G-Formula Completa ---\n")
  
  # --- ETAPA 1: SIMULAÇÃO DAS COVARIÁVEIS ---
  cat("  - Etapa 1: Simulando trajetórias das covariáveis...\n")
  covariate_sims_list = run_gformula_covariates_for_simulation(
    prepared_data, scenario_config, regime_a, regime_b
  )
  
  # --- ETAPA 2: AJUSTE DO MODELO DE DESFECHO ---
  cat("  - Etapa 2: Ajustando modelo de desfecho...\n")
  fit_data_y = na.omit(prepared_data)
  fit_data_y[, id := as.factor(id)]
  
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(fit_data_y), value = TRUE)
  base_formula_y = "s(t, k = 30)"
  
  cov_form_y = paste(
    sapply(
      cov_names,
      function(k) paste0("s(", k, ") + ti(t, ", k, ")")
    ),
    collapse = " + "
  )
  
  treat_form_y = paste(
    sapply(
      treatment_times,
      function(s) {
        # A fórmula agora usa a variável 'by' causalmente correta
        paste0("s(time_since_s", s, ", by = by_s", s, ", k = 20)")
      }
    ),
    collapse = " + "
  )
  
  model_formula_y = as.formula(
    paste("y ~", base_formula_y, "+", cov_form_y, "+", treat_form_y)
  )
  
  outcome_model = bam(
    model_formula_y,
    data = fit_data_y,
    discrete = TRUE,
    select = TRUE
  )
  
  # --- ETAPA 3: PREDIÇÃO DO DESFECHO ---
  cat("  - Etapa 3: Prevendo desfechos contrafactuais...\n")
  
  prepared_data[, id := as.factor(id)]
  prediction_data = merge(
    prepared_data,
    covariate_sims_list$simulated_data,
    by = c("id", "t")
  )
  
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  a_s_cols = paste0("A_s", treatment_times)
  
  newdata_a = copy(prediction_data)
  setnames(newdata_a, paste0("sim_a_", cov_names), cov_names)
  # 1. Força as colunas de tratamento para o regime A
  for (i in 1:length(a_s_cols)) {
    set(newdata_a, j = a_s_cols[i], value = regime_a[i])
  }
  for (s_time in treatment_times) {
    by_col = paste0("by_s", s_time)
    A_s_col = paste0("A_s", s_time)
    newdata_a[
      , 
      (by_col) := .SD[[A_s_col]] * fifelse(t >= s_time, 1, 0), 
      .SDcols = A_s_col
    ]
  }
  
  newdata_b = copy(prediction_data)
  setnames(newdata_b, paste0("sim_b_", cov_names), cov_names)
  # 1. Força as colunas de tratamento para o regime B
  for (i in 1:length(a_s_cols)) {
    set(newdata_b, j = a_s_cols[i], value = regime_b[i])
  }
  for (s_time in treatment_times) {
    by_col = paste0("by_s", s_time)
    A_s_col = paste0("A_s", s_time)
    newdata_b[
      , 
      (by_col) := .SD[[A_s_col]] * fifelse(t >= s_time, 1, 0), 
      .SDcols = A_s_col
    ]
  }
  
  Y_a_sim = predict(outcome_model, newdata = newdata_a)
  Y_b_sim = predict(outcome_model, newdata = newdata_b)
  
  # --- ETAPA 4: EMPACOTANDO OS RESULTADOS ---
  cat("  - Etapa 4: Empacotando resultados...\n")
  results_data = prediction_data
  results_data[, y_a_sim := Y_a_sim]
  results_data[, y_b_sim := Y_b_sim]
  results_data[, effect_estimate := Y_a_sim - Y_b_sim]
  
  return(list(
    results_data = results_data,
    fitted_models = list(
      covariate_models = covariate_sims_list$models,
      outcome_model = outcome_model
    ),
    formulas = list(
      outcome_model_formula = model_formula_y
    )
  ))
}