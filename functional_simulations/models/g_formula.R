# g_formula.R (CORRIGIDO)
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
  
  # --- ETAPA 1: SIMULAÇÃO DAS COVARIÁVEIS (AGORA CORRETA) ---
  cat("  - Etapa 1: Simulando trajetórias das covariáveis sequencialmente...\n")
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
      function(s) paste0("s(t, by = A_s", s, ", k = 20)")
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
  prediction_data = merge(
    prepared_data, 
    covariate_sims_list$simulated_data, 
    by = c("id", "t")
  )
  
  # Prepara newdata para predição (substituindo X pelos X simulados)
  newdata_a = copy(prediction_data); setnames(
    newdata_a, paste0("sim_a_", cov_names), cov_names
  )
  newdata_b = copy(prediction_data); setnames(
    newdata_b, paste0("sim_b_", cov_names), cov_names
  )
  
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