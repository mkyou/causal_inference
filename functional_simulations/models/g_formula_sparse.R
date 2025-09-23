library(mgcv)
library(data.table)

source("functional_simulations/models/g_formula_sparse_covariates_model.R")

is_constant_func = function(f) {
  isTRUE(all.equal(f(5), f(15)))
}

#' @title Executa a G-formula com modelos esparsos e informados pelo PGD
run_gformula_sparse = function(
    prepared_data,
    scenario_config,
    regime_a,
    regime_b
) {
  
  cat("--- Iniciando G-Formula (Estratégia ESPARSA) ---\n")
  
  # ETAPA 1: SIMULAÇÃO DAS COVARIÁVEIS (ESPARSA)
  covariate_sims_list = run_gformula_covariates_sparse(
    prepared_data, scenario_config, regime_a, regime_b
  )
  
  # ETAPA 2: AJUSTE DO MODELO DE DESFECHO (ESPARSO)
  cat("  - Ajustando modelo de desfecho (ESPARSO)...\n")
  fit_data_y = na.omit(prepared_data)
  fit_data_y[, id := as.factor(id)]
  
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(fit_data_y), value = TRUE)
  
  # Constrói a fórmula do Desfecho dinamicamente
  concurrent_funcs = scenario_config$concurrent_effect_funcs
  covariate_terms = c()
  for (k_idx in 1:length(cov_names)) {
    k_name = cov_names[k_idx]
    func = concurrent_funcs[[k_idx]]
    
    if (is_constant_func(func)) {
      term = k_name
    } else {
      term = paste0("te(t, ", k_name, ", k = c(15, 15))")
    }
    covariate_terms = c(covariate_terms, term)
  }
  cov_form_y = paste(covariate_terms, collapse = " + ")
  
  # A fórmula do tratamento é a mesma que já havíamos corrigido
  treat_form_y = paste(sapply(treatment_times, function(s) {
    paste0("s(time_since_s", s, ", by = by_s", s, ", k = 20)")
  }), collapse = " + ")
  
  base_formula_y = "s(t, k = 30)"
  model_formula_y = as.formula(
    paste("y ~", base_formula_y, "+", cov_form_y, "+", treat_form_y)
  )
  
  cat(
    "Fórmula de estimação do modelo de desfecho: ", 
    as.character(model_formula_y),
    "\n"
  )
  
  outcome_model = bam(
    model_formula_y, data = fit_data_y, discrete = TRUE, select = TRUE
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
  
  return(
    list(
      results_data = results_data,
      fitted_models = list(
        covariate_models = covariate_sims_list$models,
        outcome_model = outcome_model
      ),
      formulas = list(
        outcome_model_formula = model_formula_y
      )
    )
  )
}