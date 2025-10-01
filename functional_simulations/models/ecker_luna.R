# ===================================================================
# MODELO 1: Ecker & Luna
# ===================================================================
library(mgcv)
library(data.table)

#' @title Ajusta o modelo Ecker & Luna usando GAM
#' @description Recebe dados processados, converte-os para o formato longo e
#'   ajusta o modelo Ecker & Luna, que estratifica pela história de
#'   tratamento e controla por covariáveis basais.
#'
#' @param scenario_data A lista de dados processados.
#' @param scenario_config A configuração do cenário.
#' @param regime_a,regime_b Vetores definindo os regimes a comparar.
#'
#' @return Um data.table com id, t, e a estimativa do efeito causal.

fit_ecker_luna_mgcv = function(
    scenario_data, scenario_config, regime_a, regime_b, verbose = FALSE
) {
  if (verbose) cat("--- Ajustando modelo: Ecker & Luna (mgcv) ---\n")
  
  model_data = revert_to_long_format(scenario_data)
  treatment_times = scenario_config$treatment_times
  
  regime_a_label = paste(regime_a, collapse = "_")
  regime_b_label = paste(regime_b, collapse = "_")
  
  # Identifica a história de tratamento completa para cada indivíduo
  treatment_history_dt = model_data[
    t %in% treatment_times,
    .(treatment_history = paste(A, collapse = "_")),
    by = id
  ]
  
  # Identifica o tempo basal (último ponto antes do primeiro tratamento)
  first_treatment_time = min(treatment_times)
  all_times = sort(unique(model_data$t))
  baseline_time = max(all_times[all_times < first_treatment_time])
  
  # Extrai os valores das covariáveis no tempo basal
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  base_cov_names = paste0(cov_names, "_base")
  
  baseline_data = model_data[
    t == baseline_time, .SD, .SDcols = c("id", cov_names)
  ]
  setnames(baseline_data, cov_names, base_cov_names)
  
  # Junta a história de tratamento e os dados basais ao dataset principal
  model_data = merge(model_data, baseline_data, by = "id")
  model_data = merge(model_data, treatment_history_dt, by = "id")
  
  # --- 3. Ajuste dos Modelos por Regime ---
  
  # A fórmula usa termos de interação entre o tempo e as covariáveis basais
  smooth_terms = paste0("s(t, by = ", base_cov_names, ", k = 15)")
  model_formula = as.formula(
    paste("y ~ s(t, k = 30) +", paste(smooth_terms, collapse = " + "))
  )
  
  if (verbose) {
    cat("    > Fórmula do Desfecho (para cada estrato de tratamento):\n")
    print(model_formula)
  }
  
  # Filtra os dados para cada regime de tratamento
  data_a = model_data[treatment_history == regime_a_label]
  data_b = model_data[treatment_history == regime_b_label]
  
  if (nrow(data_a) < 20 || nrow(data_b) < 20) {
    stop("Dados insuficientes para um ou ambos os regimes de tratamento.")
  }
  
  # Ajusta um modelo GAM para cada regime
  fit_a = bam(model_formula, data = data_a, discrete = TRUE)
  fit_b = bam(model_formula, data = data_b, discrete = TRUE)
  
  # --- 4. Predição e Retorno ---
  
  # Prevê o desfecho para todos os indivíduos sob cada modelo
  pred_Y_a = predict(fit_a, newdata = model_data)
  pred_Y_b = predict(fit_b, newdata = model_data)
  
  model_data[, effect_estimate := pred_Y_a - pred_Y_b]
  
  return(model_data)
}

run_ecker_luna = function(
    scenario_object, prepared_data, scenario_config, regime_a, regime_b
) {
  # Esta função ignora 'prepared_data' e usa os dados processados originais
  estimates_dt = fit_ecker_luna_mgcv(
    scenario_data = scenario_object$data,
    scenario_config = scenario_config,
    regime_a = regime_a,
    regime_b = regime_b,
    verbose = FALSE
  )
  
  # O 'prepared_data' é usado aqui apenas como esqueleto para o merge
  full_results_data = merge(prepared_data, estimates_dt, by = c("id", "t"))
  
  return(list(results_data = full_results_data))
}

assess_ecker_luna_process = function(model_results, true_data_refund) {
  cat("--- Iniciando Avaliação Simplificada para Ecker & Luna ---\n")
  
  true_data_long = revert_to_long_format(true_data_refund)
  setnames(true_data_long, "effect_a_vs_b", "effect_a_vs_b_true")
  
  results_with_truth = merge(
    model_results$results_data,
    true_data_long[, .(id, t, effect_a_vs_b_true)],
    by = c("id", "t")
  )
  
  setnames(results_with_truth, "effect_estimate", "effect_estimate_msm")
  setnames(results_with_truth, "effect_estimate_msm", "effect_estimate")
  
  assessment_causal = assess_effect_estimation(results_with_truth)
  
  return(list(causal_effect_assessment = assessment_causal))
}
