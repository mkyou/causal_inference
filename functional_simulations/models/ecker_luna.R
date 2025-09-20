# ===================================================================
# MODELO 1: Ecker & Luna
# ===================================================================
library(mgcv)
library(data.table)

source("functional_simulations/utils/revert_to_long_format.R")


#' @title Ajusta o modelo Ecker & Luna usando GAM
#' @description Recebe dados no formato `refund`, converte-os de volta para
#'   longo, e então ajusta o modelo Ecker & Luna com `mgcv::bam`.
#'
#' @param scenario_data A lista de dados gerada por `prepare_refund_data`.
#' @param scenario_config A configuração do cenário.
#' @param regime_a,regime_b Vetores definindo os regimes a comparar.
#'
#' @return Um data.table com `id`, `t`, e `effect_estimate`.

fit_ecker_luna_mgcv = function(
    scenario_data, scenario_config, regime_a, regime_b
) {
  cat(
    "--- Ajustando modelo: Ecker & Luna (mgcv com dados pré-suavizados) ---\n"
  )
  
  # --- 1. Etapa de Conversão de Dados ---
  # Converte os dados do formato refund/largo de volta para o formato longo
  model_data = revert_to_long_format(scenario_data)
  
  # --- 2. Preparação dos Dados (lógica de Ecker & Luna) ---
  treatment_times = scenario_config$treatment_times
  
  regime_a_label = paste(regime_a, collapse = "_")
  regime_b_label = paste(regime_b, collapse = "_")
  
  treatment_history_dt = model_data[
    t %in% treatment_times, 
    .(treatment_history = paste(A, collapse = "_")), 
    by = id
  ]
  
  first_treatment_time = min(treatment_times)
  all_times = sort(unique(model_data$t))
  baseline_time = max(all_times[all_times < first_treatment_time])
  
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  base_cov_names = paste0(cov_names, "_base")
  
  baseline_data = model_data[
    t == baseline_time, .SD, .SDcols = c("id", cov_names)
  ]
  setnames(baseline_data, cov_names, base_cov_names)
  
  model_data = merge(model_data, baseline_data, by = "id")
  model_data = merge(model_data, treatment_history_dt, by = "id")
  
  # --- 3. Ajuste dos Modelos por Regime ---
  
  # A fórmula usa a sintaxe s(t, by = X_base) que você pediu
  smooth_terms = paste0("s(t, by = ", base_cov_names, ", k = 15)")
  model_formula = as.formula(
    paste("y ~ s(t, k = 30) +", paste(smooth_terms, collapse = " + "))
  )
  
  data_a = model_data[treatment_history == regime_a_label]
  data_b = model_data[treatment_history == regime_b_label]
  
  if(nrow(data_a) < 10 || nrow(data_b) < 10) {
    stop("Dados insuficientes para um ou ambos os regimes.")
  }
  
  fit_a = bam(model_formula, data = data_a, discrete = TRUE)
  fit_b = bam(model_formula, data = data_b, discrete = TRUE)
  
  # --- 4. Predição e Retorno ---
  
  pred_Y_a = predict(fit_a, newdata = model_data)
  pred_Y_b = predict(fit_b, newdata = model_data)
  
  model_data[, effect_estimate := pred_Y_a - pred_Y_b]
  
  cat("--- Estimação Ecker & Luna (mgcv) concluída ---\n")
  return(model_data[, .(id, t, effect_estimate)])
}