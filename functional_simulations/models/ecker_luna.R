# ===================================================================
# MODELO 1: Ecker & Luna
# ===================================================================
library(mgcv)
library(data.table)

#' @title Ajusta e Estima o Efeito Causal usando o Método Ecker & Luna
#' @description Versão refatorada para ser compatível com o pipeline de
#'   simulação. A preparação dos dados e o ajuste do modelo são
#'   encapsulados em uma única função.
#'
#' @param scenario_data O data.table completo gerado pelo PGD.
#' @param scenario_config A configuração do cenário.
#' @param regime_a,regime_b Vetores definindo os regimes a comparar.
#'
#' @return Um data.table com as colunas `id`, `t`, e `effect_estimate`.

fit_ecker_luna = function(
    scenario_data, scenario_config, regime_a, regime_b
) {
  cat("--- Ajustando modelo: Ecker & Luna ---\n")
  
  # --- 1. Preparação dos Dados (lógica agora interna) ---
  
  # A função agora espera o data.table diretamente
  model_data = copy(scenario_data)
  treatment_times = scenario_config$treatment_times
  
  regime_a_label = paste(regime_a, collapse = "_")
  regime_b_label = paste(regime_b, collapse = "_")
  
  # Cria o histórico de tratamento (ex: "1_1_1")
  treatment_history_dt = model_data[
    t %in% treatment_times, 
    .(treatment_history = paste(A, collapse = "_")), 
    by = id
  ]
  
  # Extrai as covariáveis de baseline (pré-primeiro tratamento)
  first_treatment_time = min(treatment_times)
  all_times = sort(unique(model_data$t))
  baseline_time = max(all_times[all_times < first_treatment_time])
  
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  base_cov_names = paste0(cov_names, "_base")
  
  baseline_data = model_data[
    t == baseline_time, .SD, .SDcols = c("id", cov_names)
  ]
  setnames(baseline_data, cov_names, base_cov_names)
  
  # Junta a história e a baseline aos dados principais
  model_data = merge(model_data, baseline_data, by = "id", all.x = TRUE)
  model_data = merge(model_data, treatment_history_dt, by = "id")
  
  # --- 2. Ajuste dos Modelos por Regime ---
  
  # Define a fórmula do modelo
  smooth_terms = paste0("s(t, by = ", base_cov_names, ", k = 8)")
  model_formula = as.formula(
    paste("y ~ s(t, k = 15) +", paste(smooth_terms, collapse = " + "))
  )
  
  # Filtra os dados para cada regime de interesse
  data_a = model_data[treatment_history == regime_a_label]
  data_b = model_data[treatment_history == regime_b_label]
  
  # Validação
  if(nrow(data_a) < 10 || nrow(data_b) < 10) {
    stop("Dados insuficientes para um ou ambos os regimes.")
  }
  
  # Ajusta os modelos
  fit_a = bam(model_formula, data = data_a, discrete = TRUE)
  fit_b = bam(model_formula, data = data_b, discrete = TRUE)
  
  # --- 3. Predição Contrafactual ---
  
  # Prevê o desfecho para TODOS os indivíduos usando cada modelo
  pred_Y_a = predict(fit_a, newdata = model_data)
  pred_Y_b = predict(fit_b, newdata = model_data)
  
  # --- 4. Retorno ---
  
  # Calcula o efeito e retorna um data.table limpo
  model_data[, effect_estimate := pred_Y_a - pred_Y_b]
  
  cat("--- Estimação Ecker & Luna concluída ---\n")
  return(model_data[, .(id, t, effect_estimate)])
}