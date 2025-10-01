library(mgcv)
library(data.table)

#' @title Ajusta um MSM "Semi-Ingênuo"
#' @description Implementa um Marginal Structural Model com um modelo de
#'   propensão linear (ingênuo) e um modelo de desfecho marginal
#'   funcional (com splines).
#'
#' @param prepared_data O data.table de entrada.
#' @param scenario_config A configuração do cenário.
#' @param regime_a,regime_b Vetores definindo os regimes a comparar.
#'
#' @return Uma lista com os resultados, modelos ajustados e fórmulas usadas.

run_fmsm = function(
    prepared_data,
    scenario_config,
    regime_a,
    regime_b,
    ...
) {
  cat(
    "--- Iniciando Semi-Naive MSM (Propensão Linear, Desfecho Funcional) ---\n"
  )
  
  # --- Setup ---
  dt = copy(prepared_data)
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  
  # --- ETAPA 1: AJUSTE DOS MODELOS DE PROPENSÃO ---
  cat("  - Etapa 1: Ajustando modelos de propensão com splines e seleção...\n")
  propensity_models = list()
  history_formula = paste(lagged_names, collapse = " + ")
  ps_formula = as.formula(paste("A ~", history_formula))
  
  cat("    > Fórmula de Propensão:\n")
  print(ps_formula)
  
  treatment_times = scenario_config$treatment_times
  for (s_time in treatment_times) {
    data_at_s = dt[t == s_time]
    
    # MUDANÇA DEFINITIVA: Criar um data.frame limpo para o modelo
    # Apenas com a resposta (A) e os preditores permitidos (lags).
    # Isso impede fisicamente qualquer vazamento de dados das colunas A_s*.
    model_df = data_at_s[, c("A", lagged_names), with = FALSE]
    model_df = na.omit(model_df)
    
    ps_model = bam(
      ps_formula, 
      data = model_df, # Usa o data.frame limpo
      family = binomial(link = "logit"),
      select = TRUE, 
      discrete = TRUE
    )
    propensity_models[[as.character(s_time)]] = ps_model
  }
  
  # --- ETAPA 2: CÁLCULO DOS PESOS ESTABILIZADOS ---
  cat("  - Etapa 2: Calculando pesos estabilizados (SW)...\n")
  weight_component_cols = paste0("sw_comp_s", treatment_times)
  for (i in 1:length(treatment_times)) {
    s_time = treatment_times[i]
    col_name = weight_component_cols[i]
    ps_model = propensity_models[[as.character(s_time)]]
    
    # A variável 'dt' já está definida e é a correta a ser usada aqui
    data_at_s = dt[t == s_time]
    valid_data_at_s = na.omit(data_at_s, cols = c("t", lagged_names))
    
    if (nrow(valid_data_at_s) > 0) {
      p_denom = predict(ps_model, newdata = valid_data_at_s, type = "response")
      prob_received_denom = ifelse(valid_data_at_s$A == 1, p_denom, 1-p_denom)
      p_num = mean(valid_data_at_s$A)
      prob_received_num = ifelse(valid_data_at_s$A == 1, p_num, 1-p_num)
      
      prob_received_denom[prob_received_denom < 1e-6] = 1e-6
      
      valid_data_at_s[, (col_name) := prob_received_num/prob_received_denom]
      dt = merge(
        dt, valid_data_at_s[, c("id", "t", col_name), with = FALSE],
        by = c("id", "t"), all.x = TRUE
      )
    }
  }
  dt[, (weight_component_cols) := lapply(.SD, zoo::na.locf, na.rm = FALSE),
     by = id, .SDcols = weight_component_cols]
  for (col in weight_component_cols) { dt[is.na(get(col)), (col) := 1.0] }
  dt[, sw := Reduce(`*`, .SD), .SDcols = weight_component_cols]
  dt[, (weight_component_cols) := NULL]
  
  quantile_trunc = quantile(dt$sw, probs = c(0.01, 0.99), na.rm = TRUE)
  dt[sw < quantile_trunc[1], sw := quantile_trunc[1]]
  dt[sw > quantile_trunc[2], sw := quantile_trunc[2]]
  
  # --- ETAPA 3: AJUSTE DO MODELO MARGINAL PONDERADO (FUNCIONAL) ---
  cat("  - Etapa 3: Ajustando o modelo de desfecho marginal ponderado...\n")
  
  # A fórmula usa as colunas 'A_s...' que já existem em 'dt'
  treat_form_y = paste(
    sapply(
      treatment_times,
      function(s) paste0("s(t, by = A_s", s, ", k = 20)")
    ),
    collapse = " + "
  )
  outcome_formula = as.formula(paste("y ~ s(t, k = 30) +", treat_form_y))
  
  cat("    > Fórmula do Desfecho:\n")
  print(outcome_formula)
  
  msm_outcome_model = bam(
    outcome_formula, data = na.omit(dt), weights = sw,
    select = TRUE, discrete = TRUE
  )
  
  # --- ETAPA 4: PREDIÇÃO DOS DESFECHOS CONTRAFRACTUAIS ---
  cat("  - Etapa 4: Prevendo desfechos contrafactuais...\n")
  
  newdata_a = copy(dt)
  newdata_b = copy(dt)
  for (i in 1:length(treatment_times)) {
    s_time = treatment_times[i]
    A_s_col = paste0("A_s", s_time)
    
    newdata_a[, (A_s_col) := regime_a[i]]
    newdata_b[, (A_s_col) := regime_b[i]]
  }
  
  dt[, y_a_msm := predict(msm_outcome_model, newdata = newdata_a)]
  dt[, y_b_msm := predict(msm_outcome_model, newdata = newdata_b)]
  dt[, effect_estimate_msm := y_a_msm - y_b_msm]
  
  # --- ETAPA 5: EMPACOTAR RESULTADOS ---
  return(
    list(
      results_data = dt,
      fitted_models = list(
        propensity_models = propensity_models,
        outcome_model = msm_outcome_model
      ),
      formulas = list(
        propensity_formula = ps_formula,
        outcome_formula = outcome_formula
      )
    )
  )
}