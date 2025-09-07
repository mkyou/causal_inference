# ===================================================================
# G-FORMULA FUNCIONAL: SUÍTE COMPLETA DE VALIDAÇÃO (VERSÃO 1.2)
#
# - Modelo do desfecho atualizado para usar te(t, X) para maior
#   flexibilidade.
# - Lógica de retorno dos modelos de covariável corrigida.
# ===================================================================

library(mgcv)
library(ggplot2)
library(data.table)

# Carrega as funções de modelagem das covariáveis
source("scripts_functional/models/g_formula_covariates_model.R")
source("scripts_functional/models/g_formula_covariates_model_lag_alternative.R")

#' @title Ajusta o Modelo GAM para o Desfecho Y
#' @description Ajusta um modelo para Y em função do tempo, tratamentos e
#'   os efeitos concorrentes e funcionais das covariáveis.
#' @param dt_long O data.table completo no formato longo.
#' @param treatment_times Vetor numérico com os tempos de tratamento.
#' @return Um objeto de modelo 'bam'.
fit_outcome_model = function(dt_long, treatment_times) {
  
  cat("Ajustando o modelo para o desfecho Y...\n")
  
  # Prepara os dados (adiciona colunas A_s)
  prepared_data = prepare_gformula_data(dt_long, treatment_times)
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  
  # Efeitos funcionais do tratamento
  treatment_effects_formula = paste(
    sapply(treatment_times, function(s) {
      paste0("s(t, by = I(A_s", s, " * (t >= ", s, ")))")
    }),
    collapse = " + "
  )
  
  # Efeitos de interação 2D (tempo, covariável)
  covariate_effects_formula = paste(
    sapply(cov_names, function(cov) {
      paste0("s(", cov, ") + ti(", cov, ", t)") 
    }),
    collapse = " + "
  )
  
  # Fórmula completa
  model_formula_str = paste(
    "y ~ s(t, k = 15) +",
    treatment_effects_formula, "+",
    covariate_effects_formula
  )
  
  outcome_model = bam(
    as.formula(model_formula_str),
    data = prepared_data,
    discrete = TRUE,
    select = TRUE
  )
  
  cat("Ajuste do modelo de desfecho concluído.\n")
  return(outcome_model)
}

#' @title Prevê Desfechos Contrafactuais
#' @description Usa o modelo de desfecho e as covariáveis simuladas para
#' prever os desfechos potenciais.
#' @return Um data.table com id, t, e o desfecho previsto (y_pred).
predict_potential_outcomes = function(
    outcome_model,
    simulated_covariates,
    treatment_times,
    regime
) {
  
  pred_data = copy(simulated_covariates)
  
  # **INÍCIO DA CORREÇÃO**
  # A lógica de renomeação foi ajustada para remover os prefixos ('a_' ou
  # 'another_a_') além do sufixo '_pred', garantindo que as colunas
  # correspondam ao que o 'outcome_model' espera (ex: 'X1', 'X2').
  cov_names_pred = grep("_pred$", names(pred_data), value = TRUE)
  
  # Remove o prefixo (ex: 'a_X1_pred' -> 'X1_pred')
  cov_names_no_prefix = gsub("^(a_|another_a_)", "", cov_names_pred)
  # Remove o sufixo (ex: 'X1_pred' -> 'X1')
  cov_names_orig = gsub("_pred$", "", cov_names_no_prefix)
  
  setnames(pred_data, cov_names_pred, cov_names_orig)
  # **FIM DA CORREÇÃO**
  
  # Adiciona as colunas do regime de tratamento fixo
  for (i in seq_along(treatment_times)) {
    s = treatment_times[i]
    regime_val = regime[i]
    A_var_name = paste0("A_s", s)
    pred_data[, (A_var_name) := regime_val]
  }
  
  # Prevê o desfecho
  predictions = predict(outcome_model, newdata = pred_data)
  
  return(data.table(
    id = pred_data$id,
    t = pred_data$t,
    y_pred = predictions
  ))
}

#' @title Executa a Suíte Completa de Validação da G-Formula
#' @description Orquestra todo o processo de ajuste e validação.
#' @param covariate_model_type O tipo de modelo para as covariáveis:
#'   "functional" ou "autoregressive".
#' @return Uma lista com métricas e um gráfico de validação.
run_full_gformula_suite = function(
    scenario_data, 
    scenario_config,
    covariate_model_type = "functional"
) {
  
  dt_long = scenario_data$outcome$long
  
  # --- PASSO 1: AJUSTAR MODELOS DAS COVARIÁVEIS E SIMULAR ---
  cat(paste0(
    "\n--- Iniciando G-Formula com modelo de covariáveis: '",
    covariate_model_type, "' ---\n"
  ))
  
  if (covariate_model_type == "functional") {
    cov_results = run_gformula_validation_suite(
      scenario_data, scenario_config
    )
  } else if (covariate_model_type == "autoregressive") {
    cov_results = run_autoregressive_validation_suite(
      scenario_data, scenario_config
    )
  } else {
    stop("Tipo de modelo de covariável desconhecido:", covariate_model_type)
  }
  
  simulated_a = cov_results$simulations$a_seq
  simulated_another_a = cov_results$simulations$another_a_seq
  
  # --- PASSO 2: AJUSTAR MODELO DO DESFECHO ---
  outcome_model = fit_outcome_model(
    dt_long, 
    scenario_config$treatment_times
  )
  
  # --- PASSO 3: PREVER DESFECHOS CONTRAFACTUAIS ---
  cat("\nPrevendo desfechos contrafactuais...\n")
  
  y_pred_a = predict_potential_outcomes(
    outcome_model, 
    simulated_a,
    scenario_config$treatment_times,
    scenario_config$a_seq
  )
  
  y_pred_another_a = predict_potential_outcomes(
    outcome_model, 
    simulated_another_a,
    scenario_config$treatment_times,
    scenario_config$another_a_seq
  )
  
  # --- PASSO 4: ESTIMAR E VALIDAR O EFEITO CAUSAL ---
  cat("Calculando e validando o efeito causal...\n")
  
  effect_data = merge(y_pred_a, y_pred_another_a, by = c("id", "t"))
  setnames(effect_data, c("y_pred.x", "y_pred.y"), c("y_a", "y_another_a"))
  
  effect_data[, effect_est := y_a - y_another_a]
  estimated_effect = effect_data[, .(estimated = mean(effect_est)), by = t]
  
  true_effect = data.table(
    t = scenario_data$outcome$wide$t_vec,
    real = colMeans(scenario_data$outcome$wide$effect_a_vs_another)
  )
  
  validation_dt = merge(estimated_effect, true_effect, by = "t")
  
  metrics = data.table(
    covariate_model = covariate_model_type,
    RMSE = sqrt(mean((validation_dt$estimated - validation_dt$real)^2)),
    Bias = mean(validation_dt$estimated - validation_dt$real),
    Correlation = cor(validation_dt$estimated, validation_dt$real)
  )
  
  plot_title = "Validação do Efeito Causal Estimado via G-Formula"
  subtitle_text = paste0(
    "Modelo Covariáveis: '", covariate_model_type, "' | ",
    "RMSE: ", round(metrics$RMSE, 3), " | ",
    "Bias: ", round(metrics$Bias, 3), " | ",
    "Corr: ", round(metrics$Correlation, 3)
  )
  
  validation_plot = ggplot(validation_dt, aes(x = t)) +
    geom_line(
      aes(y = real, color = "Efeito Causal Real"),
      linetype = "dashed", size = 1.1
    ) +
    geom_line(
      aes(y = estimated, color = "Efeito Causal Estimado"),
      size = 1.1, alpha = 0.8
    ) +
    scale_color_manual(
      values = c(
        "Efeito Causal Real" = "black",
        "Efeito Causal Estimado" = "#E69F00"
      )
    ) +
    labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Tempo", y = "Efeito Causal Médio", color = ""
    ) +
    theme_minimal(base_size = 14)
  
  cat("--- Validação da G-Formula Concluída ---\n")
  
  return(list(
    metrics = metrics,
    plot = validation_plot,
    plot_data = validation_dt
  ))
}

