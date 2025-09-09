
library(mgcv)
library(ggplot2)
library(data.table)

# Carrega as funções de modelagem das covariáveis
source("scripts_functional/models/g_formula_covariates_model.R")

#' @title Ajusta o Modelo GAM para o Desfecho Y
fit_outcome_model = function(dt_long, treatment_times, return_prepared_data = FALSE) {
  cat("Ajustando o modelo para o desfecho Y...\n")
  prepared_data = prepare_gformula_data(dt_long, treatment_times)
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  
  treatment_effects_formula = paste(
    sapply(treatment_times, function(s) paste0("s(t, by = I(A_s", s, " * (t >= ", s, ")))")),
    collapse = " + "
  )
  covariate_main_effects = paste(cov_names, collapse = " + ")
  covariate_interaction_effects = paste(
    sapply(cov_names, function(cov) paste0("s(t, by = ", cov, ")")),
    collapse = " + "
  )
  
  model_formula_str = paste(
    "y ~ s(t, k = 15) +",
    treatment_effects_formula, "+",
    covariate_main_effects, "+",
    covariate_interaction_effects
  )
  
  cat("Fórmula do modelo de desfecho:\n", model_formula_str, "\n")
  
  outcome_model = bam(
    as.formula(model_formula_str), data = prepared_data,
    discrete = TRUE, select = TRUE, method = "fREML"
  )
  
  cat("Ajuste do modelo de desfecho concluído.\n")
  
  if (return_prepared_data) {
    return(list(model = outcome_model, data = prepared_data))
  } else {
    return(outcome_model)
  }
}

#' @title Prevê Desfechos Contrafactuais (VERSÃO CORRIGIDA E SIMPLIFICADA)
predict_potential_outcomes = function(outcome_model, simulated_covariates, treatment_times, regime) {
  pred_data = copy(simulated_covariates)
  
  # Encontra todas as colunas de covariáveis simuladas (ex: 'a_X1_pred', 'another_a_X2_pred')
  sim_cov_names = grep("^(a_|another_a_)X[0-9]+_pred$", names(pred_data), value = TRUE)
  
  # Se encontrou alguma, renomeia para o nome base que o modelo de desfecho espera (ex: 'X1', 'X2')
  if (length(sim_cov_names) > 0) {
    # Remove os prefixos ('a_', 'another_a_') e o sufixo ('_pred') de uma só vez
    new_cov_names = gsub("^(a_|another_a_)|_pred$", "", sim_cov_names)
    
    data.table::setnames(pred_data, old = sim_cov_names, new = new_cov_names)
  }
  
  # Adiciona as colunas do regime de tratamento
  for (i in seq_along(treatment_times)) {
    s = treatment_times[i]
    regime_val = regime[i]
    A_var_name = paste0("A_s", s)
    pred_data[, (A_var_name) := regime_val]
  }
  
  # Prevê o desfecho
  predictions = predict(outcome_model, newdata = pred_data)
  
  return(data.table(id = pred_data$id, t = pred_data$t, y_pred = predictions))
}

#' @title Executa a Suíte Completa de Validação da G-Formula ("Super Validador")
run_full_gformula_suite = function(scenario_data, scenario_config, scenario_name) {
  
  # --- PASSO 1: G-FORMULA CORE (MODELAGEM E PREDIÇÃO) ---
  cat(paste("\n--- Iniciando G-Formula para o cenário:", scenario_name, "---\n"))
  
  cov_results = run_gformula_validation_suite(scenario_data, scenario_config)
  
  outcome_model_fit = fit_outcome_model(
    scenario_data$outcome$long, 
    scenario_config$treatment_times, 
    return_prepared_data = TRUE
  )
  outcome_model = outcome_model_fit$model
  
  y_pred_a = predict_potential_outcomes(
    outcome_model, cov_results$simulations$a_seq,
    scenario_config$treatment_times, scenario_config$a_seq
  )
  y_pred_another_a = predict_potential_outcomes(
    outcome_model, cov_results$simulations$another_a_seq,
    scenario_config$treatment_times, scenario_config$another_a_seq
  )
  
  # --- PASSO 2: CONSOLIDAR TODOS OS DADOS (VERDADEIROS E SIMULADOS) ---
  cat("Consolidando todos os dados para validação...\n")
  
  # Começa com os dados observados e ajustados de Y
  validation_master = outcome_model_fit$data[, .(id, t, y_observed = y)]
  validation_master[, y_fitted := predict(outcome_model, newdata = outcome_model_fit$data)]
  validation_master[, id := as.numeric(as.character(id))]
  
  # Adicionar predições de Y
  setnames(y_pred_a, "y_pred", "y_pred_a")
  setnames(y_pred_another_a, "y_pred", "y_pred_another_a")
  validation_master = merge(validation_master, y_pred_a, by = c("id", "t"), all = TRUE)
  validation_master = merge(validation_master, y_pred_another_a, by = c("id", "t"), all = TRUE)
  
  # Adicionar verdades contrafactuais de Y (DIRETO DO DATA.TABLE LONGO)
  true_y_long = scenario_data$outcome$long[, .(id, t, y_true_a = a_y, y_true_another_a = another_a_y)]
  validation_master = merge(validation_master, true_y_long, by = c("id", "t"), all = TRUE)
  
  # Adicionar covariáveis (preditas e verdadeiras)
  cov_names = grep("^X[0-9]+$", names(scenario_data$confounders), value = TRUE)
  
  sim_a_renamed = copy(cov_results$simulations$a_seq)
  setnames(sim_a_renamed, paste0("a_", cov_names, "_pred"), paste0(cov_names, "_pred_a"))
  
  sim_another_a_renamed = copy(cov_results$simulations$another_a_seq)
  setnames(sim_another_a_renamed, paste0("another_a_", cov_names, "_pred"), paste0(cov_names, "_pred_another_a"))
  
  validation_master = merge(validation_master, sim_a_renamed, by = c("id", "t"), all = TRUE)
  validation_master = merge(validation_master, sim_another_a_renamed, by = c("id", "t"), all = TRUE)
  
  # Adicionar verdades contrafactuais de X (DIRETO DO DATA.TABLE LONGO)
  true_covs_long = scenario_data$confounders[, .SD, .SDcols = c("id", "t", paste0("a_", cov_names), paste0("another_a_", cov_names))]
  setnames(true_covs_long, paste0("a_", cov_names), paste0(cov_names, "_true_a"))
  setnames(true_covs_long, paste0("another_a_", cov_names), paste0(cov_names, "_true_another_a"))
  validation_master = merge(validation_master, true_covs_long, by = c("id", "t"), all = TRUE)
  
  # --- PASSO 3: GERAR MÉTRICAS E GRÁFICOS DE VALIDAÇÃO ---
  cat("Gerando suíte completa de validação...\n")
  
  # Helper para criar gráficos e métricas
  create_validation_report = function(data, true_col, pred_col, group_col, title, subtitle) {
    if (!all(c(true_col, pred_col) %in% names(data))) return(list(metrics=data.table(), plot=ggplot()))
    
    metrics = data[, .(
      RMSE = sqrt(mean((get(true_col) - get(pred_col))^2, na.rm = TRUE)),
      Bias = mean(get(pred_col) - get(true_col), na.rm = TRUE)
    )]
    
    plot_data = data[, .(True = mean(get(true_col), na.rm = TRUE), Pred = mean(get(pred_col), na.rm = TRUE)), by = get(group_col)]
    setnames(plot_data, "get", group_col)
    
    plot = ggplot(plot_data, aes_string(x = group_col)) +
      geom_line(aes(y = True, color = "Verdadeiro"), linetype = "dashed", size = 1) +
      geom_line(aes(y = Pred, color = "Simulado"), size = 1) +
      scale_color_manual(values = c("Verdadeiro" = "black", "Simulado" = "orange")) +
      labs(title = title, subtitle = paste(subtitle, "| RMSE:", round(metrics$RMSE, 3), "| Bias:", round(metrics$Bias, 3)), x=group_col, y="Valor Médio") +
      theme_minimal()
    return(list(metrics = metrics, plot = plot))
  }
  
  # 3.1 Validação das Covariáveis
  covariate_reports = list()
  for (cov in cov_names) {
    covariate_reports[[cov]] = list(
      a_seq = create_validation_report(validation_master, paste0(cov, "_true_a"), paste0(cov, "_pred_a"), "t", paste("Validação da Covariável:", cov), "Regime: Sempre Tratar"),
      another_a_seq = create_validation_report(validation_master, paste0(cov, "_true_another_a"), paste0(cov, "_pred_another_a"), "t", paste("Validação da Covariável:", cov), "Regime: Nunca Tratar")
    )
  }
  
  # 3.2 Validação do Desfecho
  outcome_reports = list(
    observed_fit = create_validation_report(validation_master, "y_observed", "y_fitted", "t", "Validação do Desfecho: Ajuste vs. Observado", scenario_name),
    a_seq = create_validation_report(validation_master, "y_true_a", "y_pred_a", "t", "Validação do Desfecho: Simulado vs. Verdadeiro", "Regime: Sempre Tratar"),
    another_a_seq = create_validation_report(validation_master, "y_true_another_a", "y_pred_another_a", "t", "Validação do Desfecho: Simulado vs. Verdadeiro", "Regime: Nunca Tratar")
  )
  
  # 3.3 Validação Causal
  causal_data = validation_master[, .(
    true_effect = mean(y_true_a - y_true_another_a, na.rm = TRUE),
    pred_effect = mean(y_pred_a - y_pred_another_a, na.rm = TRUE)
  ), by = t]
  causal_report = create_validation_report(causal_data, "true_effect", "pred_effect", "t", "Validação Final: Efeito Causal", scenario_name)
  
  cat("--- Validação da G-Formula Concluída ---\n")
  
  return(list(
    covariate_validation = covariate_reports,
    outcome_validation = outcome_reports,
    causal_validation = causal_report
  ))
}