library(mgcv)
library(ggplot2)
library(data.table)

# Carrega as funções de modelagem das covariáveis.
# Certifique-se de que este script contém a versão final de fit_covariate_models
# (com bs="cr" e bs="tp") e as outras funções auxiliares.
source("scripts_functional/models/g_formula_covariates_model.R")

# ===================================================================
# FUNÇÕES AUXILIARES (VERSÕES FINAIS E CORRIGIDAS)
# ===================================================================

#' @title Ajusta o Modelo GAM para o Desfecho Y
fit_outcome_model = function(dt_long, treatment_times, return_prepared_data = FALSE) {
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
  
  outcome_model = bam(
    as.formula(model_formula_str), data = prepared_data,
    discrete = TRUE, select = TRUE, method = "fREML"
  )
  
  if (return_prepared_data) {
    return(list(model = outcome_model, data = prepared_data))
  } else {
    return(outcome_model)
  }
}

#' @title Prevê Desfechos Contrafactuais (VERSÃO FINAL E SIMPLIFICADA)
predict_potential_outcomes = function(outcome_model, simulated_covariates, treatment_times, regime) {
  pred_data = copy(simulated_covariates)
  
  # A entrada é 'simulated_covariates' com colunas como 'X1_sim', 'X2_sim'.
  # O modelo de desfecho espera 'X1', 'X2'. Esta é a única renomeação necessária.
  sim_names = names(pred_data)[grep("_sim$", names(pred_data))]
  
  if (length(sim_names) > 0) {
    new_names = gsub("_sim$", "", sim_names)
    data.table::setnames(pred_data, old = sim_names, new = new_names)
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


# ===================================================================
# FUNÇÃO PRINCIPAL COM CROSS-FITTING (FLUXO LIMPO)
# ===================================================================

#' @title Executa a Suíte Completa da G-Formula com Cross-Fitting ("Super Validador v5.0")
run_full_gformula_suite = function(
    scenario_data, 
    scenario_config, 
    scenario_name,
    k_folds = 5
) {
  
  cat(paste("\n--- Iniciando G-Formula com Cross-Fitting (K=", k_folds, ") para o cenário:", scenario_name, "---\n"))
  
  all_ids = unique(scenario_data$outcome$long$id)
  if (length(all_ids) < k_folds) stop("Número de indivíduos menor que k_folds.")
  set.seed(123) 
  id_folds = data.table(id = all_ids, fold = sample(1:k_folds, length(all_ids), replace = TRUE))
  
  all_fold_predictions = list()
  
  # --- INÍCIO DO LOOP DE CROSS-FITTING ---
  for (k in 1:k_folds) {
    cat(paste("\n--- Processando Fold", k, "/", k_folds, "---\n"))
    
    train_ids = id_folds[fold != k, id]
    predict_ids = id_folds[fold == k, id]
    
    train_data_conf = scenario_data$confounders[id %in% train_ids]
    train_data_outcome = scenario_data$outcome$long[id %in% train_ids]
    predict_data_conf = scenario_data$confounders[id %in% predict_ids]
    
    # --- PASSO 1: AJUSTAR MODELOS (SOMENTE NOS DADOS DE TREINO) ---
    
    # 1.1 Ajustar modelos de covariáveis
    prepared_train_data = prepare_gformula_data(train_data_conf, scenario_config$treatment_times)
    cov_names_to_model = grep("^X[0-9]+$", names(train_data_conf), value = TRUE)
    
    cov_models_train = fit_covariate_models(
      prepared_train_data,
      scenario_config$treatment_times,
      cov_names_to_model
    )
    
    # 1.2 Ajustar Modelo de Desfecho
    outcome_model_train = fit_outcome_model(
      train_data_outcome, 
      scenario_config$treatment_times
    )
    
    # --- PASSO 2: SIMULAR E PREVER (PARA DADOS DE PREDIÇÃO) ---
    
    # 2.1 Simular covariáveis
    sim_a_pred = simulate_counterfactual_trajectories(
      cov_models_train, predict_data_conf,
      scenario_config$treatment_times, scenario_config$a_seq
    )
    sim_another_a_pred = simulate_counterfactual_trajectories(
      cov_models_train, predict_data_conf,
      scenario_config$treatment_times, scenario_config$another_a_seq
    )
    
    # 2.2 Prever desfechos
    y_pred_a = predict_potential_outcomes(
      outcome_model_train, sim_a_pred,
      scenario_config$treatment_times, scenario_config$a_seq
    )
    y_pred_another_a = predict_potential_outcomes(
      outcome_model_train, sim_another_a_pred,
      scenario_config$treatment_times, scenario_config$another_a_seq
    )
    
    # --- PASSO 2.5: COLETAR PREDIÇÕES DO FOLD ---
    predictions_k = merge(y_pred_a, y_pred_another_a, by = c("id", "t"))
    setnames(predictions_k, c("y_pred.x", "y_pred.y"), c("y_pred_a", "y_pred_another_a"))
    
    cov_names = gsub("_sim", "", names(sim_a_pred)[-c(1,2)])
    setnames(sim_a_pred, paste0(cov_names, "_sim"), paste0(cov_names, "_pred_a"))
    setnames(sim_another_a_pred, paste0(cov_names, "_sim"), paste0(cov_names, "_pred_another_a"))
    
    predictions_k = merge(predictions_k, sim_a_pred, by = c("id", "t"))
    predictions_k = merge(predictions_k, sim_another_a_pred, by = c("id", "t"))
    
    all_fold_predictions[[k]] = predictions_k
  }
  # --- FIM DO LOOP ---
  
  # --- PASSO 3: CONSOLIDAR PREDIÇÕES E VALIDAR ---
  # Esta seção agora deve funcionar, pois lida com a saída consolidada do loop
  
  final_predictions = rbindlist(all_fold_predictions)
  
  # A fonte de verdade (observada e contrafactual) vem de scenario_data$outcome$long
  validation_master = merge(final_predictions, scenario_data$outcome$long, by = c("id", "t"))
  
  cov_names = grep("^X[0-9]+$", names(scenario_data$confounders), value = TRUE)
  setnames(validation_master, paste0("a_", cov_names), paste0(cov_names, "_true_a"), skip_absent = TRUE)
  setnames(validation_master, paste0("another_a_", cov_names), paste0(cov_names, "_true_another_a"), skip_absent = TRUE)
  setnames(validation_master, "y", "y_observed")
  setnames(validation_master, "a_y", "y_true_a")
  setnames(validation_master, "another_a_y", "y_true_another_a")
  
  # Ajusta um modelo final no dataset completo apenas para o diagnóstico 'y_fitted'
  final_outcome_model_fit = fit_outcome_model(scenario_data$outcome$long, scenario_config$treatment_times, return_prepared_data = TRUE)
  final_outcome_model_fit$data[, id := as.numeric(as.character(id))]
  validation_master = merge(validation_master, final_outcome_model_fit$data[, .(id, t, y_fitted=predict(final_outcome_model_fit$model, newdata=.SD))], by=c("id", "t"))
  
  # --- PASSO 4: GERAR MÉTRICAS E GRÁFICOS DE VALIDAÇÃO ---
  cat("Gerando suíte completa de validação...\n")
  
  create_validation_report = function(data, true_col, pred_col, group_col, title, subtitle) {
    if (!all(c(true_col, pred_col) %in% names(data))) return(list(metrics=data.table(RMSE=NA, Bias=NA), plot=ggplot() + labs(title="Dados de verdade fundamental ausentes")))
    
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
  
  covariate_reports = list()
  for (cov in cov_names) {
    renamed_pred_a = paste0(cov, "_pred_a")
    renamed_pred_another_a = paste0(cov, "_pred_another_a")
    if(!renamed_pred_a %in% names(validation_master)){
      setnames(validation_master, paste0("a_", cov, "_pred"), renamed_pred_a, skip_absent = TRUE)
      setnames(validation_master, paste0("another_a_", cov, "_pred"), renamed_pred_another_a, skip_absent = TRUE)
    }
    
    covariate_reports[[cov]] = list(
      a_seq = create_validation_report(validation_master, paste0(cov, "_true_a"), paste0(cov, "_pred_a"), "t", paste("Validação da Covariável:", cov), "Regime: Sempre Tratar"),
      another_a_seq = create_validation_report(validation_master, paste0(cov, "_true_another_a"), paste0(cov, "_pred_another_a"), "t", paste("Validação da Covariável:", cov), "Regime: Nunca Tratar")
    )
  }
  
  outcome_reports = list(
    observed_fit = create_validation_report(validation_master, "y_observed", "y_fitted", "t", "Validação do Desfecho: Ajuste vs. Observado", scenario_name),
    a_seq = create_validation_report(validation_master, "y_true_a", "y_pred_a", "t", "Validação do Desfecho: Simulado vs. Verdadeiro", "Regime: Sempre Tratar"),
    another_a_seq = create_validation_report(validation_master, "y_true_another_a", "y_pred_another_a", "t", "Validação do Desfecho: Simulado vs. Verdadeiro", "Regime: Nunca Tratar")
  )
  
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