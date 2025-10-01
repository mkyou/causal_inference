library(data.table)
library(ggplot2)
library(mgcv)

#' @title Avalia a Estimação de um Efeito Causal
assess_effect_estimation = function(results_data_with_truth) {
  eval_data = results_data_with_truth
  eval_data[, error := effect_estimate - effect_a_vs_b_true]
  
  functional_metrics = eval_data[, .(
    Bias = mean(error), RMSE = sqrt(mean(error^2))
  ), by = t]
  
  summary_metrics = eval_data[, .(
    Mean_Bias = mean(error),
    Median_Bias = median(error),
    Overall_RMSE = sqrt(mean(error^2))
  )]
  
  mean_effects = eval_data[, .(
    est = mean(effect_estimate, na.rm = TRUE),
    true = mean(effect_a_vs_b_true, na.rm = TRUE)
  ), by = t]
  
  summary_metrics[, Correlation := cor(mean_effects$est, mean_effects$true)]
  
  plot_comparison = ggplot(mean_effects, aes(x = t)) +
    geom_line(aes(y = true, color = "Verdadeiro"), lty = "dashed", lwd = 1.1) +
    geom_line(aes(y = est, color = "Estimado"), lwd = 1.1, alpha = 0.8) +
    labs(
      title = "Comparação: Efeito Causal Verdadeiro vs. Estimado",
      subtitle = paste0(
        "Corr: ", round(summary_metrics$Correlation, 3),
        ", RMSE Total: ", round(summary_metrics$Overall_RMSE, 3),
        ", Viés Médio: ", round(summary_metrics$Mean_Bias, 3)
      ),
      x = "Tempo", y = "Efeito Causal Médio", color = "Série"
    ) +
    theme_bw()
  
  return(list(
    summary = summary_metrics,
    functional = functional_metrics,
    plot = plot_comparison
  ))
}



#' @title Avalia a Simulação das Covariáveis
assess_covariate_simulation = function(results_data_with_truth) {
  cov_names_base = grep(
    "^X[0-9]+$", 
    names(results_data_with_truth), value = TRUE
  )
  all_metrics = list()
  all_plots = list()
  
  for (k in cov_names_base) {
    # MUDANÇA: Usa os nomes de coluna com sufixo '_true'
    true_a_col = paste0("a_", k, "_true")
    true_b_col = paste0("b_", k, "_true")
    sim_a_col = paste0("sim_a_", k)
    sim_b_col = paste0("sim_b_", k)
    
    metrics_a = results_data_with_truth[, .(
      RMSE = sqrt(mean((get(sim_a_col) - get(true_a_col))^2)),
      Bias = mean(get(sim_a_col) - get(true_a_col))
    )]
    metrics_b = results_data_with_truth[, .(
      RMSE = sqrt(mean((get(sim_b_col) - get(true_b_col))^2)),
      Bias = mean(get(sim_b_col) - get(true_b_col))
    )]
    
    all_metrics[[k]] = list(regime_a = metrics_a, regime_b = metrics_b)
    
    plot_data = results_data_with_truth[, .(
      True_A = mean(get(true_a_col)), Sim_A = mean(get(sim_a_col)),
      True_B = mean(get(true_b_col)), Sim_B = mean(get(sim_b_col))
    ), by = t]
    
    plot = ggplot(plot_data, aes(x = t)) +
      geom_line(aes(y = True_A, color = "Verdadeiro A"), lty = "dashed") +
      geom_line(aes(y = Sim_A, color = "Simulado A")) +
      geom_line(aes(y = True_B, color = "Verdadeiro B"), lty = "dashed") +
      geom_line(aes(y = Sim_B, color = "Simulado B")) +
      scale_color_manual(
        name = "Trajetória",
        values = c("Verdadeiro A" = "darkblue", "Simulado A" = "lightblue",
                   "Verdadeiro B" = "darkred", "Simulado B" = "pink")
      ) +
      labs(
        title = paste("Avaliação da Simulação da Covariável:", k),
        subtitle = paste0(
          "Regime A (RMSE/Bias): ", 
          round(metrics_a$RMSE, 2), "/", round(metrics_a$Bias, 2),
          " | Regime B (RMSE/Bias): ", 
          round(metrics_b$RMSE, 2), "/", round(metrics_b$Bias, 2)
        ),
        x = "Tempo", y = "Valor Médio da Covariável"
      ) +
      theme_bw()
    
    all_plots[[k]] = plot
  }
  
  return(list(metrics = all_metrics, plots = all_plots))
}

#' @title Avalia o Ajuste Factual do Modelo de Desfecho com CV
assess_factual_fit = function(
    fitted_outcome_model, 
    factual_data,
    outcome_model_formula,
    k_folds = 5
) {
  summary_model = summary(fitted_outcome_model)
  in_sample_metrics = data.table(
    Metric = c("R_squared_adj", "Deviance_explained"),
    Value = c(summary_model$r.sq, summary_model$dev.expl),
    Type = "In-Sample"
  )
  
  # Lógica de Validação Cruzada
  fit_data = na.omit(factual_data)[sample(.N)]
  folds = cut(seq(1, nrow(fit_data)), breaks = k_folds, labels = FALSE)
  cv_errors = vector("list", length = k_folds)
  
  for (i in 1:k_folds) {
    test_indices = which(folds == i)
    train_data = fit_data[-test_indices, ]
    test_data = fit_data[test_indices, ]
    
    fold_model = bam(outcome_model_formula, data = train_data, discrete = TRUE, select = TRUE)
    predictions_fold = predict(fold_model, newdata = test_data)
    cv_errors[[i]] = test_data$y - predictions_fold
  }
  
  all_cv_errors = unlist(cv_errors)
  cv_metrics = data.table(
    Metric = c("RMSE_CV", "Bias_CV"),
    Value = c(sqrt(mean(all_cv_errors^2)), mean(all_cv_errors)),
    Type = paste0(k_folds, "-Fold CV")
  )
  
  return(rbind(in_sample_metrics, cv_metrics))
}


# --- FUNÇÃO DE AVALIAÇÃO DO MODELO DE PROPENSÃO (RECONSTRUÍDA) ---
#' @title Avalia os Modelos de Propensão com Validação Cruzada
#' @description Usa validação cruzada para fornecer uma avaliação robusta
#'   da performance preditiva do modelo de propensão.
#'
#' @param propensity_models Lista de modelos de propensão ajustados no F-MSM.
#' @param prepared_data O data.table com os dados preparados.
#' @param k_folds Número de folds para a validação cruzada.
#'
#' @return Uma lista de diagnósticos para cada tempo de tratamento.

assess_propensity_models = function(
    propensity_models, 
    prepared_data, 
    k_folds = 5
) {
  cat("    - Avaliando modelos de propensão com validação cruzada...\n")
  all_assessments = list()
  
  for (s_time_char in names(propensity_models)) {
    s_time = as.numeric(s_time_char)
    model = propensity_models[[s_time_char]]
    data_at_s = na.omit(prepared_data[t == s_time])
    
    # Prepara para CV
    data_at_s = data_at_s[sample(.N)] # Embaralha os dados
    folds = cut(seq(1, nrow(data_at_s)), breaks = k_folds, labels = FALSE)
    
    cv_predictions = vector("list", length = k_folds)
    
    for (i in 1:k_folds) {
      test_indices = which(folds == i, arr.ind = TRUE)
      train_data = data_at_s[-test_indices, ]
      test_data = data_at_s[test_indices, ]
      
      # Ajusta o modelo no fold de treino
      fold_model = bam(
        formula(model), 
        data = train_data, 
        family = binomial(),
        select = TRUE,
        discrete = TRUE
        )
      
      # Prevê no fold de teste
      preds = predict(fold_model, newdata = test_data, type = "response")
      
      cv_predictions[[i]] = data.table(
        prob = preds, obs = test_data$A, fold = i
      )
    }
    
    # Consolida predições fora da amostra
    oof_preds = rbindlist(cv_predictions)
    
    # Calcula Métricas
    brier = mean((oof_preds$obs - oof_preds$prob)^2)
    
    # Métricas de classificação (com limiar 0.5)
    preds_class = ifelse(oof_preds$prob > 0.5, 1, 0)
    tp = sum(preds_class == 1 & oof_preds$obs == 1)
    fp = sum(preds_class == 1 & oof_preds$obs == 0)
    fn = sum(preds_class == 0 & oof_preds$obs == 1)
    tn = sum(preds_class == 0 & oof_preds$obs == 0)
    
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    accuracy = (tp + tn) / nrow(oof_preds)
    f1_score = 2 * (precision * recall) / (precision + recall)
    
    # AUC
    roc_dt = oof_preds[order(-prob)]
    roc_dt[, tpr := cumsum(obs) / sum(obs)]
    roc_dt[, fpr := cumsum(1 - obs) / sum(1 - obs)]
    auc = sum(
      (roc_dt$fpr[-1] - roc_dt$fpr[-nrow(roc_dt)]) * 
        roc_dt$tpr[-nrow(roc_dt)]
    )
    
    metrics_summary = data.table(
      Metrica = c("Brier_Score", "AUC", "F1_Score", "Acuracia", 
                  "Precisao", "Recall"),
      Valor = c(brier, auc, f1_score, accuracy, precision, recall)
    )
    
    # Gráfico de Calibração (com predições OOF)
    calib_plot = ggplot(oof_preds, aes(x = prob, y = obs)) +
      geom_smooth(method = "loess", color = "blue") +
      geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "red") +
      labs(
        title = paste(
          "Calibração (CV) do Modelo de Propensão em t =", 
          s_time
        ),
        x = "Probabilidade Prevista (Fora da Amostra)", 
        y = "Probabilidade Observada"
      ) +
      coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
      theme_bw()
    
    all_assessments[[s_time_char]] = list(
      metrics = metrics_summary,
      plot = calib_plot
    )
  }
  return(all_assessments)
}


#' @title Avalia o Processo Completo do MSM
assess_full_msm_process = function(model_results, true_data_refund) {
  cat("--- Iniciando Avaliação Completa do Processo MSM ---\n")
  
  # ETAPA 1: AVALIAR OS MODELOS DE PROPENSÃO
  assessment_propensity = assess_propensity_models(
    propensity_models = model_results$fitted_models$propensity_models,
    prepared_data = model_results$results_data
  )
  
  # ETAPA 2: AVALIAR A DISTRIBUIÇÃO DOS PESOS (com avisos)
  cat("    - Avaliando a distribuição dos pesos (SW)...\n")
  weights = model_results$results_data$sw
  weight_summary = summary(weights)
  if (max(weights, na.rm = TRUE) > 50) {
    warning(
      "Instabilidade de pesos detectada: Peso máximo (", 
      round(max(weights, na.rm = TRUE), 2), ") excede o limiar de 50."
    )
  }
  if (sd(weights, na.rm = TRUE) > 5) {
    warning(
      "Alta variabilidade de pesos detectada: Desvio Padrão (",
      round(sd(weights, na.rm = TRUE), 2), ") excede o limiar de 5."
    )
  }
  
  weight_plot = ggplot(data.table(sw = weights), aes(x = sw)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    labs(
      title = "Distribuição dos Pesos Estabilizados (SW)", 
      x = "Peso (SW)"
    ) +
    theme_bw()
  
  # ETAPA 3: AVALIAR O AJUSTE FACTUAL DO MODELO DE DESFECHO
  assessment_outcome_model_fit = assess_factual_fit(
    fitted_outcome_model = model_results$fitted_models$outcome_model,
    # O modelo foi ajustado com pesos, a avaliação factual pode ser otimista
    factual_data = model_results$results_data,
    outcome_model_formula = model_results$formulas$outcome_formula
  )
  
  cat("  - Preparando dados para avaliação causal...\n")
  true_data_long = revert_to_long_format(true_data_refund)
  
  cols_to_rename = names(true_data_long)[!names(true_data_long) %in% 
                                           c("id", "t")]
  new_names = paste0(cols_to_rename, "_true")
  setnames(true_data_long, cols_to_rename, new_names)
  
  results_data_with_truth = merge(
    model_results$results_data, true_data_long, by = c("id", "t")
  )
  
  setnames(results_data_with_truth, "effect_estimate_msm", "effect_estimate")
  
  cat("  - Avaliando estimação do efeito causal...\n")
  assessment_causal = assess_effect_estimation(results_data_with_truth)
  
  return(
    list(
      propensity_model_assessment = assessment_propensity,
      weight_assessment = list(summary = weight_summary, plot = weight_plot),
      outcome_model_factual_fit = assessment_outcome_model_fit,
      causal_effect_assessment = assessment_causal
    )
  )
}
