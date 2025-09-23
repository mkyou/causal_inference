library(data.table)
library(ggplot2)
library(mgcv)

#' @title Avalia a Estimação de um Efeito Causal (Corrigido)
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



#' @title Avalia a Simulação das Covariáveis (Alinhado)
assess_covariate_simulation = function(results_data_with_truth) {
  cov_names_base = grep("^X[0-9]+$", names(results_data_with_truth), value = TRUE)
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


#' @title Avalia o Ajuste Factual do Modelo de Desfecho com CV (Corrigido)
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



#' @title Avalia o Processo Completo da G-Formula (Corrigido)
assess_full_gformula_process = function(model_results, true_data_refund) {
  cat("--- Iniciando Avaliação Completa do Processo ---\n")
  assessment_outcome_model_fit = assess_factual_fit(
    fitted_outcome_model = model_results$fitted_models$outcome_model,
    factual_data = model_results$results_data,
    outcome_model_formula = model_results$formulas$outcome_model_formula
  )
  results_data = model_results$results_data
  
  cat("  - Preparando dados para avaliação...\n")
  true_data_long = revert_to_long_format(true_data_refund)
  
  cols_to_rename = names(true_data_long)[
    !names(true_data_long) %in% c("id", "t")
  ]
  new_names = paste0(cols_to_rename, "_true")
  setnames(true_data_long, cols_to_rename, new_names)
  
  results_data_with_truth = merge(
    results_data, 
    true_data_long, 
    by = c("id", "t")
  )
  
  # ETAPA 1: AVALIAR SIMULAÇÃO DAS COVARIÁVEIS
  cat("  - Etapa 1: Avaliando simulação de covariáveis...\n")
  assessment_covariates = assess_covariate_simulation(results_data_with_truth)
  
  # ETAPA 2: AVALIAR SIMULAÇÃO DO DESFECHO
  cat("  - Avaliando simulação de desfechos...\n")
  outcome_metrics = results_data_with_truth[, .(
    RMSE_A = sqrt(mean((y_a_sim - a_y_true)^2, na.rm = TRUE)),
    Bias_A = mean(y_a_sim - a_y_true, na.rm = TRUE),
    RMSE_B = sqrt(mean((y_b_sim - b_y_true)^2, na.rm = TRUE)),
    Bias_B = mean(y_b_sim - b_y_true, na.rm = TRUE)
  )]
  
  # ETAPA 3: AVALIAR EFEITO CAUSAL
  cat("  - Etapa 3: Avaliando estimação do efeito causal...\n")
  assessment_causal = assess_effect_estimation(results_data_with_truth)
  
  return(
    list(
      outcome_model_factual_fit = assessment_outcome_model_fit,
      covariate_assessment = assessment_covariates,
      outcome_assessment_metrics = outcome_metrics,
      causal_effect_assessment = assessment_causal
    )
  )
}

#' @title Avalia os Modelos de Propensão do MSM
#' @description Fornece diagnósticos básicos para os modelos de propensão ajustados.
assess_propensity_models = function(propensity_models, prepared_data) {
  
  cat("    - Avaliando modelos de propensão...\n")
  all_summaries = list()
  
  for (s_time_char in names(propensity_models)) {
    s_time = as.numeric(s_time_char)
    model = propensity_models[[s_time_char]]
    data_at_s = prepared_data[t == s_time]
    
    # Prever as propensões
    predicted_probs = predict(model, newdata = data_at_s, type = "response")
    
    # Resumo da distribuição das propensões (diagnóstico de positividade)
    prob_summary = summary(predicted_probs)
    
    # Gráfico de calibração simples
    calib_data = data.table(predicted_probs, observed = data_at_s$A)
    calib_plot = ggplot(calib_data, aes(x = predicted_probs, y = observed)) +
      geom_smooth(method = "loess", color = "blue") +
      geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "red") +
      labs(title = paste("Calibração do Modelo de Propensão em t =", s_time),
           x = "Probabilidade Prevista", y = "Probabilidade Observada") +
      coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
      theme_bw()
    
    all_summaries[[s_time_char]] = list(
      summary = prob_summary,
      plot = calib_plot
    )
  }
  return(all_summaries)
}


#' @title Avalia o Processo Completo do MSM
assess_full_msm_process = function(model_results, true_data_refund) {
  cat("--- Iniciando Avaliação Completa do Processo MSM ---\n")
  
  # ETAPA 1: AVALIAR OS MODELOS DE PROPENSÃO
  assessment_propensity = assess_propensity_models(
    propensity_models = model_results$fitted_models$propensity_models,
    prepared_data = model_results$results_data # O MSM results_data tem os dados originais
  )
  
  # ETAPA 2: AVALIAR A DISTRIBUIÇÃO DOS PESOS
  cat("    - Avaliando a distribuição dos pesos (SW)...\n")
  weight_summary = summary(model_results$results_data$sw)
  weight_plot = ggplot(model_results$results_data, aes(x = sw)) +
    geom_histogram(bins = 50, fill = "skyblue", color = "black") +
    labs(title = "Distribuição dos Pesos Estabilizados (SW)", x = "Peso (SW)") +
    theme_bw()
  
  # ETAPA 3: AVALIAR O AJUSTE FACTUAL DO MODELO DE DESFECHO
  assessment_outcome_model_fit = assess_factual_fit(
    fitted_outcome_model = model_results$fitted_models$outcome_model,
    # Atenção: o modelo foi ajustado com pesos, a avaliação factual pode ser otimista
    factual_data = model_results$results_data,
    outcome_model_formula = model_results$formulas$outcome_formula
  )
  
  # --- Preparação dos dados para avaliação causal ---
  cat("  - Preparando dados para avaliação causal...\n")
  true_data_long = revert_to_long_format(true_data_refund)
  
  # Renomeia as colunas da "verdade" para evitar conflitos
  cols_to_rename = names(true_data_long)[!names(true_data_long) %in% c("id", "t")]
  new_names = paste0(cols_to_rename, "_true")
  setnames(true_data_long, cols_to_rename, new_names)
  
  # Junta os resultados do MSM com a verdade
  results_data_with_truth = merge(
    model_results$results_data, 
    true_data_long, 
    by = c("id", "t")
  )
  
  # Renomeia a coluna de efeito para ser compatível com a função de avaliação
  setnames(results_data_with_truth, "effect_estimate_msm", "effect_estimate")
  
  # ETAPA 4: AVALIAR EFEITO CAUSAL
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