# ===================================================================
# SUÍTE DE AVALIAÇÃO DE ESTIMAÇÃO DE EFEITO
# ===================================================================
library(data.table)
library(ggplot2)

#' @title Avalia a Qualidade da Estimação de um Efeito Causal
#' @description Compara um efeito estimado com o efeito verdadeiro,
#'   calculando métricas de erro funcionais e sumarizadas.
#'
#' @param estimation_dt data.table com `id`, `t`, e `effect_estimate`.
#' @param true_data data.table do PGD, que contém `id`, `t`, e o
#'   verdadeiro efeito em `effect_a_vs_b`.
#'
#' @return Uma lista contendo um data.table de métricas sumarizadas,
#'   um data.table de métricas funcionais (ao longo do tempo), e um
#'   gráfico comparativo.

assess_effect_estimation = function(model_output, true_data_refund) {
  
  true_effect_mat = true_data_refund$functional_vars$a_y -
    true_data_refund$functional_vars$b_y
  
  true_data_long = as.data.table(true_effect_mat)
  setnames(
    true_data_long, 
    names(true_data_long),
    as.character(true_data_refund$time_grid)
  )
  
  true_data_long[, id := true_data_refund$id]
  true_data_long = melt(
    true_data_long, id.vars = "id", variable.name = "t",
    value.name = "true_effect"
  )
  true_data_long[, t := as.numeric(as.character(t))]
  
  estimation_dt = model_output
  estimation_dt[, id := as.factor(id)]
  true_data_long[, id := as.factor(id)]
  
  eval_data = merge(estimation_dt, true_data_long, by = c("id", "t"))
  
  eval_data[, error := effect_estimate - true_effect]
  
  functional_metrics = eval_data[
    , 
    .(
      Bias = mean(error),
      RMSE = sqrt(mean(error^2)),
      MAE = mean(abs(error))
    ), 
    by = t
  ]
  
  summary_metrics = eval_data[
    , 
    .(
      Mean_Bias = mean(error),
      Median_Bias = median(error),
      Mean_RMSE = mean(functional_metrics$RMSE),
      Overall_RMSE = sqrt(mean(error^2))
    )
  ]
  
  mean_effects = eval_data[
    , 
    .(
      est = mean(effect_estimate),
      true = mean(true_effect)
    ), 
    by = t
  ]
  
  summary_metrics[, Correlation := cor(mean_effects$est,
                                       mean_effects$true)]
  
  plot_comparison = ggplot(mean_effects, aes(x = t)) +
    geom_line(
      aes(y = true, color = "Verdadeiro"),
      linetype = "dashed", linewidth = 1.1
    ) +
    geom_line(
      aes(y = est, color = "Estimado"), linewidth = 1.1, alpha = 0.8
    ) +
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


#' @title Avalia a Qualidade da Simulação das Covariáveis
#' @description Compara as trajetórias de covariáveis simuladas com as
#'   trajetórias contrafactuais verdadeiras do PGD.
#' @param results_dt O data.table de saída da `run_gformula_covariates`.
#' @return Uma lista com métricas de erro e gráficos comparativos.

assess_covariate_simulation = function(model_output, true_data_refund) {
  results_dt = model_output 
  cov_names = grep("^X[0-9]+$", names(results_dt), value = TRUE)
  all_metrics = list()
  all_plots = list()
  
  for (k in cov_names) {
    # Nomes das colunas para esta covariável
    true_a_col = paste0("a_", k)
    true_b_col = paste0("b_", k)
    sim_a_col = paste0("sim_a_", k)
    sim_b_col = paste0("sim_b_", k)
    
    # Prepara os dados para análise
    eval_data = results_dt[, .(
      t, 
      true_a = get(true_a_col), sim_a = get(sim_a_col),
      true_b = get(true_b_col), sim_b = get(sim_b_col)
    )]
    
    # Métricas sumarizadas
    metrics_a = eval_data[, .(
      RMSE = sqrt(mean((sim_a - true_a)^2)),
      Bias = mean(sim_a - true_a)
    )]
    metrics_b = eval_data[, .(
      RMSE = sqrt(mean((sim_b - true_b)^2)),
      Bias = mean(sim_b - true_b)
    )]
    
    all_metrics[[k]] = list(regime_a = metrics_a, regime_b = metrics_b)
    
    # Gráfico comparativo
    plot_data = eval_data[, .(
      True_A = mean(true_a), Sim_A = mean(sim_a),
      True_B = mean(true_b), Sim_B = mean(sim_b)
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


#' @title Avalia o Ajuste de um Modelo de Desfecho Já Treinado
assess_factual_fit = function(fitted_outcome_model, k_folds = 5) {
  # Extrai métricas in-sample do sumário do modelo
  summary_model = summary(fitted_outcome_model)
  metrics = data.table(
    Metric = c("R_squared_adj", "Deviance_explained"),
    Value = c(summary_model$r.sq, summary_model$dev.expl),
    Type = "In-Sample"
  )
  # (A lógica de Validação Cruzada, se desejada, seria adicionada aqui)
  return(metrics)
}



#' @title Avalia o Processo Completo da G-Formula
#' @description Recebe um objeto rico da `run_gformula_full` e avalia
#'   cada etapa: covariáveis, desfechos e efeito causal.
assess_full_gformula_process = function(model_results, true_data_refund) {
  
  cat("--- Iniciando Avaliação Completa do Processo ---\n")
  
  # ETAPA 0: AVALIAR O AJUSTE DO MODELO DE DESFECHO NOS DADOS FACTUAIS
  cat("  - Etapa 0: Avaliando performance do modelo de desfecho...\n")
  assessment_outcome_model_fit = assess_factual_fit(
    fitted_outcome_model = model_results$fitted_models$outcome_model
  )
  
  # Extrai o data.table principal de resultados para as próximas etapas
  results_data = model_results$results_data
  
  # ETAPA 1: AVALIAR SIMULAÇÃO DAS COVARIÁVEIS
  cat("  - Etapa 1: Avaliando simulação de covariáveis...\n")
  assessment_covariates = assess_covariate_simulation(
    model_output = results_data, # Passa o data.table correto
    true_data_refund = true_data_refund
  )
  
  # 2. Avaliar Simulação do Desfecho
  cat("  - Avaliando simulação de desfechos...\n")
  outcome_eval_data = model_output$outcomes
  outcome_metrics = outcome_eval_data[, .(
    RMSE_A = sqrt(mean((y_a_sim - a_y)^2)),
    Bias_A = mean(y_a_sim - a_y),
    RMSE_B = sqrt(mean((y_b_sim - b_y)^2)),
    Bias_B = mean(y_b_sim - b_y)
  )]
  
  # ETAPA 3: AVALIAR EFEITO CAUSAL
  cat("  - Etapa 3: Avaliando estimação do efeito causal...\n")
  assessment_causal = assess_effect_estimation(
    model_output = results_data, # Passa o data.table correto
    true_data_refund = true_data_refund
  )
  
  return(list(
    outcome_model_factual_fit = assessment_outcome_model_fit,
    covariate_assessment = assessment_covariates,
    outcome_assessment_metrics = outcome_metrics,
    causal_effect_assessment = assessment_causal
  ))
}