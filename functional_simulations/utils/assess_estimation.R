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

assess_effect_estimation = function(model_output, true_data) {
  estimation_dt = model_output
  # Junta a estimativa com a verdade
  eval_data = merge(
    estimation_dt,
    true_data[, .(id, t, true_effect = effect_a_vs_b)],
    by = c("id", "t")
  )
  
  # Calcula o erro ponto a ponto
  eval_data[, error := effect_estimate - true_effect]
  
  # 1. Calcula as métricas funcionais (para cada ponto no tempo)
  functional_metrics = eval_data[, .(
    Bias = mean(error),
    RMSE = sqrt(mean(error^2)),
    MAE = mean(abs(error))
  ), by = t]
  
  # 2. Calcula as métricas sumarizadas (agregadas sobre todo o tempo)
  summary_metrics = eval_data[, .(
    Mean_Bias = mean(error),
    Median_Bias = median(error),
    Mean_RMSE = mean(functional_metrics$RMSE),
    Overall_RMSE = sqrt(mean(error^2))
  )]
  
  # 3. Calcula a correlação e métrica invariante à escala
  mean_effects = eval_data[, .(
    est = mean(effect_estimate),
    true = mean(true_effect)
  ), by = t]
  
  summary_metrics[, Correlation := cor(mean_effects$est, mean_effects$true)]
  
  # 4. Cria o gráfico comparativo
  plot_comparison = ggplot(mean_effects, aes(x = t)) +
    geom_line(
      aes(y = true, color = "Verdadeiro"), linetype = "dashed", linewidth = 1.1
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

assess_covariate_simulation = function(model_output, true_data) {
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