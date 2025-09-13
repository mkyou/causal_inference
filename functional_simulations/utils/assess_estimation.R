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

assess_effect_estimation = function(estimation_dt, true_data) {
  
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