analyze_ate_cate = function(results) {
  library(data.table)
  library(tibble)
  
  # Criando uma lista para armazenar métricas por simulação
  metrics_list = lapply(results, function(sim_result) {
    # Resgatar o id_scenario de cada simulação
    id_scenario = unique(sim_result$id_scenario)
    
    # Cálculo do ATE: Média e Viés
    ate_estimate_mean = mean(sim_result$estimate)
    ate_true_mean = mean(sim_result$true_effect)
    
    mse_ate = (ate_estimate_mean - ate_true_mean)^2
    bias_ate = ate_estimate_mean - ate_true_mean
    
    # Cálculo do CATE: MSE e Viés
    mse_cate = mean((sim_result$estimate - sim_result$true_effect)^2)
    bias_cate = mean(sim_result$estimate - sim_result$true_effect)
    
    # Retornando os resultados para a simulação
    list(
      id_scenario = id_scenario,
      mse_ate = mse_ate,
      bias_ate = bias_ate,
      mse_cate = mse_cate,
      bias_cate = bias_cate
    )
  })
  
  # Transformar a lista em data.table para sumarização
  metrics_dt = rbindlist(metrics_list, idcol = "sim_id")
  
  # Consolidar os resultados agrupando por id_scenario
  consolidated_results = metrics_dt[, .(
    avg_mse_ate = mean(mse_ate),
    avg_bias_ate = mean(bias_ate),
    avg_mse_cate = mean(mse_cate),
    avg_bias_cate = mean(bias_cate)
  ), by = id_scenario]
  
  # Convertendo para tibble
  return(as_tibble(consolidated_results))
}