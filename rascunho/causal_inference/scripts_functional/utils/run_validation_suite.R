#' @title Executa a Suíte Completa de Validação
#' @description Itera sobre todos os cenários, produzindo gráficos e
#' uma tabela de resumo das validações estatisticamente robustas.
#' @param all_scenarios_data A lista completa gerada por `generate_simulation_scenarios`.
#' @param scenario_configs A lista de configurações usada para gerar os dados.

run_validation_suite = function(all_scenarios_data, scenario_configs) {
  
  validation_results = list()
  
  # Tabela de resumo
  summary_table = rbindlist(
    lapply(names(all_scenarios_data), function(name) {
      # Passa tanto os dados gerados quanto a configuração original para a 
      # validação
      res = check_causal_structure(
        scenario_data = all_scenarios_data[[name]],
        scenario_config = scenario_configs[[name]]
      )
      res[, Scenario := name]
      return(res)
    }),
    fill = TRUE
  )
  
  print("--- Tabela de Resumo da Validação Causal ---")
  print(summary_table)
  
  validation_results$summary_table = summary_table
  
  # Geração de Gráficos (sem alteração)
  for (name in names(all_scenarios_data)) {
    cat(paste("\n--- Gerando gráficos para:", name, "---\n"))
    p_traj = plot_scenario_trajectories(all_scenarios_data[[name]])
    p_effects = plot_scenario_effects(all_scenarios_data[[name]])
    validation_results$plots[[name]] = list(
      trajectories = p_traj, 
      effects = p_effects
    )
  }
  
  return(validation_results)
}