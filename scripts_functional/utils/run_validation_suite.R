#' @title Executa a Suíte Completa de Validação
#' @description Itera sobre todos os cenários gerados, produzindo gráficos e
#' uma tabela de resumo das validações.
#' @param all_scenarios_data A lista completa gerada por `generate_simulation_scenarios`.
#' @return Uma lista contendo todos os gráficos e a tabela de validação.

run_validation_suite = function(all_scenarios_data) {
  
  validation_results = list()
  
  # Tabela de resumo
  summary_table = rbindlist(
    lapply(names(all_scenarios_data), function(name) {
      res = validate_scenario_assumptions(all_scenarios_data[[name]])
      res[, Scenario := name]
      return(res)
    }),
    fill = TRUE
  )
  
  print("--- Tabela de Resumo da Validação ---")
  print(summary_table)
  
  validation_results$summary_table = summary_table
  
  # Gráficos
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