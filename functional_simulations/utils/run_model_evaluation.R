#' @title Executa a Avaliação de um Modelo
#' @description Versão generalizada que aceita uma função de avaliação como
#'   parâmetro, permitindo avaliar diferentes aspectos (ex: covariáveis
#'   ou desfecho) com o mesmo fluxo de trabalho.
#'
#' @param all_scenarios_data Lista nomeada com os data.tables gerados.
#' @param scenario_configs A lista de configurações.
#' @param model_function A função de ajuste do modelo a ser testada.
#' @param assessment_function A função de avaliação a ser usada.
#'   Ex: `assess_covariate_simulation` ou `assess_effect_estimation`.
#' @param model_name Uma string para nomear o modelo nos resultados.
#'
#' @return Uma lista nomeada com os resultados da avaliação.

run_model_evaluation = function(
    all_scenarios_data,
    scenario_configs,
    model_function,
    assessment_function, # <-- NOVO PARÂMETRO
    model_name
) {
  
  cat(paste("\n--- Iniciando avaliação para o modelo:", model_name, "---\n"))
  
  all_results = purrr::imap(
    all_scenarios_data,
    function(scenario_data, scenario_name) {
      
      cat(paste("\n>>> Processando cenário:", scenario_name, "<<<\n"))
      
      config = scenario_configs[[scenario_name]]
      S = length(config$treatment_times)
      regime_a_dynamic = rep(1, S)
      regime_b_dynamic = rep(0, S)
      
      # 1. Executa a função do modelo para obter as estimativas
      # A saída desta função (`model_output`) será passada para a avaliação
      model_output = model_function(
        scenario_data = scenario_data,
        scenario_config = config,
        regime_a = regime_a_dynamic,
        regime_b = regime_b_dynamic
      )
      
      # 2. Executa a função de avaliação fornecida
      assessment = assessment_function(
        model_output = model_output,
        true_data = scenario_data
      )
      
      return(assessment)
    }
  )
  
  return(all_results)
}