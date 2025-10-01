library(purrr)

#' @title Orquestra a Avaliação de um Modelo Causal em Múltiplos Cenários
#' @description Itera sobre uma lista de cenários de dados, prepara os dados,
#'   executa uma função de modelo e uma função de avaliação para cada um.
#'
#' @param all_scenarios_data Lista de objetos de cenário 
#'   (cada um com $data e $config).
#' @param model_function A função que ajusta o modelo causal (ex: run_fmsm).
#' @param assessment_function A função que avalia os resultados do modelo.
#' @param model_name Uma string com o nome do modelo para os logs.
#' @param ... Argumentos adicionais a serem passados para 
#' `prepare_analysis_data`.
#'
#' @return Uma lista nomeada com os resultados da avaliação para cada cenário.

run_model_evaluation = function(
    all_scenarios_data,
    model_function,
    assessment_function,
    model_name,
    ...
) {
  cat(paste("\n--- Iniciando avaliação para o modelo:", model_name, "---\n"))
  
  all_results = purrr::imap(
    all_scenarios_data,
    function(scenario_object, scenario_name) {
      cat(paste("\n>>> Processando cenário:", scenario_name, "<<<\n"))
      
      # Desempacota o objeto para obter dados e configuração
      config = scenario_object$config
      scenario_data_processed = scenario_object$data
      
      scenario_data_long = revert_to_long_format(scenario_data_processed)
      
      prepared_data = prepare_analysis_data(
        data_long = scenario_data_long,
        treatment_times = config$treatment_times,
        ...
      )
      
      S = length(config$treatment_times)
      regime_a_to_test = config$regime_a %||% rep(1, S)
      regime_b_to_test = config$regime_b %||% rep(0, S)
      
      model_output_object = model_function(
        scenario_object = scenario_object,
        prepared_data = prepared_data,
        scenario_config = config,
        regime_a = regime_a_to_test,
        regime_b = regime_b_to_test
      )
      
      assessment = assessment_function(
        model_results = model_output_object,
        true_data_refund = scenario_data_processed
      )
      
      return(assessment)
    }
  )
  
  return(all_results)
}