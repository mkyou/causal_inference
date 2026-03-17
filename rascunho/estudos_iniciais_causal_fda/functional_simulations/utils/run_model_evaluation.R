library(purrr)

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
      
      config = scenario_object$config
      scenario_data_processed = scenario_object$data
      
      if ("functional_vars" %in% names(scenario_data_processed)) {
        data_for_analysis = revert_to_long_format(scenario_data_processed)
      } else {
        data_for_analysis = copy(scenario_data_processed)
      }
      
      prepared_data = prepare_analysis_data(
        data_long = data_for_analysis,
        treatment_times = config$treatment_times,
        ...
      )
      
      S = length(config$treatment_times)
      regime_a_to_test = config$regime_a %||% rep(1, S)
      regime_b_to_test = config$regime_b %||% rep(0, S)
      
      model_output_object = model_function(
        unprepared_data = data_for_analysis,
        prepared_data = prepared_data,
        scenario_config = config,
        regime_a = regime_a_to_test,
        regime_b = regime_b_to_test
      )
      
      assessment = assessment_function(
        model_results = model_output_object,
        true_data_refund = data_for_analysis
      )
      
      return(assessment)
    }
  )
  
  return(all_results)
}