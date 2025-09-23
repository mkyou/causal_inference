library(purrr)
source("functional_simulations/utils/revert_to_long_format.R")
source("functional_simulations/utils/prepare_analysis_data.R")

run_model_evaluation = function(
    all_scenarios_data,
    scenario_configs,
    model_function,
    assessment_function,
    model_name
) {
  
  cat(paste("\n--- Iniciando avaliação para o modelo:", model_name, "---\n"))
  
  all_results = purrr::imap(
    all_scenarios_data,
    function(scenario_data_refund, scenario_name) { # Renomeado para clareza
      
      cat(paste("\n>>> Processando cenário:", scenario_name, "<<<\n"))
      
      config = scenario_configs[[scenario_name]]
      S = length(config$treatment_times)
      regime_a_dynamic = rep(1, S)
      regime_b_dynamic = rep(0, S)
      
      # --- CORREÇÃO CRÍTICA ---
      # A chamada a 'revert_to_long_format' é a primeira etapa!
      scenario_data_long = revert_to_long_format(scenario_data_refund)
      
      # 1. Prepara os dados uma única vez
      prepared_data = prepare_analysis_data(
        scenario_data_long, 
        config$treatment_times
      )
      
      # 2. Executa o modelo, que agora recebe os dados preparados
      model_output_object = model_function(
        prepared_data = prepared_data,
        scenario_config = config,
        regime_a = regime_a_dynamic,
        regime_b = regime_b_dynamic
      )
      
      # 3. Avalia o resultado
      assessment = assessment_function(
        model_results = model_output_object,
        true_data_refund = scenario_data_refund
      )
      
      return(assessment)
    }
  )
  
  return(all_results)
}