#' @title Gera Múltiplos Cenários de Simulação
#' @description Orquestra a geração de cenários. Versão refatorada para ser
#'   compatível com o PGD com heterogeneidade individual e com a suavização 
#'   otimizada.
#'
#' @param scenario_configs Lista nomeada de configurações de cenário. Espera
#'  que cada config contenha `baseline_func_templates` e
#'  `baseline_params_list`.
#' @param global_params Lista de parâmetros globais.
#' @param process_and_smooth Booleano. Se TRUE, a saída será processada pela
#'   função `process_and_smooth_data`.
#' @param smoothing_nbasis Número de bases para a suavização.
#'
#' @return Uma lista nomeada de cenários de dados.

generate_simulation_scenarios = function(
    scenario_configs, 
    global_params = list(),
    process_and_smooth = TRUE,
    smoothing_nbasis = 25
) {
  # Validação
  if (missing(scenario_configs) || !is.list(scenario_configs) || 
      length(scenario_configs) == 0) {
    stop("scenario_configs deve ser uma lista não vazia de configurações.")
  }
  
  labels = names(scenario_configs) %||% paste0(
    "Cenario_", seq_along(scenario_configs)
  )
  
  all_scenarios_data = vector("list", length(scenario_configs))
  
  for (i in seq_along(scenario_configs)) {
    cat(paste("\n--- Gerando Cenário:", labels[i], "---\n"))
    
    # 1. Combina parâmetros globais e específicos do cenário
    params = utils::modifyList(global_params, scenario_configs[[i]])
    
    # 2. Etapa 1: Gerar as covariáveis e tratamentos
    cat("  - Etapa 1: Gerando covariáveis (generate_sequential_dgp)\n")
    
    dgp_data = generate_sequential_dgp(
      n_subjects = params$n_subjects,
      time_grid = params$time_grid,
      treatment_times = params$treatment_times,
      
      # Argumentos para heterogeneidade individual
      baseline_func_templates = params$baseline_func_templates,
      baseline_params_list = params$baseline_params_list,
      
      # Argumentos restantes
      beta_funcs_matrix = params$beta_funcs_matrix,
      gamma_matrix = params$gamma_matrix,
      theta_matrix = params$theta_matrix,
      zeta_funcs_matrix = params$zeta_funcs_matrix,
      noise_sd_vec = params$noise_sd_vec,
      regime_a = params$regime_a,
      regime_b = params$regime_b
    )
    
    # 3. Etapa 2: Gerar o desfecho
    cat("  - Etapa 2: Gerando desfecho (generate_functional_outcome)\n")
    final_dt = generate_functional_outcome(
      dgp_data = dgp_data,
      baseline_y_func = params$baseline_y_func,
      direct_effect_funcs = params$direct_effect_funcs,
      concurrent_effect_funcs = params$concurrent_effect_funcs,
      treatment_times = params$treatment_times,
      noise_sd = params$outcome_noise_sd
    )
    
    # 4. Etapa 3 (Opcional): Processar e suavizar os dados
    if (process_and_smooth) {
      cat("  - Etapa 3: Suavizando dados para formato funcional\n")
      
      # Usa o nbasis do cenário, se definido, senão usa o padrão da função
      nbasis = params$smoothing_nbasis %||% smoothing_nbasis
      
      processed_data = process_and_smooth_data(
        dgp_data = final_dt,
        treatment_times = params$treatment_times,
        nbasis = nbasis
      )
      all_scenarios_data[[i]] = processed_data
    } else {
      all_scenarios_data[[i]] = final_dt
    }
  }
  
  names(all_scenarios_data) = labels
  cat("\n--- Geração de todos os cenários concluída. ---\n")
  return(all_scenarios_data)
}

# Operador auxiliar para definir valores padrão
`%||%` = function(a, b) {
  if (is.null(a)) b else a
}