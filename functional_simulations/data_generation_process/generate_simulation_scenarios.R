#' @title Gera Múltiplos Cenários de Simulação
#' @description Orquestra a geração de múltiplos cenários de simulação,
#'   chamando `generate_sequential_dgp` e em seguida
#'   `generate_functional_outcome` para cada configuração.
#'
#' @param scenario_configs Uma lista nomeada. Cada elemento da lista é outra
#'   lista contendo todos os parâmetros necessários para gerar um cenário
#'   completo (tanto para as covariáveis quanto para o desfecho).
#' @param global_params Uma lista de parâmetros que serão aplicados a
#'   todos os cenários, a menos que sejam sobrescritos por uma
#'   configuração específica do cenário.
#'
#' @return Uma lista nomeada, onde cada elemento contém o data.table
#'   completo (covariáveis e desfecho) para um cenário.

generate_simulation_scenarios = function(
    scenario_configs, 
    global_params = list()
) {
  # Validação
  if (missing(scenario_configs) || !is.list(scenario_configs) || 
      length(scenario_configs) == 0) {
    stop("scenario_configs deve ser uma lista não vazia de configurações.")
  }
  
  labels = if (!is.null(names(scenario_configs))) {
    names(scenario_configs)
  } else {
    paste0("Cenario_", seq_along(scenario_configs))
  }
  
  all_scenarios_data = vector("list", length(scenario_configs))
  
  for (i in seq_along(scenario_configs)) {
    cat(paste("\n--- Gerando Cenário:", labels[i], "---\n"))
    
    # 1. Combina parâmetros globais e específicos do cenário
    cfg_specific = scenario_configs[[i]]
    params = utils::modifyList(global_params, cfg_specific)
    
    # 2. Primeira Etapa: Gerar as covariáveis e tratamentos
    cat("  - Etapa 1: Gerando covariáveis (generate_sequential_dgp)\n")
    dgp_data = generate_sequential_dgp(
      n_subjects = params$n_subjects,
      time_grid = params$time_grid,
      treatment_times = params$treatment_times,
      baseline_funcs = params$baseline_funcs,
      beta_funcs_matrix = params$beta_funcs_matrix,
      gamma_matrix = params$gamma_matrix,
      theta_matrix = params$theta_matrix,
      theta_integral_matrix = params$theta_integral_matrix,
      intercept_sd_vec = params$intercept_sd_vec,
      noise_sd_vec = params$noise_sd_vec,
      regime_a = params$regime_a,
      regime_b = params$regime_b
    )
    
    # 3. Segunda Etapa: Gerar o desfecho usando os dados da etapa 1
    cat("  - Etapa 2: Gerando desfecho (generate_functional_outcome)\n")
    final_dt = generate_functional_outcome(
      dgp_data = dgp_data,
      baseline_y_func = params$baseline_y_func,
      direct_effect_funcs = params$direct_effect_funcs,
      concurrent_effect_funcs = params$concurrent_effect_funcs,
      treatment_times = params$treatment_times,
      noise_sd = params$outcome_noise_sd
    )
    
    all_scenarios_data[[i]] = final_dt
  }
  
  names(all_scenarios_data) = labels
  cat("\n--- Geração de todos os cenários concluída. ---\n")
  return(all_scenarios_data)
}