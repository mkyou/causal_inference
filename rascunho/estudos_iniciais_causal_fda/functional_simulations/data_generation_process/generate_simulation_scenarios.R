generate_simulation_scenarios = function(
    scenario_configs,
    global_params = list(),
    process_and_smooth = TRUE
) {
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
    
    params = utils::modifyList(global_params, scenario_configs[[i]])
    
    cat("  - Etapa 1: Gerando covariáveis (generate_sequential_dgp)\n")
    dgp_data = generate_sequential_dgp(
      n_subjects = params$n_subjects,
      time_grid = params$time_grid,
      treatment_times = params$treatment_times,
      baseline_func_templates = params$baseline_func_templates,
      baseline_params_list = params$baseline_params_list,
      beta_funcs_matrix = params$beta_funcs_matrix,
      gamma_matrix = params$gamma_matrix,
      theta_matrix = params$theta_matrix,
      zeta_funcs_matrix = params$zeta_funcs_matrix,
      noise_sd_vec = params$noise_sd_vec,
      regime_a = params$regime_a,
      regime_b = params$regime_b,
      static_covariate_indices = params$static_covariate_indices,
      propensity_max = params$propensity_max %||% 0.999,
      propensity_min = params$propensity_min %||% 0.001
    )
    
    cat("  - Etapa 2: Gerando desfecho (generate_functional_outcome)\n")
    final_dt = generate_functional_outcome(
      dgp_data = dgp_data,
      baseline_y_func = params$baseline_y_func,
      direct_effect_funcs = params$direct_effect_funcs,
      concurrent_effect_funcs = params$concurrent_effect_funcs,
      treatment_times = params$treatment_times,
      noise_sd = params$outcome_noise_sd
    )
    
    final_data_object = if (process_and_smooth) {
      cat("  - Etapa 3: Suavizando dados para formato funcional\n")
      nbasis = params$smoothing_nbasis %||% 50
      process_and_smooth_data(
        dgp_data = final_dt,
        treatment_times = params$treatment_times,
        nbasis = nbasis
      )
    } else {
      final_dt
    }
    
    all_scenarios_data[[i]] = list(
      data = final_data_object,
      config = params
    )
  }
  
  names(all_scenarios_data) = labels
  cat("\n--- Geração de todos os cenários concluída. ---\n")
  return(all_scenarios_data)
}

`%||%` = function(a, b) {
  if (is.null(a)) b else a
}