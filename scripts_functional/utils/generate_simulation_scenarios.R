generate_simulation_scenarios = function(scenario_configs, global_params = list()) {
  # Validação
  if (missing(scenario_configs) || !is.list(scenario_configs) || length(scenario_configs) == 0) {
    stop("scenario_configs deve ser uma lista não vazia de configurações.")
  }
  
  labels = if (!is.null(names(scenario_configs))) names(scenario_configs) else paste0("Cenário ", seq_along(scenario_configs))
  
  scenarios = vector("list", length(scenario_configs))
  
  for (i in seq_along(scenario_configs)) {
    
    # --- Lógica de Parâmetros ---
    cfg_specific = scenario_configs[[i]]
    params = utils::modifyList(global_params, cfg_specific)
    
    # --- Geração das Curvas Individuais ---
    baseline_list_heterogeneous = create_heterogeneous_baselines_advanced(
      n = params$n,
      mean_funcs = params$mean_baseline_covariates,
      intercept_sd = params$intercept_sd %||% 0.2,
      scale_sd = params$scale_sd %||% 0.1,
      slope_sd = params$slope_sd %||% 0.05
    )
    
    # 1) Geração de X e tratamentos
    dt_conf = generate_confounder_pointwise(
      n = params$n,
      effects_function_list = params$effects_covariates_from_A,
      p_treat = params$p_treat,
      noise_sd = params$noise_sd,
      time_grid = params$time_grid,
      baseline_function_list = baseline_list_heterogeneous,
      propensity_function_list = params$propensity_function_list,
      treatment_times = params$treatment_times,
      a_seq = params$a_seq,
      another_a_seq = params$another_a_seq,
      return_format = "long"
    )
    
    # 2) Geração de Y 
    y_model_type = params$y_model_type %||% "concurrent"
    
    if (y_model_type == "concurrent") {
      dt_y = generate_y_concurrent(
        dt = dt_conf,
        base_function = params$base_y_fun,
        effects_function_list = params$treatment_effect,
        covariates_concurrent_effect_list = params$concurrent_effects_X,
        noise_sd = params$noise_sd,
        treatment_times = params$treatment_times, 
        return_format = "both"
      )
    } else if (y_model_type == "historical") {
      dt_y = generate_y_historical(
        dt = dt_conf,
        base_function = params$base_y_fun,
        effects_function_list = params$treatment_effect,
        covariates_historical_effect_list = params$historical_effects_X,
        noise_sd = params$noise_sd,
        treatment_times = params$treatment_times, 
        return_format = "both"
      )
    } else {
      stop("y_model_type desconhecido: ", y_model_type)
    }
    
    scenarios[[i]] = list(
      confounders = dt_conf,
      outcome = dt_y
    )
  }
  
  names(scenarios) = labels
  return(scenarios)
}

# Operador auxiliar para definir valores padrão (requer R >= 4.0)
`%||%` = function(a, b) {
  if (is.null(a)) b else a
}