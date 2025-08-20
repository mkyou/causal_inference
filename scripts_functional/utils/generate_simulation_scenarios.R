generate_simulation_scenarios = function(scenario_configs) {
  # Validação
  if (missing(scenario_configs) || !is.list(scenario_configs) || length(scenario_configs) == 0) {
    stop("scenario_configs deve ser uma lista não vazia de configurações.")
  }
  
  labels = if (!is.null(names(scenario_configs))) names(scenario_configs) else paste0("Cenário ", seq_along(scenario_configs))
  
  scenarios = vector("list", length(scenario_configs))
  
  for (i in seq_along(scenario_configs)) {
    cfg = scenario_configs[[i]]
    
    # --- Parâmetros Padrão (podem ser sobrescritos pelo cfg) ---
    params = list(
      n = 1000,
      time_grid = seq(0, 10, by = 0.2),
      base_y_fun = function(t) 1.5 + 0.5 * sin(0.3 * t),
      treatment_effect = list(lapply(1:10, function(s) function(t, s) 1.2 * exp(-0.05 * (t - s)))), # Nota: treatment_effect agora tem que ter o tamanho dos tempos de tratamento
      mean_baseline_covariates = list( # Renomeado para clareza
        function(t) 0.5 + 0.15 * t^1.2,
        function(t) 0.3 + 0.1 * sqrt(t),
        function(t) 0.6 + 0.1 * log1p(t),
        function(t) 0.4 + 0.05 * t
      ),
      noise_sd = 0.5,
      p_treat = 0.5
    )
    # Atualiza os padrões com o que foi fornecido no cenário
    params = utils::modifyList(params, cfg)
    
    # --- Ponto Crítico: Geração das Curvas Individuais ---
    baseline_list_heterogeneous = create_heterogeneous_baselines_advanced(
      n = params$n,
      mean_funcs = params$mean_baseline_covariates
    )
    
    # 1) Geração de X e tratamentos
    dt_conf = generate_confounder_pointwise(
      n = params$n,
      effects_function_list = params$effects_covariates_from_A,
      noise_sd = params$noise_sd,
      p_treat = params$p_treat,
      time_grid = params$time_grid,
      baseline_function_list = baseline_list_heterogeneous, # Passando a lista correta!
      propensity_function_list = params$prop_func,
      return_format = "long"
    )
    
    # 2) Geração de Y 
    dt_y = generate_y_concurrent(
      dt = dt_conf,
      base_function = params$base_y_fun,
      effects_function_list = params$treatment_effect[[1]], # Acessando o primeiro (e único) elemento
      covariates_concurrent_effect_list = params$concurrent_effects_X,
      noise_sd = params$noise_sd,
      return_format = "both"
    )
    
    scenarios[[i]] = list(
      # smooth = s_dt, # Ignorando por enquanto
      raw_long = dt_y$long,
      raw_wide = dt_y$wide
    )
  }
  
  names(scenarios) = labels
  return(scenarios)
}


