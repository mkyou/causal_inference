generate_test_scenarios = function(
    n_subjects = 2000, smoothing_nbasis = 50
) {
  
  scenario_1 = list(
    # Parâmetros "Globais"
    n_subjects = n_subjects,
    time_grid = seq(0, 25, by = 0.1),
    noise_sd_vec = c(0.1, 0.2, 0.1, 0.3),
    outcome_noise_sd = 2.0,
    smoothing_nbasis = smoothing_nbasis,
    
    # Parâmetros do DGP
    K = 4, S = 3,
    treatment_times = c(2, 8, 14),
    static_covariate_indices = 1,
    baseline_func_templates = list(
      function(t, p) p[1],
      function(t, p) p[1] - exp(-p[2] * t),
      function(t, p) p[1] + p[2] * t,
      function(t, p) p[1] + p[2] * t^0.5
    ),
    baseline_params_list = list(
      matrix(rnorm(n_subjects, 100, 10), ncol = 1),
      matrix(c(1, 0.4), n_subjects, 2, byrow = TRUE),
      matrix(c(18, 0.2), n_subjects, 2, byrow = TRUE),
      matrix(c(2, 0.5), n_subjects, 2, byrow = TRUE)
    ),
    gamma_matrix = {
      gamma_g = matrix(list(function(t,s) 0), 4, 3)
      gamma_g[[2, 1]] = function(t, s) ifelse(t > s, 0.5*exp(-0.5*(t-s)), 0)
      gamma_g[[4, 2]] = function(t, s) ifelse(t > s, 1.5*exp(-0.3*(t-s)), 0)
      gamma_g[[3, 3]] = function(t, s) ifelse(t > s, 0.8*exp(-0.4*(t-s)), 0)
      gamma_g
    },
    beta_funcs_matrix = {
      beta_g = matrix(list(NULL), 4, 4)
      beta_g[[3, 2]] = function(t) 0.1
      beta_g[[3, 4]] = function(t) -0.05
      beta_g
    },
    
    theta_matrix = rbind(
      # A1 (t=2): Aleatório, P=0.5
      c(0, 0.0, 0.0, 0.0, 0.0),
      # A2 (t=8): Coeficientes mais fortes
      c(1.5, -0.0, 1.0, -0.2, 0.1),
      # A3 (t=14): Coeficientes mais fortes
      c(-1.0, -0.0, 1.5, 0.1, -0.1)
    ),
    
    zeta_funcs_matrix = NULL,
    baseline_y_func = function(t) 180 / (1 + exp(-0.4 * (t - 13))),
    direct_effect_funcs = rep(list(function(t, s) 0), 3),
    concurrent_effect_funcs = list(
      function(t) 0.08, function(t) 0.015 * exp(-0.01 * (t-10)^2),
      function(t) -0.05, function(t) 0
    ),
    propensity_min = 0,
    propensity_max = 1
  )
  
  scenario_2 = list(
    # Parâmetros "Globais"
    n_subjects = n_subjects,
    time_grid = seq(0, 15, by = 0.25),
    noise_sd_vec = c(1, 0, 0.1, 0.2),
    outcome_noise_sd = 50.0,
    smoothing_nbasis = smoothing_nbasis,
    
    # Parâmetros do DGP
    K = 4, S = 2,
    treatment_times = c(3, 8),
    static_covariate_indices = 2,
    baseline_func_templates = list(
      function(t, p) p[1] + p[2]*t, function(t, p) p[1],
      function(t, p) p[1] * (1-exp(-p[2]*t)), function(t, p) p[1]*log(t + 1)
    ),
    baseline_params_list = list(
      matrix(c(5, 0.2), n_subjects, 2, byrow = TRUE),
      matrix(rnorm(n_subjects, 10, 2), ncol = 1),
      matrix(c(1, 0.3), n_subjects, 2, byrow = TRUE),
      matrix(c(1, 0), n_subjects, 2, byrow = TRUE)
    ),
    gamma_matrix = {
      gamma_s = matrix(list(function(t,s) 0), 4, 2)
      gamma_s[[1, 1]] = function(t, s) 
        ifelse(t > s, 3 * exp(-0.2 * (t-s)), 0)
      gamma_s[[3, 1]] = function(t, s) 
        ifelse(t > s, 0.5 * exp(-0.3 * (t-s)), 0)
      gamma_s[[1, 2]] = function(t, s) 
        ifelse(t > s, 5*(t-s)*exp(-0.8*(t-s)), 0)
      gamma_s
    },
    beta_funcs_matrix = {
      beta_s = matrix(list(NULL), 4, 4)
      beta_s[[1, 4]] = function(t) 0.1
      beta_s
    },
    theta_matrix = matrix(c(
      -3, 0.5, 0, 0, 0,
      -4, 0.8, 0, 0, 0
    ), nrow = 2, ncol = 5, byrow = TRUE),
    zeta_funcs_matrix = NULL,
    baseline_y_func = function(t) 1200 + 80*log(t+1),
    direct_effect_funcs = list(
      function(t,s) ifelse(t>s, 50*exp(-0.1*(t-s)), 0), function(t,s) 0
    ),
    concurrent_effect_funcs = list(
      function(t) 15, function(t) 5, function(t) 40, function(t) 20*exp(-0.1*t)
    ),
    propensity_min = 0,
    propensity_max = 1
  )
  
  # --- ORQUESTRAÇÃO FINAL ---
  scenario_configs = list(
    `CrescimentoHumano` = scenario_1,
    `AquisicaoHabilidade` = scenario_2
  )
  
  all_scenarios = generate_simulation_scenarios(
    scenario_configs = scenario_configs,
    process_and_smooth = FALSE
  )
  
  return(all_scenarios)
}