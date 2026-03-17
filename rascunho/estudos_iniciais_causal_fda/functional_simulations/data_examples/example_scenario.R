#' Gera Cenários de Teste
#'
#' Versão final que implementa corretamente o efeito X2 -> X1 para gerar
#' heterogeneidade temporal em X1, mantendo as correções de positividade
#' e de sinal-ruído da versão anterior.
#'
#' @param n_subjects Número de indivíduos a simular.
#' @param time_grid Grade de tempo para a simulação.
#' @param smoothing_nbasis Número de bases para a etapa de suavização.
#'
#' @return Um objeto de cenários pronto para análise.

generate_ecker_luna_scenarios = function(
    n_subjects = 2000,
    time_grid = seq(0, 25, by = 0.25),
    smoothing_nbasis = 50
) {
  
  K = 2
  S = 2
  
  GLOBAL_PARAMS = list(
    n_subjects = n_subjects,
    time_grid = time_grid,
    noise_sd_vec = c(1.5, 1.5),
    outcome_noise_sd = 2.0,
    regime_a = c(1, 1),
    regime_b = c(0, 0),
    smoothing_nbasis = smoothing_nbasis
  )
  
  TREATMENT_TIMES = c(5, 15)
  
  BASE_TEMPLATES = list(
    function(t, p) p[1] + p[2]*t,
    function(t, p) p[1] + p[2]*sin(t*pi/p[3])
  )
  PARAMS_LIST = list(
    matrix(c(2, 0.15), nrow = n_subjects, ncol = 2, byrow = TRUE),
    matrix(c(0, 3, 10), nrow = n_subjects, ncol = 3, byrow = TRUE)
  )
  
  GAMMA_MATRIX = matrix(list(function(t,s) 0), K, S)
  GAMMA_MATRIX[[1, 1]] = function(t, s) {
    ifelse(t > s, 6 * exp(-0.2 * (t - s)), 0)
  }
  
  THETA_CONF = matrix(c(
    0.0, 0.0, 0.0,
    -2.0, 1.5, 0.0
  ), nrow = S, ncol = K + 1, byrow = TRUE)
  
  THETA_NO_CONF = THETA_CONF
  THETA_NO_CONF[2, 2] = 0.0
  
  # MUDANÇA CRÍTICA: Definindo o efeito X2 -> X1.
  # beta_funcs_matrix[[k, j]] é o efeito de Xj em Xk.
  BETA_MATRIX = matrix(list(NULL), K, K)
  BETA_MATRIX[[1, 2]] = function(t) 0.3 # Efeito constante de X2 em X1
  
  ZETA_MATRIX = NULL
  
  BASELINE_Y = function(t) 20 + 0.5*t
  DIRECT_EFFECT_A_Y = list(function(t,s) 0, function(t,s) 0)
  CONCURRENT_EFFECT_X_Y = list(function(t) 1.5 * log(t + 2), function(t) 0)
  
  scenario_configs = list(
    `TimeVaryingConfounding` = list(
      treatment_times = TREATMENT_TIMES,
      baseline_func_templates = BASE_TEMPLATES,
      baseline_params_list = PARAMS_LIST,
      beta_funcs_matrix = BETA_MATRIX, # Matriz com o efeito X2->X1
      gamma_matrix = GAMMA_MATRIX,
      theta_matrix = THETA_CONF,
      zeta_funcs_matrix = ZETA_MATRIX,
      direct_effect_funcs = DIRECT_EFFECT_A_Y,
      concurrent_effect_funcs = CONCURRENT_EFFECT_X_Y,
      baseline_y_func = BASELINE_Y
    ),
    `NoConfounding` = list(
      treatment_times = TREATMENT_TIMES,
      baseline_func_templates = BASE_TEMPLATES,
      baseline_params_list = PARAMS_LIST,
      beta_funcs_matrix = BETA_MATRIX, # Matriz com o efeito X2->X1
      gamma_matrix = GAMMA_MATRIX,
      theta_matrix = THETA_NO_CONF,
      zeta_funcs_matrix = ZETA_MATRIX,
      direct_effect_funcs = DIRECT_EFFECT_A_Y,
      concurrent_effect_funcs = CONCURRENT_EFFECT_X_Y,
      baseline_y_func = BASELINE_Y
    )
  )
  
  all_scenarios = generate_simulation_scenarios(
    scenario_configs = scenario_configs,
    global_params = GLOBAL_PARAMS,
    process_and_smooth = TRUE
  )
  
  return(all_scenarios)
}