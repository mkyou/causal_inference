# R/scenario_runner.R
source("R/utils_dgp_functions.R")
source("R/generate_sequential_dgp.R")
source("R/generate_functional_outcomes.R")
source("R/config_ecker.R")

#' Gera um Dataset Completo para um Cenário
#' @param func_type Nome do tipo de função ("Sigmoid", "Doppler", etc)
#' @param magnitude Magnitude do efeito
#' @param N Numero de individuos
#' @param T_pts Pontos no tempo
#' @return Objeto dgp_out pronto
generate_single_scenario_dataset = function(
    func_type, magnitude, N = 2000, T_pts = 512
    ) {
  
  # Grid fixo Ecker
  t_start = 20.0
  t_end = 63.0
  time_grid = seq(t_start, t_end, length.out = T_pts)
  t_treat = 20.0 
  
  # Seleciona a função baseada no nome
  effect_func = switch(
    func_type,
    "Sigmoid" = effect_sigmoid,
    "Doppler" = effect_doppler,
    "Bumps" = effect_bumps,
    "Heavisine" = effect_heavisine
  )
  
  # Função Alvo com parâmetros fixados
  target_f = function(t, s) {
    effect_func(t, s, t_min = t_start, t_max = t_end, magnitude = magnitude)
  }
  
  # Config Ecker
  config = get_ecker_configuration(N, T_pts)
  
  # 1. História (X, A)
  dgp_covs = generate_sequential_dgp(
    n_subjects = N, time_grid = time_grid, treatment_times = c(t_treat),
    baseline_func_templates = config$baseline_func_templates,
    baseline_params_list = config$baseline_params_list,
    beta_funcs_matrix = config$beta_funcs_matrix,
    gamma_matrix = config$gamma_matrix,
    theta_matrix = config$theta_matrix,
    noise_sd_vec = rep(0, 9), 
    static_covariate_indices = config$static_indices
  )
  
  # 2. Desfecho (Y)
  dgp_out = generate_functional_outcome(
    dgp_data = dgp_covs,
    baseline_y_func = config$baseline_y_func,
    direct_effect_funcs = list(target_f),
    concurrent_effect_funcs = config$concurrent_effect_funcs,
    treatment_times = c(t_treat),
    noise_sd = 0.5
  )
  
  # Anexa a função alvo verdadeira para cálculo de erro depois
  attr(dgp_out, "true_effect_curve") = target_f(time_grid, t_treat)
  attr(dgp_out, "time_grid") = time_grid
  attr(dgp_out, "X_base") = as.matrix(
    dgp_covs[t == time_grid[1], .(X1)]
    ) # Simplificação: Pegar todas X depois se precisar
  X_cols = grep("^X[0-9]+", names(dgp_covs), value=TRUE)
  attr(dgp_out, "X_mat") = as.matrix(unique(dgp_covs[, ..X_cols]))
  attr(dgp_out, "A_vec") = unique(dgp_covs[, .(id, Trat_1)])[, Trat_1]
  
  return(dgp_out)
}