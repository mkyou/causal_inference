# =============================================================================
# R/generate_ecker_dgp.R
# DGP Heterogêneo de Ecker + Erro MGP
# =============================================================================

library(data.table)

# --- Cholesky pré-computado do MGP ---

#' Gera a decomposição de Cholesky da covariância MGP (exponencial)
#' @param time_grid Vetor de tempos
#' @param sigma Desvio padrão do processo
#' @param rho Parâmetro de decaimento (no grid normalizado [0,1])
#' @return Matriz triangular superior L tal que L'L = Sigma
make_mgp_chol = function(time_grid, sigma = 0.5, rho = 0.8) {
  tg_norm = (time_grid - min(time_grid)) / diff(range(time_grid))
  D = as.matrix(dist(tg_norm))
  Sigma = sigma^2 * exp(-D / rho)
  chol(Sigma)
}

#' Gera N realizações do MGP via Cholesky pré-computado
#' @param N Número de indivíduos
#' @param T_pts Número de pontos no tempo
#' @param L_chol Decomposição de Cholesky (output de make_mgp_chol)
#' @return Matriz N x T_pts de erros
generate_mgp_errors = function(
  N, T_pts, L_chol,
  error_dist = "gaussian", df_t = 5
) {
  Z = matrix(rnorm(N * T_pts), nrow = N, ncol = T_pts)
  eps = Z %*% L_chol
  if (error_dist == "t") {
    scale = sqrt(rchisq(N, df_t) / df_t)
    eps = sweep(eps, 1, scale, "/")
  }
  eps
}

# --- DGP principal ---

#' Gera dataset completo com DGP heterogêneo de Ecker
#'
#' Y_i(t) = mu_Y(t) + X_i^T beta^(A_i)(t) + eps_i(t)
#' onde beta^(1)(t) - beta^(0)(t) = delta(t) distribui o efeito,
#' e eps_i(t) ~ MGP(0, sigma, I_N).
#'
#' @param scen Lista com: name, func, mag, family
#' @param N Número de indivíduos
#' @param T_pts Número de pontos no tempo
#' @param time_grid Vetor de tempos
#' @param L_chol Cholesky pré-computado (ou NULL para ruído branco)
#' @param config Configuração Ecker (output de get_ecker_configuration)
#' @param t_start,t_end,t_treat Limites temporais
#' @return Lista com Y_wide, X_mat, A_vec, time_grid, true_curve
generate_ecker_dataset = function(
  scen, N, T_pts, time_grid, L_chol,
  config = NULL,
  t_start = 20.0, t_end = 80.0, t_treat = 50.0,
  error_dist = "gaussian", df_t = 5
) {
  # Config Ecker (gera covariáveis e propensity)
  if (is.null(config)) config = get_ecker_configuration(N, T_pts)

  # --- Gera X e A via generate_sequential_dgp ---
  dgp_covs = generate_sequential_dgp(
    n_subjects = N,
    time_grid = time_grid,
    treatment_times = c(t_treat),
    baseline_func_templates = config$baseline_func_templates,
    baseline_params_list = config$baseline_params_list,
    beta_funcs_matrix = config$beta_funcs_matrix,
    gamma_matrix = config$gamma_matrix,
    theta_matrix = config$theta_matrix,
    noise_sd_vec = rep(0, 9),
    static_covariate_indices = config$static_indices
  )

  X_cols = grep("^X[0-9]+", names(dgp_covs), value = TRUE)
  X_mat = as.matrix(unique(dgp_covs[, ..X_cols]))  # N x 9
  A_vec = unique(dgp_covs[, .(id, Trat_1)])[, Trat_1]
  rm(dgp_covs); gc(verbose = FALSE)

  # X com intercepto para regressão: N x 10
  X_reg = cbind(1, X_mat)
  K1 = ncol(X_reg)  # 10

  # --- Função de efeito do cenário ---
  effect_func = function(t) {
    scen$func(t, t_treat, t_min = t_start, t_max = t_end, magnitude = scen$mag)
  }

  # --- Efeitos concorrentes base: beta^(0)_k(t) ---
  concurrent_fns = config$concurrent_effect_funcs

  # Efeito heterogêneo: beta^(1) - beta^(0) = delta(t) / (K+1) p/ componente
  # Isso torna o FATE = X_bar^T * delta(t)
  beta_diff_func = function(t) rep(effect_func(t) / K1, K1)

  # --- Monta matrizes de contribuição ---
  # Como o modelo OLS assume:
  # E[Y(0) | X] = mu_Y(t) + X_mat * beta^(0)(t)
  # E[Y(1) | X] = mu_Y(t) + X_mat * beta^(1)(t) + intercepto_extra

  # 1. Efeitos concorrentes base: beta^(0)_k(t)
  beta0_mat = sapply(time_grid, function(t) {
    sapply(concurrent_fns, function(f) f(t))
  }) # 9 x T

  # 2. Delta intercepto e slope
  # beta_diff_func já tem K1 = 10 elementos.
  # delta_mat[1, ] é o intercepto extra para tratados.
  # delta_mat[2:10, ] são as diferenças de slope.
  delta_mat = sapply(time_grid, function(t) beta_diff_func(t))  # 10 x T

  # 3. Construindo as médias condicionais
  # Controle
  mu_0_mat = X_mat %*% beta0_mat  # N x T

  # Tratados (Intercepto extra + Slope extra)
  # delta_mat[1,] é reciclado pra virar N x T
  intercept_effect_1 = sweep(matrix(0, N, T_pts), 2, delta_mat[1, ], "+")
  slope_effect_1 = X_mat %*% delta_mat[2:10, ]

  mu_1_mat = mu_0_mat + intercept_effect_1 + slope_effect_1 # N x T

  # FATE exato é E[Y(1) - Y(0)] sobre toda a amostra N
  # E[Y(1) - Y(0)] = mean(mu_1_mat - mu_0_mat) ponto a ponto
  true_curve = colMeans(mu_1_mat - mu_0_mat)  # T

  # Agora monta Y_wide de acordo com o tratamento:
  # Se A_i = 1 usa mu_1_mat, se 0 usa mu_0_mat
  cov_eff = mu_0_mat * (1 - A_vec) + mu_1_mat * A_vec

  # 4. mu_Y(t) universal (cancela na diferença)
  mu_vec = 2 * log(time_grid + 0.5)

  # 5. Erros MGP
  if (!is.null(L_chol)) {
    eps_mat = generate_mgp_errors(N, T_pts, L_chol, error_dist, df_t)
  } else {
    eps_mat = matrix(rnorm(N * T_pts, 0, 0.5), nrow = N, ncol = T_pts)
  }

  # Y_wide final
  Y_wide = sweep(cov_eff + eps_mat, 2, mu_vec, "+")

  # Limpeza
  rm(
    beta0_mat, delta_mat, mu_0_mat, 
    mu_1_mat, intercept_effect_1, slope_effect_1, 
    cov_eff, eps_mat
  )
  gc(verbose = FALSE)

  list(
    Y_wide = Y_wide, # N x T
    X_mat = X_mat, # N x 9 (sem intercepto, para compatibilidade)
    X_reg = X_reg, # N x 10 (com intercepto)
    A_vec = A_vec, # N
    time_grid = time_grid, # T
    true_curve = true_curve # T
  )
}
