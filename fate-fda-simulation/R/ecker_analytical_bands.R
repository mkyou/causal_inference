# =============================================================================
# R/ecker_analytical_bands.R
# Bandas simultâneas analíticas de Ecker et al. (2024)
# Eq: theta_hat(t) +/- u_alpha(t) * sqrt(C_theta(t,t))
# =============================================================================

library(fda)

#' Estima FATE via splines (abordagem Ecker) com bandas analíticas
#'
#' @param Y_wide Matriz N x T de desfechos ruidosos
#' @param X_reg Matriz N x (K+1) com intercepto
#' @param A_vec Vetor binário de tratamento
#' @param time_grid Vetor T de tempos
#' @param nbasis Número de bases B-spline
#' @param alpha Nível de significância
#' @param n_sim Número de simulações gaussianas para calibrar u_alpha
#' @return Lista: fate_est, lower, upper, C_theta_diag
estimate_fate_ecker = function(
    Y_wide, X_reg, A_vec, time_grid,
    nbasis = 50, alpha = 0.05, n_sim = 5000,
    ridge_lambda = 1e-7
) {
  T_pts = length(time_grid)
  p = ncol(X_reg)
  x_bar = colMeans(X_reg)

  idx_1 = which(A_vec == 1)
  idx_0 = which(A_vec == 0)
  N1 = length(idx_1)
  N0 = length(idx_0)

  # --- 1. Suaviza Y com B-splines ---
  basis = create.bspline.basis(range(time_grid), nbasis = min(nbasis, T_pts / 4))
  fd_par = fdPar(basis, Lfdobj = 2, lambda = 1)  # lambda será otimizado

  # GCV para encontrar lambda ótimo (subamostra)
  sample_idx = if (nrow(Y_wide) > 50) sample(nrow(Y_wide), 50) else seq_len(nrow(Y_wide))
  log_lambdas = seq(-10, 1, length.out = 20)
  gcvs = sapply(log_lambdas, function(ll) {
    fp = fdPar(basis, 2, 10^ll)
    sum(smooth.basis(time_grid, t(Y_wide[sample_idx, ]), fp)$gcv)
  })
  best_lambda = 10^log_lambdas[which.min(gcvs)]

  # Suavização final
  final_par = fdPar(basis, 2, best_lambda)
  Y_smooth_t = smooth.basis(time_grid, t(Y_wide), final_par)$fd
  Y_smooth = t(eval.fd(time_grid, Y_smooth_t))  # N x T

  # --- 2. OLS ponto a ponto nos dados suavizados ---
  X1 = X_reg[idx_1, , drop = FALSE]
  X0 = X_reg[idx_0, , drop = FALSE]

  XtX1 = crossprod(X1); diag(XtX1) = diag(XtX1) + ridge_lambda
  XtX0 = crossprod(X0); diag(XtX0) = diag(XtX0) + ridge_lambda
  XtX1_inv = solve(XtX1)
  XtX0_inv = solve(XtX0)

  H1 = XtX1_inv %*% t(X1)  # p x N1
  H0 = XtX0_inv %*% t(X0)  # p x N0

  Y1_smooth = Y_smooth[idx_1, , drop = FALSE]  # N1 x T
  Y0_smooth = Y_smooth[idx_0, , drop = FALSE]  # N0 x T

  beta1 = H1 %*% Y1_smooth  # p x T
  beta0 = H0 %*% Y0_smooth  # p x T

  fate_est = as.numeric(x_bar %*% (beta1 - beta0))  # vetor T

  res1 = Y1_smooth - X1 %*% beta1  # N1 x T (X1: N1 x p, beta1: p x T)
  res0 = Y0_smooth - X0 %*% beta0  # N0 x T

  # sigma_a(t,t) = var dos resíduos em cada ponto t
  # Estimativa: (1/(N_a - p)) * sum(res_i(t)^2)
  sigma1_tt = colSums(res1^2) / (N1 - p)  # vetor T
  sigma0_tt = colSums(res0^2) / (N0 - p)  # vetor T

  # --- 4. C_theta(t,t) ---
  # C_theta(t,t) = sigma_1(t,t) * x_bar' (X1'X1)^{-1} x_bar
  #              + sigma_0(t,t) * x_bar' (X0'X0)^{-1} x_bar
  vf1 = as.numeric(t(x_bar) %*% XtX1_inv %*% x_bar)
  vf0 = as.numeric(t(x_bar) %*% XtX0_inv %*% x_bar)

  C_theta_diag = sigma1_tt * vf1 + sigma0_tt * vf0  # vetor T

  # --- 5. Calibrar u_alpha via simulação gaussiana ---
  # Simular Z(t) ~ N(0, C_theta_normalized) e obter quantil de max|Z(t)|
  # C_theta_normalized(s,t) = C_theta(s,t) / sqrt(C_theta(s,s) * C_theta(t,t))
  # Para bandas simultâneas, precisamos da covariância completa.
  # Aproximação: usamos a covariância dos resíduos empíricos.

  # Estima a covariância temporal completa dos resíduos
  # Sigma_res^(a)(s,t) = (1/(N_a-p)) * sum(res_i(s) * res_i(t))
  Sigma_res1 = crossprod(res1) / (N1 - p)  # T x T
  Sigma_res0 = crossprod(res0) / (N0 - p)  # T x T

  # C_theta(s,t) = vf1 * Sigma_res1(s,t) + vf0 * Sigma_res0(s,t)
  C_theta = vf1 * Sigma_res1 + vf0 * Sigma_res0  # T x T

  # Normalizada para calibração de u_alpha
  C_theta_sd = sqrt(pmax(diag(C_theta), 1e-12))
  C_theta_cor = C_theta / outer(C_theta_sd, C_theta_sd)

  # Simulação: gerar Z ~ N(0, C_theta_cor), calcular max|Z/1| = max|Z|
  L_cor = tryCatch(chol(C_theta_cor), error = function(e) {
    # Fallback: adiciona jitter para positividade
    C_jit = C_theta_cor + diag(1e-6, T_pts)
    chol(C_jit)
  })

  max_abs = numeric(n_sim)
  for (s in seq_len(n_sim)) {
    z = rnorm(T_pts)
    Z_t = as.numeric(L_cor %*% z)  # vetor T ~ N(0, C_cor)
    # Excurso normalizado: Z(t) / sqrt(C_theta(t,t)) já é 1 por construção
    max_abs[s] = max(abs(Z_t))
  }

  u_alpha = quantile(max_abs, 1 - alpha)

  # --- 6. Bandas finais ---
  se = sqrt(pmax(C_theta_diag, 0))
  lower = fate_est - u_alpha * se
  upper = fate_est + u_alpha * se

  list(
    fate_est = fate_est,
    lower = lower,
    upper = upper,
    C_theta_diag = C_theta_diag,
    u_alpha = u_alpha,
    Y_smooth = Y_smooth
  )
}
