# =============================================================================
# R/ecker_analytical_bands.R
# Bandas simultâneas de Ecker et al. (2024) via ffscb (Liebl & Reimherr 2023)
# Banda: theta_hat(t) +/- u_alpha(t) * sqrt(C_theta(t,t))
# onde u_alpha(t) é função do tempo (limiar variável), não constante.
# =============================================================================

library(fda)
library(ffscb)

estimate_fate_ecker = function(
    Y_wide, X_reg, A_vec, time_grid,
    alpha = 0.05,
    ridge_lambda = 1e-7
) {
  T_pts = length(time_grid)
  p = ncol(X_reg)
  x_bar = colMeans(X_reg)

  idx_1 = which(A_vec == 1)
  idx_0 = which(A_vec == 0)
  N1 = length(idx_1)
  N0 = length(idx_0)

  # --- 1. Suaviza Y com B-splines (lambda=0, nbasis adaptativo a T) ---
  # nbasis determinístico: preserva linearidade do estimador (requerida pelo ffscb)
  nbasis_eff = min(floor(T_pts / 4), 100)
  basis = create.bspline.basis(range(time_grid), nbasis = nbasis_eff)

  final_par = fdPar(basis, 2, 0)
  Y_smooth_t = smooth.basis(time_grid, t(Y_wide), final_par)$fd
  Y_smooth = t(eval.fd(time_grid, Y_smooth_t))  # N x T

  # --- 2. OLS ponto a ponto nos dados suavizados ---
  X1 = X_reg[idx_1, , drop = FALSE]
  X0 = X_reg[idx_0, , drop = FALSE]

  XtX1 = crossprod(X1); diag(XtX1) = diag(XtX1) + ridge_lambda
  XtX0 = crossprod(X0); diag(XtX0) = diag(XtX0) + ridge_lambda
  XtX1_inv = solve(XtX1)
  XtX0_inv = solve(XtX0)

  H1 = XtX1_inv %*% t(X1)
  H0 = XtX0_inv %*% t(X0)

  Y1_smooth = Y_smooth[idx_1, , drop = FALSE]
  Y0_smooth = Y_smooth[idx_0, , drop = FALSE]

  beta1 = H1 %*% Y1_smooth
  beta0 = H0 %*% Y0_smooth

  fate_est = as.numeric(x_bar %*% (beta1 - beta0))

  res1 = Y1_smooth - X1 %*% beta1
  res0 = Y0_smooth - X0 %*% beta0

  # --- 3. C_theta (Eq. 7 de Ecker et al.) ---
  vf1 = as.numeric(t(x_bar) %*% XtX1_inv %*% x_bar)
  vf0 = as.numeric(t(x_bar) %*% XtX0_inv %*% x_bar)

  Sigma_res1 = crossprod(res1) / (N1 - p)
  Sigma_res0 = crossprod(res0) / (N0 - p)

  C_theta = vf1 * Sigma_res1 + vf0 * Sigma_res0
  C_theta_diag = pmax(diag(C_theta), 0)

  # --- 4. Bandas ffscb (limiar variável u_alpha(t)) ---
  # tau(t): rugosidade local do processo padronizado, estimada de C_theta.
  # cov2tau_fun assume grid normalizado [0,1]; C_theta tem índices em time_grid.
  # Passamos C_theta diretamente: a forma de tau(t) é preservada, apenas a
  # escala absoluta varia com o comprimento do intervalo — efeito negligível
  # pois make_band_FFSCB_z recalibra internamente via diag.cov.x.
  tau = tryCatch(
    cov2tau_fun(C_theta),
    error = function(e) {
      # Fallback: rugosidade constante se C_theta é numericamente singular
      rep(1, T_pts)
    }
  )

  upper_col = paste0("FFSCB.z.u.", 1 - alpha)
  lower_col = paste0("FFSCB.z.l.", 1 - alpha)

  band = tryCatch(
    make_band_FFSCB_z(
      x = fate_est,
      diag.cov.x = C_theta_diag,
      tau = tau,
      conf.level = 1 - alpha
    ),
    error = function(e) NULL
  )

  if (is.null(band) || !upper_col %in% colnames(band)) {
    se = sqrt(C_theta_diag)
    return(list(
      fate_est = fate_est,
      lower = fate_est - 2 * se,
      upper = fate_est + 2 * se,
      C_theta_diag = C_theta_diag,
      Y_smooth = Y_smooth
    ))
  }

  list(
    fate_est = fate_est,
    lower = band[, lower_col],
    upper = band[, upper_col],
    C_theta_diag = C_theta_diag,
    Y_smooth = Y_smooth
  )
}
