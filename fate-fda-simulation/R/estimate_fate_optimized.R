# =============================================================================
# R/estimate_fate_optimized.R
# Estimação vetorizada do FATE com projeções pré-computadas
# =============================================================================

library(data.table)

#' Pré-computa matrizes de projeção para ambos os grupos
#'
#' H^(a) = (X_a'X_a + lambda I)^{-1} X_a'
#' var_factor^(a) = x_bar' (X_a'X_a + lambda I)^{-1} x_bar
#'
#' @param X_reg Matriz N x (K+1) com intercepto
#' @param A_vec Vetor binário de tratamento
#' @param ridge_lambda Penalidade Ridge
#' @return Lista com H1, H0, x_bar, var_factor_1, var_factor_0, idx_1, idx_0
precompute_projection = function(X_reg, A_vec, ridge_lambda = 1e-7) {
  x_bar = colMeans(X_reg)
  idx_1 = which(A_vec == 1)
  idx_0 = which(A_vec == 0)
  p = ncol(X_reg)

  calc = function(X) {
    XtX = crossprod(X)
    diag(XtX) = diag(XtX) + ridge_lambda
    XtX_inv = solve(XtX)
    H = XtX_inv %*% t(X)
    var_factor = as.numeric(t(x_bar) %*% XtX_inv %*% x_bar)
    list(H = H, var_factor = var_factor)
  }

  m1 = calc(X_reg[idx_1, , drop = FALSE])
  m0 = calc(X_reg[idx_0, , drop = FALSE])

  list(
    H1 = m1$H, H0 = m0$H,
    var_factor_1 = m1$var_factor, var_factor_0 = m0$var_factor,
    x_bar = x_bar,
    idx_1 = idx_1, idx_0 = idx_0,
    N1 = length(idx_1), N0 = length(idx_0),
    p = p
  )
}

#' Estima FATE vetorizado usando projeções pré-computadas
#'
#' Para cada t_j: theta(t_j) = x_bar' (H1 y1_j - H0 y0_j)
#'
#' @param Y_wide Matriz N x T de desfechos
#' @param proj Output de precompute_projection
#' @return Lista com fate_est (vetor T)
estimate_fate_fast = function(Y_wide, proj) {
  Y1 = Y_wide[proj$idx_1, , drop = FALSE]  # N1 x T
  Y0 = Y_wide[proj$idx_0, , drop = FALSE]  # N0 x T

  beta1 = proj$H1 %*% Y1  # p x T
  beta0 = proj$H0 %*% Y0  # p x T

  fate_est = as.numeric(proj$x_bar %*% (beta1 - beta0))  # vetor T

  list(fate_est = fate_est)
}

#' Estima FATE para uma amostra bootstrap (dados brutos, sem suavização)
#'
#' Reamostra indivíduos, re-estima via OLS ponto a ponto.
#' Usa a mesma lógica de projeção mas com índices reamostrados.
#'
#' @param Y_wide Matriz N x T de desfechos RUIDOSOS
#' @param X_reg Matriz N x (K+1) com intercepto
#' @param A_vec Vetor binário de tratamento
#' @param n_boot Número de amostras bootstrap
#' @param ridge_lambda Penalidade Ridge
#' @return Matriz n_boot x T de curvas bootstrap
run_bootstrap_raw = function(Y_wide, X_reg, A_vec, n_boot, ridge_lambda = 1e-7) {
  N = nrow(Y_wide)
  T_pts = ncol(Y_wide)
  x_bar = colMeans(X_reg)
  p = ncol(X_reg)

  boot_curves = matrix(NA_real_, nrow = n_boot, ncol = T_pts)

  for (b in seq_len(n_boot)) {
    idx = sample(N, N, replace = TRUE)

    A_b = A_vec[idx]
    i1 = which(A_b == 1)
    i0 = which(A_b == 0)

    if (length(i1) < p + 5 || length(i0) < p + 5) next

    X1 = X_reg[idx[i1], , drop = FALSE]
    X0 = X_reg[idx[i0], , drop = FALSE]
    Y1 = Y_wide[idx[i1], , drop = FALSE]
    Y0 = Y_wide[idx[i0], , drop = FALSE]

    # Projeções para esta amostra
    XtX1 = crossprod(X1); diag(XtX1) = diag(XtX1) + ridge_lambda
    XtX0 = crossprod(X0); diag(XtX0) = diag(XtX0) + ridge_lambda

    H1 = tryCatch(solve(XtX1) %*% t(X1), error = function(e) NULL)
    H0 = tryCatch(solve(XtX0) %*% t(X0), error = function(e) NULL)
    if (is.null(H1) || is.null(H0)) next

    beta1 = H1 %*% Y1  # p x T
    beta0 = H0 %*% Y0  # p x T

    boot_curves[b, ] = as.numeric(x_bar %*% (beta1 - beta0))
  }

  boot_curves
}

#' Constrói bandas simultâneas bootstrap
#'
#' D_b = max_t |theta_b(t) - theta_hat(t)|
#' q_{1-alpha} = quantil de {D_b}
#' Banda: theta_hat(t) +/- q_{1-alpha}
#'
#' @param theta_hat Vetor T (estimativa pontual)
#' @param boot_curves Matriz B x T de curvas bootstrap
#' @param alpha Nível de significância
#' @return Lista com lower, upper (vetores T)
build_simultaneous_band = function(theta_hat, boot_curves, alpha = 0.05) {
  valid = complete.cases(boot_curves)
  if (sum(valid) < 50) {
    return(list(
      lower = rep(-Inf, length(theta_hat)),
      upper = rep(Inf, length(theta_hat))
    ))
  }

  devs = apply(boot_curves[valid, , drop = FALSE], 1,
                function(r) max(abs(r - theta_hat), na.rm = TRUE))
  q = quantile(devs, 1 - alpha, na.rm = TRUE)

  list(lower = theta_hat - q, upper = theta_hat + q)
}
