library(data.table)

#' Estima FATE e Intervalos de Confiança Ponto-a-Ponto
estimate_fate_pointwise = function(
    smoothed_y, covs_matrix, treatment_vec, 
    time_grid, ridge_lambda = 1e-7
    ) {
  
  n_t = length(time_grid)
  X_design = cbind(1, covs_matrix) 
  x_bar = colMeans(X_design) # Perfil médio para o efeito marginal
  
  idx_1 = which(treatment_vec == 1)
  idx_0 = which(treatment_vec == 0)
  N1 = length(idx_1)
  N0 = length(idx_0)
  p = ncol(X_design)
  
  if(N1 <= p || N0 <= p) return(NULL) # Graus de liberdade insuficientes
  
  X1 = X_design[idx_1, , drop=FALSE]
  X0 = X_design[idx_0, , drop=FALSE]
  
  # Pré-cálculo de Matrizes de Variância
  # Precisamos da diagonal de x_bar * (X'X)^-1 * x_bar' 
  # para o erro padrão do efeito médio
  
  calc_matrices = function(X) {
    XtX = crossprod(X)
    diag(XtX) = diag(XtX) + ridge_lambda
    # Inversa de X'X
    XtX_inv = tryCatch(solve(XtX), error = function(e) NULL)
    if(is.null(XtX_inv)) return(NULL)
    
    # Hat Matrix H = (X'X)^-1 X'
    H = XtX_inv %*% t(X)
    
    # Fator de variância para a média x_bar: 
    # var_factor = x_bar * (X'X)^-1 * x_bar'
    # Isso é um escalar que multiplica o sigma^2 depois
    var_factor = as.numeric(t(x_bar) %*% XtX_inv %*% x_bar)
    
    list(H = H, var_factor = var_factor)
  }
  
  m1 = calc_matrices(X1)
  m0 = calc_matrices(X0)
  
  if(is.null(m1) || is.null(m0)) 
    return(data.table(t = time_grid, fate_est = NA))
  
  # Vetores de saída
  theta_est = numeric(n_t)
  ci_lower = numeric(n_t)
  ci_upper = numeric(n_t)
  
  # Valor crítico t-student (95%)
  # Graus de liberdade conservadores (minimo dos dois grupos)
  t_crit = qt(0.975, df = min(N1, N0) - p)
  
  # Loop Ponto-a-Ponto
  for(j in 1:n_t) {
    y_t = smoothed_y[, j]
    
    # Grupo 1
    y1 = y_t[idx_1]
    beta1 = m1$H %*% y1
    res1 = y1 - (X1 %*% beta1)
    sigma2_1 = sum(res1^2) / (N1 - p)
    mu1 = sum(x_bar * beta1)
    se1_sq = sigma2_1 * m1$var_factor
    
    # Grupo 0
    y0 = y_t[idx_0]
    beta0 = m0$H %*% y0
    res0 = y0 - (X0 %*% beta0)
    sigma2_0 = sum(res0^2) / (N0 - p)
    mu0 = sum(x_bar * beta0)
    se0_sq = sigma2_0 * m0$var_factor
    
    # Efeito e IC
    theta_est[j] = mu1 - mu0
    se_diff = sqrt(se1_sq + se0_sq) # Erro padrão da diferença
    
    ci_lower[j] = theta_est[j] - t_crit * se_diff
    ci_upper[j] = theta_est[j] + t_crit * se_diff
  }
  
  return(
    data.table(
      t = time_grid, fate_est = theta_est, 
      ci_lower = ci_lower, ci_upper = ci_upper
      )
    )
}