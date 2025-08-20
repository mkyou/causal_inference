generate_y_historical = function(
    dt,
    base_function,
    effects_function_list,
    covariates_historical_effect_list,
    noise_sd = 1,
    return_format = c("long", "wide", "both")
) {
  requireNamespace("data.table")
  return_format = match.arg(return_format)
  
  required_cols = c("id", "t", "A", "A_a", "A_another")
  if (!all(required_cols %in% names(dt))) {
    stop("dt deve conter colunas: id, t, A, A_a, A_another")
  }
  
  # Nomes das covariáveis
  Xnames = grep("^X\\d+", names(dt), value = TRUE)
  a_Xnames = grep("^a_X\\d+", names(dt), value = TRUE)
  another_Xnames = grep("^another_a_X\\d+", names(dt), value = TRUE)
  num_covariates = length(Xnames)
  
  if (num_covariates != length(covariates_historical_effect_list)) {
    stop("Número de covariáveis incompatível com 
         covariates_historical_effect_list")
  }
  
  # Preparação
  dt = data.table::copy(data.table::setorder(dt, id, t))
  ids = unique(dt$id)
  times = sort(unique(dt$t))
  n_ids = length(ids)
  n_t = length(times)
  
  # Tempos de tratamento
  S = length(effects_function_list)
  treatment_times = seq_len(S)
  treatment_col_indices = sapply(
    treatment_times, function(s) which.min(abs(times - s))
  )
  
  dist_warn = abs(times[treatment_col_indices] - treatment_times)
  if (any(dist_warn > (min(diff(times))/2 + .Machine$double.eps))) {
    warning("Algum tempo de tratamento não coincide exatamente com time_grid; 
            usando índice mais próximo.")
  }
  
  # Formar matrizes e arrays (idêntico à função anterior)
  A_mat = matrix(dt$A, nrow = n_ids, ncol = n_t, byrow = TRUE)
  A_a_mat = matrix(dt$A_a, nrow = n_ids, ncol = n_t, byrow = TRUE)
  A_another_mat = matrix(dt$A_another, nrow = n_ids, ncol = n_t, byrow = TRUE)
  
  X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  a_X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  another_X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  
  for (k in seq_len(num_covariates)) {
    X_arr[,,k] = matrix(dt[[Xnames[k]]], nrow = n_ids, ncol = n_t, byrow = TRUE)
    a_X_arr[,,k] = matrix(
      dt[[a_Xnames[k]]], 
      nrow = n_ids, 
      ncol = n_t, 
      byrow = TRUE
    )
    another_X_arr[,,k] = matrix(
      dt[[another_Xnames[k]]], 
      nrow = n_ids, 
      ncol = n_t, 
      byrow = TRUE
    )
  }
  
  # Base e efeito direto do tratamento (idêntico à função anterior)
  base_vals = matrix(base_function(times), n_ids, n_t, byrow = TRUE)
  beta_mat = matrix(0, nrow = S, ncol = n_t)
  for (s_idx in seq_len(S)) {
    s_val = treatment_times[s_idx]
    for (t_idx in seq_len(n_t)) {
      t_val = times[t_idx]
      if (t_val >= s_val) {
        beta_mat[s_idx, t_idx] = effects_function_list[[s_idx]](t_val, s_val)
      }
    }
  }
  
  # Para cada covariável k, criamos uma matriz Gamma_k (n_t x n_t)
  # Onde Gamma_k[j, i] = gamma_k(t_i, t_j) para j <= i, e 0 caso contrário.
  Gamma_mat_list = lapply(
    covariates_historical_effect_list, 
    function(gamma_func) {
      gamma_mat = matrix(0, n_t, n_t)
      for (i in 1:n_t) {
        for (j in 1:i) { # Linha u_j
          # Influência de X(u_j) em Y(t_i)
          gamma_mat[j, i] = gamma_func(times[i], times[j])
        }
      }
    return(gamma_mat)
    }
  )
  
  # Ruído 
  noise_common = matrix(
    stats::rnorm(n_ids * n_t, 0, noise_sd), 
    nrow = n_ids, 
    ncol = n_t
  )
  
  # Função auxiliar que calcula Y | A, X
  compute_y = function(A_full, X_array) {
    A_treated = A_full[, treatment_col_indices, drop = FALSE]
    direct_effects = A_treated %*% beta_mat
    
    # --- Cálculo do efeito das covariáveis ---
    # Agora, usamos produto matricial para integrar a história
    covar_effects = matrix(0, n_ids, n_t)
    for (k in seq_len(num_covariates)) {
      X_k_mat = X_array[,,k, drop=FALSE]
      Gamma_k_mat = Gamma_mat_list[[k]]
      # (n_ids x n_t) %*% (n_t x n_t) -> (n_ids x n_t)
      historical_k_effects = X_k_mat %*% Gamma_k_mat
      covar_effects = covar_effects + historical_k_effects
    }
    base_vals + direct_effects + covar_effects + noise_common
  }
  
  # Gerar os três Ys (lógica idêntica à função anterior)
  y_mat = compute_y(A_mat, X_arr)
  a_y_mat = compute_y(A_a_mat, a_X_arr)
  another_a_y_mat = compute_y(A_another_mat, another_X_arr)
  effect_a_vs_another = a_y_mat - another_a_y_mat
  
  # Acrescentar colunas ao dt (long)
  dt[, y := as.vector(t(y_mat))]
  dt[, a_y := as.vector(t(a_y_mat))]
  dt[, another_a_y := as.vector(t(another_a_y_mat))]
  dt[, effect_a_vs_another := as.vector(t(effect_a_vs_another))]
  
  # Retorno (lógica idêntica à função anterior)
  switch(
    return_format,
    long = dt,
    wide = list(
      y_mat = y_mat, a_y_mat = a_y_mat, another_a_y_mat = another_a_y_mat,
      effect_a_vs_another = effect_a_vs_another, t_vec = times
    ),
    both = list(
      long = dt,
      wide = list(
        y_mat = y_mat, a_y_mat = a_y_mat, another_a_y_mat = another_a_y_mat,
        effect_a_vs_another = effect_a_vs_another, t_vec = times
      )
    )
  )
}