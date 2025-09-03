#' @title Gera Covariáveis Longitudinais com Confundimento Variante no Tempo
#' @description Simula as trajetórias de múltiplas covariáveis funcionais ao
#'   longo do tempo para 'n' indivíduos. O processo permite a inclusão de
#'   feedback do tratamento (tratamentos passados afetam covariáveis futuras) e
#'   a atribuição de tratamento dependente da história das covariáveis,
#'   induzindo o confundimento.
#'
#' @param n Número de indivíduos.
#' @param effects_function_list Lista de listas de funções que definem o efeito
#'   do tratamento A(s) nas covariáveis em tempos t > s.
#' @param p_treat Probabilidade de tratamento em cada ponto de tempo, usada
#'   quando \code{propensity_function_list} é NULL.
#' @param noise_sd Desvio padrão do ruído branco adicionado à trajetória de
#'   cada covariável em cada ponto de tempo.
#' @param time_grid Vetor numérico definindo a grade de tempo para a simulação.
#' @param baseline_function_list Lista de listas de funções de base individuais,
#'   tipicamente gerada por \code{create_baseline_heterogeneity}.
#' @param propensity_function_list Lista de funções que calculam a
#'   probabilidade de tratamento com base na história das covariáveis.
#' @param a_seq,another_a_seq Vetores definindo regimes de tratamento
#'   contrafactuais a serem simulados.
#' @param return_format Formato do output: "long", "wide", ou "both".
#' @param treatment_times Vetor especificando os tempos exatos onde o
#'   tratamento pode ser administrado.
#'
#' @return Um objeto contendo os dados simulados no formato especificado.
#'   Inclui as trajetórias observadas e contrafactuais para tratamentos (A) e
#'   covariáveis (X).
#'
generate_confounder_pointwise = function(
    n, 
    effects_function_list, 
    p_treat, 
    noise_sd = 1,
    time_grid = NULL, 
    baseline_function_list = NULL,
    propensity_function_list = NULL, 
    a_seq = NULL,
    another_a_seq = NULL, 
    return_format = c("long","wide","both"),
    treatment_times = NULL
) {
  
  # --- Input Validations ---
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) 
    stop("n deve ser um inteiro positivo")
  if (!is.list(effects_function_list) || length(effects_function_list)==0) 
    stop("effects_function_list deve ser uma lista não vazia")
  if (!is.null(baseline_function_list)) {
    if (
      !is.list(baseline_function_list) || 
      any(
        sapply(
          baseline_function_list, 
          function(sub) !is.list(sub) || length(sub) != n
        )
      )
    )
      stop("baseline_function_list deve ser uma lista de n listas de funções")
  }
  
  # --- Parameter Setup ---
  return_format = match.arg(return_format)
  num_covariates = length(effects_function_list)
  total_treatment_times = length(effects_function_list[[1]])
  
  # Logic to handle the new treatment_times parameter
  if (is.null(treatment_times)) {
    treatment_times = 1:total_treatment_times
  }
  if (length(treatment_times) != total_treatment_times) {
    stop("O comprimento de 'treatment_times' deve ser igual ao número de 
         funções de efeito/propensão.")
  }
  
  if (is.null(time_grid)) time_grid = 0:ceiling(max(treatment_times) + 1)
  total_times = length(time_grid)
  if (is.null(a_seq)) a_seq = rep(1, total_treatment_times)
  if (is.null(another_a_seq)) another_a_seq = rep(0, total_treatment_times)
  
  # --- Data Structure Initialization (no changes here) ---
  Z_list = replicate(
    num_covariates, 
    matrix(rnorm(n * total_times, 0, noise_sd), n, total_times), 
    simplify = FALSE
  )
  A = A_a = A_another = matrix(0, n, total_times)
  X_list = a_X_list = another_a_X_list = lapply(
    1:num_covariates, 
    function(k) matrix(0, n, total_times)
  )
  
  # Pre-computation of effects (no changes here)
  effects_mat_list = lapply(1:num_covariates, function(k){
    effects_mat = matrix(0, total_treatment_times, total_times)
    for (s_idx in 1:total_treatment_times){
      s_val = treatment_times[s_idx]
      for (t_idx in 1:total_times){
        t_val = time_grid[t_idx]
        if (t_val >= s_val) {
          # Pass the actual treatment time s_val to the effect function
          effects_mat[s_idx, t_idx] = effects_function_list[[k]][[s_idx]](
            t_val, 
            s_val
          )
        }
      }
    }
    effects_mat
  })
  
  # --- Initialization at t=1 (vectorized) ---
  for (k in 1:num_covariates){
    base_vals_t0 = sapply(
      baseline_function_list[[k]], 
      function(f) f(time_grid[1])
    )
    X_list[[k]][, 1] = base_vals_t0 + Z_list[[k]][, 1]
    a_X_list[[k]][, 1] = X_list[[k]][, 1]
    another_a_X_list[[k]][, 1] = X_list[[k]][, 1]
  }
  
  # --- Dynamic Process (t>=2) ---
  for (t_index in 2:total_times) {
    current_time = time_grid[t_index]
    
    # --- Treatment Assignment ---
    if (current_time %in% treatment_times) {
      s_index = which(treatment_times == current_time) # Gets the position
      prop_func = if (!is.null(propensity_function_list)) 
        propensity_function_list[[s_index]] else NULL
      
      A[, t_index] = if (!is.null(prop_func)) {
        rbinom(n, 1, prop_func(X_list, t_index))
      } else {
        rbinom(n, 1, p_treat)
      }
      
      A_a[, t_index] = rep(a_seq[s_index], n)
      A_another[, t_index] = rep(another_a_seq[s_index], n)
    }
    
    # --- Covariate Evolution ---
    treated_cols_nat = which(
      time_grid < current_time & 
        time_grid %in% treatment_times
    )
    
    for (k in 1:num_covariates) {
      base_vals_t = sapply(
        baseline_function_list[[k]], 
        function(f) f(current_time)
      )
      
      if (length(treated_cols_nat) > 0) {
        # The index of the treatment in our list (e.g., 1st treatment, ...)
        s_indices = sapply(
          time_grid[treated_cols_nat], 
          function(tt) which(treatment_times == tt)
        )
        
        effects_mat = effects_mat_list[[k]]
        effects_vec_t = effects_mat[s_indices, t_index, drop = FALSE]
        
        # Vectorized calculation remains the same and is highly efficient
        effects_obs_nat = A[, treated_cols_nat, drop = FALSE] %*% effects_vec_t
        effects_a = A_a[, treated_cols_nat, drop = FALSE] %*% effects_vec_t
        effects_another = A_another[, treated_cols_nat, drop = FALSE] %*% 
          effects_vec_t
        
      } else {
        effects_obs_nat = effects_a = effects_another = rep(0, n)
      }
      
      X_list[[k]][, t_index] = base_vals_t + effects_obs_nat + 
        Z_list[[k]][, t_index]
      
      a_X_list[[k]][, t_index] = base_vals_t + effects_a + 
        Z_list[[k]][, t_index]
      
      another_a_X_list[[k]][, t_index] = base_vals_t + effects_another + 
        Z_list[[k]][, t_index]
    }
  }
  
  # --- Data Formatting and Return (no changes here) ---
  dt = data.table(
    id = rep(1:n,each=total_times),
    t = rep(time_grid,n),
    A = as.vector(t(A)),
    A_a = as.vector(t(A_a)),
    A_another = as.vector(t(A_another))
  )
  
  for (k in 1:num_covariates){
    dt[,paste0("X",k):=as.vector(t(X_list[[k]]))]
    dt[,paste0("a_X",k):=as.vector(t(a_X_list[[k]]))]
    dt[,paste0("another_a_X",k):=as.vector(t(another_a_X_list[[k]]))]
  }
  
  switch(
    return_format,
    "long" = dt,
    "wide" = list(
      A_mat = A, A_a_mat = A_a, A_another_mat = A_another,
      X_list = X_list, a_X_list = a_X_list,
      another_a_X_list = another_a_X_list,
      t_vec = time_grid
    ),
    "both" = list(
      long = dt,
      wide = list(
        A_mat = A, A_a_mat = A_a, A_another_mat = A_another,
        X_list = X_list, a_X_list = a_X_list,
        another_a_X_list = another_a_X_list,
        t_vec = time_grid
      )
    )
  )
}
