#' @title Gera Desfecho Funcional com Efeitos Concorrentes das Covariáveis
#' @description A partir de um conjunto de dados simulado (geralmente de
#'   \code{generate_confounder_pointwise}), calcula a trajetória do desfecho
#'   Y(t). Este modelo assume que o efeito das covariáveis X(t) sobre o
#'   desfecho Y(t) é **concorrente**, ou seja, depende apenas do valor da
#'   covariável no mesmo ponto do tempo t.
#'
#' @param dt O data.table no formato "long" gerado por
#'   \code{generate_confounder_pointwise}.
#' @param base_function Uma função que define a trajetória de base do desfecho,
#'   \eqn{\mu_Y(t)}.
#' @param effects_function_list Lista de funções para o efeito direto do
#'   tratamento A(s) em Y(t).
#' @param covariates_concurrent_effect_list Lista de funções que definem o
#'   efeito concorrente \eqn{\delta_k(t)} de cada covariável X_k(t) em Y(t).
#' @param noise_sd Desvio padrão do ruído branco adicionado ao desfecho.
#' @param return_format Formato do output: "long", "wide", ou "both".
#' @param treatment_times Vetor com os tempos exatos de tratamento.
#'
#' @return Um objeto com os dados de entrada e as colunas adicionais para os
#'   desfechos observados (y) e contrafactuais (a_y, another_a_y).
#'
generate_y_concurrent = function(
    dt, 
    base_function,
    effects_function_list,
    covariates_concurrent_effect_list, 
    noise_sd = 1, 
    return_format = c("long", "wide", "both"),
    treatment_times = NULL
) {
  requireNamespace("data.table")
  return_format = match.arg(return_format)
  
  required_cols = c("id", "t", "A", "A_a", "A_another")
  if (!all(required_cols %in% names(dt))) {
    stop("dt deve conter colunas: id, t, A, A_a, A_another")
  }
  
  # nomes das covariáveis
  Xnames = grep("^X\\d+", names(dt), value = TRUE)
  a_Xnames = grep("^a_X\\d+", names(dt), value = TRUE)
  another_Xnames = grep("^another_a_X\\d+", names(dt), value = TRUE)
  num_covariates = length(Xnames)
  
  if (num_covariates != length(covariates_concurrent_effect_list)) {
    stop("Número de covariáveis incompatível com 
         covariates_concurrent_effect_list")
  }
  
  # preparação
  dt = data.table::copy(data.table::setorder(dt, id, t))
  ids = unique(dt$id)
  times = sort(unique(dt$t))
  n_ids = length(ids)
  n_t = length(times)
  
  # tempos de tratamento assumidos: s = 1,2,...,S
  S = length(effects_function_list)
  
  # Se os tempos não forem fornecidos, assume o comportamento antigo 1, 2, ...
  if (is.null(treatment_times)) {
    treatment_times = 1:S
  }
  
  if (length(treatment_times) != S) {
    stop("O comprimento de 'treatment_times' deve ser igual ao número de 
         funções de efeito.")
  }
  
  # mapa: para cada s encontre o índice da coluna de 'times' mais próxima de s
  treatment_col_indices = sapply(
    treatment_times, function(s) which.min(abs(times - s))
  )
  
  # alerta se o s não coincide bem com o grid:
  dist_warn = abs(times[treatment_col_indices] - treatment_times)
  if (any(dist_warn > (min(diff(times))/2 + .Machine$double.eps))) {
    warning("Algum tempo de tratamento não coincide exatamente com time_grid; 
            usando índice mais próximo.")
  }
  
  # formar matrizes (assume dt ordenado por id,t e completo)
  A_mat = matrix(dt$A, nrow = n_ids, ncol = n_t, byrow = TRUE)
  A_a_mat = matrix(dt$A_a, nrow = n_ids, ncol = n_t, byrow = TRUE)
  A_another_mat = matrix(
    dt$A_another, 
    nrow = n_ids, 
    ncol = n_t, 
    byrow = TRUE
  )
  
  X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  a_X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  another_X_arr = array(NA, dim = c(n_ids, n_t, num_covariates))
  
  for (k in seq_len(num_covariates)) {
    X_arr[,,k] = matrix(
      dt[[Xnames[k]]], 
      nrow = n_ids, 
      ncol = n_t, 
      byrow = TRUE
    )
    
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
  
  # base
  base_vals = matrix(base_function(times), n_ids, n_t, byrow = TRUE)
  
  # pré-computar beta_mat (S x T) com zeros para t < s
  beta_mat = matrix(0, nrow = S, ncol = n_t)
  for (s_idx in seq_len(S)) {
    # Usa o valor real do tempo de tratamento, não o índice
    s_val = treatment_times[s_idx] 
    for (t_idx in seq_len(n_t)) {
      t_val = times[t_idx]
      if (t_val >= s_val) {
        beta_mat[s_idx, t_idx] = effects_function_list[[s_idx]](t_val, s_val)
      }
    }
  }
  
  # pré-computar gammas para covariáveis (vetor por tempo)
  gamma_list = lapply(
    covariates_concurrent_effect_list, 
    function(gf) sapply(times, gf)
  )
  
  # gerar RUÍDO uma vez e reutilizar (mesmo ruído para todos os cenários)
  noise_common = matrix(
    stats::rnorm(n_ids * n_t, 0, noise_sd), 
    nrow = n_ids, 
    ncol = n_t
  )
  
  # função que calcula Y|A (n x T) e o array X (n x T x K)
  compute_y = function(A_full, X_array) {
    # extrair colunas de A correspondentes aos tempos de tratamento s=1..S
    A_treated = A_full[, treatment_col_indices, drop = FALSE] # n x S
    # efeitos diretos: A_treated %*% beta_mat  -> n x T
    direct_effects = A_treated %*% beta_mat
    # efeitos de covariáveis concorrentes
    covar_effects = matrix(0, n_ids, n_t)
    for (k in seq_len(num_covariates)) {
      gamma_k = matrix(gamma_list[[k]], n_ids, n_t, byrow = TRUE)
      covar_effects = covar_effects + X_array[,,k] * gamma_k
    }
    base_vals + direct_effects + covar_effects + noise_common
  }
  
  # gerar os três Ys
  y_mat = compute_y(A_mat, X_arr)
  a_y_mat = compute_y(A_a_mat, a_X_arr)
  another_a_y_mat = compute_y(A_another_mat, another_X_arr)
  
  # efeito pedido: diferença entre os dois regimes contrafactuais
  effect_a_vs_another = a_y_mat - another_a_y_mat
  
  # calcular efeitos diretos "teóricos" (lista de vetores 1xT, 0 antes de s)
  true_effect_list = lapply(seq_len(S), function(s_idx) {
    s_val = treatment_times[s_idx]
    vec = numeric(n_t)
    idxs = which(times >= s_val)
    if (length(idxs) > 0) {
      vec[idxs] = sapply(
        times[idxs], 
        function(t) effects_function_list[[s_idx]](t, s_val)
      )
    }
    vec
  })
  
  # acrescentar colunas ao dt (long)
  dt[, y := as.vector(t(y_mat))]
  dt[, a_y := as.vector(t(a_y_mat))]
  dt[, another_a_y := as.vector(t(another_a_y_mat))]
  dt[, effect_a_vs_another := as.vector(t(effect_a_vs_another))]
  
  # retorno
  switch(
    return_format,
    long = dt,
    wide = list(
      y_mat = y_mat,
      a_y_mat = a_y_mat,
      another_a_y_mat = another_a_y_mat,
      effect_a_vs_another = effect_a_vs_another,
      t_vec = times,
      true_effect_list = true_effect_list
    ),
    both = list(
      long = dt,
      wide = list(
        y_mat = y_mat,
        a_y_mat = a_y_mat,
        another_a_y_mat = another_a_y_mat,
        effect_a_vs_another = effect_a_vs_another,
        t_vec = times,
        true_effect_list = true_effect_list
      )
    )
  )
}
