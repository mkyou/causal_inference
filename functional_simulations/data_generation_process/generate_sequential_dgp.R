library(data.table)

#' @title Gera Dados Longitudinais com PGD Sequencial e Funcional
#' @description Gera dados de covariáveis funcionais. Esta versão utiliza
#'   parâmetros de base específicos por indivíduo para gerar trajetórias
#'   de base únicas, permitindo a simulação de confundimento não medido.
#'
#' @param n_subjects Número de indivíduos a simular.
#' @param time_grid Vetor numérico definindo a grade de tempo.
#' @param treatment_times Vetor numérico dos tempos onde o tratamento ocorre.
#'
#' @param baseline_func_templates Lista de K "templates" de função, onde cada
#'   função aceita `t` e um vetor de `params`. 
#'   Ex: `function(t, p) p[1] + p[2]*t`.
#' @param baseline_params_list Lista de K matrizes de parâmetros. Cada matriz
#'   tem dimensão `n_subjects` x `n_params` para a função correspondente.
#'
#' @param beta_funcs_matrix Matriz KxK de funções beta_kj(t) para o efeito
#'   de Xj(t-1) em Xk(t).
#' @param gamma_matrix Matriz KxS de funções gamma_ks(t) para o efeito do
#'   tratamento A(s) em Xk(t) (feedback).
#'
#' @param theta_matrix Matriz Sx(K+1) de coeficientes para o modelo de
#'   propensão (baseado em X(s-1)).
#' @param zeta_funcs_matrix Matriz SxK opcional de funções de ponderação
#'   zeta_sk(u) para a integral histórica.
#'
#' @param noise_sd_vec Vetor de K desvios padrão para o ruído estocástico.
#'
#' @param regime_a,regime_b Vetores definindo regimes contrafactuais.
#'
#' @return Um data.table longo com trajetórias factuais e contrafactuais.

generate_sequential_dgp = function(
    n_subjects,
    time_grid,
    treatment_times,
    baseline_func_templates,
    baseline_params_list,
    beta_funcs_matrix,
    gamma_matrix,
    theta_matrix,
    noise_sd_vec,
    zeta_funcs_matrix = NULL,
    regime_a = NULL,
    regime_b = NULL
) {
  
  # --- Setup e Validações ---
  K = length(baseline_funcs)
  S = length(treatment_times)
  total_times = length(time_grid)
  dt = time_grid[2] - time_grid[1]
  
  # Validações de input robustas
  if (nrow(beta_funcs_matrix) != K || ncol(beta_funcs_matrix) != K)
    stop("beta_funcs_matrix deve ter dimensões K x K.")
  if (nrow(gamma_matrix) != K || ncol(gamma_matrix) != S)
    stop("gamma_matrix deve ter dimensões K x S.")
  if (nrow(theta_matrix) != S || ncol(theta_matrix) != (K + 1))
    stop("theta_matrix deve ter dimensões S x (K+1).")
  if (length(noise_sd_vec) != K) stop("noise_sd_vec deve ter comprimento K.")
  if (length(baseline_params_list) != K)
    stop("baseline_params_list deve ter um elemento para cada covariável.")
  if (!all(sapply(baseline_params_list, nrow) == n_subjects))
    stop("As matrizes em baseline_params_list devem ter n_subjects linhas.")
  
  if (is.null(regime_a)) regime_a = rep(1, S)
  if (is.null(regime_b)) regime_b = rep(0, S)
  
  # --- Inicialização ---
  X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  a_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  b_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  A = A_a = A_b = matrix(0, n_subjects, S)
  
  # --- Estágio 1: Inicialização em t=0 ---
  t0 = time_grid[1]
  for (k in 1:K) {
    noise_t0 = rnorm(n_subjects, 0, noise_sd_vec[k])
    
    func_template = baseline_func_templates[[k]]
    params_matrix = baseline_params_list[[k]]
    
    subject_specific_baselines = sapply(1:n_subjects, function(i) {
      func_template(t0, params_matrix[i, ])
    })
    
    base_val = subject_specific_baselines + noise_t0
    
    X_list[[k]][, 1] = base_val
    a_X_list[[k]][, 1] = base_val
    b_X_list[[k]][, 1] = base_val
  }
  
  # --- Estágio 2: Loop de Simulação Sequencial ---
  for (t_idx in 2:total_times) {
    current_time = time_grid[t_idx]
    
    if (current_time %in% treatment_times) {
      s_idx = which(treatment_times == current_time)
      prev_t_idx = t_idx - 1
      
      X_prev = sapply(X_list, function(X) X[, prev_t_idx])
      linear_predictor = cbind(1, X_prev) %*% theta_matrix[s_idx, ]
      
      integral_predictor = 0
      all_integrals = sapply(1:K, function(k) {
        history_matrix = X_list[[k]][, 1:prev_t_idx, drop = FALSE]
        if (
          !is.null(zeta_funcs_matrix) &&
          !is.null(zeta_funcs_matrix[[s_idx, k]])
          ) {
          zeta_func = zeta_funcs_matrix[[s_idx, k]]
          time_weights = zeta_func(time_grid[1:prev_t_idx])
          history_matrix = history_matrix * matrix(
            time_weights,
            nrow = n_subjects,
            ncol = prev_t_idx,
            byrow = TRUE
          )
        }
        if (ncol(history_matrix) < 2) return(rep(0, n_subjects))
        y_i = history_matrix[, -1, drop = FALSE]
        y_im1 = history_matrix[, -ncol(history_matrix), drop = FALSE]
        return(rowSums(0.5 * (y_i + y_im1)) * dt)
      })
      
      integral_predictor = rowSums(
        all_integrals *
          matrix(
            theta_integral_matrix[s_idx, ],
            nrow = n_subjects, ncol = K, byrow = TRUE
          )
      )
      
      linear_predictor = linear_predictor + integral_predictor
      
      propensity = plogis(linear_predictor)
      A[, s_idx] = rbinom(n_subjects, 1, propensity)
      A_a[, s_idx] = regime_a[s_idx]
      A_b[, s_idx] = regime_b[s_idx]
    }
    
    for (k in 1:K) {
      func_template = baseline_func_templates[[k]]
      params_matrix = baseline_params_list[[k]]
      subject_specific_baselines = sapply(1:n_subjects, function(i) {
        func_template(current_time, params_matrix[i, ])
      })
      new_val = subject_specific_baselines
      
      cross_effect_factual = 0; cross_effect_a = 0; cross_effect_b = 0
      for (j in 1:K) {
        if (k == j || is.null(beta_funcs_matrix[[k, j]])) next
        beta_val = beta_funcs_matrix[[k, j]](current_time)
        cross_effect_factual = cross_effect_factual + beta_val * 
          X_list[[j]][, t_idx-1]
        cross_effect_a = cross_effect_a + beta_val * a_X_list[[j]][, t_idx-1]
        cross_effect_b = cross_effect_b + beta_val * b_X_list[[j]][, t_idx-1]
      }
      feedback_effect = 0; feedback_a = 0; feedback_b = 0
      for (s_idx in 1:S) {
        s_time = treatment_times[s_idx]
        if (current_time >= s_time) {
          gamma_val = gamma_matrix[[k, s_idx]](current_time, s_time)
          feedback_effect = feedback_effect + gamma_val * A[, s_idx]
          feedback_a = feedback_a + gamma_val * A_a[, s_idx]
          feedback_b = feedback_b + gamma_val * A_b[, s_idx]
        }
      }
      
      noise = rnorm(n_subjects, 0, noise_sd_vec[k])
      
      X_list[[k]][, t_idx] = new_val + cross_effect_factual + 
        feedback_effect + noise
      a_X_list[[k]][, t_idx] = new_val + cross_effect_a + feedback_a + noise
      b_X_list[[k]][, t_idx] = new_val + cross_effect_b + feedback_b + noise
    }
  }
  
  dt_long = data.table(
    id = rep(1:n_subjects, each = total_times), 
    t = rep(time_grid, n_subjects)
  )
  
  A_long = as.data.table(A); A_long[, id := .I]
  A_long = melt(
    A_long, id.vars = "id", 
    variable.name = "s_idx_char", value.name = "A"
  )
  A_long[, s_idx := as.numeric(gsub("V", "", s_idx_char))][
    , t := treatment_times[s_idx]
  ]
  dt_long[A_long, on=.(id,t), A := i.A]; dt_long[is.na(A), A := 0]
  
  for (k in 1:K) {
    dt_long[, paste0("X", k) := as.vector(t(X_list[[k]]))]
    dt_long[, paste0("a_X", k) := as.vector(t(a_X_list[[k]]))]
    dt_long[, paste0("b_X", k) := as.vector(t(b_X_list[[k]]))]
  }
  
  for (k in 1:K) {
    params_matrix = baseline_params_list[[k]]
    for (j in 1:ncol(params_matrix)) {
      param_col_name = paste0("X", k, "_b", j - 1)
      dt_long[
        , 
        (param_col_name) := rep(params_matrix[, j], each = total_times)
      ]
    }
  }
  
  dt_long[, A_a := 0]; dt_long[, A_b := 0]
  for(s_idx in 1:S) {
    s_time = treatment_times[s_idx]
    dt_long[t == s_time, A_a := regime_a[s_idx]]
    dt_long[t == s_time, A_b := regime_b[s_idx]]
  }
  
  return(dt_long)
}
