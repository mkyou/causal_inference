#' @title Gera Dados Longitudinais com PGD Sequencial e Funcional
#' @description Gera dados de covariáveis funcionais.
#'
#' @param n_subjects Número de indivíduos a simular.
#' @param time_grid Vetor numérico definindo a grade de tempo.
#' @param treatment_times Vetor numérico dos tempos onde o tratamento ocorre.
#'
#' @param baseline_func_templates Lista de K "templates" de função.
#' @param baseline_params_list Lista de K matrizes de parâmetros para as bases.
#'
#' @param beta_funcs_matrix Matriz KxK de funções para o efeito de Xj(t-1) em
#' Xk(t).
#' @param gamma_matrix Matriz KxS de funções para o efeito do tratamento A(s)
#' em Xk(t).
#'
#' @param theta_matrix Matriz Sx(K+1) de coeficientes para o modelo de
#' propensão.
#' @param zeta_funcs_matrix Matriz SxK opcional para o termo histórico no
#' modelo de propensão.
#'
#' @param noise_sd_vec Vetor de K desvios padrão para o ruído estocástico.
#'
#' @param regime_a,regime_b Vetores definindo regimes contrafactuais.
#' @param static_covariate_indices Vetor numérico opcional indicando os índices
#'   (de 1 a K) das covariáveis que devem ser estáticas no tempo (exceto por
#'   feedback ou efeitos cruzados).
#'
#' @return Um data.table longo com trajetórias factuais e contrafactuais.

library(data.table)

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
    regime_b = NULL,
    static_covariate_indices = NULL,
    propensity_min = 0.001,
    propensity_max = 0.999
) {

  K = length(baseline_func_templates)
  S = length(treatment_times)
  total_times = length(time_grid)

  if (is.null(regime_a)) regime_a = rep(1, S)
  if (is.null(regime_b)) regime_b = rep(0, S)

  # Encontra o índice exato no grid para cada tempo de tratamento
  treatment_indices = integer(S)
  for(i in 1:S) {
    treatment_indices[i] = which.min(abs(time_grid - treatment_times[i]))
  }

  # Inicialização
  X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  a_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  b_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))

  A_mat = matrix(0, n_subjects, S)
  A_a_mat = matrix(0, n_subjects, S)
  A_b_mat = matrix(0, n_subjects, S)

  # Pré-cálculos
  baseline_matrices = lapply(1:K, function(k) {
    sapply(time_grid, function(t) {
      sapply(1:n_subjects, function(i) {
        baseline_func_templates[[k]](t, baseline_params_list[[k]][i, ])
      })
    })
  })

  noise_matrices = lapply(1:K, function(k) {
    if (k %in% static_covariate_indices) {
      noise_t0 = rnorm(n_subjects, 0, noise_sd_vec[k])
      matrix(noise_t0, nrow = n_subjects, ncol = total_times)
    } else {
      matrix(rnorm(n_subjects * total_times, 0, noise_sd_vec[k]),
             nrow = n_subjects, ncol = total_times)
    }
  })

  # Inicialização t=1
  t_idx = 1
  current_time = time_grid[t_idx]

  for (k in 1:K) {
    base_val = baseline_matrices[[k]][, 1] + noise_matrices[[k]][, 1]
    X_list[[k]][, 1] = base_val
    a_X_list[[k]][, 1] = base_val
    b_X_list[[k]][, 1] = base_val
  }

  # Tratamento no t=1 (Usando Índice)
  if (t_idx %in% treatment_indices) {
    # Recupera qual tratamento é (pode haver mais de um, mas assumimos ordem)
    # Qual s_idx corresponde a este t_idx?
    s_idx_hits = which(treatment_indices == t_idx)

    for(s_idx in s_idx_hits) {
      X_curr = sapply(X_list, function(X) X[, 1])
      linear_predictor = cbind(1, X_curr) %*% theta_matrix[s_idx, ]

      prop = plogis(linear_predictor)
      prop = pmax(propensity_min, pmin(propensity_max, prop))

      A_mat[, s_idx] = rbinom(n_subjects, 1, prop)
      A_a_mat[, s_idx] = regime_a[s_idx]
      A_b_mat[, s_idx] = regime_b[s_idx]
    }
  }

  # --- Loop Sequencial ---
  for (t_idx in 2:total_times) {
    current_time = time_grid[t_idx]

    # Tratamento t > 1 (Usando Índice)
    if (t_idx %in% treatment_indices) {
      s_idx_hits = which(treatment_indices == t_idx)

      for(s_idx in s_idx_hits) {
        prev_t_idx = t_idx - 1
        X_prev = sapply(X_list, function(X) X[, prev_t_idx])
        linear_predictor = cbind(1, X_prev) %*% theta_matrix[s_idx, ]

        prop = plogis(linear_predictor)
        prop = pmax(propensity_min, pmin(propensity_max, prop))

        A_mat[, s_idx] = rbinom(n_subjects, 1, prop)
        A_a_mat[, s_idx] = regime_a[s_idx]
        A_b_mat[, s_idx] = regime_b[s_idx]
      }
    }

    # Evolução Covariáveis
    for (k in 1:K) {
      new_val = baseline_matrices[[k]][, t_idx]

      cross_ef = 0; cross_a = 0; cross_b = 0
      for (j in 1:K) {
        if (k == j || is.null(beta_funcs_matrix[[k, j]])) next
        beta_val = beta_funcs_matrix[[k, j]](current_time)
        cross_ef = cross_ef + beta_val * X_list[[j]][, t_idx-1]
        cross_a = cross_a + beta_val * a_X_list[[j]][, t_idx-1]
        cross_b = cross_b + beta_val * b_X_list[[j]][, t_idx-1]
      }

      fb_ef = 0; fb_a = 0; fb_b = 0
      for (s_idx in 1:S) {
        # Feedback só ocorre APÓS o tratamento ter sido atribuído
        # Comparamos índices para segurança absoluta
        treat_idx_s = treatment_indices[s_idx]
        if (t_idx > treat_idx_s) {
          s_time = treatment_times[s_idx]
          gamma_val = gamma_matrix[[k, s_idx]](current_time, s_time)
          fb_ef = fb_ef + gamma_val * A_mat[, s_idx]
          fb_a = fb_a  + gamma_val * A_a_mat[, s_idx]
          fb_b = fb_b  + gamma_val * A_b_mat[, s_idx]
        }
      }

      noise = noise_matrices[[k]][, t_idx]

      X_list[[k]][, t_idx] = new_val + cross_ef + fb_ef + noise
      a_X_list[[k]][, t_idx] = new_val + cross_a + fb_a + noise
      b_X_list[[k]][, t_idx] = new_val + cross_b + fb_b + noise
    }
  }

  # Formatação Output
  dt_long = data.table(
    id = rep(1:n_subjects, each = total_times),
    t = rep(time_grid, n_subjects)
  )

  for(s in 1:S) {
    col_A = A_mat[, s]
    col_A_a = A_a_mat[, s]
    col_A_b = A_b_mat[, s]

    dt_long[, paste0("Trat_", s) := rep(col_A, each = total_times)]
    dt_long[, paste0("Trat_", s, "_a") := rep(col_A_a, each = total_times)]
    dt_long[, paste0("Trat_", s, "_b") := rep(col_A_b, each = total_times)]
  }

  for (k in 1:K) {
    dt_long[, paste0("X", k) := as.vector(t(X_list[[k]]))]
    dt_long[, paste0("a_X", k) := as.vector(t(a_X_list[[k]]))]
    dt_long[, paste0("b_X", k) := as.vector(t(b_X_list[[k]]))]
  }

  return(dt_long)
}
