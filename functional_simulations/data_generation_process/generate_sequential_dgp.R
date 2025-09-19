library(data.table)

#' @title Gera Dados Longitudinais com PGD Sequencial e Funcional
#' @description Gera dados de covariáveis funcionais para múltiplos
#'   indivíduos com base em um PGD explicitamente sequencial (Markoviano).
#'   Este PGD é ideal para testar métodos como a G-formula, onde as
#'   suposições de trocabilidade sequencial são válidas por desenho.
#'
#' @param n_subjects Número de indivíduos a simular.
#' @param time_grid Vetor numérico definindo a grade de tempo.
#' @param treatment_times Vetor numérico dos tempos onde o tratamento ocorre.
#'
#' @param baseline_funcs Lista de K funções, uma para cada covariável,
#'   definindo a trajetória média da população mu_k(t).
#' @param beta_funcs_matrix Matriz KxK de funções. A célula [k, j] contém a
#'   função beta_kj(t) que descreve o efeito de Xj(t-1) em Xk(t).
#'   A diagonal é ignorada. Use NULL para especificar ausência de efeito.
#' @param gamma_matrix Matriz KxS de funções. A célula [k, s] contém a
#'   função gamma_ks(t, s_time) que descreve o efeito do s-ésimo tratamento
#'   (ocorrido em s_time) na covariável k no tempo t.
#'
#' @param theta_matrix Matriz Sx(K+1) de coeficientes numéricos para o
#'   modelo de propensão logística (baseado em X(s-1)). A 1ª coluna é o
#'   intercepto.
#' @param theta_integral_matrix Matriz SxK opcional (default NULL). Se
#'   fornecida, adiciona um termo baseado na integral da história da
#'   covariável ao modelo de propensão.
#'
#' @param intercept_sd_vec Vetor de K desvios padrão para os interceptos
#'   aleatórios b_ik (heterogeneidade persistente).
#' @param noise_sd_vec Vetor de K desvios padrão para o ruído estocástico
#'   epsilon_ik(t).
#'
#' @param regime_a,regime_b Vetores definindo regimes de tratamento
#'   contrafactuais a serem simulados.
#'
#' @return Um data.table em formato longo contendo as trajetórias
#'   factuais e contrafactuais.

generate_sequential_dgp = function(
    n_subjects,
    time_grid,
    treatment_times,
    baseline_funcs,
    beta_funcs_matrix,
    gamma_matrix,
    theta_matrix,
    theta_integral_matrix = NULL,
    intercept_sd_vec,
    noise_sd_vec,
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
  if (!is.null(theta_integral_matrix)) {
    if (nrow(theta_integral_matrix) != S || ncol(theta_integral_matrix) != K)
      stop("theta_integral_matrix deve ter dimensões S x K.")
  }
  if (length(intercept_sd_vec) != K || length(noise_sd_vec) != K)
    stop("intercept_sd_vec e noise_sd_vec devem ter comprimento K.")
  
  if (is.null(regime_a)) regime_a = rep(1, S)
  if (is.null(regime_b)) regime_b = rep(0, S)
  
  # --- Inicialização das Estruturas de Dados ---
  X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  a_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  b_X_list = lapply(1:K, function(k) matrix(0, n_subjects, total_times))
  
  A = A_a = A_b = matrix(0, n_subjects, S)
  colnames(A) = colnames(A_a) = colnames(A_b) = paste0("A", 1:S)
  
  # --- Estágio 1: Gerar Heterogeneidade Invariante no Tempo ---
  b_ik_matrix = sapply(intercept_sd_vec, function(sd) {
    rbinom(n_subjects, 1, 1)
  })
  
  # --- Estágio 2: Inicialização em t=0 ---
  t0 = time_grid[1]
  for (k in 1:K) {
    noise_t0 = rnorm(n_subjects, 0, noise_sd_vec[k])
    base_val = baseline_funcs[[k]](t0) + b_ik_matrix[, k] + noise_t0
    
    X_list[[k]][, 1] = base_val
    a_X_list[[k]][, 1] = base_val
    b_X_list[[k]][, 1] = base_val
  }
  
  # --- Estágio 3: Loop de Simulação Sequencial ---
  for (t_idx in 2:total_times) {
    current_time = time_grid[t_idx]
    
    # --- Atribuição de Tratamento ---
    if (current_time %in% treatment_times) {
      s_idx = which(treatment_times == current_time)
      prev_t_idx = t_idx - 1
      
      X_prev = sapply(X_list, function(X) X[, prev_t_idx])
      linear_predictor = cbind(1, X_prev) %*% theta_matrix[s_idx, ]
      
      if (!is.null(theta_integral_matrix)) {
        X_integral = sapply(X_list, function(X) {
          rowSums(X[, 1:prev_t_idx, drop = FALSE]) * dt
        })
        integral_predictor = rowSums(
          X_integral * matrix(
          theta_integral_matrix[s_idx, ],
          nrow = n_subjects, ncol = K, byrow = TRUE
          )
        )
        linear_predictor = linear_predictor + integral_predictor
      }
      
      propensity = plogis(linear_predictor)
      A[, s_idx] = rbinom(n_subjects, 1, propensity)
      A_a[, s_idx] = regime_a[s_idx]
      A_b[, s_idx] = regime_b[s_idx]
    }
    
    # --- Evolução das Covariáveis ---
    for (k in 1:K) {
      new_val = baseline_funcs[[k]](current_time) + b_ik_matrix[, k]
      
      cross_effect_factual = 0
      cross_effect_a = 0
      cross_effect_b = 0
      
      for (j in 1:K) {
        if (k == j || is.null(beta_funcs_matrix[[k, j]])) next
        
        beta_val = beta_funcs_matrix[[k, j]](current_time)
        
        cross_effect_factual = cross_effect_factual + 
          beta_val * X_list[[j]][, t_idx-1]
        
        cross_effect_a = cross_effect_a + 
          beta_val * a_X_list[[j]][, t_idx-1]
        
        cross_effect_b = cross_effect_b + 
          beta_val * b_X_list[[j]][, t_idx-1]
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
      a_X_list[[k]][, t_idx] = new_val + cross_effect_a + 
        feedback_a + noise
      b_X_list[[k]][, t_idx] = new_val + cross_effect_b + 
        feedback_b + noise
    }
  }
  
  # --- Formatação do Output para data.table longo ---
  dt_long = data.table(
    id = rep(1:n_subjects, each = total_times),
    t = rep(time_grid, n_subjects)
  )
  
  # Adiciona tratamentos factuais e contrafactuais ao data.table
  A_long = as.data.table(A)
  A_long[, id := .I]
  A_long = melt(
    A_long, id.vars = "id", variable.name = "A_idx", value.name = "A"
  )
  # (Esta parte pode ser otimizada, mas é clara)
  # Cria uma coluna A geral para o tempo de tratamento correspondente
  A_map = data.table(
    A_idx = colnames(A), treatment_time = treatment_times
  )
  A_long[A_map, on = .(A_idx), t := i.treatment_time]
  
  dt_long[A_long, on = .(id, t), A := i.A]
  dt_long[is.na(A), A := 0] # Assume 0 se não for um tempo de tratamento
  
  # Adiciona as covariáveis
  for (k in 1:K) {
    dt_long[, paste0("X", k) := as.vector(t(X_list[[k]]))]
    dt_long[, paste0("a_X", k) := as.vector(t(a_X_list[[k]]))]
    dt_long[, paste0("b_X", k) := as.vector(t(b_X_list[[k]]))]
  }
  
  # Adiciona as colunas de regime para referência
  # (Isso é um pouco redundante, mas útil para validação)
  dt_long[, A_a := 0]; dt_long[, A_b := 0]
  for(s_idx in 1:S) {
    s_time = treatment_times[s_idx]
    dt_long[t == s_time, A_a := regime_a[s_idx]]
    dt_long[t == s_time, A_b := regime_b[s_idx]]
  }
  
  return(dt_long)
}
