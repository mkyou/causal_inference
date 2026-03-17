library(data.table)

#' @title Gera dados funcionais para um problema causal (X -> A -> Y)
#' @description Esta função é o "motor" do DGP. Ela recebe
#'   as funções de parâmetros (alpha, betas) como argumentos.
#'
#' @param N (Integer) Número de indivíduos (tamanho da amostra).
#' @param T (Integer) Número de pontos de tempo (discretização).
#' @param t_min (Numeric) Valor mínimo do grid de tempo.
#' @param t_max (Numeric) Valor máximo do grid de tempo.
#' @param t_A_time (Numeric) O ponto no tempo em que o
#'   tratamento A é aplicado.
#' @param propensity_intercept (Numeric) Intercepto da regressão
#'   logística para P(A=1 | X).
#' @param noise_sd (Numeric) Desvio padrão do ruído branco.
#' @param seed (Integer) Semente aleatória.
#' @param fun_alpha (Function) Uma função `f(t)` que retorna o
#'   intercepto funcional alpha(t).
#' @param fun_betas (list) Uma lista de funções.
#'   - O primeiro elemento, `fun_betas[[1]]`, DEVE ser a função
#'     beta_0(t) para o tratamento A. Ela DEVE aceitar
#'     dois argumentos: `(t, t_A)`.
#'   - Os elementos seguintes, `fun_betas[[2]]`, ...,
#'     `fun_betas[[k+1]]`, são as funções beta_j(t) para as
#'     covariáveis X_j. Elas DEVEM aceitar um argumento: `(t)`.
#'
#' @return (list) Uma lista contendo:
#'   $static_data: Um data.table (N linhas) com id, A, X1...Xk.
#'   $curve_data: Um data.table (N*T linhas) com id, t, Y_obs,
#'     Y_cf_A1, Y_cf_A0, e as versões "true" (sem ruído).
#'   $params_curves: Um data.table (T linhas) com os
#'     parâmetros funcionais usados (alpha_t, beta0_t, ...).

generate_fdg_causal = function(
    N = 2000,
    T = 2000,
    t_min = 0.0,
    t_max = 1.0,
    t_A_time = 0.5,
    propensity_intercept = 0.0,
    noise_sd = 0.5,
    seed = 20251111,
    fun_alpha,
    fun_betas
) {
  
  set.seed(seed)
  
  # Validação de inputs e definição de parâmetros
  
  # Inferir k (número de covariáveis)
  k = length(fun_betas) - 1
  if (k < 0) {
    stop("`fun_betas` must contain at least one function (beta0).")
  }
  
  # Definição de grid temporal
  t_grid = seq(t_min, t_max, length.out = T)
  
  if (t_A_time < t_min | t_A_time > t_max) {
    warning("t_A_time está fora do grid [t_min, t_max].")
  }
  
  # Geração de dados estáticos (Covariáveis X e Tratamento A)
  
  dt_static = data.table(id = 1:N)
  X_cols = c()
  
  if (k > 0) {
    # Geração de k covariáveis X (N x k)
    X_matrix = matrix(rnorm(N * k), ncol = k)
    X_cols = paste0("X", 1:k)
    colnames(X_matrix) = X_cols
    
    dt_static = cbind(dt_static, as.data.table(X_matrix))
    
    # Geração de Tratamento A | X (Modelo Logístico)
    # Definição de pesos aleatórios para o efeito de X em A
    gamma_slopes = runif(k, -0.5, 0.5)
    
    # Cálculo do preditor linear (eta)
    eta = propensity_intercept +
      as.matrix(dt_static[, ..X_cols]) %*% gamma_slopes
  } else {
    # Se k=0, A depende apenas do intercepto
    eta = propensity_intercept
  }
  
  # Cálculo da propensão p(x) = P(A=1 | X)
  propensity = 1 / (1 + exp(-eta))
  
  # Geração de A
  dt_static[, A := rbinom(.N, 1, propensity)]
  
  # Geração de Funções de Parâmetros (alpha(t), beta(t))
  
  # Avaliação das funções fornecidas pelo usuário no grid
  alpha_t = fun_alpha(t_grid)
  
  # A função beta0 (a primeira da lista) recebe t_A_time
  beta0_t = fun_betas[[1]](t_grid, t_A_time)
  
  dt_params = data.table(t = t_grid, alpha_t, beta0_t)
  beta_cols = c()
  
  if (k > 0) {
    beta_matrix = matrix(0, nrow = T, ncol = k)
    for (j in 1:k) {
      # As funções beta_j (da segunda em diante) recebem apenas t
      beta_matrix[, j] = fun_betas[[j + 1]](t_grid)
    }
    beta_cols = paste0("beta", 1:k)
    colnames(beta_matrix) = beta_cols
    
    dt_params = cbind(dt_params, as.data.table(beta_matrix))
  }
  
  # Criação de tabela no formato longo (N * T linhas)
  
  # Cross-join de indivíduos e tempo
  dt_full = CJ(id = 1:N, t = t_grid)
  
  # Merge dos dados estáticos (A, X)
  dt_full = dt_static[dt_full, on = "id"]
  
  # Merge dos parâmetros funcionais (alpha, betas)
  dt_full = dt_params[dt_full, on = "t"]
  
  # Cálculo de efeitos e curvas (sem ruído)
  
  # Cálculo do efeito combinado das covariáveis: sum(beta_j(t) * X_j)
  cov_effect = 0
  if (k > 0) {
    X_matrix_long = as.matrix(dt_full[, ..X_cols])
    beta_matrix_long = as.matrix(dt_full[, ..beta_cols])
    
    # rowSums é vetorizado e rápido
    cov_effect = rowSums(X_matrix_long * beta_matrix_long)
  }
  
  # Cálculo dos desfechos reais (sem ruído)
  
  # true para Y Observado
  dt_full[, true_Y_obs := alpha_t + (beta0_t * A) + cov_effect]
  
  # true para Contrafactual Y | A = 1
  dt_full[, true_Y_cf_A1 := alpha_t + (beta0_t * 1) + cov_effect]
  
  # true para Contrafactual Y | A = 0
  dt_full[, true_Y_cf_A0 := alpha_t + (beta0_t * 0) + cov_effect]
  
  # Adição de ruído branco
  
  # Gerar UMA realização de ruído por (id, t)
  dt_full[, noise := rnorm(.N, 0, noise_sd)]
  
  dt_full[, Y_obs := true_Y_obs + noise]
  dt_full[, Y_cf_A1 := true_Y_cf_A1 + noise]
  dt_full[, Y_cf_A0 := true_Y_cf_A0 + noise]
  
  # Preparação de saída
  
  static_data = dt_static
  
  curve_data = dt_full[, .(
    id,
    t,
    Y_obs,
    Y_cf_A1,
    Y_cf_A0,
    true_Y_obs,
    true_Y_cf_A1,
    true_Y_cf_A0
  )]
  
  return(
    list(
      static_data = static_data,
      curve_data = curve_data,
      params_curves = dt_params
    )
  )
}


#' @title Gera UMA OBSERVAÇÃO de dados funcionais para um tempo 't'
#' @description Esta é uma versão "em tempo real" que gera os dados
#'   para um único ponto de tempo, dado um conjunto estático de
#'   indivíduos.
#'
#' @param t_scalar (Numeric) O ponto de tempo único a ser avaliado.
#' @param static_data (data.table) O data.table $static_data
#'   retornado por `generate_fdg_causal`.
#' @param fun_alpha (Function) Função alpha(t).
#' @param fun_betas (list) Lista de funções beta.
#' @param t_A_time (Numeric) Tempo de aplicação do tratamento.
#' @param noise_sd (Numeric) Desvio padrão do ruído.
#'
#' @return (data.table) Um data.table (N linhas) com id, t, Y_obs, ...
#'
generate_fdg_causal_step = function(
    t_scalar,
    static_data,
    fun_alpha,
    fun_betas,
    t_A_time,
    noise_sd = 0.5
) {
  
  N = nrow(static_data)
  k = length(fun_betas) - 1
  
  # Avaliação dos parâmetros em t_scalar
  alpha_t = fun_alpha(t_scalar)
  beta0_t = fun_betas[[1]](t_scalar, t_A_time)
  
  # Cálculo do efeito das covariáveis em t_scalar ---
  cov_effect = 0
  X_cols = c()
  
  if (k > 0) {
    X_cols = paste0("X", 1:k)
    # Avalaiação dos betas das covariáveis em t_scalar
    betas_j_t = vapply(
      2:(k + 1),
      function(i) fun_betas[[i]](t_scalar),
      numeric(1)
    )
    
    X_matrix = as.matrix(static_data[, ..X_cols])
    cov_effect = as.vector(X_matrix %*% betas_j_t)
  }
  
  # Cálculo das curvas em t_scalar
  # Usar `copy()` é vital para não modificar 'static_data'
  dt_step = copy(static_data)
  
  # Cálculo dos valores reais
  dt_step[, true_Y_obs := alpha_t + (beta0_t * A) + cov_effect]
  dt_step[, true_Y_cf_A1 := alpha_t + (beta0_t * 1) + cov_effect]
  dt_step[, true_Y_cf_A0 := alpha_t + (beta0_t * 0) + cov_effect]
  
  # Adição de ruído
  dt_step[, noise := rnorm(N, 0, noise_sd)]
  dt_step[, Y_obs := true_Y_obs + noise]
  dt_step[, Y_cf_A1 := true_Y_cf_A1 + noise]
  dt_step[, Y_cf_A0 := true_Y_cf_A0 + noise]
  
  # Adição de tempo à DT
  dt_step[, t := t_scalar]
  
  # Saída
  return(dt_step[, .(
    id,
    t,
    A,
    Y_obs,
    Y_cf_A1,
    Y_cf_A0,
    true_Y_obs,
    true_Y_cf_A1,
    true_Y_cf_A0
  )])
}