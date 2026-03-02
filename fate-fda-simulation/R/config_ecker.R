library(data.table)

#' Configura os parâmetros baseado em Ecker et al.
#' @param N Número de indivíduos
#' @param T_pts Número de pontos no tempo
#' @return Lista com todos os argumentos necessários para
#'  generate_sequential_dgp
get_ecker_configuration = function(N, T_pts) {
  # 1. Definição das 9 Covariáveis (Tabela D.1 - Men)
  # K = 9
  # 1. year_first_inc (Normal)
  # 2. prev_income (Normal)
  # 3. children (Binomial)
  # 4. sec_edu (Binomial)
  # 5. siblings (Poisson)
  # 6. parents_inc (Normal)
  # 7. foreign_par (Binomial)
  # 8. par_sec_edu (Binomial)
  # 9. par_tert_edu (Binomial)

  # Matriz de Parâmetros Basais (N x Parametros)
  # Como a função base é constante (p[1]), a matriz é N x 1 (o próprio valor)

  params_list = list(
    matrix(rnorm(N, -2.3, 1.2), ncol=1), # 1
    matrix(rnorm(N, 8.3, 4.1), ncol=1), # 2
    matrix(rbinom(N, 1, 0.019), ncol=1), # 3
    matrix(rbinom(N, 1, 0.710), ncol=1), # 4
    matrix(rpois(N, 2.3), ncol=1), # 5
    matrix(rnorm(N, 11.8, 0.9), ncol=1), # 6
    matrix(rbinom(N, 1, 0.023), ncol=1), # 7
    matrix(rbinom(N, 1, 0.297), ncol=1), # 8
    matrix(rbinom(N, 1, 0.112), ncol=1)  # 9
  )

  # Templates: funções constantes no tempo (X é estático)
  func_templates = rep(list(function(t, p) p[1]), 9)

  # 2. Propensão / Confundimento (X -> A)
  # Tabela E.1 (Men)
  # Ordem: Intercepto, X1, ..., X9
  theta_vec = c(
    -6.194, -0.052, -0.001, 0.020, -0.097, -0.056,
    0.510, 0.568, 0.674, 0.851
    )

  # Matriz S x (K+1). Como S=1 (um tratamento), é 1 x 10.
  theta_matrix = matrix(theta_vec, nrow = 1, ncol = 10)

  # 3. Efeitos no Desfecho (X -> Y(t))
  # Funções beta_X(t). Ecker não dá a fórmula, mas sabemos que variam no tempo.
  # Criamos funções suaves razoáveis para simular o efeito na renda.

  concurrent_funcs = list(
    function(t) -0.02, # 1. year (constante)
    function(t) 0.5 * exp(-2 * t), # 2. prev_inc (decai)
    function(t) -0.1 + 0.2*(t-0.5)^2, # 3. children (U-shape)
    function(t) 0.2 * t, # 4. sec_edu (cresce)
    function(t) 0.01 * sin(10*t), # 5. siblings (ruído cíclico)
    function(t) 0.05, # 6. par_inc (constante)
    function(t) -0.1, # 7. foreign (constante)
    function(t) 0.1 * exp(-t), # 8. par_sec (decai)
    function(t) 0.15 * exp(-t) # 9. par_tert (decai)
  )

  # 4. Desligando Dinâmica Sequencial
  K = 9
  S = 1
  # X -> X (zero)
  beta_funcs_matrix = matrix(list(), K, K)
  # A -> X (zero feedback)
  gamma_matrix = matrix(list(function(t,s) 0), K, S)

  # Intercepto do Desfecho (Tendência global de renda)
  baseline_y = function(t) 2 * log(t + 0.5)

  return(
    list(
      baseline_func_templates = func_templates,
      baseline_params_list = params_list,
      theta_matrix = theta_matrix,
      concurrent_effect_funcs = concurrent_funcs,
      beta_funcs_matrix = beta_funcs_matrix,
      gamma_matrix = gamma_matrix,
      baseline_y_func = baseline_y,
      static_indices = 1:9 # Todos os Xs são estáticos
      )
    )
}
