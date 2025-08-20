#' @title Cria Trajetórias de Base Heterogêneas (Avançado)
#' @description Gera curvas de base individuais a partir de funções médias,
#' permitindo heterogeneidade no intercepto, escala e inclinação.
#'
#' @param n Número de indivíduos.
#' @param mean_funcs Lista de funções de trajetória média.
#' @param intercept_sd Desvio padrão do intercepto aleatório (variação no nível).
#' @param scale_sd Desvio padrão da escala aleatória (variação na amplitude). A média é 1.
#' @param slope_sd Desvio padrão da inclinação aleatória (variação na tendência).
#' @return Uma lista de listas de funções, pronta para generate_confounder_pointwise.

create_heterogeneous_baselines_advanced = function(
    n, 
    mean_funcs, 
    intercept_sd = 0.2, 
    scale_sd = 0.1, 
    slope_sd = 0.05
) {
  num_covariates = length(mean_funcs)
  
  # Gera n coeficientes aleatórios para cada covariável
  random_intercepts = replicate(num_covariates, rnorm(n, 0, intercept_sd), simplify = FALSE)
  random_scales = replicate(num_covariates, rnorm(n, 1, scale_sd), simplify = FALSE)
  random_slopes = replicate(num_covariates, rnorm(n, 0, slope_sd), simplify = FALSE)
  
  # Para cada covariável k...
  lapply(1:num_covariates, function(k) {
    # ...cria n funções, uma para cada indivíduo i
    lapply(1:n, function(i) {
      
      # Armazena os parâmetros do indivíduo para não recalculá-los
      b = random_intercepts[[k]][i]
      a = random_scales[[k]][i]
      c = random_slopes[[k]][i]
      mean_f = mean_funcs[[k]]
      
      # A função do indivíduo i é uma combinação dos efeitos aleatórios
      function(t) {
        # escala * média(t) + intercepto + inclinação * t
        a * mean_f(t) + b + c * t
      }
    })
  })
}