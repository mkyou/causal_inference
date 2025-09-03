#' @title Cria Heterogeneidade de Intercepto
#' @description Gera curvas de base individuais com desvios aleatórios
#' no intercepto (nível) de cada indivíduo.
#' @param n Número de indivíduos.
#' @param mean_funcs Lista de funções de trajetória média.
#' @param intercept_sd Desvio padrão do intercepto aleatório.
#' @param slope_sd Desvio padrão da inclinação aleatória.
#' @return Uma lista de listas de funções para generate_confounder_pointwise.

create_intercept_heterogeneity = function(
    n, 
    mean_funcs,
    slope_sd = 0.2,
    intercept_sd = 0.2
) {
  # Esta função não precisa de scale_sd, removemos para clareza
  num_covariates = length(mean_funcs)
  
  random_intercepts = replicate(
    num_covariates, 
    rnorm(n, 0, intercept_sd), 
    simplify = FALSE
  )
  
  random_slope = replicate(
    num_covariates, 
    rnorm(n, 1, intercept_sd), 
    simplify = FALSE
  )

  lapply(1:num_covariates, function(k) {
    lapply(1:n, function(i) {
      
      a = random_slope[[k]][i]
      b = random_intercepts[[k]][i]
      mean_f = mean_funcs[[k]]
      
      # PGD: f_media(t) + intercepto_i
      function(t) {
        a * mean_f(t) + b
      }
    })
  })
}