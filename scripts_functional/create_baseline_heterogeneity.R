#' @title Cria Heterogeneidade de Base para Trajetórias Individuais
#' @description Gera uma lista de funções de trajetória de base individuais.
#' Cada função incorpora heterogeneidade ao adicionar desvios aleatórios tanto
#' no intercepto (nível) quanto na inclinação (amplitude) de uma trajetória
#' média populacional.
#'
#' @param n Número de indivíduos a serem simulados.
#' @param mean_funcs Uma lista de funções, onde cada função representa a
#'   trajetória média populacional para uma covariável.
#' @param intercept_sd Desvio padrão do ruído gaussiano adicionado ao
#'   intercepto de cada indivíduo. Controla a variabilidade de nível entre os
#'   indivíduos.
#' @param slope_sd Desvio padrão do fator multiplicativo gaussiano aplicado à
#'   trajetória média. Controla a variabilidade na amplitude/forma entre os
#'   indivíduos.
#'
#' @return Uma lista (uma por covariável), contendo listas (uma por
#'   indivíduo) de funções que definem as trajetórias de base individuais.
#'   Esta estrutura é projetada para ser usada diretamente pelo argumento
#'   \code{baseline_function_list} em \code{generate_confounder_pointwise}.
#'

create_baseline_heterogeneity = function(
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
    rnorm(n, 1, slope_sd), 
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