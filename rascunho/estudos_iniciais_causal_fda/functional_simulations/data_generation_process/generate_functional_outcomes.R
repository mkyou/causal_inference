#' @title Gera Desfecho Funcional
#' @description Calcula a trajetória do desfecho Y(t) a partir dos dados do DGP.
#'
#' @param dgp_data O data.table "long" gerado por `generate_sequential_dgp`.
#' @param baseline_y_func Função para a trajetória de base do desfecho, mu_Y(t).
#' @param direct_effect_funcs Lista de S funções para o efeito direto de A(s).
#' @param concurrent_effect_funcs Lista de K funções para o efeito de X_k(t).
#' @param treatment_times Vetor com os tempos de tratamento.
#' @param noise_sd Desvio padrão do ruído do desfecho.
#'
#' @return O data.table de entrada com colunas de desfecho adicionadas.

generate_functional_outcome = function(
    dgp_data,
    baseline_y_func,
    direct_effect_funcs,
    concurrent_effect_funcs,
    treatment_times,
    noise_sd = 1
) {
  dt = copy(dgp_data)
  K = length(grep("^X[0-9]+$", names(dt), value = TRUE))
  S = length(treatment_times)
  
  if (length(direct_effect_funcs) != S)
    stop("direct_effect_funcs deve ter um elemento para cada tratamento.")
  if (length(concurrent_effect_funcs) != K)
    stop("concurrent_effect_funcs deve ter um elemento para cada covariável.")
  
  dt[, y_base := baseline_y_func(t)]
  dt[, y_noise := rnorm(.N, 0, noise_sd)]
  
  factual_x = paste0("X", 1:K)
  regime_a_x = paste0("a_X", 1:K)
  regime_b_x = paste0("b_X", 1:K)
  
  dt[, cov_effect_factual := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, factual_x)
  ), .SDcols = factual_x]
  
  dt[, cov_effect_a := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, regime_a_x)
  ), .SDcols = regime_a_x]
  
  dt[, cov_effect_b := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, regime_b_x)
  ), .SDcols = regime_b_x]
  
  dt[, direct_effect_factual := 0.0]
  dt[, direct_effect_a := 0.0]
  dt[, direct_effect_b := 0.0]
  
  for (s_idx in 1:S) {
    s_time = treatment_times[s_idx]
    
    ids_treated_factually = dt[t == s_time & A == 1, id]
    if (length(ids_treated_factually) > 0) {
      dt[
        id %in% ids_treated_factually & t >= s_time,
        direct_effect_factual := direct_effect_factual +
          direct_effect_funcs[[s_idx]](t, s_time)
      ]
    }
    if (unique(dt[t == s_time, A_a]) == 1) {
      dt[
        t >= s_time,
        direct_effect_a := direct_effect_a +
          direct_effect_funcs[[s_idx]](t, s_time)
      ]
    }
    if (unique(dt[t == s_time, A_b]) == 1) {
      dt[
        t >= s_time,
        direct_effect_b := direct_effect_b +
          direct_effect_funcs[[s_idx]](t, s_time)
      ]
    }
  }
  
  dt[, y := y_base + cov_effect_factual + direct_effect_factual + y_noise]
  dt[, a_y := y_base + cov_effect_a + direct_effect_a + y_noise]
  dt[, b_y := y_base + cov_effect_b + direct_effect_b + y_noise]
  
  dt[, effect_a_vs_b := a_y - b_y]
  
  cols_to_remove = grep("(_base|_noise|_effect)", names(dt), value = TRUE)
  dt[, (cols_to_remove) := NULL]
  
  return(dt)
}