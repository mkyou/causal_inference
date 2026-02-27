#' @title Gera Desfecho Funcional
#' @description Calcula a trajetória do desfecho Y(t) a partir dos dados do DGP.
#' Utiliza lógica de "tempo mais próximo" para evitar erros de ponto flutuante.
#'
#' @param dgp_data O data.table "long" gerado por `generate_sequential_dgp`.
#' @param baseline_y_func Função para a trajetória de base do desfecho, mu_Y(t).
#' @param direct_effect_funcs Lista de S funções para o efeito direto de A(s).
#' @param concurrent_effect_funcs Lista de K funções para o efeito de X_k(t).
#' @param treatment_times Vetor com os tempos de tratamento.
#' @param noise_sd Desvio padrão do ruído do desfecho.
#'
#' @return O data.table de entrada com colunas de desfecho adicionadas.

library(data.table)
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
  
  # Extrai o grid de tempos que REALMENTE existe nos dados
  actual_time_grid = sort(unique(dt$t))
  
  if (length(direct_effect_funcs) != S) 
    stop("Erro no tamanho de direct_effect_funcs")
  
  # 1. Base e Ruído
  dt[, y_base := baseline_y_func(t)]
  dt[, y_noise := rnorm(.N, 0, noise_sd)]
  
  # 2. Efeitos Concorrentes (X -> Y)
  factual_x = paste0("X", 1:K)
  regime_a_x = paste0("a_X", 1:K)
  regime_b_x = paste0("b_X", 1:K)
  
  dt[, cov_effect_factual := rowSums(
    mapply(
      function(f, x) f(t) * .SD[[x]], 
      concurrent_effect_funcs, factual_x
      )
    ), .SDcols = factual_x]
  
  dt[, cov_effect_a := rowSums(
    mapply(
      function(f, x) f(t) * .SD[[x]], 
      concurrent_effect_funcs, regime_a_x
      )
    ), .SDcols = regime_a_x]
  
  dt[, cov_effect_b := rowSums(
    mapply(
      function(f, x) f(t) * .SD[[x]], 
      concurrent_effect_funcs, 
      regime_b_x
      )
    ), .SDcols = regime_b_x]
  
  # 3. Efeitos Diretos (A -> Y)
  dt[, direct_effect_factual := 0.0]
  dt[, direct_effect_a := 0.0]
  dt[, direct_effect_b := 0.0]
  
  for (s_idx in 1:S) {
    target_s_time = treatment_times[s_idx]
    
    # Encontra o valor no grid mais próximo do tempo de tratamento teórico 
    # Isso evita falhas de comparação como 20.0 != 20.00000001
    idx_closest = which.min(abs(actual_time_grid - target_s_time))
    s_time_data = actual_time_grid[idx_closest]
  
    col_trat_factual = paste0("Trat_", s_idx)
    col_trat_a = paste0("Trat_", s_idx, "_a")
    col_trat_b = paste0("Trat_", s_idx, "_b")
    
    # Aplica efeito factual (Quem tem A=1 no tempo de tratamento encontrado)
    # Usamos s_time_data para filtrar, 
    # mas passamos target_s_time para a função matemática
    ids_treated = dt[t == s_time_data & get(col_trat_factual) == 1, id]
    
    if (length(ids_treated) > 0) {
      # Aplica para todo o futuro (t >= s_time_data)
      dt[
        id %in% ids_treated,
        direct_effect_factual := direct_effect_factual + 
          direct_effect_funcs[[s_idx]](t, target_s_time)
        ]
    }
    
    # Aplica efeito Regime A
    # Verifica se o regime impõe tratamento
    is_treated_a = unique(dt[t == s_time_data, get(col_trat_a)])
    if (length(is_treated_a) > 0 && is_treated_a == 1) {
      dt[
        , 
        direct_effect_a := direct_effect_a + 
          direct_effect_funcs[[s_idx]](t, target_s_time)
        ]
    }
    
    # Aplica efeito Regime B
    is_treated_b = unique(dt[t == s_time_data, get(col_trat_b)])
    if (length(is_treated_b) > 0 && is_treated_b == 1) {
      dt[
        , 
        direct_effect_b := direct_effect_b + 
          direct_effect_funcs[[s_idx]](t, target_s_time)
        ]
    }
  }
  
  # 4. Soma Final
  dt[, y := y_base + cov_effect_factual + direct_effect_factual + y_noise]
  dt[, a_y := y_base + cov_effect_a + direct_effect_a + y_noise]
  dt[, b_y := y_base + cov_effect_b + direct_effect_b + y_noise]
  
  # Efeito Verdadeiro (Fundamental para validação)
  dt[, effect_a_vs_b := a_y - b_y]
  
  # Limpeza
  cols_to_remove = grep("(_base|_noise|_effect)", names(dt), value = TRUE)
  cols_to_remove = setdiff(cols_to_remove, "effect_a_vs_b")
  dt[, (cols_to_remove) := NULL]
  
  return(dt)
}