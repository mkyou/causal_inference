#' @title Gera Desfecho Funcional
#' @description A partir dos dados gerados por `generate_sequential_dgp`,
#'   calcula a trajetória do desfecho Y(t). A função foi refatorada para
#'   ser consistente com o novo PGD, usando inputs estruturados e operações
#'   eficientes de data.table.
#'
#' @param dgp_data O data.table no formato "long" gerado por
#'   `generate_sequential_dgp`.
#' @param baseline_y_func Uma função que define a trajetória de base do
#'   desfecho, mu_Y(t).
#' @param direct_effect_funcs Lista de S funções para o efeito direto do
#'   s-ésimo tratamento A(s) em Y(t).
#' @param concurrent_effect_funcs Lista de K funções que definem o
#'   efeito concorrente delta_k(t) de cada covariável X_k(t) em Y(t).
#' @param treatment_times Vetor com os tempos exatos de tratamento.
#' @param noise_sd Desvio padrão do ruído branco adicionado ao desfecho.
#'
#' @return O data.table de entrada com as colunas adicionais para os
#'   desfechos factuais (y) e contrafactuais (a_y, b_y).

generate_functional_outcome = function(
    dgp_data,
    baseline_y_func,
    direct_effect_funcs,
    concurrent_effect_funcs,
    treatment_times,
    noise_sd = 1
) {
  
  # --- Setup e Validações ---
  dt = copy(dgp_data)
  K = length(grep("^X[0-9]+$", names(dt), value = TRUE))
  S = length(treatment_times)
  
  if (length(direct_effect_funcs) != S)
    stop("direct_effect_funcs deve ter um elemento para cada tratamento.")
  if (length(concurrent_effect_funcs) != K)
    stop("concurrent_effect_funcs deve ter um elemento para cada covariável.")
  
  # --- Cálculos Vetorizados Diretamente no data.table ---
  
  # 1. Adiciona o termo de base e o ruído comum
  # O ruído é o mesmo para os cenários factual e contrafactual
  dt[, y_base := baseline_y_func(t)]
  dt[, y_noise := rnorm(.N, 0, noise_sd)]
  
  # 2. Calcula os efeitos concorrentes totais das covariáveis
  # Para cada cenário (factual, regime_a, regime_b)
  
  # Nomes das colunas para cada cenário
  factual_x = paste0("X", 1:K)
  regime_a_x = paste0("a_X", 1:K)
  regime_b_x = paste0("b_X", 1:K)
  
  # Calcula a soma dos efeitos: sum_k( delta_k(t) * X_k(t) )
  dt[, cov_effect_factual := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, factual_x)
  ), .SDcols = factual_x]
  
  dt[, cov_effect_a := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, regime_a_x)
  ), .SDcols = regime_a_x]
  
  dt[, cov_effect_b := rowSums(
    mapply(function(f, x) f(t) * .SD[[x]], concurrent_effect_funcs, regime_b_x)
  ), .SDcols = regime_b_x]
  
  # 3. Calcula os efeitos diretos totais dos tratamentos
  # Inicializa as colunas de efeito
  dt[, direct_effect_factual := 0.0]
  dt[, direct_effect_a := 0.0]
  dt[, direct_effect_b := 0.0]
  
  # Itera sobre cada um dos S tratamentos
  for (s_idx in 1:S) {
    s_time = treatment_times[s_idx]
    
    # Efeito Factual (depende da coluna A gerada no PGD)
    # Apenas indivíduos tratados em s_time recebem o efeito
    ids_treated_factually = dt[t == s_time & A == 1, id]
    if (length(ids_treated_factually) > 0) {
      dt[id %in% ids_treated_factually & t >= s_time, 
         direct_effect_factual := direct_effect_factual + 
           direct_effect_funcs[[s_idx]](t, s_time)]
    }
    
    # Efeito do Regime A (usa a coluna A_a gerada no PGD)
    if (unique(dt[t == s_time, A_a]) == 1) {
      dt[t >= s_time, direct_effect_a := direct_effect_a +
           direct_effect_funcs[[s_idx]](t, s_time)]
    }
    
    # Efeito do Regime B (usa a coluna A_b gerada no PGD)
    if (unique(dt[t == s_time, A_b]) == 1) {
      dt[t >= s_time, direct_effect_b := direct_effect_b +
           direct_effect_funcs[[s_idx]](t, s_time)]
    }
  }
  
  # 4. Soma todos os componentes para criar o desfecho final
  dt[, y := y_base + cov_effect_factual + direct_effect_factual + y_noise]
  dt[, a_y := y_base + cov_effect_a + direct_effect_a + y_noise]
  dt[, b_y := y_base + cov_effect_b + direct_effect_b + y_noise]
  
  # 5. Calcula o efeito causal total e limpa as colunas intermediárias
  dt[, effect_a_vs_b := a_y - b_y]
  
  cols_to_remove = grep("(_base|_noise|_effect)", names(dt), value = TRUE)
  dt[, (cols_to_remove) := NULL]
  
  return(dt)
}