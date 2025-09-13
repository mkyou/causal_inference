# ===================================================================
# FUNÇÃO GENERALIZADA PARA AVALIAÇÃO DE MODELOS
# ===================================================================
library(purrr)

#' @title Executa a Avaliação de um Modelo em Múltiplos Cenários
#' @description Uma função de alta ordem que recebe uma função de modelo
#'   e a aplica a uma lista de cenários, avaliando os resultados.
#'
#' @param all_scenarios_data Lista nomeada com os data.tables gerados.
#' @param scenario_configs A lista de configurações usada para gerar os dados.
#' @param model_function A função de ajuste do modelo a ser testada.
#'   Ex: `fit_ecker_luna`. Deve aceitar `scenario_data`,
#'   `scenario_config`, `regime_a`, `regime_b`.
#' @param model_name Uma string para nomear o modelo nos resultados.
#'
#' @return Uma lista nomeada com os resultados da avaliação para cada cenário.

run_model_evaluation = function(
    all_scenarios_data,
    scenario_configs,
    model_function,
    model_name
) {
  
  cat(paste("\n--- Iniciando avaliação para o modelo:", model_name, "---\n"))
  
  # `imap` itera sobre os dados (.x) e os nomes (.y) dos cenários
  all_results = purrr::imap(
    all_scenarios_data,
    function(scenario_data, scenario_name) {
      
      cat(paste("\n>>> Processando cenário:", scenario_name, "<<<\n"))
      
      # Extrai a configuração específica para este cenário
      config = scenario_configs[[scenario_name]]
      
      S = length(config$treatment_times)
      
      # Regime A é "tratar sempre" (um vetor de 1s de comprimento S)
      # Regime B é "nunca tratar" (um vetor de 0s de comprimento S)
      regime_a_dynamic = rep(1, S)
      regime_b_dynamic = rep(0, S)
      
      # 1. Estima o efeito usando a função de modelo fornecida
      effect_estimate_dt = model_function(
        scenario_data = scenario_data,
        scenario_config = config,
        regime_a = config$regime_a %||% regime_a_dynamic,
        regime_b = config$regime_b %||% regime_b_dynamic
      )
      
      # 2. Avalia a estimativa contra a verdade
      assessment = assess_effect_estimation(
        estimation_dt = effect_estimate_dt,
        true_data = scenario_data
      )
      
      return(assessment)
    }
  )
  
  return(all_results)
}

# Operador auxiliar para definir valores padrão
`%||%` = function(a, b) {
  if (is.null(a)) b else a
}