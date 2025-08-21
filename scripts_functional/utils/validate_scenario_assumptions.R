#' @title Valida Numericamente as Hipóteses de um Cenário
#' @description Calcula métricas que verificam se o cenário foi gerado como
#' esperado (feedback, dinamismo das covariáveis, etc.).
#' @param scenario_data Um único elemento da lista `all_scenarios_data`.
#' @return Uma data.table com um resumo da validação.

validate_scenario_assumptions = function(scenario_data) {
  require(data.table)
  dt = scenario_data$outcome$long
  
  # 1. Checagem de Feedback (A -> X)
  # Se o tratamento afeta X, a_X e another_a_X devem ser diferentes.
  feedback_check = max(abs(dt$a_X1 - dt$another_a_X1)) > 1e-9
  
  # 2. Checagem de Dinamismo das Covariáveis
  # Se as covariáveis são dinâmicas, X(t) deve ser diferente de X(0).
  static_check_val = dt[, max(abs(X1 - first(X1))), by = id][, V1]
  dynamic_check = any(static_check_val > 1e-9)
  
  # 3. Checagem de Confundimento (X -> A)
  # Se X confunde A, o valor médio de X no momento do tratamento deve
  # ser diferente entre os grupos tratado e controle.
  treatment_time = min(dt[A == 1, t])
  conf_check = dt[t == treatment_time, .(mean_X1 = mean(X1)), by = A][, diff(mean_X1)]
  confounding_check = abs(conf_check) > 1e-9 # Se a diferença é não-nula
  
  data.table(
    Feedback_Ativo = feedback_check,
    Covariaveis_Dinamicas = dynamic_check,
    Confundimento_Ativo = confounding_check
  )
}