#' @title (VERSÃO FINAL E MAIS ROBUSTA) Valida a Estrutura Causal de um Cenário
#' @description Verifica feedback e confundimento, aplicando uma correção para
#' comparações múltiplas (Bonferroni) para evitar falsos positivos.
#' @param scenario_data Um cenário da lista all_scenarios_data.
#' @param scenario_config A configuração correspondente.
#' @return Uma data.table com um resumo descritivo da validação.

check_causal_structure = function(scenario_data, scenario_config) {
  require(data.table)
  dt = scenario_data$outcome$long
  
  # --- 1. Checagem de Feedback (sem alterações) ---
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  feedback_active_on = c()
  for (cov in cov_names) {
    a_cov = paste0("a_", cov)
    another_a_cov = paste0("another_a_", cov)
    if (max(abs(dt[[a_cov]] - dt[[another_a_cov]])) > 1e-9) {
      feedback_active_on = c(feedback_active_on, cov)
    }
  }
  feedback_summary = if (length(feedback_active_on) > 0) paste(feedback_active_on, collapse = ", ") else "Nenhum"
  
  # --- 2. Checagem de Confundimento com Correção de Bonferroni ---
  treatment_times = scenario_config$treatment_times
  confounding_results = list()
  
  if (!is.null(treatment_times)) {
    all_time_points = sort(unique(dt$t))
    num_tests = length(cov_names) * length(treatment_times)
    alpha_corrected = 0.05 / num_tests # <-- CORREÇÃO DE BONFERRONI
    
    for (tt in treatment_times) {
      # ... (lógica pré-tratamento, sem mudanças)
      pre_treatment_time_index = which(all_time_points == tt) - 1
      if (pre_treatment_time_index < 1) next
      pre_treatment_time = all_time_points[pre_treatment_time_index]
      
      pre_treatment_covs = dt[t == pre_treatment_time, .SD, .SDcols = c("id", cov_names)]
      treatment_assignment = dt[t == tt, .(id, A)]
      test_data = merge(pre_treatment_covs, treatment_assignment, by = "id")
      
      conf_messages = c()
      for (cov in cov_names) {
        test_result = tryCatch({
          if (test_data[, uniqueN(A)] == 2) t.test(test_data[[cov]] ~ test_data$A) else NULL
        }, error = function(e) NULL)
        
        # Usa o alfa corrigido para determinar a significância
        if (!is.null(test_result) && test_result$p.value < alpha_corrected) { # <-- USANDO ALFA CORRIGIDO
          conf_messages = c(conf_messages, sprintf("%s (p=%.2g)", cov, test_result$p.value))
        }
      }
      
      if (length(conf_messages) > 0) {
        confounding_results[[as.character(tt)]] = paste(conf_messages, collapse = "; ")
      }
    }
  }
  
  confounding_summary = if (length(confounding_results) > 0) {
    paste(sapply(names(confounding_results), function(tt) {
      sprintf("@ t=%s: [%s]", tt, confounding_results[[tt]])
    }), collapse = " | ")
  } else {
    "Nenhum detectado (pós-correção)"
  }
  
  return(data.table(
    Feedback_Ativo = feedback_summary,
    Confundimento_Detectado = confounding_summary
  ))
}