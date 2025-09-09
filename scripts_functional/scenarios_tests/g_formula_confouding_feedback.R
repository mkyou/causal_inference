# ===================================================================
# EXEMPLO DE USO: EXECUTANDO A SUÍTE COMPLETA DA G-FORMULA
# ===================================================================

# --- 1. CARREGAR TODAS AS FUNÇÕES NECESSÁRIAS ---
# Carrega a suíte principal (que por sua vez carrega os modelos de covariável)
source("scripts_functional/models/g_formula_model.R")

# Carrega um dos arquivos de cenário para o teste
source("scripts_functional/data_generation/confouding_feedback_scenario.R")

# --- 2. EXECUTAR A G-FORMULA PARA UM CENÁRIO ESPECÍFICO ---

all_results_comb_g_formula = lapply(
  names(combined_scenarios),
  function(name) {
    cat(paste("\n>>> Analisando o cenário:", name, "<<<\n"))
    run_full_gformula_suite(
      scenario_data = combined_scenarios[[name]],
      scenario_config = scenario_configs_comb[[name]],
      scenario_name = name
    )
  }
)
names(all_results_comb_g_formula) = names(combined_scenarios)

# --- Verificação Robusta e Visualmente Aprimorada ---
purrr::iwalk(all_results_comb_g_formula, ~{
  
  # Cabeçalho principal para o cenário
  cat("\n\n===========================================================\n")
  cat(">>> ANÁLISE DO CENÁRIO:", .y, "\n")
  cat("===========================================================\n")
  
  # 1. VALIDAÇÃO DAS COVARIÁVEIS (Iteração dinâmica)
  # -----------------------------------------------------------
  cat("\n\n--- 1. AVALIAÇÃO DAS COVARIÁVEIS ---\n")
  
  # Pega os nomes das covariáveis dinamicamente dos resultados
  cov_names = names(.x$covariate_validation)
  
  if (length(cov_names) > 0) {
    for (cov in cov_names) {
      cat("\n---------------------------------------------------\n")
      cat("--> Covariável:", cov, "\n")
      cat("---------------------------------------------------\n")
      
      # Regime A ('sempre tratar')
      cat("\n  - Regime 'Sempre Tratar' (a_seq):\n")
      print(.x$covariate_validation[[cov]]$a_seq$metrics)
      print(.x$covariate_validation[[cov]]$a_seq$plot)
      
      # Regime B ('nunca tratar')
      cat("\n  - Regime 'Nunca Tratar' (another_a_seq):\n")
      print(.x$covariate_validation[[cov]]$another_a_seq$metrics)
      print(.x$covariate_validation[[cov]]$another_a_seq$plot)
    }
  } else {
    cat("\nNenhuma validação de covariável encontrada.\n")
  }
  
  # 2. VALIDAÇÃO DO DESFECHO (Y)
  # -----------------------------------------------------------
  cat("\n\n--- 2. AVALIAÇÃO DO DESFECHO (Y) ---\n")
  
  # Ajuste nos dados observados
  cat("\n  - Ajuste vs. Observado:\n")
  print(.x$outcome_validation$observed_fit$metrics)
  print(.x$outcome_validation$observed_fit$plot)
  
  # Simulação vs. Verdadeiro (Regime A)
  cat("\n  - Simulado vs. Verdadeiro (Regime 'Sempre Tratar'):\n")
  print(.x$outcome_validation$a_seq$metrics)
  print(.x$outcome_validation$a_seq$plot)
  
  # Simulação vs. Verdadeiro (Regime B)
  cat("\n  - Simulado vs. Verdadeiro (Regime 'Nunca Tratar'):\n")
  # Corrigido o erro de digitação de 'outcole_validation' para 'outcome_validation'
  print(.x$outcome_validation$another_a_seq$metrics)
  print(.x$outcome_validation$another_a_seq$plot)
  
  # 3. VALIDAÇÃO FINAL DO EFEITO CAUSAL
  # -----------------------------------------------------------
  cat("\n\n--- 3. AVALIAÇÃO FINAL: EFEITO CAUSAL ---\n\n")
  
  # Corrigido o caminho para as métricas e plot causal
  print(.x$causal_validation$metrics)
  print(.x$causal_validation$plot)
  
  # Pausa para o próximo cenário
  readline(prompt = "\n\n>>> Pressione [enter] para ver o próximo cenário...")
})


rm(list = ls())
graphics.off()
