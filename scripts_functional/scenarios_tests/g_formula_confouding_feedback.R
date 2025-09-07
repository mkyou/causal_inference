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
      covariate_model_type = "functional"
    )
  }
)
names(all_results_comb_g_formula) = names(combined_scenarios)


# --- Verificação com Estilo Tidyverse (purrr::iwalk) ---
iwalk(all_results_comb_g_formula, ~{
  
  cat("\n=====================================================\n")
  cat(">>> Análise do Cenário:", .y, "\n")
  cat("=====================================================\n\n")
  
  print(.x$metrics)
  print(.x$plot)
  
  readline(prompt = "Pressione [enter] para ver o próximo cenário...")
})


rm(list = ls())
graphics.off()
