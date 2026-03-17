# --- 1. CARREGAR TODAS AS FUNÇÕES NECESSÁRIAS ---
# Funções do PGD
source("scripts_functional/generate_y_concurrent.R")
source("scripts_functional/generate_confounder_pointwise.R")
source("scripts_functional/generate_simulation_scenarios.R")
source("scripts_functional/create_baseline_heterogeneity.R") 

# Funções de modelagem e de geração de dados
source("scripts_functional/models/g_formula_covariates_model.R")
source("scripts_functional/data_generation/confouding_feedback_scenario.R")

teste = run_gformula_validation_suite(
  scenario_data = combined_scenarios$`C3.2_A1->X->A2_Dyna`,
  scenario_config = scenario_configs_comb$`C3.2_A1->X->A2_Dyna`
)

teste$metrics
teste$plots$X1
teste$plots$X3

# --- 2. Ajustar modelos --- #
all_results_comb_g_formula_covariates = lapply(
  names(combined_scenarios),
  function(name) {
    run_gformula_validation_suite(
      scenario_data = combined_scenarios[[name]],
      scenario_config = scenario_configs_comb[[name]]
    )
  }
)
names(all_results_comb_g_formula_covariates) = names(combined_scenarios)

purrr::iwalk(all_results_comb_g_formula_covariates, ~{
  
  cat("\n=====================================================\n")
  cat(">>> Análise do Cenário:", .y, "\n")
  cat("=====================================================\n\n")
  
  print(.x$metrics)
  print(.x$plots$X1)
  print(.x$plots$X2)
  print(.x$plots$X3)
  
  readline(prompt = "Pressione [enter] para ver o próximo cenário...")
})

rm(list = ls())
graphics.off()
