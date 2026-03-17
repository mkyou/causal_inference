# --- 1. CARREGAR TODAS AS FUNÇÕES NECESSÁRIAS ---
# Funções do PGD
source("scripts_functional/generate_y_concurrent.R")
source("scripts_functional/generate_confounder_pointwise.R")
source("scripts_functional/generate_simulation_scenarios.R")
source("scripts_functional/create_baseline_heterogeneity.R") 

# Funções de modelagem e de geração de dados
source("scripts_functional/models/g_formula_covariates_model.R")
source("scripts_functional/data_generation/feedback_scenario.R")


# --- 2. Ajustar modelos --- #
all_results_fb_g_formula_covariates = lapply(
  names(feedback_scenarios),
  function(name) {
      run_gformula_validation_suite(
        scenario_data = feedback_scenarios[[name]],
        scenario_config = scenario_configs_fb[[name]]
      )
  }
)
names(all_results_fb_g_formula_covariates) = names(feedback_scenarios)

purrr::iwalk(all_results_fb_g_formula_covariates, ~{
  
  # .y é o nome do cenário (ex: "C1.1_S1_NoConf")
  # .x é o conteúdo do cenário (a lista com $metrics e $plot)
  
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
