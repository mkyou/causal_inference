source("scripts_functional/data_generation/confouding_feedback_scenario.R")
source("scripts_functional/models/ecker_luna_model.R")

validation_output_comb$summary_table

# --- ESTIMAÇÃO PARA MODELO DE ECKER E LUNA ---
all_results_comb_ecker_luna = lapply(
  names(combined_scenarios_data),
  function(name) {
    cat(paste("\n>>> Analisando o cenário:", name, "<<<\n"))
    run_sequential_eckerluna_suite(
      scenario_data = combined_scenarios_data[[name]],
      scenario_config = scenario_configs_comb[[name]],
      regime_a = scenario_configs_comb[[name]]$a_seq,
      regime_b = scenario_configs_comb[[name]]$another_a_seq
    )
  }
)
names(all_results_comb_ecker_luna) = names(combined_scenarios_data)


# --- Verificação com Estilo Tidyverse (purrr::iwalk) ---

library(purrr)

iwalk(all_results_comb_ecker_luna, ~{
  
  cat("\n=====================================================\n")
  cat(">>> Análise do Cenário:", .y, "\n")
  cat("=====================================================\n\n")
  
  print(.x$metrics)
  print(.x$plot)
  
  readline(prompt = "Pressione [enter] para ver o próximo cenário...")
})


rm(list = ls())
graphics.off()
