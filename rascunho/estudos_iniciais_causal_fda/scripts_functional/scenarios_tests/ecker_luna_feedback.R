source("scripts_functional/data_generation/feedback_scenario.R")
source("scripts_functional/models/ecker_luna_model.R")

validation_output_fb$summary_table

# --- ESTIMAÇÃO PARA MODELO DE ECKER E LUNA ---
all_results_fb_ecker_luna = lapply(
  names(feedback_scenarios),
  function(name) {
    cat(paste("\n>>> Analisando o cenário:", name, "<<<\n"))
    run_sequential_eckerluna_suite(
      scenario_data = feedback_scenarios[[name]],
      scenario_config = scenario_configs_fb[[name]],
      regime_a = scenario_configs_fb[[name]]$a_seq,
      regime_b = scenario_configs_fb[[name]]$another_a_seq
    )
  }
)
names(all_results_fb_ecker_luna) = names(feedback_scenarios)


# --- Verificação com Estilo Tidyverse (purrr::iwalk) ---

library(purrr)

iwalk(all_results_fb_ecker_luna, ~{
  
  cat("\n=====================================================\n")
  cat(">>> Análise do Cenário:", .y, "\n")
  cat("=====================================================\n\n")
  
  print(.x$metrics)
  print(.x$plot)
  
  readline(prompt = "Pressione [enter] para ver o próximo cenário...")
})


rm(list = ls())
graphics.off()
