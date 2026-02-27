source("scripts_functional/data_generation/confouding_scenario.R")
source("scripts_functional/models/ecker_luna_model.R")

# --- ESTIMAÇÃO PARA MODELO DE ECKER E LUNA ---
all_results_conf_ecker_luna = lapply(
  names(confouding_scenarios),
  function(name) {
    cat(paste("\n>>> Analisando o cenário:", name, "<<<\n"))
    run_sequential_eckerluna_suite(
      scenario_data = confouding_scenarios[[name]],
      scenario_config = scenario_configs_conf[[name]],
      regime_a = scenario_configs_conf[[name]]$a_seq,
      regime_b = scenario_configs_conf[[name]]$another_a_seq
    )
  }
)
names(all_results_conf_ecker_luna) = names(confouding_scenarios)


# --- Verificação com Estilo Tidyverse (purrr::iwalk) ---

library(purrr)

iwalk(all_results_conf_ecker_luna, ~{
  
  # .y é o nome do cenário (ex: "C1.1_S1_NoConf")
  # .x é o conteúdo do cenário (a lista com $metrics e $plot)
  
  cat("\n=====================================================\n")
  cat(">>> Análise do Cenário:", .y, "\n")
  cat("=====================================================\n\n")
  
  print(.x$metrics)
  print(.x$plot)
  
  readline(prompt = "Pressione [enter] para ver o próximo cenário...")
})



rm(list = ls())
graphics.off()
