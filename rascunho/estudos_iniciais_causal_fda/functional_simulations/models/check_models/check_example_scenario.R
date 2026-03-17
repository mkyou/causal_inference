rm(list = ls()); graphics.off()

library(purrr); library(data.table); library(ggplot2); 
library(mgcv); library(zoo)

# --- 1. CARREGAMENTO DE FUNÇÕES AUXILIARES ---

dgp_path = "functional_simulations/data_generation_process/"
utils_path = "functional_simulations/utils/"
models_path = "functional_simulations/models/"
data_path = "functional_simulations/data_examples/"

# Funções do DGP e Orquestradores
source(paste0(dgp_path, "generate_sequential_dgp.R"))
source(paste0(dgp_path, "generate_functional_outcomes.R"))
source(paste0(dgp_path, "process_functional_data.R"))
source(paste0(dgp_path, "generate_simulation_scenarios.R"))
source(paste0(data_path, "example_scenario.R"))

# Funções de Modelos e Wrappers
source(paste0(models_path, "ecker_luna.R"))
source(paste0(models_path, "marginal_structural_model.R"))

# Funções Utilitárias e de Avaliação
source(paste0(utils_path, "revert_to_long_format.R"))
source(paste0(utils_path, "prepare_analysis_data.R"))
source(paste0(utils_path, "assess_estimation.R"))
source(paste0(utils_path, "run_model_evaluation.R"))

# --- 2. GERAÇÃO DOS DADOS DE TESTE ---
all_data = generate_ecker_luna_scenarios(n_subjects = 2000)

# --- 3. DEFINIÇÃO DOS MODELOS A SEREM AVALIADOS ---
models_to_run = list(
  Ecker_Luna = list(
    model_function = run_ecker_luna,
    assessment_function = assess_ecker_luna_process
  ),
  MSM = list(
    model_function = run_fmsm,
    assessment_function = assess_full_msm_process
  )
)

# --- 4. EXECUÇÃO UNIFICADA ---
all_results = purrr::imap(
  models_to_run,
  function(model, model_name) { # 'model_name' agora é um argumento
    cat(paste("\n\n===== AVALIANDO MODELO:", model_name, "=====\n"))
    
    # A linha que causava o erro foi removida.
    
    run_model_evaluation(
      all_scenarios_data = all_data,
      model_function = model$model_function,
      assessment_function = model$assessment_function,
      model_name = model_name # Usamos o nome diretamente
    )
  }
)

# --- 5. APRESENTAÇÃO COMPARATIVA DOS RESULTADOS ---
summary_table = purrr::map_dfr(
  all_results,
  function(model_results) {
    purrr::map_dfr(
      model_results, ~ .x$causal_effect_assessment$summary, .id = "Cenario"
    )
  },
  .id = "Modelo"
)

cat("\n\n--- TABELA COMPARATIVA DE PERFORMANCE FINAL ---\n")
print(summary_table)

all_results$MSM$TimeVaryingConfounding$propensity_model_assessment$`5`$metrics
all_results$MSM$TimeVaryingConfounding$propensity_model_assessment$`5`$plot

all_results$MSM$TimeVaryingConfounding$propensity_model_assessment$`15`$metrics
all_results$MSM$TimeVaryingConfounding$propensity_model_assessment$`15`$plot

all_results$MSM$NoConfounding$propensity_model_assessment$`5`$metrics
all_results$MSM$NoConfounding$propensity_model_assessment$`5`$plot

all_results$MSM$NoConfounding$propensity_model_assessment$`15`$metrics
all_results$MSM$NoConfounding$propensity_model_assessment$`15`$plot

