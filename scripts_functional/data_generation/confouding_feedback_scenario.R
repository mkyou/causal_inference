# --- 1. SETUP E CARREGAMENTO DE FUNÇÕES ---
source("scripts_functional/generate_y_historical.R")
source("scripts_functional/generate_y_concurrent.R")
source("scripts_functional/generate_confounder_pointwise.R")
source("scripts_functional/utils/generate_simulation_scenarios.R")
source("scripts_functional/utils/create_heterogeneous_baselines_advanced.R")


source("scripts_functional/utils/run_validation_suite.R")
source("scripts_functional/utils/plot_scenario_effects.R")
source("scripts_functional/utils/check_causal_structure.R")
source("scripts_functional/utils/plot_scenario_trajectories.R")

source("scripts_functional/models/ecker_luna_model.R")

library(data.table)
library(mgcv)

TIME_GRID = seq(0, 10, by = 0.1)
GLOBAL_PARAMS = list(
  time_grid = TIME_GRID, noise_sd = 0.2, p_treat = 0.5,
  base_y_fun = function(t) 2 + 0.5*t
)

# Funções de trajetória média para 3 covariáveis
MEAN_FUNCS = list(
  X1 = function(t) 1.0 + 0.2*t, 
  X2 = function(t) 0.5 + 0.3*sqrt(t),
  X3 = function(t) 1.5 * cos(0.8*t)
)

# Efeitos (médios) de X em Y
EFFECTS_X_ON_Y = list(
  X1 = function(t) 1.5, 
  X2 = function(t) 1.0, 
  X3 = function(t) -0.8
)

# Efeitos (pequenos) de A em Y para S=3
EFFECTS_A_ON_Y_S3 = list(
  function(t, s) ifelse(t >= s, 0.5 * exp(-0.8 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.4 * exp(-0.6 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.6 * exp(-1.0 * (t - s)), 0)
)

# Função de propensão de força média
PROP_MEDIUM = function(X, ti) 1/
  (1+exp(-(-1.0 + 1.5*X[[1]][,ti-1] + 2.2*X[[2]][,ti-1] - 1.0*X[[3]][,ti-1])))

# Funções de Feedback (Médio) para A1, A2
FEEDBACK_A1_MED_S3 = list(
  X1 = list(
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0), 
    function(t,s) 0, 
    function(t,s) 0),
  X2 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X3 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0)
)
FEEDBACK_A2_MED_S3 = list(
  X1 = list(
    function(t,s) 0,
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0),
    function(t,s) 0
    ),
  X2 = list(
    function(t,s) 0, 
    function(t,s)ifelse(t>s,1.0*exp(-0.4*(t-s)),0), 
    function(t,s) 0
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) 0
  )
)

# Função de feedback completo para o cenário final
FEEDBACK_ALL_MED_S3 = list(
  X1 = list(
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) ifelse(t>s,1.0*exp(-0.5*(t-s)), 0), 
    function(t,s)0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s,1.0*exp(-0.4*(t-s)), 0), 
    function(t,s) ifelse(t>s,0.8*exp(-0.4*(t-s)), 0)
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s)ifelse(t>s,-1.2,0)
  ) # Efeito constante de A3 em X3
)


# --- 3. DEFINIÇÃO DOS CENÁRIOS COMBINADOS ---
scenario_configs_comb = list(
  # PGD: A1 é aleatório. A1 causa X via feedback. X então confunde A2. 
  # A3 é aleatório.
  # Violação: O confundidor de A2, X(1), é pós-tratamento em relação a A1.
  `C3.1_A1->X->A2` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL),
    effects_covariates_from_A = FEEDBACK_A1_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # PGD: A1 aleatório. A2 é confundido. A2 causa X via feedback. 
  # X então confunde A3.
  # Violação: O confundidor de A3, X(2), é afetado por A2 e confundido por X(1).
  `C3.2_A2->X->A3` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_A2_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # PGD: A1 causa X. A2 é aleatório. X (influenciado por A1) confunde A3.
  # Violação: Similar ao C3.1, mas a cadeia causal é mais longa.
  `C3.3_A1->X->A3` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0),
    propensity_function_list = list(NULL, NULL, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_A1_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS, 
    concurrent_effects_X = EFFECTS_X_ON_Y, treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # PGD: "Full G-Formula World". Tudo está conectado.
  # A1 é confundido por X(0). A1->X. X(1) confunde A2. A2->X. X(2) confunde A3.
  # Violação: Todas as hipóteses de E&L são violadas da forma mais complexa
  # possível.
  `C3.4_Full_Loop` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0),
    propensity_function_list = list(PROP_MEDIUM, PROP_MEDIUM, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_ALL_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

# --- 4. EXECUÇÃO E ANÁLISE ---
# (Recomenda-se rodar em uma sessão R limpa)

# Gera os dados
combined_scenarios_data = generate_simulation_scenarios(
  scenario_configs = scenario_configs_comb,
  global_params = GLOBAL_PARAMS
)

# Roda a validação
validation_output_comb = run_validation_suite(
  all_scenarios_data = combined_scenarios_data,
  scenario_configs = scenario_configs_comb
)

validation_output_comb$summary_table
validation_output_comb$plots

# Roda a estimação
all_results_comb = lapply(
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
names(all_results_comb) = names(combined_scenarios_data)

print(all_results_comb$`C3.1_A1->X->A2`$metrics)
print(all_results_comb$`C3.1_A1->X->A2`$plot)

print(all_results_comb$`C3.2_A2->X->A3`$metrics)
print(all_results_comb$`C3.2_A2->X->A3`$plot)

print(all_results_comb$`C3.3_A1->X->A3`$metrics)
print(all_results_comb$`C3.3_A1->X->A3`$plot)

print(all_results_comb$C3.4_Full_Loop$metrics)
print(all_results_comb$C3.4_Full_Loop$plot)

rm(
  all_results_comb, combined_scenarios_data, EFFECTS_A_ON_Y_S3,
  EFFECTS_X_ON_Y, FEEDBACK_A1_MED_S3, FEEDBACK_A2_MED_S3, FEEDBACK_ALL_MED_S3,
  GLOBAL_PARAMS, MEAN_FUNCS, scenario_configs_comb, validation_output_comb,
  TIME_GRID, PROP_MEDIUM
)
