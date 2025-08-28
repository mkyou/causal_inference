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

# ===================================================================
# SCRIPT 1: TESTE DE ROBUSTEZ AO FEEDBACK (A -> X)
# ===================================================================

# --- 1. PARÂMETROS GLOBAIS E FUNÇÕES BASE ---
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

# Efeitos (pequenos) de A em Y
EFFECTS_A_ON_Y_S1 = list(
  `1` = function(t, s) ifelse(t >= s, 0.5 * exp(-0.8 * (t - s)), 0)
)

EFFECTS_A_ON_Y_S3 = list(
  function(t, s) ifelse(t >= s, 0.5 * exp(-0.8 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.4 * exp(-0.6 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.6 * exp(-1.0 * (t - s)), 0)
)

# Funções de propensão com força variável nulos (sem confundimento)
PROP_NONE_S1 = NULL
PROP_NONE_S3 = list(NULL, NULL, NULL)

# Efeito de feedback com força variável
FEEDBACK_NONE_S1 = list(
  X1 = list(`1` = function(t,s) 0), 
  X2 = list(`1` = function(t,s) 0), 
  X3 = list(`1` = function(t,s) 0)
)

FEEDBACK_NONE_S3 = list(
  X1 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X2 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X3 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0)
)

FEEDBACK_MEDIUM_S1 = list(
  X1 = list(`1` = function(t,s) ifelse(t > s, 1.5 * exp( -0.5 * (t - s)), 0)), 
  X2 = list(`1` = function(t,s) 0), 
  X3 = list(`1` = function(t,s) 0)
)
FEEDBACK_HIGH_S1 = list(
  X1 = list(`1` = function(t,s) ifelse(t > s, 3.0 * exp( -0.5 * (t - s)), 0)), 
  X2 = list(`1` = function(t,s) ifelse(t > s, -2.0, 0)), 
  X3 = list(`1` =function(t,s) 0)
)

FEEDBACK_A1_MED_S3 = list(
  X1 = list(
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0),
    function(t,s) 0, 
    function(t,s) 0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) 0
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) 0
  )
)

FEEDBACK_A2_MED_S3 = list(
  X1 = list(
    function(t,s) 0,
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) 0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) 0
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) 0
  )
)

FEEDBACK_A3_MED_S3 = list(
  X1 = list(
    function(t,s) 0,
    function(t,s) 0, 
    function(t,s) 0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0)
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0)
  )
)

FEEDBACK_A2_A3_MED_S3 = list(
  X1 = list(
    function(t,s) 0,
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) ifelse(t>s, 2.5*exp(-0.5*(t-s)), 0)
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) 0
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s, 2.5*exp(-0.5*(t-s)), 0), 
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0)
  )
)



# --- 2. DEFINIÇÃO DOS CENÁRIOS DE CONFUNDIMENTO ---
scenario_configs_fb = list(
  # --- Testes com 1 tratamento (S=1) ---
  `C2.1_S1_NoFB` = list(
    n = 1000, treatment_times = 1, a_seq = 1,
    another_a_seq = 0,
    propensity_function_list = PROP_NONE_S1, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C2.2_S1_MedFB` = list(
    n = 1000, treatment_times = 1, a_seq = 1,
    another_a_seq = 0,
    propensity_function_list = PROP_NONE_S1, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_MEDIUM_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C2.3_S1_HighFB` = list(
    n = 1000, treatment_times = 1, a_seq = 1,
    another_a_seq = 0,
    propensity_function_list = PROP_NONE_S1, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_HIGH_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- Testes com 3 tratamentos (S=3) ---
  `C2.4_S3_NoFB` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), propensity_function_list = PROP_NONE_S3, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.5_S3_FB_A1` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), propensity_function_list = PROP_NONE_S3, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_A1_MED_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.6_S3_FB_A2` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), propensity_function_list = PROP_NONE_S3, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_A2_MED_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.7_S3_FB_A2` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), propensity_function_list = PROP_NONE_S3, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_A3_MED_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.8_S3_FB_A2_A3` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), propensity_function_list = PROP_NONE_S3, 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_A2_A3_MED_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

fb_scenarios = generate_simulation_scenarios(
  scenario_configs_fb,
  global_params = GLOBAL_PARAMS
)

# --- 3. EXECUÇÃO E ANÁLISE ---
validation_output = run_validation_suite(
  fb_scenarios, 
  scenario_configs_fb
)

print(validation_output$summary_table)
validation_output$plots


# --- ESTIMAÇÃO PARA MODELO DE ECKER E LUNA ---
all_results_fb_ecker_luna = lapply(
  names(fb_scenarios),
  function(name) {
    cat(paste("\n>>> Analisando o cenário:", name, "<<<\n"))
    run_sequential_eckerluna_suite(
      scenario_data = fb_scenarios[[name]],
      scenario_config = scenario_configs_fb[[name]],
      regime_a = scenario_configs_fb[[name]]$a_seq,
      regime_b = scenario_configs_fb[[name]]$another_a_seq
    )
  }
)
names(all_results_fb_ecker_luna) = names(fb_scenarios)

# --- 4. Verificando os Resultados ---
print(all_results_fb_ecker_luna$C2.1_S1_NoFB$metrics)
print(all_results_fb_ecker_luna$C2.1_S1_NoFB$plot)

print(all_results_fb_ecker_luna$C2.2_S1_MedFB$metrics)
print(all_results_fb_ecker_luna$C2.2_S1_MedFB$plot)

print(all_results_fb_ecker_luna$C2.3_S1_HighFB$metrics)
print(all_results_fb_ecker_luna$C2.3_S1_HighFB$plot)

print(all_results_fb_ecker_luna$C2.4_S3_NoFB$metrics)
print(all_results_fb_ecker_luna$C2.4_S3_NoFB$plot)

print(all_results_fb_ecker_luna$C2.5_S3_FB_A1$metrics)
print(all_results_fb_ecker_luna$C2.5_S3_FB_A1$plot)

print(all_results_fb_ecker_luna$C2.6_S3_FB_A2$metrics)
print(all_results_fb_ecker_luna$C2.6_S3_FB_A2$plot)

print(all_results_fb_ecker_luna$C2.7_S3_FB_A2$metrics)
print(all_results_fb_ecker_luna$C2.7_S3_FB_A2$plot)

print(all_results_fb_ecker_luna$C2.8_S3_FB_A2$metrics)
print(all_results_fb_ecker_luna$C2.8_S3_FB_A2$plot)


rm(
  all_results_fb_ecker_luna, current_config, current_data, 
  EFFECTS_A_ON_Y_S1, EFFECTS_A_ON_Y_S3, EFFECTS_X_ON_Y,
  fb_scenarios, FEEDBACK_A1_MED_S3, FEEDBACK_A2_A3_MED_S3,
  FEEDBACK_A2_MED_S3, FEEDBACK_A3_MED_S3, FEEDBACK_HIGH_S1,
  FEEDBACK_MEDIUM_S1, FEEDBACK_NONE_S1, FEEDBACK_NONE_S3,
  GLOBAL_PARAMS, MEAN_FUNCS, PROP_NONE_S3, PROP_NONE_S1,
  scenario_configs_fb, validation_output, scenario_name,
  TIME_GRID
)
