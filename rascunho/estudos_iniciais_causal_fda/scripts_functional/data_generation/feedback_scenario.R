# ===================================================================
# SCRIPT 2 (REFATORADO): TESTE DE ROBUSTEZ AO FEEDBACK (A -> X)
# Esta versão mantém uma estrutura de cenários completa e incorpora
# as formas funcionais "implacáveis" para maior clareza visual.
# ===================================================================

# --- 1. SETUP E CARREGAMENTO DE FUNÇÕES ---
source("scripts_functional/generate_y_historical.R")
source("scripts_functional/generate_y_concurrent.R")
source("scripts_functional/generate_confounder_pointwise.R")
source("scripts_functional/generate_simulation_scenarios.R")
source("scripts_functional/create_baseline_heterogeneity.R")

source("scripts_functional/utils/plot_scenario_effects.R")
source("scripts_functional/utils/plot_scenario_trajectories.R")
source("scripts_functional/utils/check_causal_structure.R")
source("scripts_functional/utils/run_validation_suite.R")

library(mgcv)
library(ggplot2)
library(data.table)

# ===================================================================
# --- 2. DEFINIÇÃO DOS COMPONENTES DA SIMULAÇÃO ---
# ===================================================================

# -- Componentes Padrão (reutilizados para consistência) --
MEAN_FUNCS_STATIC = list(
  X1 = function(t) 1.0, 
  X2 = function(t) 0.5, 
  X3 = function(t) 1.5
)

MEAN_FUNCS_DYNAMIC = list(
  X1 = function(t) 3.2*t, 
  X2 = function(t) 0.3*sin(t), 
  X3 = function(t) 1.5*log(1 + 0.8*t)
)

TIME_GRID = seq(0, 10, by = 0.1)

GLOBAL_PARAMS = list(
  time_grid = TIME_GRID, 
  noise_sd = 0.2, 
  p_treat = 0.5,
  base_y_fun = function(t) 2 + 0.5 * t
)

# -- Efeitos de A em Y (Versões "Killer") --
EFFECTS_A_ON_Y_S1 = list(
  `1` = function(t, s) ifelse(t >= s, 0.5 * (t - s) * exp(-.65 * (t - s)), 0)
)
EFFECTS_A_ON_Y_S3 = list(
  function(t, s) ifelse(t >= s, 0.65 * (t - s) * exp(-1.2 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.55 * (t - s) * exp(-1.5 * (t - s)), 0),
  function(t, s) ifelse(t >= s, -1.0 * (1 - exp(-1.0 * (t - s))), 0)
)

# -- Efeitos de X em Y (Versão "Killer") --
EFFECTS_X_ON_Y = list(
  X1 = function(t) 3 * log(1 + (t / 5) * pi),
  X2 = function(t) 1.5 + 0.3 * sqrt(t),
  X3 = function(t) 2.5 - cos(0.8 * t)
)

# -- Funções de Feedback (A -> X), o foco deste script --
FEEDBACK_NONE_S1 = list(
  X1 = list( `1` = function(t,s) 0), 
  X2 = list( `1` = function(t,s) 0), 
  X3 = list( `1` = function(t,s) 0)
)

FEEDBACK_NONE_S3 = list(
  X1 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X2 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X3 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0)
)
FEEDBACK_MEDIUM_S1 = list(
  X1 = list(`1` = function(t, s) ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)),
  X2 = FEEDBACK_NONE_S1$X2,
  X3 = FEEDBACK_NONE_S1$X3
)
FEEDBACK_HIGH_S1 = list(
  X1 = list(`1` = function(t, s) ifelse(t > s, 3.0 * exp(-0.5 * (t - s)), 0)),
  X2 = list(`1` = function(t, s) ifelse(t > s, -2.0, 0)),
  X3 = FEEDBACK_NONE_S1$X3
)
FEEDBACK_S3_A1_MED = list(
  X1 = list(
    function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) 0, 
    function(t,s) 0
  ),
  X2 = FEEDBACK_NONE_S3$X2,
  X3 = FEEDBACK_NONE_S3$X3
)
FEEDBACK_S3_A2_A3_MED = list(
  X1 = list(
    function(t, s) 0, 
    function(t, s) ifelse(t > s, 1.5*exp(-0.5*(t-s)), 0), 
    function(t, s) ifelse(t > s, 2.5*exp(-0.5*(t-s)),0)
  ),
  X2 = FEEDBACK_NONE_S3$X2,
  X3 = list(
    function(t, s) 0, 
    function(t, s) ifelse(t > s, 2.5*exp(-0.5*(t-s)), 0), 
    function(t, s) ifelse(t > s, 1.5*exp(-0.5*(t-s)), 0)
  )
)

# ===================================================================
# --- 3. DEFINIÇÃO DOS CENÁRIOS DE FEEDBACK ---
# ===================================================================
scenario_configs_fb = list(
  
  # --- S=1, Covariáveis Estáticas ---
  `C2.1_S1_NoFB_Static` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C2.2_S1_MedFB_Static` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_MEDIUM_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C2.3_S1_HighFB_Static` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_HIGH_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- S=3, Covariáveis Estáticas ---
  `C2.4_S3_NoFB_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.5_S3_FB_A1_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_S3_A1_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.6_S3_FB_A2_A3_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_S3_A2_A3_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # --- S=1, Covariáveis Dinâmicas ---
  `C2.7_S1_NoFB_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C2.8_S1_MedFB_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_MEDIUM_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- S=3, Covariáveis Dinâmicas ---
  `C2.9_S3_NoFB_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.10_S3_FB_A1_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_S3_A1_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

# ===================================================================
# --- 4. EXECUÇÃO E ANÁLISE ---
# ===================================================================

feedback_scenarios = generate_simulation_scenarios(
  scenario_configs_fb,
  global_params = GLOBAL_PARAMS
)

validation_output_fb = run_validation_suite(
  feedback_scenarios,
  scenario_configs_fb
)