# ===================================================================
# SCRIPT 2 (CORRIGIDO): TESTE DE ROBUSTEZ AO FEEDBACK (A -> X)
# Execute este bloco em uma sessão R limpa para garantir a reprodutibilidade.
# ===================================================================

# Carregando funções e pacotes necessários
source("scripts_functional/generate_y_historical.R")
source("scripts_functional/generate_y_concurrent.R")
source("scripts_functional/generate_confounder_pointwise.R")
source("scripts_functional/generate_simulation_scenarios.R")
source("scripts_functional/create_intercept_heterogeneity.R")

source("scripts_functional/utils/run_validation_suite.R")
source("scripts_functional/utils/plot_scenario_effects.R")
source("scripts_functional/utils/check_causal_structure.R")
source("scripts_functional/utils/plot_scenario_trajectories.R")

library(mgcv)
library(ggplot2)
library(data.table)

# ===================================================================
# --- 1. DEFINIÇÃO DOS COMPONENTES DA SIMULAÇÃO ---
# Reutilizando os mesmos componentes do script de confundimento para
# garantir uma comparação justa.
# ===================================================================

# -- Trajetórias médias das covariáveis --
MEAN_FUNCS_STATIC = list(
  X1 = function(t) 1.0, X2 = function(t) 0.5, X3 = function(t) 1.5
)
MEAN_FUNCS_DYNAMIC = list(
  X1 = function(t) 3.2*t, X2 = function(t) 0.3*sin(t), 
  X3 = function(t) 1.5*log(1 + 0.8*t)
)

# -- Parâmetros globais --
TIME_GRID = seq(0, 10, by = 0.1)
GLOBAL_PARAMS = list(
  time_grid = TIME_GRID, noise_sd = 0.2, p_treat = 0.5,
  base_y_fun = function(t) 2 + 0.5*t
)

# -- Efeitos de X em Y e de A em Y --
EFFECTS_X_ON_Y = list(
  X1 = function(t) 1.0 + 0.2*t, X2 = function(t) 0.5 + 0.3*sqrt(t),
  X3 = function(t) 1.5 * cos(0.8*t)
)
EFFECTS_A_ON_Y_S1 = list(
  `1` = function(t, s) ifelse(t >= s, 0.5 * exp(-0.8 * (t - s)), 0)
)
EFFECTS_A_ON_Y_S3 = list(
  function(t, s) ifelse(t >= s, 0.5 * exp(-0.8 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.4 * exp(-0.6 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 0.6 * exp(-1.0 * (t - s)), 0)
)

# --- Definições de Feedback (A -> X), o foco deste script ---
FEEDBACK_NONE_S1 = list(
  X1 = list(`1`=function(t,s) 0), X2 = list(`1`=function(t,s) 0),
  X3 = list(`1`=function(t,s) 0)
)
FEEDBACK_NONE_S3 = list(
  X1 = list(function(t,s) 0,function(t,s) 0,function(t,s) 0),
  X2 = list(function(t,s) 0,function(t,s) 0,function(t,s) 0),
  X3 = list(function(t,s) 0,function(t,s) 0,function(t,s) 0)
)
FEEDBACK_MEDIUM_S1 = list(
  X1 = list(`1`=function(t,s) ifelse(t>s, 1.5*exp(-0.5*(t-s)), 0)),
  X2 = FEEDBACK_NONE_S1$X2, X3 = FEEDBACK_NONE_S1$X3
)
FEEDBACK_HIGH_S1 = list(
  X1 = list(`1`=function(t,s) ifelse(t>s, 3.0*exp(-0.5*(t-s)), 0)),
  X2 = list(`1`=function(t,s) ifelse(t>s, -2.0, 0)),
  X3 = FEEDBACK_NONE_S1$X3
)
FEEDBACK_S3_A1_MED = list(
  X1 = list(function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0), 
            function(t,s) 0, function(t,s) 0),
  X2 = FEEDBACK_NONE_S3$X2, X3 = FEEDBACK_NONE_S3$X3
)
FEEDBACK_S3_A2_A3_MED = list(
  X1 = list(function(t,s) 0, 
            function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0), 
            function(t,s) ifelse(t>s,2.5*exp(-0.5*(t-s)),0)),
  X2 = FEEDBACK_NONE_S3$X2,
  X3 = list(function(t,s) 0, 
            function(t,s) ifelse(t>s,2.5*exp(-0.5*(t-s)),0), 
            function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0))
)

# ===================================================================
# --- 2. DEFINIÇÃO DOS CENÁRIOS DE FEEDBACK ---
# Estrutura explícita, análoga ao script de confundimento
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
    n = 10000, treatment_times = c(1,2,3), a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.5_S3_FB_A1_Static` = list(
    n = 10000, treatment_times = c(1,2,3), a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_S3_A1_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C2.6_S3_FB_A2_A3_Static` = list(
    n = 10000, treatment_times = c(1,2,3), a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_S3_A2_A3_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # --- S=1, Covariáveis Dinâmicas ---
  `C2.7_S1_MedFB_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_MEDIUM_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- S=3, Covariáveis Dinâmicas ---
  `C2.8_S3_FB_A1_Dyna` = list(
    n = 10000, treatment_times = c(1,2,3), a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_S3_A1_MED,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

# ===================================================================
# --- 3. EXECUÇÃO E ANÁLISE ---
# ===================================================================

fb_scenarios = generate_simulation_scenarios(
  scenario_configs_fb,
  global_params = GLOBAL_PARAMS
)

validation_output_fb = run_validation_suite(
  fb_scenarios, 
  scenario_configs_fb
)

# Para ver os resultados
# validation_output_fb$summary_table
# validation_output_fb$plots