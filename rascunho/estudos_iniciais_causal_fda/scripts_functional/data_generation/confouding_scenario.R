# ===================================================================
# SCRIPT 1 (REFATORADO): TESTE DE ROBUSTEZ AO CONFUNDIMENTO (X -> A)
# Versão com formatação melhorada para maior legibilidade.
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

# -- Fábrica de função de propensão normalizada --
create_normalized_propensity_fun = function(intercept, coeffs_x) {
  normalized_prop_fun = function(X, ti) {
    num_covariates = length(coeffs_x)
    X_prev_list = lapply(1:num_covariates, function(k) X[[k]][, ti - 1])
    X_norm_list = lapply(X_prev_list, function(x_vec) {
      mu = mean(x_vec); sigma = sd(x_vec)
      if (sigma < .Machine$double.eps) {
        return(rep(0, length(x_vec)))
      } else {
        return((x_vec - mu) / sigma)
      }
    })
    linear_predictor = intercept + Reduce(
      `+`,
      lapply(1:num_covariates, function(k) coeffs_x[k] * X_norm_list[[k]])
    )
    1 / (1 + exp(-linear_predictor))
  }
  return(normalized_prop_fun)
}

# -- Componentes Padrão --
MEAN_FUNCS_STATIC = list(
  X1 = function(t) 1.0, 
  X2 = function(t) 0.5, 
  X3 = function(t) 1.5
)

MEAN_FUNCS_DYNAMIC = list(
  X1 = function(t) 3.2 * t, 
  X2 = function(t) 0.3 * sin(t), 
  X3 = function(t) 1.5 * log(1 + 0.8*t)
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
  function(t, s) ifelse(t >= s, -0.2 * (1 - exp(-1.0 * (t - s))), 0)
)

# -- Efeitos de X em Y (Versão "Killer") --
EFFECTS_X_ON_Y = list(
  X1 = function(t) 1.5 + 0.5 * t,
  X2 = function(t) 2.0,
  X3 = function(t) -1.0
)

# -- Funções de Propensão --
PROP_LOW = create_normalized_propensity_fun(
  intercept = -0.5, coeffs_x = c(0.8, 0.6, -0.5)
)

PROP_MEDIUM = create_normalized_propensity_fun(
  intercept = -1.0, coeffs_x = c(1.5, 1.2, -1.0)
)

PROP_HIGH = create_normalized_propensity_fun(
  intercept = -1.5, 
  coeffs_x = c(2.5, 2.0, -2.0)
)

# -- Efeito de Feedback Nulo --
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

# ===================================================================
# --- 3. DEFINIÇÃO DOS 18 CENÁRIOS DE CONFUNDIMENTO ---
# ===================================================================
scenario_configs_conf = list(
  
  # --- S=1, Covariáveis Estáticas ---
  `C1.1_S1_NoConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.2_S1_LowConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_LOW),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.3_S1_MedConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.4_S1_HighConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_HIGH),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- S=3, Covariáveis Estáticas ---
  `C1.5_S3_NoConf` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.6_S3_Conf_A1` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(PROP_MEDIUM, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.7_S3_Conf_A2` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.8_S3_Conf_A3` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.9_S3_Conf_A2_A3` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # --- S=1, Covariáveis Dinâmicas ---
  `C1.10_S1_NoConf_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.11_S1_LowConf_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_LOW),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.12_S1_MedConf_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.13_S1_HighConf_Dyna` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_HIGH),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- S=3, Covariáveis Dinâmicas ---
  `C1.14_S3_NoConf_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.15_S3_Conf_A1_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(PROP_MEDIUM, NULL, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.16_S3_Conf_A2_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.17_S3_Conf_A3_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, NULL, PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.18_S3_Conf_A2_A3_Dyna` = list(
    n = 5000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, PROP_MEDIUM),
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    effects_covariates_from_A = FEEDBACK_NONE_S3,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

# ===================================================================
# --- 4. EXECUÇÃO E ANÁLISE ---
# ===================================================================

confouding_scenarios = generate_simulation_scenarios(
  scenario_configs_conf,
  global_params = GLOBAL_PARAMS
)

validation_output = run_validation_suite(
  confouding_scenarios,
  scenario_configs_conf
)