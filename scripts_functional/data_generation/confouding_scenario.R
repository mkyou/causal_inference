# ===================================================================
# SCRIPT DE VERIFICAÇÃO FINAL
# Execute este bloco em uma sessão R limpa para garantir a reprodutibilidade.
# ===================================================================
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
source("scripts_functional/models/g_formula_model.R")

library(data.table)
library(mgcv)

# ===================================================================
# SCRIPT 1: TESTE DE ROBUSTEZ AO CONFUNDIMENTO (X -> A)
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
  X1 = function(t) 3.5, 
  X2 = function(t) 2.0, 
  X3 = function(t) -1.8
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

# Funções de propensão com força variável
PROP_LOW = function(X, ti) 1/
  (1+exp(-(-0.5 + 0.8*X[[1]][,ti-1] + 0.6*X[[2]][,ti-1] - 0.5*X[[3]][,ti-1])))
PROP_MEDIUM = function(X, ti) 1/
  (1+exp(-(-1.0 + 1.5*X[[1]][,ti-1] + 1.2*X[[2]][,ti-1] - 1.0*X[[3]][,ti-1])))
PROP_HIGH = function(X, ti) 1/
  (1+exp(-(-2.0 + 4.5*X[[1]][,ti-1] + 4.0*X[[2]][,ti-1] - 3.0*X[[3]][,ti-1])))

# Efeito de feedback NULO (usado em todos os cenários deste script)
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

# --- 2. DEFINIÇÃO DOS CENÁRIOS DE CONFUNDIMENTO ---
scenario_configs_conf = list(
  # --- Testes com 1 tratamento (S=1) ---
  `C1.1_S1_NoConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = NULL,
    mean_baseline_covariates = MEAN_FUNCS,
    effects_covariates_from_A = FEEDBACK_NONE_S1,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.2_S1_LowConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_LOW), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.3_S1_MedConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_MEDIUM), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  `C1.4_S1_HighConf` = list(
    n = 1000, treatment_times = 1, a_seq = 1, another_a_seq = 0,
    propensity_function_list = list(PROP_HIGH), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S1, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S1
  ),
  
  # --- Testes com 3 tratamentos (S=3) ---
  `C1.5_S3_NoConf` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), 
    propensity_function_list = list(NULL, NULL, NULL), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.6_S3_Conf_A1` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), 
    propensity_function_list = list(PROP_MEDIUM, NULL, NULL), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.7_S3_Conf_A2` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), 
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.8_S3_Conf_A3` = list(
    n = 5000,  treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), 
    propensity_function_list = list(NULL, NULL, PROP_MEDIUM), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  `C1.9_S3_Conf_A2_A3` = list(
    n = 5000, treatment_times = c(1,2,3), a_seq = c(1,0,1), 
    another_a_seq = c(0,1,0), 
    propensity_function_list = list(NULL, PROP_MEDIUM, PROP_MEDIUM), 
    mean_baseline_covariates = MEAN_FUNCS, 
    effects_covariates_from_A = FEEDBACK_NONE_S3, 
    concurrent_effects_X = EFFECTS_X_ON_Y, 
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

confouding_scenarios = generate_simulation_scenarios(
  scenario_configs_conf,
  global_params = GLOBAL_PARAMS
)

# --- 3. EXECUÇÃO E ANÁLISE ---
validation_output = run_validation_suite(
  confouding_scenarios, 
  scenario_configs_conf
)

validation_output$summary_table
validation_output$plots

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

# --- ESTIMAÇÃO PARA G-FORMULA
all_results_conf_g_formula = lapply(
  names(confouding_scenarios),
  function(name) {
    cat(paste("\n>>> Analisando o cenário com G-formula:", name, "<<<\n"))
    run_gformula_suite(
      scenario_data = confouding_scenarios[[name]],
      scenario_config = scenario_configs_conf[[name]]
    )
  }
)
names(all_results_conf_g_formula) = names(confouding_scenarios)

# ---------  Verificando os Resultados ECKER E LUNA ---------
print(all_results_conf_ecker_luna$C1.1_S1_NoConf$metrics)
print(all_results_conf_ecker_luna$C1.1_S1_NoConf$plot)

print(all_results_conf_ecker_luna$C1.2_S1_LowConf$metrics)
print(all_results_conf_ecker_luna$C1.2_S1_LowConf$plot)

print(all_results_conf_ecker_luna$C1.3_S1_MedConf$metrics)
print(all_results_conf_ecker_luna$C1.3_S1_MedConf$plot)

print(all_results_conf_ecker_luna$C1.4_S1_HighConf$metrics)
print(all_results_conf_ecker_luna$C1.4_S1_HighConf$plot)

print(all_results_conf_ecker_luna$C1.5_S3_NoConf$metrics)
print(all_results_conf_ecker_luna$C1.5_S3_NoConf$plot)

print(all_results_conf_ecker_luna$C1.6_S3_Conf_A1$metrics)
print(all_results_conf_ecker_luna$C1.6_S3_Conf_A1$plot)

print(all_results_conf_ecker_luna$C1.7_S3_Conf_A2$metrics)
print(all_results_conf_ecker_luna$C1.7_S3_Conf_A2$plot)

print(all_results_conf_ecker_luna$C1.8_S3_Conf_A3$metrics)
print(all_results_conf_ecker_luna$C1.8_S3_Conf_A3$plot)

print(all_results_conf_ecker_luna$C1.9_S3_Conf_A2_A3$metrics)
print(all_results_conf_ecker_luna$C1.9_S3_Conf_A2_A3$plot)

# ---------  Verificando os Resultados G-FORMULA ---------
print(all_results_conf_g_formula$C1.1_S1_NoConf$metrics)
print(all_results_conf_g_formula$C1.1_S1_NoConf$plot)

print(all_results_conf_g_formula$C1.2_S1_LowConf$metrics)
print(all_results_conf_g_formula$C1.2_S1_LowConf$plot)

print(all_results_conf_g_formula$C1.3_S1_MedConf$metrics)
print(all_results_conf_g_formula$C1.3_S1_MedConf$plot)

print(all_results_conf_g_formula$C1.4_S1_HighConf$metrics)
print(all_results_conf_g_formula$C1.4_S1_HighConf$plot)

print(all_results_conf_g_formula$C1.5_S3_NoConf$metrics)
print(all_results_conf_g_formula$C1.5_S3_NoConf$plot)

print(all_results_conf_g_formula$C1.6_S3_Conf_A1$metrics)
print(all_results_conf_g_formula$C1.6_S3_Conf_A1$plot)

print(all_results_conf_g_formula$C1.7_S3_Conf_A2$metrics)
print(all_results_conf_g_formula$C1.7_S3_Conf_A2$plot)

print(all_results_conf_g_formula$C1.8_S3_Conf_A3$metrics)
print(all_results_conf_g_formula$C1.8_S3_Conf_A3$plot)

print(all_results_conf_g_formula$C1.9_S3_Conf_A2_A3$metrics)
print(all_results_conf_g_formula$C1.9_S3_Conf_A2_A3$plot)



rm(
  all_results_conf_ecker_luna, all_results_conf_g_formula,
  cfg1, cfg2, cfg3, cfg4, cfg5, concurrent_effects_on_Y, 
  confouding_scenarios, current_config, current_data, 
  EFFECTS_A_ON_Y_S1, EFFECTS_A_ON_Y_S3, EFFECTS_X_ON_Y, 
  FEEDBACK_NONE_S1, FEEDBACK_NONE_S3, global_params, 
  GLOBAL_PARAMS, MEAN_FUNCS, mean_funcs_dynamic, 
  mean_funcs_static, scenario_configs, scenario_configs_conf,
  validation_output, scenario_name, TIME_GRID,
  PROP_HIGH, PROP_LOW, PROP_MEDIUM
)
