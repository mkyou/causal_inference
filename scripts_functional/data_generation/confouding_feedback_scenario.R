# ===================================================================
# SCRIPT 3 (PADRONIZADO): TESTE COM CONFUNDIMENTO E FEEDBACK (A->X->A)
# Execute este bloco em uma sessão R limpa para garantir a reprodutibilidade.
# ===================================================================

# --- 1. SETUP E CARREGAMENTO DE FUNÇÕES ---
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
# --- 2. DEFINIÇÃO DOS COMPONENTES DA SIMULAÇÃO ---
# Componentes padronizados e reutilizados dos scripts anteriores.
# ===================================================================

# -- Fábrica de função de propensão normalizada --
create_normalized_propensity_fun = function(intercept, coeffs_x) {
  normalized_prop_fun = function(X, ti) {
    num_covariates = length(coeffs_x)
    X_prev_list = lapply(1:num_covariates, function(k) X[[k]][, ti - 1])
    X_norm_list = lapply(X_prev_list, function(x_vec) {
      mu = mean(x_vec); sigma = sd(x_vec)
      if (sigma < .Machine$double.eps) { return(rep(0, length(x_vec)))
      } else { return((x_vec - mu) / sigma) }
    })
    linear_predictor = intercept + Reduce(`+`,
                                          lapply(1:num_covariates, function(k) coeffs_x[k] * X_norm_list[[k]])
    )
    1 / (1 + exp(-linear_predictor))
  }
  return(normalized_prop_fun)
}

# -- Trajetórias médias, Parâmetros Globais, Efeitos --
MEAN_FUNCS_STATIC = list(
  X1 = function(t) 1, 
  X2 = function(t) 0.5, 
  X3 = function(t) 1.5
)

MEAN_FUNCS_DYNAMIC = list(
  X1 = function(t) 3.2*t, 
  X2 = function(t) 0.3*sin(t), 
  X3 = function(t) 1.5*log(1+0.8*t)
)

TIME_GRID = seq(0, 10, by = 0.1)
GLOBAL_PARAMS = list(
  time_grid = TIME_GRID, 
  noise_sd = 0.2, 
  p_treat = 0.5, 
  base_y_fun = function(t) 2 + 0.5*t
)

EFFECTS_X_ON_Y = list(
  X1 = function(t) 1 + 0.2*t, 
  X2 = function(t) 0.5 + 0.3*sqrt(t), 
  X3 = function(t) 1.5*cos(0.8*t)
)
EFFECTS_A_ON_Y_S3 = list(
  function(t,s) ifelse(t>=s, 0.5*exp(-0.8*(t-s)),0),
  function(t,s) ifelse(t>=s, 0.4*exp(-0.6*(t-s)),0),
  function(t,s) ifelse(t>=s, 0.6*exp(-1.0*(t-s)),0)
)

# -- Função de Propensão (Média, Normalizada) --
PROP_MEDIUM = create_normalized_propensity_fun(
  intercept = -1.0, coeffs_x = c(1.5, 1.2, -1.0)
)

# -- Funções de Feedback (Médio) --
FEEDBACK_NONE_S3 = list(
  X1 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X2 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0),
  X3 = list(function(t,s) 0, function(t,s) 0, function(t,s) 0)
)

FEEDBACK_A1_MED_S3 = list(
  X1 = list(
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)),0), 
    function(t,s) 0, 
    function(t,s) 0
  ),
  X2 = FEEDBACK_NONE_S3$X2, 
  X3 = FEEDBACK_NONE_S3$X3
)
FEEDBACK_A2_MED_S3 = list(
  X1 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) 0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s,1.0*exp(-0.4*(t-s)), 0), 
    function(t,s) 0
  ),
  X3 = FEEDBACK_NONE_S3$X3
)

FEEDBACK_ALL_MED_S3 = list(
  X1 = list(
    function(t,s) ifelse(t>s,1.5*exp(-0.5*(t-s)), 0), 
    function(t,s) ifelse(t>s,1.0*exp(-0.5*(t-s)),0), 
    function(t,s) 0
  ),
  X2 = list(
    function(t,s) 0, 
    function(t,s) ifelse(t>s,1.0*exp(-0.4*(t-s)),0), 
    function(t,s) ifelse(t>s,0.8*exp(-0.4*(t-s)),0)
  ),
  X3 = list(
    function(t,s) 0, 
    function(t,s) 0, 
    function(t,s) ifelse(t>s,-1.2,0)
  )
)

# ===================================================================
# --- 3. DEFINIÇÃO DOS CENÁRIOS COMBINADOS ---
# ===================================================================

scenario_configs_comb = list(
  # --- Cenários Combinados: Covariáveis Estáticas ---
  
  # PGD: A1->X->A2. O confundidor de A2, X(1), é pós-tratamento em relação a A1.
  `C3.1_A1->X->A2_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL),
    effects_covariates_from_A = FEEDBACK_A1_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # PGD: A2->X->A3. O confundidor de A3, X(2), é afetado por A2.
  `C3.2_A2->X->A3_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_A2_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # PGD: "Full G-Formula World". A1->X, X->A2, A2->X, X->A3.
  `C3.3_Full_Loop_Static` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(PROP_MEDIUM, PROP_MEDIUM, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_ALL_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS_STATIC,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  # --- Cenários Combinados: Covariáveis Dinâmicas ---
  
  `C3.4_A1->X->A2_Dyna` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(NULL, PROP_MEDIUM, NULL),
    effects_covariates_from_A = FEEDBACK_A1_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  ),
  
  `C3.5_Full_Loop_Dyna` = list(
    n = 10000, treatment_times = c(1,2,3), 
    a_seq = c(1,1,1), another_a_seq = c(0,0,0),
    propensity_function_list = list(PROP_MEDIUM, PROP_MEDIUM, PROP_MEDIUM),
    effects_covariates_from_A = FEEDBACK_ALL_MED_S3,
    mean_baseline_covariates = MEAN_FUNCS_DYNAMIC,
    concurrent_effects_X = EFFECTS_X_ON_Y,
    treatment_effect = EFFECTS_A_ON_Y_S3
  )
)

# ===================================================================
# --- 4. EXECUÇÃO E ANÁLISE ---
# ===================================================================

combined_scenarios_data = generate_simulation_scenarios(
  scenario_configs = scenario_configs_comb,
  global_params = GLOBAL_PARAMS
)

validation_output_comb = run_validation_suite(
  all_scenarios_data = combined_scenarios_data,
  scenario_configs = scenario_configs_comb
)

# Para ver os resultados:
# validation_output_comb$summary_table
# validation_output_comb$plots