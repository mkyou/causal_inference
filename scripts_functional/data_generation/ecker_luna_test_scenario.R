# ===================================================================
# SCRIPT DE SIMULAÇÃO
# Objetivo: Testar os limites do modelo de Ecker & Luna (2024) com
# cenários inequívocos e estrutura de código modular.
# ===================================================================

# --- 1. Carregando o Laboratório de Simulação ---
# source("utils/generate_simulation_scenarios.R") # O orquestrador
# source("utils/create_heterogeneous_baselines_advanced.R")
# source("generate_confundador_pointwise.R")
# source("generate_y_concurrent.R")
# source("generate_y_historical.R")

# ===================================================================
# Bloco 2: Definição dos Cenários de Simulação (Revisado)
# ===================================================================

library(data.table)

# --- Parâmetros Globais e Funções Base "Inequívocas" ---
TIME_GRID = seq(0, 5, by = 0.1)
N_SAMPLE = 500
NOISE = 0.2
TREATMENT_TIME = 1 # Tratamento único em t=1

# Efeito do tratamento em Y(t): Positivo e decaindo
common_treatment_effect = list(
  `1` = function(t, s) ifelse(t >= s, 2.0 * exp(-0.8 * (t - s)), 0)
)

# Funções médias para as covariáveis
mean_funcs_static = list(X1 = function(t) 1.0, X2 = function(t) 0.5)
mean_funcs_dynamic = list(
  X1 = function(t) 1.0 + 0.2*t, 
  X2 = function(t) 0.5 + 0.3*sqrt(t)
)

# Efeito das covariáveis em Y(t): Ambos positivos
concurrent_effects_on_Y = list(X1 = function(t) 1.5, X2 = function(t) 1.0)

# Função de propensão: X1 e X2 aumentam a chance de tratamento
propensity_func_confounding = list(`1` = function(X_list, t_idx) {
  score = -1.0 + 1.5 * X_list[[1]][,t_idx] + 1.2 * X_list[[2]][,t_idx]
  return(1 / (1 + exp(-score)))
})

# Efeito de A em X (Feedback): Tratamento causa aumento em X1
feedback_effect_A_on_X1 = list(
  `1` = function(t, s) ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)
)


# --- CENÁRIO 1: O MUNDO IDEAL PARA ECKER & LUNA ---
# PGD: Covariáveis X são estáticas (não mudam após t=0). O tratamento A em t=1
# é confundido por X. X e A afetam a trajetória de Y.
# Hipóteses E&L: Todas são satisfeitas. Covariáveis são pré-tratamento e
# estáticas. Ignorabilidade condicional em X basal é válida.
# Violação Esperada: Nenhuma. O modelo deve acertar.
cfg1 = list(
  mean_baseline_covariates = mean_funcs_static,
  effects_covariates_from_A = list(
    X1 = list(`1`=function(t,s) 0), 
    X2 = list(`1`=function(t,s) 0)),
  propensity_function_list = propensity_func_confounding,
  concurrent_effects_X = concurrent_effects_on_Y,
  treatment_effect = common_treatment_effect
)

# --- CENÁRIO 2: COVARIÁVEIS DINÂMICAS, MAS EXÓGENAS ---
# PGD: Covariáveis X evoluem com o tempo. O tratamento A em t=1 é confundido
# pelo valor de X em t=1. A trajetória de X NÃO é afetada por A.
# Hipóteses E&L: VIOLADA a hipótese de covariáveis estáticas/pré-tratamento.
# O confundimento agora é uma função variante no tempo.
# Violação Esperada: Viés. O modelo E&L ajusta por X(1) (ou X(0)), mas
# isso não é suficiente para controlar o confundimento ao longo de toda
# a trajetória de Y, pois a associação X-Y persiste para t>1.
cfg2 = list(
  mean_baseline_covariates = mean_funcs_dynamic, # <-- A MUDANÇA
  effects_covariates_from_A = list(
    X1 = list(`1`=function(t,s) 0), 
    X2 = list(`1`=function(t,s) 0)),
  propensity_function_list = propensity_func_confounding,
  concurrent_effects_X = concurrent_effects_on_Y,
  treatment_effect = common_treatment_effect
)

# --- CENÁRIO 3: TREATMENT-CONFOUNDER FEEDBACK ---
# PGD: Covariáveis X evoluem com o tempo. A em t=1 é confundido por X(1).
# Crucialmente, A em t=1 agora AFETA a trajetória futura de X1(t) para t>1.
# Hipóteses E&L: VIOLADAS as hipóteses de covariáveis pré-tratamento e de
# ignorabilidade condicional em X basal. X(t) para t>1 é agora um
# mediador, não um confundidor puro.
# Violação Esperada: Viés severo. O modelo não tem como saber que X1(t>1)
# é, em parte, consequência do tratamento. Ao tentar "ajustar"
# por X, ele estará, na verdade, ajustando para parte do efeito
# do próprio tratamento, o que atenua artificialmente o FATE estimado.
cfg3 = list(
  mean_baseline_covariates = mean_funcs_dynamic,
  effects_covariates_from_A = list(
    X1 = feedback_effect_A_on_X1, 
    X2 = list(`1`=function(t,s) 0)
  ), # <-- A MUDANÇA
  propensity_function_list = propensity_func_confounding,
  concurrent_effects_X = concurrent_effects_on_Y,
  treatment_effect = common_treatment_effect
)

# --- Agrupando todos os cenários para execução ---
# Cada cenário herda os parâmetros globais e só especifica o que muda
scenario_configs = list(
  `Cenario 1 (Ideal E&L)` = cfg1,
  `Cenario 2 (Cov Dinamica)` = cfg2,
  `Cenario 3 (Feedback)` = cfg3
)

# --- Exemplo de Execução ---
# (Supondo que a função orquestradora foi adaptada para usar um parâmetro
# y_model_type, que por padrão é "concurrent")
#
all_scenarios_data = generate_simulation_scenarios(
    scenario_configs = scenario_configs,
    global_params = list(
      n = N_SAMPLE, 
      time_grid = TIME_GRID, 
      noise_sd = NOISE,
      p_treat = 0.5,
      base_y_fun = function(t) 2 + 0.5*t
    )
)

rm(
  TIME_GRID, N_SAMPLE, NOISE, TREATMENT_TIME, common_treatment_effect,
  mean_funcs_static, mean_funcs_dynamic, concurrent_effects_on_Y,
  propensity_func_confounding, feedback_effect_A_on_X1, cfg1, cfg2, cfg3,
  scenario_configs
)

source('scripts_functional/utils/plot_scenario_effects.R')
source('scripts_functional/utils/plot_scenario_trajectories.R')
source('scripts_functional/utils/validate_scenario_assumptions.R')
source('scripts_functional/utils/run_validation_suite.R')

validation_output = run_validation_suite(all_scenarios_data)
validation_output$plots

