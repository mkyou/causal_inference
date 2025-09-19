# ===================================================================
# SCRIPT DE GERAÇÃO DE CENÁRIOS DE FEEDBACK (A -> X) (V2.0)
#
# Refatorado para espelhar a estrutura sistemática do 'confouding.R'.
# Foco: Testar a robustez do modelo ao feedback, tanto de forma
# isolada quanto em conjunto com efeitos cruzados entre covariáveis.
# ===================================================================

# --- 1. CARREGAMENTO DAS FUNÇÕES DE GERAÇÃO ---
source(
  "functional_simulations/data_generation_process/generate_sequential_dgp.R"
)
source(
  "functional_simulations/data_generation_process/generate_functional_outcomes.R"
)
source(
  "functional_simulations/data_generation_process/generate_simulation_scenarios.R"
)
source(
  "functional_simulations/data_generation_process/prepare_refund_data.R"
)

# ===================================================================
# --- 2. DEFINIÇÃO DOS COMPONENTES GLOBAIS DA SIMULAÇÃO ---
# ===================================================================

# -- Parâmetros Globais (consistente com 'confouding.R') --
GLOBAL_PARAMS = list(
  n_subjects = 3000,
  time_grid = seq(0, 20, by = 0.1),
  intercept_sd_vec = c(1.0, 1.0, 1.0),
  noise_sd_vec = c(2, 2, 2),
  outcome_noise_sd = 5,
  regime_a = c(1, 1, 1), # Preenchimento para S=3
  regime_b = c(0, 0, 0)
)
K = 3 # Número de covariáveis

# --- Componentes Funcionais Globais (reutilizados do 'confouding.R') ---

DIRECT_EFFECT_A_ON_Y = list(
  function(t, s) ifelse(t >= s, 1.5 * exp(-0.5 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 1.2 * exp(-0.4 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 1.0 * exp(-0.3 * (t - s)), 0)
)
CONCURRENT_EFFECT_X_ON_Y = list(
  function(t) 1.5 + 0.05 * t,
  function(t) 0.8,
  function(t) -1.0
)
BASELINE_Y_FUNC = function(t) 10 + 0.5 * t
BASELINE_STATIC = list(function(t) 2, function(t) 1, function(t) 3)
BASELINE_DYNAMIC = list(
  function(t) 0.5 * t, function(t) 1, function(t) 10 - 0.5 * t
)

# --- MATRIZES DE EFEITO PARA OS CENÁRIOS DE FEEDBACK ---

# THETA (X -> A): Sem confundimento em todos os cenários
THETA_ZERO_S1 = matrix(c(0, 0, 0, 0), nrow = 1)
THETA_ZERO_S3 = matrix(0, nrow = 3, ncol = K + 1)

# BETA (X -> X): Reutilizado do 'confouding.R' para testar efeitos cruzados
BETA_ZERO = matrix(list(NULL), K, K)
beta_low_dynamic = matrix(list(NULL), K, K)
beta_low_dynamic[[2, 1]] = function(t) 0.15 * exp(-0.02 * t)
beta_high_dynamic = matrix(list(NULL), K, K)
beta_high_dynamic[[2, 1]] = function(t) 0.5 * exp(-0.02 * t)
beta_high_dynamic[[3, 1]] = function(t) -0.3
beta_high_dynamic[[3, 2]] = function(t) 0.2 * sin(t * pi / 10)


# GAMMA (A -> X): O foco deste script, configurado para feedback médio
GAMMA_ZERO_S1 = matrix(list(function(t, s) 0), K, 1)
GAMMA_ZERO_S3 = matrix(list(function(t, s) 0), K, 3)

GAMMA_MEDIUM_S1 = matrix(list(function(t, s) 0), K, 1)
GAMMA_MEDIUM_S1[[1, 1]] = function(t, s) { # A1 afeta X1
  ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)
}

GAMMA_MEDIUM_S3 = matrix(list(function(t, s) 0), K, 3)
GAMMA_MEDIUM_S3[[1, 1]] = function(t, s) { # A1 afeta X1
  ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)
}
GAMMA_MEDIUM_S3[[2, 2]] = function(t, s) { # A2 afeta X2
  ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)
}
GAMMA_MEDIUM_S3[[3, 3]] = function(t, s) { # A3 afeta X3
  ifelse(t > s, 1.5 * exp(-0.5 * (t - s)), 0)
}


# ===================================================================
# --- 3. DEFINIÇÃO DOS CENÁRIOS DE FEEDBACK ---
# ===================================================================

scenario_configs_feedback = list()

# --- PUXADA 1: FEEDBACK BÁSICO (S=1, X ESTÁTICO, SEM EFEITO CRUZADO) ---
# O objetivo é validar o ajuste ao feedback na sua forma mais pura.

scenario_configs_feedback = c(scenario_configs_feedback, list(
  `F1_Base_NoFB` = list(
    treatment_times = 5, baseline_funcs = BASELINE_STATIC,
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO_S1, # <-- Sem Feedback
    theta_matrix = THETA_ZERO_S1, # <-- Sem Confundimento
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `F2_Base_MedFB` = list(
    treatment_times = 5, baseline_funcs = BASELINE_STATIC,
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_MEDIUM_S1, # <-- Feedback Médio
    theta_matrix = THETA_ZERO_S1, # <-- Sem Confundimento
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))

# --- PUXADA 2: FEEDBACK COM EFEITOS CRUZADOS (S=3, X DINÂMICO) ---
# O objetivo é testar o ajuste quando o feedback ocorre
# simultaneamente com a influência mútua entre as covariáveis.

scenario_configs_feedback = c(scenario_configs_feedback, list(
  `F3_Cross_No` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = BETA_ZERO, # <-- Sem Efeito Cruzado
    gamma_matrix = GAMMA_MEDIUM_S3, # <-- Feedback Médio
    theta_matrix = THETA_ZERO_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `F4_Cross_Low` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_low_dynamic, # <-- Efeito Cruzado Baixo
    gamma_matrix = GAMMA_MEDIUM_S3, # <-- Feedback Médio
    theta_matrix = THETA_ZERO_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `F5_Cross_High` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high_dynamic, # <-- Efeito Cruzado Alto
    gamma_matrix = GAMMA_MEDIUM_S3, # <-- Feedback Médio
    theta_matrix = THETA_ZERO_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))

# ===================================================================
# --- 4. EXECUÇÃO DA SIMULAÇÃO ---
# ===================================================================

# Gera todos os cenários de feedback definidos acima
cat("--- Iniciando geração de dados para cenários de FEEDBACK... ---\n")
all_feedback_scenarios = generate_simulation_scenarios(
  scenario_configs = scenario_configs_feedback,
  global_params = GLOBAL_PARAMS,
  prepare_for_refund = TRUE
)
cat("--- Geração de dados de feedback concluída. ---\n")