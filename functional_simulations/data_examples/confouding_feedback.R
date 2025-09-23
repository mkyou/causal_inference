# ===================================================================
# SCRIPT DE GERAÇÃO DE CENÁRIOS COMBINADOS (V1.0)
#
# Foco: Testar a robustez do modelo a cenários complexos que
# incluem simultaneamente:
#   1. Confundimento (X -> A)
#   2. Feedback (A -> X)
#   3. Efeitos cruzados entre covariáveis (X -> X)
#
# Este é um teste de estresse fundamental para a g-formula funcional.
# ===================================================================

# --- 1. CARREGAMENTO DAS FUNÇÕES DE GERAÇÃO ---
# Assumindo que estes scripts estão no seu ambiente de trabalho
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

# -- Parâmetros Globais (consistente com outros scripts) --
GLOBAL_PARAMS = list(
  n_subjects = 3000,
  time_grid = seq(0, 20, by = 0.1),
  intercept_sd_vec = c(1.0, 1.0, 1.0),
  noise_sd_vec = c(2, 2, 2),
  outcome_noise_sd = 5,
  regime_a = c(1, 1, 1),
  regime_b = c(0, 0, 0)
)
K = 3 # Número de covariáveis

# --- Componentes Funcionais Globais (reutilizados) ---
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
BASELINE_DYNAMIC = list(
  function(t) 0.5 * t, function(t) 1, function(t) 10 - 0.5 * t
)

# --- MATRIZES DE EFEITO PARA OS CENÁRIOS COMBINADOS ---

# THETA (X -> A): Confounding Médio, S=3 (de confounding.R)
THETA_MED_S3_DYNAMIC = matrix(
  c(
    -1.0, 2.0, 0.8, -0.5,
    -1.0, 0.5, 1.5, -1.2,
    -0.8, 0.8, 0.8, -0.8
  ),
  nrow = 3, byrow = TRUE
)

# GAMMA (A -> X): Feedback Médio, S=3 (de feedback.R)
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

# BETA (X -> X): Variaremos este componente (de confounding.R)
BETA_ZERO = matrix(list(NULL), K, K)

beta_low_dynamic = matrix(list(NULL), K, K)
beta_low_dynamic[[2, 1]] = function(t) 0.15 * exp(-0.02 * t)

beta_high_dynamic = matrix(list(NULL), K, K)
beta_high_dynamic[[2, 1]] = function(t) 0.5 * exp(-0.02 * t)
beta_high_dynamic[[3, 1]] = function(t) -0.3
beta_high_dynamic[[3, 2]] = function(t) 0.2 * sin(t * pi / 10)


# ===================================================================
# --- 3. DEFINIÇÃO DOS CENÁRIOS COMBINADOS ---
# ===================================================================

scenario_configs_combined = list()

# PUXADA ÚNICA: Confundimento Médio + Feedback Médio, com variação
# nos efeitos cruzados (X -> X). Todos os cenários são S=3 com
# covariáveis dinâmicas.

scenario_configs_combined = c(scenario_configs_combined, list(
  `CF1_Base` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = BETA_ZERO,          # <-- Sem Efeito Cruzado
    gamma_matrix = GAMMA_MEDIUM_S3,         # <-- Feedback Médio
    theta_matrix = THETA_MED_S3_DYNAMIC,    # <-- Confundimento Médio
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `CF2_LowCross` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_low_dynamic,   # <-- Efeito Cruzado Baixo
    gamma_matrix = GAMMA_MEDIUM_S3,         # <-- Feedback Médio
    theta_matrix = THETA_MED_S3_DYNAMIC,    # <-- Confundimento Médio
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `CF3_HighCross` = list(
    treatment_times = c(5, 10, 15), baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high_dynamic,  # <-- Efeito Cruzado Alto
    gamma_matrix = GAMMA_MEDIUM_S3,         # <-- Feedback Médio
    theta_matrix = THETA_MED_S3_DYNAMIC,    # <-- Confundimento Médio
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))

# ===================================================================
# --- 4. EXECUÇÃO DA SIMULAÇÃO ---
# ===================================================================

# Gera todos os cenários combinados definidos acima
cat("--- Iniciando geração de dados para cenários COMBINADOS (Conf+FB)... ---\n")
all_combined_scenarios = generate_simulation_scenarios(
  scenario_configs = scenario_configs_combined,
  global_params = GLOBAL_PARAMS,
  prepare_for_refund = TRUE # Mantendo a preparação para análise funcional
)
cat("--- Geração de dados combinados concluída. ---\n")