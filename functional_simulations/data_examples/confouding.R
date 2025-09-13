# EXECUTIVE_SUMMARY_DRAFT
library(data.table)

# ===================================================================
# SCRIPT DE GERAÇÃO DE CENÁRIOS DE SIMULAÇÃO (V2.0)
#
# Refatorado para usar a nova arquitetura de PGD sequencial.
# Gera um conjunto sistemático de cenários para testar a robustez
# de modelos de g-formula funcional.
# ===================================================================

# --- 1. CARREGAMENTO DAS FUNÇÕES DE GERAÇÃO E ANÁLISE ---
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

source("functional_simulations/data_generation_process/analyze_scenario.R")

# ===================================================================
# --- 2. DEFINIÇÃO DOS COMPONENTES GLOBAIS DA SIMULAÇÃO ---
# ===================================================================

# -- Parâmetros Globais --
GLOBAL_PARAMS = list(
  n_subjects = 3000,
  time_grid = seq(0, 20, by = 0.1),
  intercept_sd_vec = c(1.0, 1.0, 1.0), # Heterogeneidade de base
  noise_sd_vec = c(2, 2, 2), # Ruído nas covariáveis
  outcome_noise_sd = 5, # Ruído no desfecho
  regime_a = c(1, 1, 1), # Preenchimento para S=3
  regime_b = c(0, 0, 0)
)

# -- Componentes Funcionais Globais --
K = 3 # Número de covariáveis

# Efeito A -> Y (Direto, no Desfecho)
# Efeito que decai com o tempo desde a aplicação
DIRECT_EFFECT_A_ON_Y = list(
  function(t, s) ifelse(t >= s, 1.5 * exp(-0.5 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 1.2 * exp(-0.4 * (t - s)), 0),
  function(t, s) ifelse(t >= s, 1.0 * exp(-0.3 * (t - s)), 0)
)

# Efeito X -> Y (Concorrente, no Desfecho)
CONCURRENT_EFFECT_X_ON_Y = list(
  function(t) 1.5 + 0.5 * t,    # Efeito de X1 variante no tempo
  function(t) 0.8,    # Efeito de X2
  function(t) -1.0    # Efeito de X3
)

# Trajetória Base do Desfecho
BASELINE_Y_FUNC = function(t) 10 + 0.5 * t

# ===================================================================
# --- 3. DEFINIÇÃO DOS CENÁRIOS ---
# ===================================================================

scenario_configs = list()

# --- PUXADA 1: CENÁRIOS CANÔNICOS (S=1, X ESTÁTICO, SEM FEEDBACK) ---
# O objetivo aqui é validar o controle do confundimento na sua forma mais pura.

BETA_ZERO = matrix(list(NULL), K, K)
GAMMA_ZERO = matrix(list(function(t, s) 0), K, 1)

scenario_configs = c(scenario_configs, list(
  `C1_Can_NoConf` = list(
    treatment_times = 5,
    baseline_funcs = list(function(t) 2, function(t) 1, function(t) 3),
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO,
    theta_matrix = matrix(c(0, 0, 0, 0), nrow = 1),
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C2_Can_LowConf` = list(
    treatment_times = 5,
    baseline_funcs = list(function(t) 2, function(t) 1, function(t) 3),
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO,
    theta_matrix = matrix(c(-0.5, 0.8, 0.6, -0.5), nrow = 1),
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C3_Can_MedConf` = list(
    treatment_times = 5,
    baseline_funcs = list(function(t) 2, function(t) 1, function(t) 3),
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO,
    theta_matrix = matrix(c(-1.0, 1.5, 1.2, -1.0), nrow = 1),
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C4_Can_HighConf` = list(
    treatment_times = 5,
    baseline_funcs = list(function(t) 2, function(t) 1, function(t) 3),
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO,
    theta_matrix = matrix(c(-1.5, 2.5, 2.0, -2.0), nrow = 1),
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y[1],
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))


# --- PUXADA 2: EFEITOS CRUZADOS (S=3, X DINÂMICO, CONF. MÉDIO) ---
# O objetivo é testar o ajuste quando as covariáveis se influenciam.

BASELINE_DYNAMIC = list(
  function(t) 0.5*t, function(t) 1, function(t) 10 - 0.5*t
)
THETA_MED_S3 = matrix(
  c(-1, 1.5, 1.2, -1.0, -1, 1.5, 1.2, -1.0, -1, 1.5, 1.2, -1.0),
  nrow = 3, byrow = TRUE
)
GAMMA_ZERO_S3 = matrix(list(function(t,s) 0), K, 3)

beta_low = matrix(list(NULL), K, K); beta_low[[2,1]] = function(t) 0.05
beta_high = matrix(list(NULL), K, K); beta_high[[2,1]] = function(t) 0.25

scenario_configs = c(scenario_configs, list(
  `C5_Cross_No` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = BETA_ZERO,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = THETA_MED_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C6_Cross_Low` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_low,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = THETA_MED_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C7_Cross_High` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = THETA_MED_S3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))

# --- PUXADA 3: CONFUNDIMENTO VARIANTE NO TEMPO (S=3, EFEITO CRUZADO ALTO) ---
# O objetivo é testar o ajuste quando o confundimento varia no tempo.

theta_a1 = THETA_MED_S3; theta_a1[2:3, 2:4] = 0
theta_a2 = THETA_MED_S3; theta_a2[c(1,3), 2:4] = 0
theta_a1a3 = THETA_MED_S3; theta_a1a3[2, 2:4] = 0

scenario_configs = c(scenario_configs, list(
  `C8_TimeVarying_A1` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = theta_a1,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C9_TimeVarying_A2` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = theta_a2,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  ),
  `C10_TimeVarying_A1A3` = list(
    treatment_times = c(5, 10, 15),
    baseline_funcs = BASELINE_DYNAMIC,
    beta_funcs_matrix = beta_high,
    gamma_matrix = GAMMA_ZERO_S3,
    theta_matrix = theta_a1a3,
    direct_effect_funcs = DIRECT_EFFECT_A_ON_Y,
    concurrent_effect_funcs = CONCURRENT_EFFECT_X_ON_Y,
    baseline_y_func = BASELINE_Y_FUNC
  )
))


# ===================================================================
# --- 4. EXECUÇÃO DA SIMULAÇÃO E ANÁLISE ---
# ===================================================================

# Gera todos os cenários definidos acima
cat("--- Iniciando a geração de dados para todos os cenários... ---\n")
all_scenarios = generate_simulation_scenarios(
  scenario_configs = scenario_configs,
  global_params = GLOBAL_PARAMS
)
cat("--- Geração de dados concluída. ---\n")

# ===================================================================
# --- 5. EXECUÇÃO DA ANÁLISE (USANDO PURRR) ---
# ===================================================================

# Usamos `imap` para iterar sobre a lista de cenários e seus nomes,
# executando a análise e ARMAZENANDO os resultados.
cat("\n--- Iniciando a análise de todos os cenários... ---\n")
all_analysis_results = purrr::imap(
  all_scenarios,
  function(scenario_data, scenario_name) {
    
    full_config = utils::modifyList(
      GLOBAL_PARAMS,
      scenario_configs[[scenario_name]]
    )
    
    # A função `analyze_scenario` faz todo o trabalho pesado
    analysis = analyze_scenario(
      scenario_data = scenario_data,
      scenario_config = full_config,
      scenario_name = scenario_name
    )
    
    return(analysis)
  }
)
cat("--- Análise de todos os cenários concluída. ---\n")

# ===================================================================
# --- 6. APRESENTAÇÃO DOS RESULTADOS (USANDO IWALK) ---
# ===================================================================

# Agora, usamos `iwalk` para iterar sobre os RESULTADOS da análise
# apenas para o efeito colateral de imprimir no console.
cat("\n\n--- APRESENTANDO RESULTADOS DA SUÍTE DE VALIDAÇÃO ---\n")
purrr::iwalk(
  all_analysis_results,
  function(results, scenario_name) {
    
    # Imprime um cabeçalho claro para cada cenário
    cat(paste("\n\n=======================================================\n"))
    cat(paste("  Resultados para o Cenário:", scenario_name))
    cat("\n=======================================================\n\n")
    
    # Imprime a tabela de resumo da estrutura causal
    print(results$summary_table)
    
    # Imprime o gráfico do efeito causal total
    print(results$plots$causal_effect)
    
    # Opcional: descomente para imprimir os outros gráficos
    print(results$plots$trajectories)
    print(results$plots$confounding_structure)
  }
)
