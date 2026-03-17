library(data.table)
library(ggplot2)
library(patchwork)

source("scripts/dgp/generator_engine.R")
source("scripts/dgp/engine_checks.R")

# Parâmetros globais da simulação

N_SAMP = 2000
T_POINTS = 2000
T_MIN_GRID = 0.0 # Grid de tempo mínimo
T_MAX_GRID = 5.0 # Grid de tempo máximo
T_TREATMENT = 0.5 # Ponto de tratamento (agora em t=0.5 de 5)
NOISE = 0.5
PROP_INT = 0.0
SEED_VAL = 20251111

# Parâmetros para a função check_sanity
params_base = list(
  N = N_SAMP,
  T = T_POINTS,
  t_min = T_MIN_GRID,
  t_max = T_MAX_GRID,
  t_A_time = T_TREATMENT,
  propensity_intercept = PROP_INT,
  noise_sd = NOISE,
  seed = SEED_VAL
)

# Definição das Funções de Efeito Base

# Intercepto basal (comum a todos)
fun_alpha_base = function(t) {
  return(asinh(t))
}

# Efeito da covariável (comum a todos)
fun_beta1_base = function(t) {
  # Efeito em 25% do grid [0, 5] = t=1.25
  mean_point = 0.25 * (T_MAX_GRID - T_MIN_GRID) + T_MIN_GRID
  # Escalar o desvio padrão com o grid
  sd_scaled = 0.05 * (T_MAX_GRID - T_MIN_GRID)
  return(3 * dnorm(t, mean = mean_point, sd = sd_scaled))
}

#------------------------------------------------------------------------
# Cenário 1: HEAVYSINE 

beta0_heavysine = function(t, t_A) {
  # Normalizar t para [0, 1] para manter o padrão
  t_norm = (t - T_MIN_GRID) / (T_MAX_GRID - T_MIN_GRID)
  effect = 4 * sin(4 * pi * t_norm)
  effect[t < t_A] = 0
  return(effect)
}

fun_betas_1 = list(beta0_heavysine, fun_beta1_base)
params_1 = c(params_base, k = length(fun_betas_1) - 1)

dgp_data_1 = generate_fdg_causal(
  N = N_SAMP, T = T_POINTS, t_min = T_MIN_GRID, t_max = T_MAX_GRID,
  t_A_time = T_TREATMENT, propensity_intercept = PROP_INT,
  noise_sd = NOISE, seed = SEED_VAL,
  fun_alpha = fun_alpha_base, fun_betas = fun_betas_1
)

check_sanity(dgp_data_1, params_1)
print(check_propensity_scores(dgp_data_1))
print(check_ate(dgp_data_1, T_TREATMENT))

#------------------------------------------------------------------------
#Cenário 2: DOPPLER

beta0_doppler = function(t, t_A) {
  t_norm = (t - T_MIN_GRID) / (T_MAX_GRID - T_MIN_GRID)
  t_norm[t_norm == 0] = 1e-6 # Evita divisão por zero
  
  effect = sqrt(t_norm * (1 - t_norm)) *
    sin((2.1 * pi) / (t_norm + 0.05))
  effect[t < t_A] = 0
  return(effect * 5)
}

fun_betas_2 = list(beta0_doppler, fun_beta1_base)
params_2 = c(params_base, k = length(fun_betas_2) - 1)

dgp_data_2 = generate_fdg_causal(
  N = N_SAMP, T = T_POINTS, t_min = T_MIN_GRID, t_max = T_MAX_GRID,
  t_A_time = T_TREATMENT, propensity_intercept = PROP_INT,
  noise_sd = NOISE, seed = SEED_VAL,
  fun_alpha = fun_alpha_base, fun_betas = fun_betas_2
)

check_sanity(dgp_data_2, params_2)
print(check_propensity_scores(dgp_data_2))
print(check_ate(dgp_data_2, T_TREATMENT))

#------------------------------------------------------------------------
# Cenário 3: BUMP

beta0_bump = function(t, t_A) {
  # Pulso em 75% do grid [0, 5] = t=3.75
  mean_point = 0.75 * (T_MAX_GRID - T_MIN_GRID) + T_MIN_GRID
  sd_scaled = 0.03 * (T_MAX_GRID - T_MIN_GRID)
  
  effect = 10 * dnorm(t, mean = mean_point, sd = sd_scaled)
  effect[t < t_A] = 0
  return(effect)
}

fun_betas_3 = list(beta0_bump, fun_beta1_base)
params_3 = c(params_base, k = length(fun_betas_3) - 1)

dgp_data_3 = generate_fdg_causal(
  N = N_SAMP, T = T_POINTS, t_min = T_MIN_GRID, t_max = T_MAX_GRID,
  t_A_time = T_TREATMENT, propensity_intercept = PROP_INT,
  noise_sd = NOISE, seed = SEED_VAL,
  fun_alpha = fun_alpha_base, fun_betas = fun_betas_3
)

check_sanity(dgp_data_3, params_3)
print(check_propensity_scores(dgp_data_3))
print(check_ate(dgp_data_3, T_TREATMENT))

#------------------------------------------------------------------------
# Cenário 4: TEMPO REAL (BUMP)

# Usando as mesmas funções do Cenário 3
fun_betas_s4 = fun_betas_3
k_s4 = length(fun_betas_s4) - 1

# Geração dos dados estáticos (necessário para a função step)
# [Este bloco é uma réplica da Seção 2 de generator_engine.R]
set.seed(SEED_VAL)
dt_static_s4 = data.table(id = 1:N_SAMP)
X_cols_s4 = c()

if (k_s4 > 0) {
  X_matrix = matrix(rnorm(N_SAMP * k_s4), ncol = k_s4)
  X_cols_s4 = paste0("X", 1:k_s4)
  colnames(X_matrix) = X_cols_s4
  dt_static_s4 = cbind(dt_static_s4, as.data.table(X_matrix))
  
  gamma_slopes = runif(k_s4, -0.5, 0.5)
  eta = PROP_INT +
    as.matrix(dt_static_s4[, ..X_cols_s4]) %*% gamma_slopes
} else {
  eta = PROP_INT
}
propensity = 1 / (1 + exp(-eta))
dt_static_s4[, A := rbinom(.N, 1, propensity)]

# Geração dos dados funcionais iterativamente (em tempo real)
t_grid_s4 = seq(T_MIN_GRID, T_MAX_GRID, length.out = T_POINTS)
step_results = list() # Lista para armazenar os data.tables

for (i in 1:length(t_grid_s4)) {
  t_step = t_grid_s4[i]
  dt_step = generate_fdg_causal_step(
    t_scalar = t_step,
    static_data = dt_static_s4,
    fun_alpha = fun_alpha_base,
    fun_betas = fun_betas_s4,
    t_A_time = T_TREATMENT,
    noise_sd = NOISE
  )
  step_results[[i]] = dt_step
}

# Combinação de todos os data.tables da iteração
dt_curve_s4 = rbindlist(step_results)

# Reconstrução do objeto dgp_data para verificação
dgp_data_4 = list(
  static_data = dt_static_s4,
  curve_data = dt_curve_s4
)

# Parâmetros para verificação
params_4 = c(params_base, k = k_s4)

# Executar verificações
check_sanity(dgp_data_4, params_4)
print(check_propensity_scores(dgp_data_4))
print(check_ate(dgp_data_4, T_TREATMENT))
