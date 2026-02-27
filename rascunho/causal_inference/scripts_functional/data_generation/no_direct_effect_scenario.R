source("scripts_functional/scenarios_utils.R")
library(gridExtra)

n = 1000
# --- PROCESSO GERADOR DE DADOS 2 --- #
# Ausência de efeitos diretos; fortes efeitos indiretos e de covariáveis

# 1) Efeito direto de A em Y (zero)
treatment_effect_zero = list(function(t, s) 0)

# 2) Efeitos indiretos FORTES (A -> X), 4 covariáveis
effects_covariates_from_A_strong = lapply(1:4, function(k) {
  list(function(t, s) {
    if (t >= s) {
      amplitude = 0.5 + 0.2 * k
      decay_rate = 2.0
      shape_param = 0.5
      time_diff = t - s
      effect = amplitude * 
        (1 - exp(-decay_rate * time_diff)) * 
        exp(-shape_param * time_diff)
      return(effect)
    } else {
      return(0)
    }
  })
})

# 2a) Versão "zero" (sem A -> X)
effects_covariates_from_A_zero = lapply(1:4, function(k) list(function(t, s) 0))

# 3) Forte efeito das covariáveis em Y (X -> Y)
concurrent_effects_X_strong = list(
  function(t) 1.5 * log1p(t),
  function(t) 1.2 * log1p(t),
  function(t) 1.8 * log1p(t),
  function(t) 1.0 * log1p(t)
)

# 3a) Efeito das covariáveis em Y (X -> Y) fixo
concurrent_effects_X_fixed = list(
  function(t) 1.5 * log1p(t),
  function(t) 1.2 * log1p(t),
  function(t) 1.8 * log1p(t),
  function(t) 1.0 * log1p(t)
)

# 4) Propensão com confundimento baseline (X -> A)
prop_func_conf = list(function(X_list, t_index) {
  plogis(-2 + 2 * X_list[[1]][, t_index] + 1.5 * X_list[[2]][, t_index])
})

# 5) Baseline das covariáveis (mantido)
baseline_covs = generate_individual_baselines(n)

# 6) Base de Y (mantido)
base_y_fun = function(t) 1.5 + 0.5 * sin(0.3 * t)

## --- Montagem dos cenários ---

scenario_configs = list(
  "Stress 1" = list(
    concurrent_effects_X = concurrent_effects_X_fixed, # X -> Y constante
    effects_covariates_from_A = effects_covariates_from_A_zero, # sem A->X
    prop_func = NULL # sem X->A
  ),
  "Stress 2" = list(
    concurrent_effects_X = concurrent_effects_X_strong, # X -> Y muda com t
    effects_covariates_from_A = effects_covariates_from_A_zero, # sem A->X
    prop_func = NULL # sem X->A
  ),
  "Stress 3" = list(
    concurrent_effects_X = concurrent_effects_X_strong,
    effects_covariates_from_A = effects_covariates_from_A_zero, # sem A->X
    prop_func = prop_func_conf # com X->A (confundimento)
  ),
  "Stress 4" = list(
    concurrent_effects_X = concurrent_effects_X_strong,
    effects_covariates_from_A = effects_covariates_from_A_strong, # A->X forte
    prop_func = NULL # sem X->A
  ),
  "Stress 5" = list(
    concurrent_effects_X = concurrent_effects_X_strong,
    effects_covariates_from_A = effects_covariates_from_A_strong, # A->X forte
    prop_func = prop_func_conf # com X->A (confundimento)
  )
)

## --- Chamada final ---
stress_test = generate_simulation_scenarios(
  n = n,
  time_grid = seq(0, 10, by = 0.1),
  base_y_fun = base_y_fun,
  treatment_effect = treatment_effect_zero,   # direto = 0
  baseline_covariates = baseline_covs,
  scenario_configs = scenario_configs,
  noise_sd = 0.5,
  p_treat = 0.5
)

rm(
  treatment_effect_zero, 
  effects_covariates_from_A_strong,
  effects_covariates_from_A_zero,
  concurrent_effects_X_fixed,
  concurrent_effects_X_strong,
  prop_func_conf,
  baseline_covs,
  base_y_fun,
  scenario_configs
)

plots_list = unlist(lapply(stress_test, function(cenario) {
  # plot_effects_raw_from_dt() retorna uma lista com 2 plots
  plots = plot_effects_raw_from_dt(cenario$outcome)
  plots
}), recursive = FALSE)

# Organizar num grid (exemplo 3 colunas, ajuste conforme quiser)
grid.arrange(grobs = plots_list, nrow = 5)

rm(plots_list, cenario)


