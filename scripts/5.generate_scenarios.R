library(AIPW)


# Parameters for generating data
data_params = list(
  n = 1000,
  effect_formula = "3.5 + X1",
  propensity_score_formula = ".5*X1",
  outcome_formula = "2 + 0.8 * X1 - 0.5 * X2 + X3",
  features_dist = list(rnorm, runif, rnorm, rnorm),
  noise_level = 1
)

# Run simulation pipeline
results_dm = simulate_pipeline(
  n_sim = 1000, 
  estimator = difference_in_means,
  data_params = data_params,
  id_scenario = "linear_outcome_scenario_means_diff"
)

results_dm[[1]]

results_lmm = simulate_pipeline(
  n_sim = 1000, 
  estimator = linear_regression_main_effects,
  data_params = data_params,
  id_scenario = "linear_outcome_scenario_lm_main"
)

results_lmi = simulate_pipeline(
  n_sim = 1000, 
  estimator = linear_regression_interactions,
  data_params = data_params,
  id_scenario = "linear_outcome_scenario_lm_int"
)

results_naive_ipw = simulate_pipeline(
  n_sim = 1000,
  estimator = naive_ipw,
  data_params = data_params,
  id_scenario = "naive_propensity_scenario"
)

results_oracle_ipw = simulate_pipeline(
  n_sim = 1000,
  estimator = oracle_ipw,
  data_params = data_params,
  id_scenario = "true_propensity_scenario"
)

results_aipw_double_naive = simulate_pipeline(
  n_sim = 1000,
  estimator = aipw_double_naive,
  data_params = data_params,
  id_scenario = 'true_propensity_scenario'
)

results_aipw_outcome_naive = simulate_pipeline(
  n_sim = 1000,
  estimator = aipw_outcome_naive,
  data_params = data_params,
  id_scenario = 'naive_propensity_scenario'
)


results_aipw_propensity_naive = simulate_pipeline(
  n_sim = 1000,
  estimator = aipw_propensity_naive,
  data_params = data_params,
  id_scenario = "true_propensity_scenario"
)

results_aipw_double_robust = simulate_pipeline(
  n_sim = 1000,
  estimator = aipw_double_robust,
  data_params = data_params,
  id_scenario = "true_propensity_scenario"
)
