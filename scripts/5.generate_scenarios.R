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

results_model_ipw = simulate_pipeline(
  n_sim = 1000,
  estimator = model_ipw,
  data_params = data_params,
  id_scenario = "model_propensity_scenario"
)

results_oracle_ipw = simulate_pipeline(
  n_sim = 1000,
  estimator = oracle_ipw,
  data_params = data_params,
  id_scenario = "true_propensity_scenario"
)
