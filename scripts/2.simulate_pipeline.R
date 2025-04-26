#' Simulate Pipeline
#'
#' @description Runs a simulation pipeline for causal effect estimation, 
#' generating multiple datasets based on fixed parameters and applying an 
#' estimator to each dataset in parallel. 
#' Tracks execution progress with a real-time progress bar.
#'
#' @param n_sim Integer. The number of simulations to run.
#' @param estimator A list containing the estimator's metadata:
#' - `name`: Character. The name of the estimator (e.g., "difference_in_means").
#' - `estimates_heterogeneous_effects`: Logical. Whether the estimator supports
#'  heterogeneous effects.
#' - `estimator`: Function. Function for applying the estimator to a dataset.
#' @param data_params A list of fixed parameters passed to `generate_data`. 
#' Should include all the required arguments for the `generate_data` function.
#' @param id_scenario Character. A string identifier for the scenario being 
#' simulated.
#'
#' @return A list of `data.table` objects, where each `data.table` contains:
#'  - `id_scenario`: The scenario identifier (repeated for all rows).
#'  - `sim_id`: The simulation ID (unique for each simulation).
#'  - `effect_formula`: The formula representing the true causal effect
#'    (repeated for all rows).
#'  - `estimate`: The estimated causal effect (average or heterogeneous).
#'  - `true_effect`: The true causal effect for each observation.
#'  - `true_propensity_score`: The true propensity score for each observation 
#'    (included for reference).
#'  - `relevant_covariates`: Covariates affecting the causal effect or 
#'  propensity score.
#'
#' @details 
#' - Uses parallel execution (`future.apply::future_lapply`) to speed up 
#'  simulations.
#' - Implements real-time progress tracking using the `progressr` package.
#' - Applies user-defined estimators, which compute propensity scores 
#' internally as needed.
#'
#' @examples
#' # Define the difference in means estimator
#' difference_in_means = list(
#'   name = "difference_in_means",
#'   estimates_heterogeneous_effects = FALSE,
#'   estimator = function(data) {
#'     treated_mean = mean(data[data$treatment == 1, Y])
#'     control_mean = mean(data[data$treatment == 0, Y])
#'     treated_mean - control_mean
#'   }
#' )
#'
#' # Parameters for generating data
#' data_params = list(
#'   n = 10000,
#'   effect_formula = "3.5",
#'   propensity_score_formula = ".5 * X1",
#'   outcome_formula = "2 + 0.8 * X1 - 0.5 * X2",
#'   features_dist = list(rnorm, runif),
#'   noise_level = 1
#' )
#'
#' # Run the simulation pipeline with progress tracking
#' results = simulate_pipeline(
#'   n_sim = 10, # Run 10 simulations
#'   estimator = difference_in_means,
#'   data_params = data_params,
#'   id_scenario = "linear_outcome_scenario"
#' )
#'
#' # Inspect the structure of the first simulation result
#' str(results[[1]])
#'
#' # Example analysis: Compare estimates to true effects
#' head(results[[1]])
#' 

simulate_pipeline = function(n_sim, estimator, data_params, id_scenario) {
  library(future.apply)
  library(dplyr)
  library(data.table)
  library(progressr)
  
  # Initialize progress tracking
  progressr::handlers(global = TRUE)  
  p = progressr::progressor(along = seq_len(n_sim))  
  
  # Verifications
  if (is.null(estimator$name) || is.null(estimator$estimator) || 
      is.null(estimator$estimates_heterogeneous_effects)) {
    stop("Estimator must be a list containing 'name', 'estimator', 
         and 'estimates_heterogeneous_effects'.")
  }
  
  if (missing(id_scenario) || !is.character(id_scenario)) {
    stop("The parameter 'id_scenario' must be a valid string 
         identifying the scenario.")
  }
  
  # Enable parallelization
  future::plan(future::multisession)
  
  # Run simulations in parallel with progress tracking
  results = future_lapply(seq_len(n_sim), function(sim_id) {
    p()  
    
    # Generate dataset
    data = do.call(generate_data, data_params)
    
    # Apply estimator
    estimator_result = estimator$estimator(data)
    
    # Extract estimates
    estimates = if (
      is.list(estimator_result) && 
      "estimates" %in% names(estimator_result)
      ) {
      estimator_result$estimates
    } else {
      estimator_result
    }
    
    # Extract true values
    true_effects = data$true_effect
    true_propensity_scores = data$propensity_score  # useful for reference
    
    # Identify relevant covariates
    relevant_vars = unique(c(
      all.vars(parse(text = data_params$effect_formula)),
      all.vars(parse(text = data_params$propensity_score_formula)),
      all.vars(parse(text = data_params$outcome_formula))
    ))
    
    relevant_covariates = if (length(relevant_vars) > 0) {
      data[, ..relevant_vars]
    } else {
      data.table()  
    }
    
    # Structure result as a data.table
    result_dt = as.data.table(cbind(
      id_scenario = rep(id_scenario, nrow(data)),
      sim_id = rep(sim_id, nrow(data)),
      outcome_formula = rep(data_params$outcome_formula, nrow(data)),
      effect_formula = rep(data_params$effect_formula, nrow(data)),
      estimate = if (is.numeric(estimates) && length(estimates) == 1) {
        rep(estimates, nrow(data))
      } else {
        estimates
      },
      true_effect = true_effects,
      true_propensity_score = true_propensity_scores,
      relevant_covariates
    ))
    
    return(result_dt)
  }, future.seed = TRUE)  # Ensures proper random number generation in parallel
  
  future::plan(future::sequential)  
  return(results)
}

