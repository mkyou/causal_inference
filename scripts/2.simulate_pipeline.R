#' Simulate Pipeline
#'
#' @description A function to run a simulation pipeline for causal effect 
#' estimation, including optional estimation of propensity scores. Generates 
#' multiple datasets based on fixed parameters and applies an estimator to each 
#' dataset in parallel. Returns structured results with key components for 
#' validation.
#'
#' @param n_sim Integer. The number of simulations to run.
#' @param estimator A list containing the estimator's metadata:
#'  - `name`: Character. The name of the estimator (e.g., 
#'    "difference_in_means").
#'  - `estimates_heterogeneous_effects`: Logical. Whether the estimator supports 
#'    heterogeneous effects.
#'  - `estimator`: Function. The function for applying the estimator to a 
#'    dataset.
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
#'  - `true_propensity_score`: The true propensity score for each observation.
#'  - `estimated_propensity_score`: The estimated propensity score, if the 
#'    estimator computes it, or the proportion of treated individuals otherwise.
#'  - `relevant_covariates`: Covariates affecting the causal effect or 
#'    propensity score.
#'
#' @examples
#' # Define the difference in means estimator
#' difference_in_means = list(
#'   name = "difference_in_means",
#'   estimates_heterogeneous_effects = FALSE,
#'   estimates_propensity_scores = FALSE,
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
#' # Run the simulation pipeline
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

simulate_pipeline = function(n_sim, estimator, data_params, id_scenario) {
  library(future.apply)
  library(dplyr)
  library(data.table)
  
  # Verificações preliminares
  if (is.null(estimator$name) || is.null(estimator$estimator) || 
      is.null(estimator$estimates_heterogeneous_effects)) {
    stop("Estimator must be a list containing 'name', 'estimator', and 'estimates_heterogeneous_effects'.")
  }
  
  if (missing(id_scenario) || !is.character(id_scenario)) {
    stop("The parameter 'id_scenario' must be a valid string identifying the scenario.")
  }
  
  # Planejamento de paralelização
  future::plan(future::multisession) # Multithreading compatível com Windows e Unix
  
  # Simulações paralelas
  results = future_lapply(seq_len(n_sim), function(sim_id) {
    # Gerar dataset
    data = do.call(generate_data, data_params)
    
    # Aplicar estimador para estimar efeitos causais
    estimator_result = estimator$estimator(data)
    
    # Extração das estimativas
    estimates = if (is.list(estimator_result) && "estimates" %in% names(estimator_result)) {
      estimator_result$estimates
    } else {
      estimator_result
    }
    
    # Escore de propensão estimado
    estimated_propensity_score = if (is.list(estimator_result) && "propensity_scores" %in% names(estimator_result)) {
      estimator_result$propensity_scores
    } else {
      rep(sum(data$treatment) / nrow(data), nrow(data)) # Proporção de tratados como fallback
    }
    
    # Valores verdadeiros
    true_effects = data$true_effect
    true_propensity_scores = data$propensity_score
    
    # Variáveis relevantes
    relevant_vars = unique(c(
      all.vars(parse(text = data_params$effect_formula)),
      all.vars(parse(text = data_params$propensity_score_formula)),
      all.vars(parse(text = data_params$outcome_formula))
    ))
    
    relevant_covariates = if (length(relevant_vars) > 0) {
      data[, ..relevant_vars]
    } else {
      data.table() # Retorna vazio se não houver variáveis relevantes
    }
    
    # Estruturar resultado como data.table
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
      estimated_propensity_score = estimated_propensity_score,
      relevant_covariates
    ))
    
    # Retorna o resultado estruturado
    return(result_dt)
  }, future.seed = TRUE) # Assegura geração correta de números aleatórios em paralelo
  
  # Resetar plano de paralelização para evitar interferências futuras
  future::plan(future::sequential)
  
  # Retorna os resultados como uma lista de data.tables
  return(results)
}





