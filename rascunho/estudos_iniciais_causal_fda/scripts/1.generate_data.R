#' Generate Simulated Data
#'
#' @description Function to generate simulated data with control over
#' the distributions for all features, including relevant and nuisance features.
#' Allows modeling heterogeneous causal effects, propensity scores, 
#' and observed outcomes.
#'
#' @param n Number of sample units.
#' @param effect_formula Formula's right side for the true causal effect. 
#' @param propensity_score_formula Formula's right side for calculating # MUDAR
#' the true propensity score. The true propensity score will be the logistic
#' function of the formula.
#' @param outcome_formula Formula's right side for calculating the baseline 
#' outcome (before adding the causal effect).
#' @param features_dist Vector or list of distribution functions to generate 
#' all features.
#' If the number of features exceeds the length of `features_dist`, 
#' standard normal distribution (`rnorm`) will be used to generate the 
#' remaining features, with a warning.
#' @param noise_level Amount of noise added to the observed outcome. Default: 0.
#' @param seed Seed for reproducibility. Default: NULL (no seed).
#'
#' @return A data.table containing relevant features, nuisance features, 
#' propensity scores, treatment assignments, true causal effects, 
#' and observed outcomes.
#'
#' @examples
#' generate_data(
#'   n = 100,
#'   effect_formula = "0.5 * X1 + 0.3 * X2^2",
#'   propensity_score_formula = "0.3 * X1 - 0.2 * X2",
#'   outcome_formula = "2 + 0.8 * X1 - 0.5 * X2",
#'   features_dist = list(rnorm, runif, rpois(n = ., lambda = 2)),
#'   noise_level = 1,
#'   seed = 42
#' )
generate_data = function(n, effect_formula, propensity_score_formula, 
                         outcome_formula, features_dist, noise_level = 0, 
                         seed = NULL) {
  library(data.table)
  
  # Set seed for reproducibility
  if (!is.null(seed)) set.seed(seed)
  
  # Helper function to extract variables from an expression
  extract_variables = function(expression_text) {
    expr = parse(text = expression_text)[[1]]
    all.vars(expr)
  }
  
  # Determine relevant variables from the formulas
  relevant_variables = unique(c(
    extract_variables(effect_formula),
    extract_variables(propensity_score_formula),
    extract_variables(outcome_formula)
  ))
  
  n_relevant_features = length(relevant_variables)
  n_total_features = length(features_dist)
  
  # Warning if there are not enough distributions provided
  if (n_total_features < n_relevant_features) {
    warning("The number of distributions in 'features_dist' is less than 
            the number of features. Using 'rnorm' for the remainder.")
    features_dist = c(
      features_dist,
      rep(list(rnorm),
          n_relevant_features - n_total_features)
      )
  }
  
  # Initialize relevant features data.table with n rows
  # Initialize with an ID column for row tracking
  relevant_features = data.table(id = 1:n) 
  
  # Generate relevant features
  for (i in seq_along(relevant_variables)) {
    dist_func = features_dist[[i]]
    relevant_features[, (relevant_variables[i]) := dist_func(n)]
  }
  
  # Drop the ID column after initialization
  relevant_features[, id := NULL]
  
  # Generate nuisance features
  # Initialize with an ID column for row tracking
  nuisance_features = data.table(id = 1:n) 
  if (n_total_features > n_relevant_features) {
    for (i in seq((n_relevant_features + 1), n_total_features)) {
      dist_func = features_dist[[i]]
      nuisance_features[, paste0("Z", i - n_relevant_features) := dist_func(n)]
    }
  }
  nuisance_features[, id := NULL]
  
  # Calculate propensity score
  relevant_features[, propensity_score := eval(parse(
    text = propensity_score_formula
    ), envir = .SD)]
  relevant_features[, propensity_score := plogis(propensity_score)]
  
  # Assign treatment
  relevant_features[, treatment := rbinom(.N, 1, propensity_score)]
  
  # Calculate true causal effect
  relevant_features[, true_effect := eval(parse(text = effect_formula), 
                                          envir = .SD)]
  
  # Generate observed outcome
  relevant_features[, Y := eval(parse(text = outcome_formula), envir = .SD) +
                      true_effect * treatment + rnorm(.N, sd = noise_level)] 
  
  # Combine the final data
  data = cbind(relevant_features, nuisance_features)
  
  return(data)
}

