difference_in_means = list(
  name = "difference_in_means",
  estimates_heterogeneous_effects = FALSE,
  estimator = function(data) {
    treated_mean = mean(data[data$treatment == 1, Y])
    control_mean = mean(data[data$treatment == 0, Y])
    estimates = treated_mean - control_mean
    return(estimates)
  }
)


linear_regression_main_effects = list(
  name = "linear_regression_main_effects",
  estimates_heterogeneous_effects = FALSE,
  estimator = function(data) {
    # Identifica todas as covariáveis disponíveis (excluindo "treatment" e "Y")
    covariates = setdiff(names(data), c("treatment", "Y"))
    # Fórmula da regressão com efeitos principais
    formula = as.formula(paste("Y ~ treatment +", paste(covariates, 
                                                        collapse = " + ")))
    # Ajusta o modelo de regressão linear
    model = lm(formula, data = data)
    # Retorna o coeficiente do tratamento (tau)
    estimates = coef(model)["treatment"]
    
    return(estimates)
  }
)


linear_regression_interactions = list(
  name = "linear_regression_interactions",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    covariates = setdiff(names(data), c("treatment", "Y"))
    formula = as.formula(paste("Y ~ treatment * (", paste(covariates, collapse = " + "), ")"))
    model = lm(formula, data = data)
    
    # Predição considerando todos tratados
    data_treated = data
    data_treated$treatment = 1
    pred_treated = predict(model, newdata = data_treated)
    
    # Predição considerando todos controles
    data_control = data
    data_control$treatment = 0
    pred_control = predict(model, newdata = data_control)
    
    # Retorna a diferença entre os cenários
    estimates = pred_treated - pred_control
    
    return(estimates)
  }
)

naive_ipw = list(
  name = "propensity_ipw_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    covariates = setdiff(names(data), c("treatment", "Y", "true_effects", 
                                        "propensity_score"))
  
    propensities = rep(sum(data$treatment) / nrow(data), nrow(data))
    
    # Calcula as estimativas IPW individuais
    ipw_estimates = mean(
      (data$treatment * data$Y /
         propensities) - ((1 - data$treatment) * data$Y /
                            (1 - propensities)))
    
    # Retorna tanto as estimativas quanto os escores de propensão
    list(estimates = ipw_estimates, propensity_scores = propensities)
  }
)


model_ipw = list(
  name = "propensity_ipw_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    covariates = setdiff(names(data), c("treatment", "Y", "true_effects", 
                                        "propensity_score"))
    
    formula = as.formula(paste("treatment ~",
                               paste(covariates, collapse = " + ")))
    
    propensity_model = gbm::gbm(
      formula,
      data = data,
      distribution = "bernoulli",
      n.trees = 100,
      interaction.depth = 3,
      shrinkage = 0.01,
      cv.folds = 5,
      verbose = FALSE
    )
    
    best_trees = gbm::gbm.perf(propensity_model, method = "cv", 
                               plot.it = FALSE)
    
    propensities = predict(propensity_model,
                           newdata = data, 
                           n.trees = best_trees, 
                           type = "response")
    
    # Calcula as estimativas IPW individuais
    ipw_estimates = mean(
      (data$treatment * data$Y /
         propensities) - ((1 - data$treatment) * data$Y /
                            (1 - propensities)))
    
    # Retorna tanto as estimativas quanto os escores de propensão
    list(estimates = ipw_estimates, propensity_scores = propensities)
  }
)


oracle_ipw = list(
  name = "propensity_ipw_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    covariates = setdiff(names(data), c("treatment", "Y", "true_effects", 
                                        "propensity_score"))
    propensities = data$propensity_score
    
    # Calcula as estimativas IPW individuais
    ipw_estimates = mean(
      (data$treatment * data$Y /
         propensities) - ((1 - data$treatment) * data$Y /
                            (1 - propensities)))
    
    # Retorna tanto as estimativas quanto os escores de propensão
    list(estimates = ipw_estimates, propensity_scores = propensities)
  }
)