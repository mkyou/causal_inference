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
    covariates = setdiff(
      names(data), 
      c("treatment", "Y", "true_effect", "propensity_score")
      )
    formula = as.formula(paste("Y ~ treatment +", paste(covariates, 
                                                        collapse = " + ")))

    model = lm(formula, data = data)
    estimates = coef(model)["treatment"]
    
    return(estimates)
  }
)

linear_regression_interactions = list(
  name = "linear_regression_interactions",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    covariates = setdiff(
      names(data), 
      c("treatment", "Y", "true_effect", "propensity_score")
    )
    
    formula = as.formula(
      paste("Y ~ treatment * (", paste(covariates, collapse = " + "), ")")
    )
    
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

oracle_ipw = list(
  name = "propensity_ipw_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    propensities = data$propensity_score
    
    # Calcula estimativa de efeito IPW
    ipw_estimates = mean(
      (
        data$treatment * data$Y / propensities) - (
          (1 - data$treatment) * data$Y / (1 - propensities)
          )
      )
    
    # Retorna tanto as estimativas quanto os escores de propensão
    list(estimates = ipw_estimates, propensity_scores = propensities)
  }
)

aipw_double_naive = list(
  name = "aipw_double_naive_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    library(SuperLearner)
    propensities = rep(sum(data$treatment) / nrow(data), nrow(data))
    
    covariates = setdiff(
      names(data), 
      c("treatment", "Y", "true_effect", "propensity_score")
    )
    
    aipw_model = AIPW::AIPW$new(
      Y = data$Y,
      A = data$treatment,
      W = data[, ..covariates],
      Q.SL.library = "SL.glm",
      g.SL.library = "SL.glm",
      k_split = 5,
      verbose = FALSE
    )
    
    aipw_model$obs_est$p_score = propensities  
    
    aipw_model$fit()
    
    aipw_estimates_raw = aipw_model$summary()$estimates
    aipw_estimates = mean(unlist(aipw_estimates_raw), na.rm = TRUE)
    
    list(estimates = aipw_estimates, propensity_scores = propensities)
  }
)

aipw_outcome_naive = list(
  name = "aipw_outcome_naive_estimator",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    library(SuperLearner)
    propensities = data$propensity_score
    
    covariates = setdiff(
      names(data), c("treatment", "Y", "true_effect", "propensity_score")
    )
    
    X = data[, ..covariates]
    
    n = nrow(data)
    fold_idx = sample(rep(1:5, length.out = n)) 
    
    mu_1 = numeric(n)  
    mu_0 = numeric(n)  
    
    for (fold in 1:5) {
      train_idx = which(fold_idx != fold)
      test_idx = which(fold_idx == fold)
      
      # Modelo de desfecho para tratados
      outcome_model_treated = lm(
        y ~ ., data = select(data, covariates)[train_idx, ]
      )
      
      mu_1[test_idx] = predict(
        outcome_model_treated, X[test_idx, ]
      )$predictions
      
      # Modelo de desfecho para controles
      outcome_model_control = lm(
        y ~ ., data = select(data, covariates)[train_idx, ]
      )
      
      mu_0[test_idx] = predict(
        outcome_model_control, X[test_idx, ]
      )$predictions
    }
    
    aipw_estimates = (
      mu_1 - mu_0 +
        data$treatment * (data$Y - mu_1) / propensities -
        (1 - data$treatment) * (data$Y - mu_0) / (1 - propensities)
    )
    
    list(estimates = aipw_estimates, propensity_scores = propensities)
  }
)

aipw_propensity_naive = list(
  name = "aipw_propensity_naive",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    library(grf)
    propensities = rep(sum(data$treatment) / nrow(data), nrow(data))
    covariates = setdiff(
      names(data), c("treatment", "Y", "true_effect", "propensity_score")
    )
    
    X = data[, ..covariates]
    
    n = nrow(data)
    fold_idx = sample(rep(1:5, length.out = n)) 
    
    mu_1 = numeric(n)  
    mu_0 = numeric(n)  
    
    for (fold in 1:5) {
      train_idx = which(fold_idx != fold)
      test_idx = which(fold_idx == fold)
      
      # Modelo de desfecho para tratados
      outcome_model_treated = regression_forest(
        X[train_idx, ], data$Y[train_idx]
      )
      
      mu_1[test_idx] = predict(
        outcome_model_treated, X[test_idx, ]
      )$predictions
      
      # Modelo de desfecho para controles
      outcome_model_control = regression_forest(
        X[train_idx, ], data$Y[train_idx]
      )
      
      mu_0[test_idx] = predict(
        outcome_model_control, X[test_idx, ]
      )$predictions
    }
    
    aipw_estimates = mean(
      mu_1 - mu_0 +
        data$treatment * (data$Y - mu_1) / propensities -
        (1 - data$treatment) * (data$Y - mu_0) / (1 - propensities)
    )
    
    list(estimates = aipw_estimates, propensity_scores = propensities)
  }
)

aipw_double_robust = list(
  name = "aipw_double_robust",
  estimates_heterogeneous_effects = TRUE,
  estimator = function(data) {
    library(grf)
    propensities = data$propensity_score
    covariates = setdiff(
      names(data), c("treatment", "Y", "true_effect", "propensity_score")
    )
    
    X = data[, ..covariates]
    
    n = nrow(data)
    fold_idx = sample(rep(1:5, length.out = n)) 
    
    mu_1 = numeric(n)  
    mu_0 = numeric(n)  
    
    for (fold in 1:5) {
      train_idx = which(fold_idx != fold)
      test_idx = which(fold_idx == fold)
      
      # Modelo de desfecho para tratados
      outcome_model_treated = regression_forest(
        X[train_idx, ], data$Y[train_idx]
      )
      
      mu_1[test_idx] = predict(
        outcome_model_treated, X[test_idx, ]
      )$predictions
      
      # Modelo de desfecho para controles
      outcome_model_control = regression_forest(
        X[train_idx, ], data$Y[train_idx]
      )
      
      mu_0[test_idx] = predict(
        outcome_model_control, X[test_idx, ]
      )$predictions
    }
    
    aipw_estimates = (
      mu_1 - mu_0 +
        data$treatment * (data$Y - mu_1) / propensities -
        (1 - data$treatment) * (data$Y - mu_0) / (1 - propensities)
    )
    
    list(estimates = aipw_estimates, propensity_scores = propensities)
  }
)
