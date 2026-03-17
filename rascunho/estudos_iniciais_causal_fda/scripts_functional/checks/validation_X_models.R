# ===================================================================
# SCRIPT DE VALIDAÇÃO FINAL DOS MODELOS DE COVARIÁVEIS (X)
# Objetivo: Comparar sistematicamente duas especificações de modelo para X
# (com e sem lags) contra a verdade simulada do PGD.
# ===================================================================

# --- 1. SETUP ---
library(data.table)
library(mgcv)
library(ggplot2)

`%||%` = function(a, b) if (is.null(a)) b else a

# --- 2. FUNÇÕES DO PROCESSO GERADOR DE DADOS (PGD) ---

create_intercept_slope_heterogeneity = function(n,
                                                mean_funcs,
                                                intercept_sd = 0.2,
                                                slope_sd = 0.05) {
  num_covariates = length(mean_funcs)
  random_intercepts = replicate(num_covariates,
                                rnorm(n, 0, intercept_sd), FALSE)
  random_slopes = replicate(num_covariates,
                            rnorm(n, 0, slope_sd), FALSE)
  lapply(1:num_covariates, function(k) {
    lapply(1:n, function(i) {
      b = random_intercepts[[k]][i]
      c = random_slopes[[k]][i]
      mean_f = mean_funcs[[k]]
      function(t) { mean_f(t) + b + c * t }
    })
  })
}

generate_confounder_pointwise = function(n,
                                         effects_function_list,
                                         p_treat,
                                         noise_sd = 1,
                                         time_grid = NULL,
                                         baseline_function_list = NULL,
                                         propensity_function_list = NULL,
                                         a_seq = NULL,
                                         another_a_seq = NULL,
                                         return_format = c("long", "wide", "both"),
                                         treatment_times = NULL) {
  
  return_format = match.arg(return_format)
  num_covariates = length(effects_function_list)
  total_treatment_times = length(effects_function_list[[1]])
  if (is.null(treatment_times)) {
    treatment_times = 1:total_treatment_times
  }
  if (is.null(time_grid)) {
    time_grid = 0:ceiling(max(treatment_times) + 1)
  }
  total_times = length(time_grid)
  if (is.null(a_seq)) a_seq = rep(1, total_treatment_times)
  if (is.null(another_a_seq)) another_a_seq = rep(0, total_treatment_times)
  
  Z_list = replicate(num_covariates,
                     matrix(rnorm(n * total_times, 0, noise_sd), n, total_times),
                     simplify = FALSE)
  A = A_a = A_another = matrix(0, n, total_times)
  X_list = a_X_list = another_a_X_list =
    lapply(1:num_covariates, function(k) matrix(0, n, total_times))
  
  effects_mat_list = lapply(1:num_covariates, function(k) {
    effects_mat = matrix(0, total_treatment_times, total_times)
    for (s_idx in 1:total_treatment_times) {
      s_val = treatment_times[s_idx]
      for (t_idx in 1:total_times) {
        t_val = time_grid[t_idx]
        if (t_val >= s_val) {
          effects_mat[s_idx, t_idx] =
            effects_function_list[[k]][[s_idx]](t_val, s_val)
        }
      }
    }
    effects_mat
  })
  
  for (k in 1:num_covariates) {
    base_vals_t0 = sapply(baseline_function_list[[k]],
                          function(f) f(time_grid[1]))
    X_list[[k]][, 1] = base_vals_t0 + Z_list[[k]][, 1]
    a_X_list[[k]][, 1] = X_list[[k]][, 1]
    another_a_X_list[[k]][, 1] = X_list[[k]][, 1]
  }
  
  for (t_index in 2:total_times) {
    current_time = time_grid[t_index]
    if (current_time %in% treatment_times) {
      s_index = which(treatment_times == current_time)
      prop_func = if (!is.null(propensity_function_list)) {
        propensity_function_list[[s_index]]
      } else { NULL }
      A[, t_index] = if (!is.null(prop_func)) {
        rbinom(n, 1, prop_func(X_list, t_index))
      } else { rbinom(n, 1, p_treat) }
      A_a[, t_index] = rep(a_seq[s_index], n)
      A_another[, t_index] = rep(another_a_seq[s_index], n)
    }
    
    treated_cols = which(time_grid < current_time &
                           time_grid %in% treatment_times)
    
    for (k in 1:num_covariates) {
      base_vals_t = sapply(baseline_function_list[[k]],
                           function(f) f(current_time))
      
      effects_obs = effects_a_val = effects_another_val = rep(0, n)
      if (length(treated_cols) > 0) {
        s_indices = sapply(time_grid[treated_cols],
                           function(tt) which(treatment_times == tt))
        effects_vec_t = effects_mat_list[[k]][s_indices, t_index,
                                              drop = FALSE]
        effects_obs = A[, treated_cols, drop = FALSE] %*% effects_vec_t
        effects_a_val = A_a[, treated_cols, drop = FALSE] %*% effects_vec_t
        effects_another_val = A_another[, treated_cols,
                                        drop = FALSE] %*% effects_vec_t
      }
      
      X_list[[k]][, t_index] = base_vals_t + effects_obs +
        Z_list[[k]][, t_index]
      a_X_list[[k]][, t_index] = base_vals_t + effects_a_val +
        Z_list[[k]][, t_index]
      another_a_X_list[[k]][, t_index] = base_vals_t +
        effects_another_val +
        Z_list[[k]][, t_index]
    }
  }
  
  wide_output = list(
    A_mat = A, A_a_mat = A_a, A_another_mat = A_another,
    X_list = X_list, a_X_list = a_X_list,
    another_a_X_list = another_a_X_list, t_vec = time_grid
  )
  
  return(list(wide = wide_output))
}


# --- 3. SUÍTE DE VALIDAÇÃO PARA OS MODELOS DE X ---

fit_gformula_models = function(observed_data, scenario_config) {
  require(mgcv)
  
  # --- 1. Preparar dados ---
  model_data = as.data.table(observed_data)
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  treatment_times = scenario_config$treatment_times
  
  if (!is.null(treatment_times)) {
    for (s in treatment_times) {
      treatment_status_s = model_data[t == s, .(id, A_s = A)]
      setnames(treatment_status_s, "A_s", paste0("A_", s))
      model_data = merge(model_data, treatment_status_s, by = "id", all.x = TRUE)
      # Preenche com 0 os NAs para indivíduos que saíram mais cedo do estudo
      model_data[is.na(get(paste0("A_", s))), (paste0("A_", s)) := 0]
    }
  }
  
  # --- 2. Ajustar um modelo GAM para cada covariável ---
  confounder_models = list()
  
  # Cria os termos para os efeitos funcionais do tratamento (feedback)
  treatment_varying_effects = if (!is.null(treatment_times)) {
    paste0("s(t, by = A_", treatment_times, ", k = 10)") # Reduzi k para estabilidade
  } else {
    "1"
  }
  
  for (cov in cov_names) {
    cat(paste("Ajustando modelo para", cov, "\n"))
    
    # --- FÓRMULA ATUALIZADA (SEM s(id, bs='re')) ---
    formula_text = as.formula(
      paste(cov, "~ s(t) +", paste(treatment_varying_effects, collapse=" + "))
    )
    
    confounder_models[[cov]] = bam(formula_text, data = model_data, 
                                   discrete = TRUE, select = TRUE)
  }
  
  # --- 3. Ajustar um modelo GAM para o desfecho Y ---
  cat("Ajustando modelo para o desfecho Y\n")
  
  cov_smooth_terms = paste0("s(", cov_names, ")")
  
  # A fórmula para Y também usa os efeitos funcionais do tratamento
  formula_y = as.formula(
    paste("y ~ s(t) +", 
          paste(cov_smooth_terms, collapse=" + "), "+", 
          paste(treatment_varying_effects, collapse=" + "))
  )
  
  outcome_model = bam(formula_y, data = model_data, discrete = TRUE, select = TRUE)
  
  # --- 4. Retornar os objetos necessários ---
  return(list(
    confounder_models = confounder_models,
    outcome_model = outcome_model,
    noise_sds = list(
      y = sd(outcome_model$residuals),
      X = sapply(confounder_models, function(m) sd(m$residuals))
    )
  ))
}

simulate_x_trajectories = function(fitted_x_models,
                                   observed_data,
                                   scenario_config,
                                   intervention_regime,
                                   include_lags) {
  all_times = sort(unique(observed_data$t))
  ids = unique(observed_data$id)
  n_ids = length(ids)
  cov_names = names(fitted_x_models$models)
  treatment_times = scenario_config$treatment_times
  
  sim_data = observed_data[t == min(all_times),
                           .SD, .SDcols = c("id", "t", cov_names)]
  
  for (t_idx in 2:length(all_times)) {
    current_t = all_times[t_idx]
    newdata = data.table(id = ids, t = current_t)
    
    previous_A = 0
    previous_tt = all_times[t_idx - 1]
    s_index = which(treatment_times == previous_tt)
    if (length(s_index) > 0) {
      previous_A = intervention_regime[s_index]
    }
    newdata[, A_lag1 := previous_A]
    
    if (include_lags) {
      previous_data = sim_data[t == all_times[t_idx - 1]]
      lags_to_add = previous_data[, .SD, .SDcols = c("id", cov_names)]
      setnames(lags_to_add, cov_names, paste0(cov_names, "_lag1"))
      newdata = merge(newdata, lags_to_add, by = "id")
    }
    
    sim_row = data.table(id = ids, t = current_t)
    for (cov in cov_names) {
      pred = predict(fitted_x_models$models[[cov]], newdata = newdata)
      sim_row[, (cov) := pred + rnorm(n_ids, 0,
                                      fitted_x_models$noise_sds[cov])]
    }
    sim_data = rbind(sim_data, sim_row)
  }
  return(sim_data)
}

run_x_model_validation = function(scenario_data,
                                  scenario_config,
                                  regime_to_validate) {
  
  wide_X_list = scenario_data$wide$X_list
  wide_A_mat = scenario_data$wide$A_mat
  time_vec = scenario_data$wide$t_vec
  cov_names = paste0("X", 1:length(wide_X_list))
  
  # Etapa A: Converter lista de matrizes X em um data.table
  observed_long_list = lapply(seq_along(wide_X_list), function(k){
    dt = as.data.table(wide_X_list[[k]])
    dt[, id := .I]
    melt(dt, id.vars = "id", variable.name = "t_idx",
         value.name = cov_names[k], variable.factor = FALSE)
  })
  observed_long_x = Reduce(
    function(x, y) merge(x, y, by = c("id", "t_idx")), observed_long_list
  )
  
  # Etapa B: Converter matriz A
  A_long = as.data.table(wide_A_mat)
  A_long[, id := .I]
  A_long = melt(A_long, id.vars="id", variable.name="t_idx",
                value.name="A", variable.factor=FALSE)
  
  # Etapa C: Juntar tudo
  observed_data = merge(observed_long_x, A_long, by=c("id", "t_idx"))
  observed_data[, t := time_vec[as.numeric(gsub("V", "", t_idx))]][, t_idx := NULL]
  setcolorder(observed_data, c("id", "t", "A", cov_names))
  
  # Fit both models
  cat("--- Ajustando Modelo A (Sem Lags) ---\n")
  models_no_lags = fit_x_models_for_gformula(observed_data, scenario_config,
                                             include_lags = FALSE)
  cat("--- Ajustando Modelo B (Com Lags) ---\n")
  models_with_lags = fit_x_models_for_gformula(observed_data, scenario_config,
                                               include_lags = TRUE)
  
  # Simulate with both models
  cat("--- Simulando Trajetórias Contrafactuais ---\n")
  sim_x_no_lags = simulate_x_trajectories(models_no_lags, observed_data,
                                          scenario_config, regime_to_validate,
                                          include_lags = FALSE)
  sim_x_with_lags = simulate_x_trajectories(models_with_lags, observed_data,
                                            scenario_config, regime_to_validate,
                                            include_lags = TRUE)
  
  # Get the ground truth
  is_a_seq = all(regime_to_validate == (scenario_config$a_seq %||% c(1)))
  true_x_wide = if(is_a_seq) {
    scenario_data$wide$a_X_list
  } else {
    scenario_data$wide$another_a_X_list
  }
  
  true_long_list = lapply(seq_along(true_x_wide), function(k) {
    dt = as.data.table(true_x_wide[[k]])
    dt[, id := .I]
    melt(dt, id.vars="id", value.name=cov_names[k],
         variable.name="t_idx", variable.factor=FALSE)
  })
  true_x_long = Reduce(
    function(x, y) merge(x, y, by = c("id", "t_idx")), true_long_list
  )
  true_x_long[, t := time_vec[as.numeric(gsub("V", "", t_idx))]][, t_idx := NULL]
  
  # Combine for plotting and metrics
  sim_x_no_lags[, Model := "Sem Lags (Sua Hipótese)"]
  sim_x_with_lags[, Model := "Com Lags (Minha Hipótese)"]
  true_x_long[, Model := "Real (PGD)"]
  
  combined_data = rbindlist(list(
    melt(sim_x_no_lags, id.vars=c("id","t","Model"),
         measure.vars=cov_names, variable.name="Covariável"),
    melt(sim_x_with_lags, id.vars=c("id","t","Model"),
         measure.vars=cov_names, variable.name="Covariável"),
    melt(true_x_long, id.vars=c("id","t","Model"),
         measure.vars=cov_names, variable.name="Covariável")
  ))
  
  metrics_data = merge(combined_data[Model != "Real (PGD)"],
                       combined_data[Model == "Real (PGD)"],
                       by = c("id", "t", "Covariável"),
                       suffixes = c(".sim", ".real"))
  metrics = metrics_data[, .(RMSE = sqrt(mean((value.sim - value.real)^2))),
                         by = .(Model, Covariável)]
  
  return(list(plot_data = combined_data, metrics = metrics))
}

plot_x_validation = function(results) {
  avg_trajectories = results$plot_data[, .(Mean_Value = mean(value)),
                                       by = .(t, Model, Covariável)]
  
  ggplot(avg_trajectories, aes(x = t, y = Mean_Value, color = Model,
                               linetype = Model)) +
    geom_line(linewidth = 1.1) +
    facet_wrap(~Covariável, scales = "free_y") +
    labs(
      title = "Validação dos Modelos de X: Real vs. Simulado",
      subtitle = "Comparando as trajetórias médias contrafactuais de X",
      x = "Tempo", y = "Valor Médio da Covariável"
    ) +
    scale_color_manual(values = c("Real (PGD)" = "black",
                                  "Com Lags (Minha Hipótese)" = "blue",
                                  "Sem Lags (Sua Hipótese)" = "red")) +
    scale_linetype_manual(values = c("Real (PGD)" = "dashed",
                                     "Com Lags (Minha Hipótese)" = "solid",
                                     "Sem Lags (Sua Hipótese)" = "dotted")) +
    theme_light() +
    theme(legend.position = "bottom")
}


# --- 5. EXECUÇÃO DO EXPERIMENTO ---
cat("--- Definindo o Cenário de Teste ---\n")
TIME_GRID = seq(0, 10, by = 0.1)
MEAN_FUNCS = list(X1=function(t)1+0.2*t, X2=function(t)0.5+0.3*sqrt(t),
                  X3=function(t)1.5*cos(0.8*t))
PROP_MEDIUM = function(X, ti) 1/(1+exp(-(-1+1.5*X[[1]][,ti-1]+1.2*X[[2]][,ti-1]-1*X[[3]][,ti-1])))
FEEDBACK_ALL = list(
  X1 = list(function(t,s)ifelse(t>s,1.5*exp(-0.5*(t-s)),0),
            function(t,s)ifelse(t>s,1*exp(-0.5*(t-s)),0), function(t,s)0),
  X2 = list(function(t,s)0, function(t,s)ifelse(t>s,1*exp(-0.4*(t-s)),0),
            function(t,s)ifelse(t>s,0.8*exp(-0.4*(t-s)),0)),
  X3 = list(function(t,s)0, function(t,s)0, function(t,s)ifelse(t>s,-1.2,0))
)
test_config = list(
  n = 1000, 
  treatment_times = c(1, 2, 3), 
  a_seq = c(1, 0, 1), 
  another_a_seq = c(0, 1, 0),
  propensity_function_list = list(PROP_MEDIUM, PROP_MEDIUM, PROP_MEDIUM),
  effects_covariates_from_A = FEEDBACK_ALL,
  mean_baseline_covariates = MEAN_FUNCS
)

cat("--- Gerando Dados ---\n")
baselines = create_intercept_slope_heterogeneity(
  n = test_config$n, mean_funcs = test_config$mean_baseline_covariates
)
generated_data = generate_confounder_pointwise(
  n = test_config$n,
  time_grid = TIME_GRID,
  noise_sd = 0.2,
  p_treat = 0.5,
  baseline_function_list = baselines,
  treatment_times = test_config$treatment_times,
  a_seq = test_config$a_seq,
  another_a_seq = test_config$another_a_seq,
  effects_function_list = test_config$effects_covariates_from_A,
  propensity_function_list = test_config$propensity_function_list,
  return_format = "wide"
)

regime_a_validar = test_config$a_seq

cat("--- Iniciando o Duelo de Modelos para X ---\n")
validation_results_x = run_x_model_validation(
  scenario_data = list(wide = generated_data$wide),
  scenario_config = test_config, 
  regime_to_validate = regime_a_validar
)

fit_gformula_models(generated_data$wide, test_config)

cat("\n--- Resultados do Duelo (RMSE) ---\n")
print(validation_results_x$metrics)

cat("\n--- Gerando Gráfico de Comparação ---\n")
print(plot_x_validation(validation_results_x))