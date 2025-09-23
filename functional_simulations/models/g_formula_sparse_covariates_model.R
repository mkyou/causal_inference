library(mgcv)
library(data.table)

#' @title Ajusta modelos de covariáveis esparsos baseados no PGD
fit_covariate_models_sparse = function(
    fit_data,
    treatment_times,
    scenario_config
) {
  fitted_models = list()
  cov_names = grep("^X[0-9]+$", names(fit_data), value = TRUE)
  
  beta_matrix = scenario_config$beta_funcs_matrix
  gamma_matrix = scenario_config$gamma_matrix
  
  for (k_idx in 1:length(cov_names)) {
    k_name = cov_names[k_idx]
    
    # Constrói a fórmula do Histórico dinamicamente
    history_terms = c()
    for (j_idx in 1:length(cov_names)) {
      if (k_idx != j_idx && !is.null(beta_matrix[[k_idx, j_idx]])) {
        lag_name = paste0(cov_names[j_idx], "_lag1")
        term = paste0(
          "s(", lag_name, ", k = 15) + ti(t, ", lag_name, ", k = c(12, 12))"
        )
        history_terms = c(history_terms, term)
      }
    }
    history_formula = if (length(history_terms) > 0) {
      paste(history_terms, collapse = " + ")
    } else { "1" }
    
    # Constrói a fórmula do Tratamento dinamicamente
    treatment_terms = c()
    for (s_idx in 1:length(treatment_times)) {
      s_time = treatment_times[s_idx]
      
      # --- INÍCIO DA CORREÇÃO ---
      # A verificação agora inspeciona o corpo da função, em vez de avaliá-la
      # em um único ponto. Isso diferencia corretamente as funções nulas.
      func_body_is_zero = isTRUE(
        deparse(body(gamma_matrix[[k_idx, s_idx]])) == "0"
      )
      
      if (!is.null(gamma_matrix[[k_idx, s_idx]]) && !func_body_is_zero) {
        # --- FIM DA CORREÇÃO ---
        term = paste0(
          "s(time_since_s", s_time, ", by = by_s", s_time, ", k = 15)"
        )
        treatment_terms = c(treatment_terms, term)
      }
    }
    treatment_formula = if (length(treatment_terms) > 0) {
      paste(treatment_terms, collapse = " + ")
    } else { "1" }
    
    base_formula = "s(t, k = 30)"
    final_formula = paste(
      k_name, "~", base_formula, "+", treatment_formula, "+", history_formula
    )
    final_formula = gsub("\\+ 1 \\+ 1", "+ 1", final_formula)
    final_formula = gsub("\\+ 1$", "", final_formula)
    
    cat(
      "Fórmula de estimação do modelo das covariáveis: ", 
      final_formula,
      "\n"
    )
    
    model = bam(
      as.formula(final_formula), data = fit_data, discrete = TRUE, select = TRUE
    )
    
    fitted_models[[k_name]] = list(
      model = model, residual_sd = summary(model)$scale
    )
  }
  return(fitted_models)
}


#' @title Roda a Simulação de Covariáveis com o modelo esparso (Corrigido)
run_gformula_covariates_sparse = function(
    prepared_data,
    scenario_config,
    regime_a,
    regime_b,
    add_noise = TRUE,
    chunk_size = 1500
) {
  
  # A ESTRUTURA DESTA FUNÇÃO AGORA É IDÊNTICA À SUA VERSÃO FUNCIONAL
  
  cat("--- Iniciando G-Formula (Etapa de Covariáveis ESPARSA) ---\n")
  
  fit_data = na.omit(prepared_data)
  fit_data[, id := as.factor(id)]
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
  
  cat("  - Ajustando modelos de covariáveis (ESPARSOS)...\n")
  # MUDANÇA: Chama a nova função de ajuste esparsa
  fitted_models_list = fit_covariate_models_sparse(
    fit_data,
    treatment_times,
    scenario_config
  )
  
  simulate_one_regime = function(regime_to_sim, regime_name) {
    cat(
      paste0(
        "  -- Iniciando simulação para o Regime: ", regime_name, " --\n"
      )
    )
    
    ids = unique(fit_data$id)
    n_subjects = length(ids)
    time_grid = sort(unique(fit_data$t))
    
    cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
    
    first_treatment_time = min(treatment_times)
    
    # Inicializa a simulação com uma cópia completa dos dados factuais
    sim_dt = copy(
      prepared_data[
        , .SD, 
        .SDcols = c(
          "id", "t", 
          cov_names, grep("^time_since_s", names(prepared_data), value = TRUE)
        )
      ]
    )
    
    # Adiciona TODAS as variáveis necessárias para a predição ao sim_dt
    
    # 1. Adiciona colunas de tratamento do regime
    for (j in 1:length(treatment_times)) {
      s_time = treatment_times[j]
      sim_dt[, paste0("A_s", s_time) := regime_to_sim[j]]
    }
    
    # 2. Adiciona colunas de tempo e de 'by' correspondentes
    for (s_time in treatment_times) {
      # Variável de tempo desde o tratamento
      time_since_col = paste0("time_since_s", s_time)
      sim_dt[, (time_since_col) := fifelse(t >= s_time, t - s_time, 0)]
      
      # Variável 'by' para o spline, baseada no TRATAMENTO DO REGIME
      by_col = paste0("by_s", s_time)
      A_s_col = paste0("A_s", s_time)
      sim_dt[
        , 
        (by_col) := .SD[[A_s_col]] * fifelse(t >= s_time, 1, 0), 
        .SDcols = A_s_col
      ]
    }
    
    # Encontra o índice de tempo para começar a simulação iterativa
    start_sim_idx = which(time_grid == first_treatment_time)
    
    # Loop de simulação principal
    for (i in start_sim_idx:length(time_grid)) {
      current_t = time_grid[i]
      previous_t = time_grid[i-1]
      
      history = sim_dt[t == previous_t, .SD, .SDcols = c("id", cov_names)]
      setnames(history, cov_names, paste0(cov_names, "_lag1"))
      
      pred_data_t = sim_dt[t == current_t]
      pred_data_t = merge(pred_data_t, history, by = "id")
      
      # Prediz e atualiza os valores para o tempo atual
      for (k in cov_names) {
        model_info = fitted_models_list[[k]]
        all_predictions = numeric(n_subjects)
        
        id_chunks = split(
          1:n_subjects,
          ceiling((1:n_subjects) / chunk_size)
        )
        
        for (chunk_indices in id_chunks) {
          
          # Prediz a média para o bloco atual
          chunk_predictions = predict(
            model_info$model,
            newdata = pred_data_t[chunk_indices, ]
          )
          
          # Adiciona ruído opcionalmente
          if (add_noise) {
            # MUDANÇA: Usa rnorm com o desvio padrão armazenado
            noise = rnorm(
              length(chunk_predictions),
              mean = 0,
              sd = model_info$residual_sd
            )
            chunk_predictions = chunk_predictions + noise
          }
          
          all_predictions[chunk_indices] = chunk_predictions
        }
        
        # Atualiza o data.table principal com todas as predições
        sim_dt[t == current_t, (k) := all_predictions]
      }
    }
    
    return(sim_dt[, .SD, .SDcols = c("id", "t", cov_names)])
  }
  
  sim_a = simulate_one_regime(regime_a, "Regime A")
  sim_b = simulate_one_regime(regime_b, "Regime B")
  
  setnames(sim_a, cov_names, paste0("sim_a_", cov_names))
  setnames(sim_b, cov_names, paste0("sim_b_", cov_names))
  simulation_results = merge(sim_a, sim_b, by = c("id", "t"))
  
  cat("--- Estágio de Covariáveis (Esparso) Concluído ---\n")
  return(list(
    models = lapply(fitted_models_list, `[[`, "model"),
    simulated_data = simulation_results
  ))
}