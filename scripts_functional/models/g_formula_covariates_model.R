# ===================================================================
# G-FORMULA FUNCIONAL: MODELAGEM E VALIDAÇÃO DAS COVARIÁVEIS (VERSÃO 6.0)
#
# Abordagem com Ajuste para Confundidores Variantes no Tempo:
# - Adiciona os valores das covariáveis pré-tratamento como preditores
#   no modelo GAM para controlar o viés de seleção.
# - A função de simulação foi reescrita para ser sequencial, refletindo
#   a dependência temporal do novo modelo.
# ===================================================================

library(mgcv)
library(ggplot2)
library(data.table)

#' @title Prepara Dados para Modelagem da G-Formula (VERSÃO 8.1)
#' @description Corrige o bug de tipos incompatíveis no merge, garantindo
#'  que a coluna 'id' seja do tipo 'factor' em ambas as tabelas
#'  antes da junção.
prepare_gformula_data = function(dt_long, treatment_times) {
  dt_prepared = copy(dt_long)
  dt_prepared[, id := as.factor(id)]
  
  setkey(dt_prepared, id, t)
  
  all_cov_names = grep("^X[0-9]+$", names(dt_prepared), value = TRUE)
  
  lagged_cov_names = paste0(all_cov_names, "_lag1")
  
  dt_prepared[, (lagged_cov_names) := shift(.SD, n = 1, type = "lag"), by = id, .SDcols = all_cov_names]
  
  for (s in treatment_times) {
    time_idx = which.min(abs(unique(dt_prepared$t) - s))
    actual_time = unique(dt_prepared$t)[time_idx]
    
    # Cria a tabela de tratamento a partir do dt_long original
    treatment_at_s = dt_long[t == actual_time, .SD, .SDcols = c("id", "A")]
    
    # --- INÍCIO DA CORREÇÃO ---
    # Garante que o tipo da coluna de junção 'id' seja o mesmo (factor)
    treatment_at_s[, id := as.factor(id)]
    # --- FIM DA CORREÇÃO ---
    
    setnames(treatment_at_s, "A", paste0("A_s", s))
    dt_prepared = merge(dt_prepared, treatment_at_s, by = "id", all.x = TRUE)
  }
  
  # dt_prepared = na.omit(dt_prepared, cols = lagged_cov_names)
  
  return(dt_prepared)
}

#' @title Ajusta Modelos GAM para Covariáveis (VERSÃO 10.0 - Modelo Misto)
#' @description Aborda o confundimento por heterogeneidade individual.
#'  O modelo agora é um GAM Misto (GAMM), que inclui um smooth aleatório
#'  para cada indivíduo (s(t, id, bs="fs")). Isso modela explicitamente
#'  os interceptos e slopes aleatórios (a_ik, b_ik) do PGD.
fit_covariate_models = function(
    prepared_data,
    treatment_times,
    cov_to_model
) {
    fitted_models = list()
    all_cov_names = grep("^X[0-9]+$", names(prepared_data), value = TRUE)
    lagged_cov_names = paste0(all_cov_names, "_lag1")
    
    cat(paste0("Ajustando modelos GAM (V10.0 - GAMM com Efeitos Aleatórios) para: ",
             paste(cov_to_model, collapse=", "), "\n"))

    for (cov_name in cov_to_model) {
    
        # 1. Termo de base para a MÉDIA da população (mu_k(t))
        base_formula = 's(t, bs = "cr", k = 20)'
        
        treatment_effects_formula = paste(
            sapply(treatment_times, function(s) {
                paste0('s(t, bs = "tp", k = 10, by = I(A_s', s, ' * (t >= ', s, ')))')
            }), collapse = " + ")

        # 4. Termos de ajuste para o histórico lagado (sem alteração)
        history_terms = paste(
            sapply(lagged_cov_names, function(lag_name) {
                paste0('s(', lag_name, ', bs = "cr", k = 8)')
            }), collapse = " + ")

        # 5. Montagem da fórmula final
        model_formula_str = paste(
            cov_name, "~",
            base_formula, "+",
            treatment_effects_formula, "+",
            history_terms
        )
        
        cat("\n--- Fórmula para a covariável:", cov_name, "---\n")
        cat(model_formula_str, "\n")
        cat("---------------------------------------------\n\n")

        # É crucial que 'id' seja um fator para os efeitos aleatórios
        prepared_data[, id := as.factor(id)]
        
        fitted_models[[cov_name]] = bam(
            as.formula(model_formula_str), data = prepared_data, discrete = TRUE, select = TRUE, method = "fREML"
        )
    }
    cat("Ajuste de modelos GAM concluído.\n")
    return(fitted_models)
}


#' @title Simula Trajetórias Contrafactuais (VERSÃO 8.4 - CORREÇÃO DE BUG)
#' @description Corrige o bug na inicialização dos dados. O 'update join'
#'  agora é feito usando tanto 'id' quanto 't', garantindo que apenas o
#'  primeiro ponto de tempo de cada indivíduo seja inicializado com os
#'  valores observados.
simulate_counterfactual_trajectories = function(
    fitted_models,
    original_data,
    treatment_times,
    regime,
    simple_covs_data = NULL
) {
  # 1. PREPARAÇÃO (sem alterações)
  covs_gam = names(fitted_models)
  covs_simple = if (!is.null(simple_covs_data)) setdiff(names(simple_covs_data), "id") else character(0)
  all_covs = c(covs_gam, covs_simple)
  
  ids = unique(original_data$id)
  time_grid = sort(unique(original_data$t))
  
  simulated_dt = CJ(id = ids, t = time_grid)
  for (cov_name in all_covs) {
    simulated_dt[, (paste0(cov_name, "_sim")) := NA_real_]
  }
  
  # 2. PREENCHIMENTO DAS COVARIÁVEIS ESTÁTICAS (sem alterações)
  if (length(covs_simple) > 0) {
    simple_sim_names = paste0(covs_simple, "_sim")
    setnames(simple_covs_data, covs_simple, simple_sim_names)
    simulated_dt[simple_covs_data, on = .(id), (simple_sim_names) := mget(paste0("i.", simple_sim_names))]
  }
  
  # 3. SIMULAÇÃO DAS COVARIÁVEIS DINÂMICAS (SEQUENCIAL)
  
  # --- INÍCIO DA CORREÇÃO DO BUG ---
  # 3.1. INICIALIZAÇÃO (t = t_1) para covariáveis GAM
  if (length(covs_gam) > 0) {
    sim_names_gam = paste0(covs_gam, "_sim")
    
    # Pega os dados do primeiro ponto do tempo
    first_time_data = original_data[t == time_grid[1], .SD, .SDcols = c("id", covs_gam)]
    
    # Adiciona a coluna 't' para a junção correta
    first_time_data[, t := time_grid[1]]
    
    # Renomeia as colunas de dados para corresponder aos nomes em simulated_dt
    setnames(first_time_data, covs_gam, sim_names_gam)
    
    # Realiza o update join usando ID e T, atualizando apenas as linhas de t=0
    simulated_dt[first_time_data, on = .(id, t), (sim_names_gam) := mget(paste0("i.", sim_names_gam))]
  }
  # --- FIM DA CORREÇÃO DO BUG ---
  
  # 3.2. LOOP DE SIMULAÇÃO (t = t_2, ..., t_T) para covariáveis GAM (sem alterações na lógica)
  if (length(covs_gam) > 0) {
    for (i in 2:length(time_grid)) {
      current_t = time_grid[i]
      previous_t = time_grid[i-1]
      
      history_sim_names = paste0(all_covs, "_sim")
      lagged_names = paste0(all_covs, "_lag1")
      history_data = simulated_dt[t == previous_t, .SD, .SDcols = c("id", history_sim_names)]
      setnames(history_data, history_sim_names, lagged_names)
      
      pred_data_t = data.table(id = ids, t = current_t)
      for (j in seq_along(treatment_times)) {
        pred_data_t[, (paste0("A_s", treatment_times[j])) := regime[j]]
      }
      pred_data_t = merge(pred_data_t, history_data, by = "id")
      
      for (cov_name in covs_gam) {
        predictions = predict(fitted_models[[cov_name]], newdata = pred_data_t)
        simulated_dt[t == current_t, (paste0(cov_name, "_sim")) := predictions]
      }
    }
  }
  
  # 4. FINALIZAÇÃO
  simulated_dt[, id := as.numeric(as.character(id))]
  return(simulated_dt)
}


#' @title Executa a Suíte de Validação da G-Formula (VERSÃO 8.2 - HÍBRIDA)
#' @description Reintroduz a lógica híbrida para modelar covariáveis.
#'  Covariáveis dinâmicas ou com feedback usam GAMs (V8.1), enquanto
#'  covariáveis estáticas são modeladas pela média individual.
run_gformula_validation_suite = function(
    scenario_data, 
    scenario_config
) {
  dt_conf = scenario_data$confounders
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(dt_conf), value = TRUE)
  
  # --- LÓGICA DE DECISÃO HÍBRIDA (RESTAURADA) ---
  model_strategies = list()
  for (i in seq_along(cov_names)) {
    cov_name = cov_names[i]
    # Checa se a função de base da covariável é constante no tempo
    is_static = isTRUE(all.equal(
      scenario_config$mean_baseline_covariates[[i]](1),
      scenario_config$mean_baseline_covariates[[i]](10)
    ))
    # Checa se há algum efeito de feedback do tratamento na covariável
    has_feedback = !all(sapply(
      scenario_config$effects_covariates_from_A[[i]],
      function(f) f(5, 1) == 0
    ))
    
    # Decide a estratégia
    if (!is_static || has_feedback) {
      model_strategies[[cov_name]] = "gam"
    } else {
      model_strategies[[cov_name]] = "mean_trajectory"
    }
  }
  cat("Estratégias de modelagem definidas:\n")
  print(unlist(model_strategies))
  
  # --- AJUSTE DOS MODELOS GAM (APENAS PARA COVARIÁVEIS DINÂMICAS) ---
  covs_requiring_gam = names(model_strategies)[model_strategies == "gam"]
  if (length(covs_requiring_gam) > 0) {
    # Prepara os dados com lags e colunas de tratamento
    prepared_data = prepare_gformula_data(dt_conf, treatment_times)
    
    fitted_gam_models = fit_covariate_models(
      prepared_data, treatment_times, covs_requiring_gam
    )
  } else {
    fitted_gam_models = list()
  }
  
  # --- PRÉ-CÁLCULO PARA AS COVARIÁVEIS ESTÁTICAS ---
  covs_requiring_mean = names(model_strategies)[model_strategies == "mean_trajectory"]
  mean_trajectories_data = NULL
  if (length(covs_requiring_mean) > 0) {
    cat("\nCalculando trajetórias médias para covariáveis estáticas:", paste(covs_requiring_mean, collapse=", "), "\n")
    # Calcula a média POR INDIVÍDUO para as covariáveis estáticas
    mean_trajectories_data = dt_conf[
      ,
      lapply(.SD, mean), by = id, .SDcols = covs_requiring_mean
    ]
  }
  
  # --- SIMULAÇÃO HÍBRIDA ---
  cat("\nSimulando para o regime 'a_seq'...\n")
  simulated_a = simulate_counterfactual_trajectories(
    fitted_models = fitted_gam_models,
    original_data = dt_conf,
    treatment_times = treatment_times,
    regime = scenario_config$a_seq,
    simple_covs_data = mean_trajectories_data # Passando os dados das médias
  )
  
  cat("\nSimulando para o regime 'another_a_seq'...\n")
  simulated_another_a = simulate_counterfactual_trajectories(
    fitted_models = fitted_gam_models,
    original_data = dt_conf,
    treatment_times = treatment_times,
    regime = scenario_config$another_a_seq,
    simple_covs_data = mean_trajectories_data # Passando os dados das médias
  )
  
  
  # --- SEÇÃO DE VALIDAÇÃO ---
  
  # 1. Renomeia as colunas de predição ANTES de qualquer merge
  # Para simulated_a: X1_sim -> a_X1_pred
  sim_names_a = names(simulated_a)[grep("_sim$", names(simulated_a))]
  pred_names_a = paste0("a_", gsub("_sim$", "", sim_names_a), "_pred")
  setnames(simulated_a, sim_names_a, pred_names_a)
  
  # Para simulated_another_a: X1_sim -> another_a_X1_pred
  sim_names_another_a = names(simulated_another_a)[grep("_sim$", names(simulated_another_a))]
  pred_names_another_a = paste0("another_a_", gsub("_sim$", "", sim_names_another_a), "_pred")
  setnames(simulated_another_a, sim_names_another_a, pred_names_another_a)
  
  # 2. Faz o merge das duas tabelas de predição
  sim_data = merge(simulated_a, simulated_another_a, by = c("id", "t"))
  
  # 3. Faz o merge final com os dados verdadeiros
  validation_data = merge(dt_conf, sim_data, by = c("id", "t"))
  
  # --- FIM DA CORREÇÃO ---
  
  plot_data = validation_data[, lapply(.SD, mean, na.rm=TRUE), by = t]
  
  all_metrics = list()
  all_plots = list()
  
  for (cov_name in cov_names) {
    true_a_col = paste0("a_", cov_name)
    pred_a_col = paste0("a_", cov_name, "_pred")
    true_another_a_col = paste0("another_a_", cov_name)
    pred_another_a_col = paste0("another_a_", cov_name, "_pred")
    
    # Adicionando uma verificação para garantir que as colunas existem antes de calcular
    if (!all(c(true_a_col, pred_a_col, true_another_a_col, pred_another_a_col) %in% names(validation_data))) {
      warning(paste("Colunas de validação ausentes para a covariável:", cov_name))
      next # Pula para a próxima covariável
    }
    
    model_used_label = model_strategies[[cov_name]]
    if(model_used_label == "gam") {
      model_used_label = "GAM-simplified"
    }
    
    metrics = data.table(
      covariate = cov_name, model_used = model_used_label,
      RMSE_a = sqrt(mean( (validation_data[[true_a_col]] - validation_data[[pred_a_col]])^2, na.rm = TRUE )),
      Bias_a = mean( validation_data[[pred_a_col]] - validation_data[[true_a_col]], na.rm = TRUE ),
      Cor_a = tryCatch(cor( plot_data[[true_a_col]], plot_data[[pred_a_col]], use = "complete.obs" ), error = function(e) NA_real_),
      RMSE_another_a = sqrt(mean( (validation_data[[true_another_a_col]] - validation_data[[pred_another_a_col]])^2, na.rm = TRUE )),
      Bias_another_a = mean( validation_data[[pred_another_a_col]] - validation_data[[true_another_a_col]], na.rm = TRUE ),
      Cor_another_a = tryCatch(cor( plot_data[[true_another_a_col]], plot_data[[pred_another_a_col]], use = "complete.obs" ), error = function(e) NA_real_)
    )
    all_metrics[[cov_name]] = metrics
    
    plot_title = paste("Validação da G-Formula para", cov_name, "| Modelo:", model_used_label)
    subtitle_text = paste0("RMSE: ", round(metrics$RMSE_a, 3), " (Sempre) / ", round(metrics$RMSE_another_a, 3), " (Nunca) | ", "Corr: ", round(metrics$Cor_a, 3), " (Sempre) / ", round(metrics$Cor_another_a, 3), " (Nunca)")
    
    plot = ggplot(plot_data, aes(x = t)) +
      geom_line(aes(y = .data[[true_a_col]], color = "Verdadeiro (Sempre Tratar)"), linetype = "dashed", size = 1) +
      geom_line(aes(y = .data[[pred_a_col]], color = "Simulado (Sempre Tratar)"), size = 1) +
      geom_line(aes(y = .data[[true_another_a_col]], color = "Verdadeiro (Nunca Tratar)"), linetype = "dashed", size = 1) +
      geom_line(aes(y = .data[[pred_another_a_col]], color = "Simulado (Nunca Tratar)"), size = 1) +
      scale_color_manual(values = c("Verdadeiro (Sempre Tratar)"="blue", "Simulado (Sempre Tratar)"="lightblue", "Verdadeiro (Nunca Tratar)"="red", "Simulado (Nunca Tratar)"="pink")) +
      labs(title = plot_title, subtitle = subtitle_text, x = "Tempo", y = "Valor da Covariável", color = "Trajetória") +
      theme_minimal()
    
    all_plots[[cov_name]] = plot
  }
  
  return(
    list(
      metrics = rbindlist(all_metrics), 
      plots = all_plots,
      simulations = list(
        a_seq = sim_data[, .SD, .SDcols = c("id", "t", pred_names_a)],
        another_a_seq = sim_data[, .SD, .SDcols = c("id", "t", pred_names_another_a)]
      )
    )
  )
}


