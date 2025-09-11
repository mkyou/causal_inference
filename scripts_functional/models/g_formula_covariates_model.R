# ===================================================================
# G-FORMULA FUNCIONAL: MODELAGEM E VALIDAÇÃO DAS COVARIÁVEIS (VERSÃO 5.2)
#
# Abordagem Híbrida com Detecção Granular por Covariável:
# - Foco exclusivo no modelo 'simplified' para maior velocidade.
# - Lógica de deteção de feedback e tipos de dados corrigidos para robustez.
# - Código formatado para respeitar o limite de 75 caracteres por linha.
# ===================================================================

library(mgcv)
library(ggplot2)
library(data.table)

#' @title Prepara Dados para Modelagem da G-Formula Funcional
#' @description Garante que 'id' seja um fator e adiciona as colunas de
#'   tratamento A_s para uso posterior.
#' @param dt_long O data.table no formato longo.
#' @param treatment_times Vetor numérico com os tempos de tratamento.
#' @return Um data.table modificado, pronto para a modelagem.
prepare_gformula_data = function(dt_long, treatment_times) {
  
  dt_prepared = copy(dt_long)
  dt_prepared[, id := as.factor(id)]
  
  for (s in treatment_times) {
    time_idx = which.min(abs(unique(dt_prepared$t) - s))
    actual_time = unique(dt_prepared$t)[time_idx]
    
    treatment_at_s = dt_prepared[
      t == actual_time, .SD, .SDcols = c("id", "A")
    ]
    A_var_name = paste0("A_s", s)
    setnames(treatment_at_s, "A", A_var_name)
    
    dt_prepared = merge(dt_prepared, treatment_at_s, by = "id")
  }
  return(dt_prepared)
}

#' @title Ajusta Modelos GAM para Covariáveis (VERSÃO DEFINITIVA E ROBUSTA)
#' @description Usa bases de suavização distintas para a tendência de base (bs="cr")
#'   e os efeitos do tratamento (bs="tp"), garantindo estabilidade numérica e
#'   um ajuste preciso para ambos os regimes contrafactuais.
fit_covariate_models = function(
    prepared_data, 
    treatment_times, 
    cov_to_model
) {
  
  fitted_models = list()
  
  cat(paste0(
    "Ajustando modelos GAM (com bases distintas) para: ", 
    paste(cov_to_model, collapse=", "), "\n"
  ))
  
  for (cov_name in cov_to_model) {
    
    # --- INÍCIO DA MUDANÇA PRINCIPAL ---
    
    # Define a fórmula da linha de base com uma base de spline cúbica,
    # que é ótima para tendências gerais. Aumentamos 'k' para mais flexibilidade.
    base_formula = 's(t, bs = "cr", k = 30)'
    
    # Define os efeitos do tratamento com uma base de thin plate spline,
    # que é um bom padrão para desvios. 'k' pode ser menor aqui.
    treatment_effects_formula = paste(
      sapply(treatment_times, function(s) {
        paste0('s(t, bs = "tp", k = 10, by = I(A_s', s, ' * (t >= ', s, ')))')
      }),
      collapse = " + "
    )
    
    # A fórmula combina a linha de base estável com os desvios bem definidos.
    model_formula_str = paste(
      cov_name, "~", base_formula, "+", treatment_effects_formula
    )
    
    # --- FIM DA MUDANÇA PRINCIPAL ---
    
    model_formula = as.formula(model_formula_str)
    
    fitted_models[[cov_name]] = bam(
      model_formula,
      data = prepared_data,
      discrete = TRUE,
      # 'select = TRUE' continua sendo importante para regularizar os termos
      select = TRUE,
      method = "fREML"
    )
  }
  
  cat("Ajuste de modelos GAM concluído.\n")
  return(fitted_models)
}


#' @title Simula Trajetórias Contrafactuais (VERSÃO ATUALIZADA)
simulate_counterfactual_trajectories = function(
    fitted_models, 
    original_data,
    treatment_times, 
    regime
) {
  
  original_ids = unique(original_data$id)
  pred_data = CJ(id = original_ids, t = unique(original_data$t))
  
  # Adiciona as colunas do regime para construir a história
  A_vars = paste0("A_s", treatment_times)
  for (i in seq_along(treatment_times)) {
    pred_data[, (A_vars[i]) := regime[i]]
  }
  
  # --- CRIA A MESMA COLUNA 'history' NOS DADOS DE PREDIÇÃO ---
  pred_data[, history := "H_Pre"]
  sorted_times = sort(treatment_times)
  for (i in seq_along(sorted_times)) {
    s_current = sorted_times[i]
    s_next = if (i < length(sorted_times)) sorted_times[i+1] else Inf
    
    history_values_for_regime = paste("H", paste(regime[1:i], collapse="_"), sep="_")
    
    pred_data[t >= s_current & t < s_next, history := history_values_for_regime]
  }
  
  pred_data[, history := as.factor(history)]
  
  # Simulação
  simulated_dt = pred_data[, .(id, t)]
  
  for (cov_name in names(fitted_models)) {
    predictions = predict(
      fitted_models[[cov_name]],
      newdata = pred_data
    )
    simulated_dt[, (paste0(cov_name, "_sim")) := predictions]
  }
  
  simulated_dt[, id := as.numeric(as.character(id))]
  return(simulated_dt)
}

#' @title Executa a Suíte de Validação da G-Formula (Detecção Granular)
run_gformula_validation_suite = function(
    scenario_data, 
    scenario_config, 
    static_model_type = "mean_trajectory"
) {
  
  dt_conf = scenario_data$confounders
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(dt_conf), value = TRUE)
  
  # --- LÓGICA DE DECISÃO HÍBRIDA POR COVARIÁVEL ---
  model_strategies = list()
  for (i in seq_along(cov_names)) {
    cov_name = cov_names[i]
    
    is_static = isTRUE(all.equal(
      scenario_config$mean_baseline_covariates[[i]](1),
      scenario_config$mean_baseline_covariates[[i]](10)
    ))
    
    # CORRIGIDO: Usa '==' para uma verificação lógica robusta
    has_feedback = !all(sapply(
      scenario_config$effects_covariates_from_A[[i]],
      function(f) f(5, 1) == 0
    ))
    
    if (!is_static || has_feedback) {
      model_strategies[[cov_name]] = "gam"
    } else {
      model_strategies[[cov_name]] = static_model_type
    }
  }
  
  # --- AJUSTE DOS MODELOS NECESSÁRIOS ---
  covs_requiring_gam = names(model_strategies)[model_strategies == "gam"]
  
  if (length(covs_requiring_gam) > 0) {
    prepared_data = prepare_gformula_data(dt_conf, treatment_times)
    fitted_gam_models = fit_covariate_models(
      prepared_data, treatment_times, covs_requiring_gam
    )
  } else {
    fitted_gam_models = list()
  }
  
  # --- SIMULAÇÃO HÍBRIDA ---
  simulate_regime = function(regime) {
    full_simulation = CJ(id = unique(dt_conf$id), t = unique(dt_conf$t))
    
    # Simula covariáveis que usam GAM
    if (length(covs_requiring_gam) > 0) {
      gam_sims = simulate_counterfactual_trajectories(
        fitted_gam_models, dt_conf, treatment_times, regime
      )
      setnames(
        gam_sims, 
        paste0(covs_requiring_gam, "_sim"), 
        covs_requiring_gam
      )
      full_simulation = merge(full_simulation, gam_sims, by = c("id", "t"))
    }
    
    # Simula covariáveis que usam média
    covs_requiring_mean = names(model_strategies)[model_strategies != "gam"]
    if (length(covs_requiring_mean) > 0) {
      mean_trajectories = dt_conf[, 
                                  lapply(.SD, mean), by = id, .SDcols = covs_requiring_mean
      ]
      mean_sims = merge(
        CJ(id = unique(dt_conf$id), t = unique(dt_conf$t)),
        mean_trajectories, by = "id"
      )
      full_simulation = merge(full_simulation, mean_sims, by = c("id", "t"))
    }
    
    return(full_simulation[, .SD, .SDcols = c("id", "t", cov_names)])
  }
  
  cat("\nSimulando para o regime 'a_seq'...\n")
  simulated_a = simulate_regime(scenario_config$a_seq)
  
  cat("\nSimulando para o regime 'another_a_seq'...\n")
  simulated_another_a = simulate_regime(scenario_config$another_a_seq)
  
  # --- VALIDAÇÃO ---
  setnames(simulated_a, cov_names, paste0("a_", cov_names, "_pred"))
  setnames(
    simulated_another_a, 
    cov_names, 
    paste0("another_a_", cov_names, "_pred")
  )
  
  sim_data = merge(simulated_a, simulated_another_a, by = c("id", "t"))
  validation_data = merge(dt_conf, sim_data, by = c("id", "t"))
  plot_data = validation_data[, lapply(.SD, mean), by = t]
  
  all_metrics = list()
  all_plots = list()
  
  for (cov_name in cov_names) {
    true_a_col = paste0("a_", cov_name)
    pred_a_col = paste0("a_", cov_name, "_pred")
    true_another_a_col = paste0("another_a_", cov_name)
    pred_another_a_col = paste0("another_a_", cov_name, "_pred")
    
    model_used_label = model_strategies[[cov_name]]
    if(model_used_label == "gam") {
      model_used_label = "GAM-simplified"
    }
    
    metrics = data.table(
      covariate = cov_name, model_used = model_used_label,
      RMSE_a = sqrt(mean(
        (validation_data[[true_a_col]]-validation_data[[pred_a_col]])^2
      )),
      Bias_a = mean(
        validation_data[[pred_a_col]] - validation_data[[true_a_col]]
      ),
      Cor_a = cor(
        plot_data[[true_a_col]], 
        plot_data[[pred_a_col]]
      ),
      RMSE_another_a = sqrt(mean(
        (validation_data[[true_another_a_col]] - 
           validation_data[[pred_another_a_col]])^2
      )),
      Bias_another_a = mean(
        validation_data[[pred_another_a_col]] - 
          validation_data[[true_another_a_col]]
      ),
      Cor_another_a = cor(
        plot_data[[true_another_a_col]], 
        plot_data[[pred_another_a_col]]
      )
    )
    all_metrics[[cov_name]] = metrics
    
    plot_title = paste(
      "Validação da G-Formula para", cov_name, "| Modelo:", model_used_label
    )
    subtitle_text = paste0(
      "RMSE: ", round(metrics$RMSE_a, 3), " (Sempre) / ", 
      round(metrics$RMSE_another_a, 3), " (Nunca) | ",
      "Corr: ", round(metrics$Cor_a, 3), " (Sempre) / ",
      round(metrics$Cor_another_a, 3), " (Nunca)"
    )
    
    plot = ggplot(plot_data, aes(x = t)) +
      geom_line(
        aes(y = .data[[true_a_col]], color = "Verdadeiro (Sempre Tratar)"), 
        linetype = "dashed", size = 1
      ) +
      geom_line(
        aes(y = .data[[pred_a_col]], color = "Simulado (Sempre Tratar)"), 
        size = 1
      ) +
      geom_line(
        aes(y = .data[[true_another_a_col]], color = "Verdadeiro (Nunca Tratar)"), 
        linetype = "dashed", size = 1
      ) +
      geom_line(
        aes(y = .data[[pred_another_a_col]], color = "Simulado (Nunca Tratar)"),
        size = 1
      ) +
      scale_color_manual(
        values = c(
          "Verdadeiro (Sempre Tratar)"="blue", 
          "Simulado (Sempre Tratar)"="lightblue", 
          "Verdadeiro (Nunca Tratar)"="red", 
          "Simulado (Nunca Tratar)"="pink"
        )
      ) +
      labs(
        title = plot_title, 
        subtitle = subtitle_text, 
        x = "Tempo", y = "Valor da Covariável", color = "Trajetória"
      ) +
      theme_minimal()
    
    all_plots[[cov_name]] = plot
  }
  
  return(
    list(
      metrics = rbindlist(all_metrics), 
      plots = all_plots,
      simulations = list(
        a_seq = simulated_a,
        another_a_seq = simulated_another_a
      )
    )
  )
}


