# ===================================================================
# G-FORMULA COM MODELO AUTORREGRESSIVO (VERSÃO 6.0 - EFEITO SUAVIZADO)
#
# Abordagem que utiliza uma estratégia de modelagem adaptativa:
# - Para covariáveis dinâmicas ou com feedback, ajusta um modelo
#   híbrido com s(t), s(lag) e EFEITOS DE TRATAMENTO SUAVIZADOS.
# - Para covariáveis estáticas e sem feedback, utiliza um modelo
#   de média individual.
# ===================================================================

library(mgcv)
library(ggplot2)
library(data.table)

#' @title Prepara Dados para Modelagem Autorregressiva
#' @description Cria colunas lagadas para todas as covariáveis e as colunas
#'   de tratamento A_s. Remove linhas com NAs resultantes do lag.
#' @param dt_long O data.table no formato longo.
#' @param treatment_times Vetor numérico com os tempos de tratamento.
#' @return Um data.table modificado, pronto para a modelagem.
prepare_autoregressive_data = function(dt_long, treatment_times) {
  
  dt_prepared = copy(dt_long)
  cov_names = grep("^X[0-9]+$", names(dt_prepared), value = TRUE)
  lag_names = paste0(cov_names, "_lag1")
  
  # Cria as colunas lagadas
  dt_prepared[
    , 
    (lag_names) := lapply(.SD, shift, n=1, type="lag"), 
    .SDcols = cov_names, 
    by = id
  ]
  
  # Cria as colunas de tratamento
  for (s in treatment_times) {
    time_idx = which.min(abs(unique(dt_prepared$t) - s))
    actual_time = unique(dt_prepared$t)[time_idx]
    treatment_at_s = dt_prepared[
      t == actual_time, 
      .SD, 
      .SDcols = c("id", "A")
    ]
    A_var_name = paste0("A_s", s)
    setnames(treatment_at_s, "A", A_var_name)
    dt_prepared = merge(dt_prepared, treatment_at_s, by = "id")
  }
  
  # Remove a primeira observação de cada indivíduo (onde o lag é NA)
  dt_prepared = na.omit(dt_prepared)
  return(dt_prepared)
}

#' @title Ajusta Modelos GAM Autorregressivos Híbridos
#' @description Ajusta um modelo para um subconjunto específico de 
#'   covariáveis com efeitos de tratamento suavizados.
#' @param prepared_data O data.table da função `prepare_autoregressive_data`.
#' @param target_covariates Vetor de nomes das covariáveis a serem modeladas.
#' @return Uma lista nomeada contendo um modelo 'bam' para cada covariável.
fit_autoregressive_models = function(
    prepared_data, 
    treatment_times,
    target_covariates
) {
  
  fitted_models = list()
  
  cat(
    "Ajustando modelos autorregressivos (efeito suavizado) para:",
    paste(target_covariates, collapse = ", "), "\n"
  )
  
  for (cov_name in target_covariates) {
    own_lag = paste0(cov_name, "_lag1")
    
    # Constrói a parte da fórmula para os efeitos suavizados do tratamento
    treatment_effects_formula = paste(
      sapply(treatment_times, function(s) {
        paste0("s(t, by = A_s", s, ")")
      }),
      collapse = " + "
    )
    
    # Fórmula Híbrida: Xk ~ s(t) + s(Xk_lag1) + s(t, by=A_s...)
    formula_str = paste(
      cov_name, "~",
      "s(t) +",
      paste0("s(", own_lag, ")"),
      "+",
      treatment_effects_formula
    )
    
    fitted_models[[cov_name]] = bam(
      as.formula(formula_str),
      data = prepared_data,
      discrete = TRUE,
      select = TRUE
    )
  }
  
  cat("Ajuste de modelos concluído.\n")
  return(fitted_models)
}

#' @title Simula Trajetórias Sequencialmente (Abordagem Híbrida)
#' @description Usa os modelos GAM para covariáveis dinâmicas e a média
#'   individual para covariáveis estáticas.
#' @return Um data.table com as trajetórias simuladas.
simulate_autoregressive_trajectories = function(
    fitted_models,
    mean_covariates,
    original_data, 
    treatment_times, 
    regime
) {
  
  all_times = sort(unique(original_data$t))
  gam_covariates = names(fitted_models)
  all_cov_names = c(gam_covariates, mean_covariates)
  lag_names = paste0(all_cov_names, "_lag1")
  treatment_vars = paste0("A_s", treatment_times)
  
  # Calcula as médias para as covariáveis estáticas
  mean_trajectories = NULL
  if (length(mean_covariates) > 0) {
    mean_trajectories = original_data[
      , 
      lapply(.SD, mean), 
      by = id, 
      .SDcols = mean_covariates
    ]
  }
  
  # Inicializa o data.table com os valores de linha de base (t=0)
  sim_dt = original_data[
    t == all_times[1], 
    .SD, 
    .SDcols = c("id", "t", all_cov_names)
  ]
  
  cat("  Iniciando simulação sequencial híbrida...\n")
  
  # Loop sequencial sobre os pontos de tempo futuros
  for (i in 2:length(all_times)) {
    current_t = all_times[i]
    prev_t = all_times[i-1]
    
    last_sim_values = sim_dt[t == prev_t, .SD, .SDcols = c("id", all_cov_names)]
    setnames(last_sim_values, all_cov_names, lag_names)
    
    pred_data_step = data.table(id = last_sim_values$id, t = current_t)
    pred_data_step = cbind(pred_data_step, last_sim_values[, ..lag_names])
    
    for (j in seq_along(treatment_times)) {
      pred_data_step[, (treatment_vars[j]) := regime[j]]
    }
    
    current_sim_values = data.table(
      id = unique(original_data$id), 
      t = current_t
    )
    
    if (length(gam_covariates) > 0) {
      for (cov_name in gam_covariates) {
        predictions = predict(
          fitted_models[[cov_name]], 
          newdata = pred_data_step
        )
        current_sim_values[, (cov_name) := predictions]
      }
    }
    
    if (length(mean_covariates) > 0) {
      mean_vals = mean_trajectories[, .SD, .SDcols = c("id", mean_covariates)]
      if (length(gam_covariates) == 0) {
        current_sim_values = mean_vals
        current_sim_values[, t := current_t]
      } else {
        current_sim_values = merge(current_sim_values, mean_vals, by = "id")
      }
    }
    
    sim_dt = rbind(sim_dt, current_sim_values)
  }
  
  cat("  Simulação sequencial concluída.\n")
  return(sim_dt)
}


#' @title Executa a Suíte de Validação para o Modelo Autorregressivo
#' @description Orquestra o ajuste e simulação, escolhendo o modelo
#'   (GAM vs. Média) para cada covariável com base em sua natureza.
#' @return Uma lista contendo métricas e gráficos.
run_autoregressive_validation_suite = function(
    scenario_data, 
    scenario_config
) {
  
  dt_conf = scenario_data$confounders
  treatment_times = scenario_config$treatment_times
  cov_names = grep("^X[0-9]+$", names(dt_conf), value = TRUE)
  
  # --- LÓGICA DE DETECÇÃO POR COVARIÁVEL ---
  gam_covariates = c()
  mean_covariates = c()
  
  for (i in seq_along(cov_names)) {
    cov_name = cov_names[i]
    is_static = isTRUE(all.equal(
      scenario_config$mean_baseline_covariates[[i]](1),
      scenario_config$mean_baseline_covariates[[i]](10)
    ))
    
    has_feedback = !isTRUE(all(sapply(
      scenario_config$effects_covariates_from_A[[i]],
      function(f) f(5, 1) == 0
    )))
    
    if (!is_static || has_feedback) {
      gam_covariates = c(gam_covariates, cov_name)
    } else {
      mean_covariates = c(mean_covariates, cov_name)
    }
  }
  
  # --- AJUSTE DOS MODELOS ---
  prepared_data = prepare_autoregressive_data(dt_conf, treatment_times)
  covariate_models = NULL
  if (length(gam_covariates) > 0) {
    covariate_models = fit_autoregressive_models(
      prepared_data, 
      treatment_times,
      target_covariates = gam_covariates
    )
  }
  
  # --- SIMULAÇÃO ---
  cat("\nSimulando para o regime 'a_seq'...\n")
  simulated_a = simulate_autoregressive_trajectories(
    covariate_models, mean_covariates, dt_conf, treatment_times, 
    scenario_config$a_seq
  )
  
  cat("\nSimulando para o regime 'another_a_seq'...\n")
  simulated_another_a = simulate_autoregressive_trajectories(
    covariate_models, mean_covariates, dt_conf, treatment_times, 
    scenario_config$another_a_seq
  )
  
  # --- VALIDAÇÃO ---
  all_cov_names = c(gam_covariates, mean_covariates)
  setnames(
    simulated_a, all_cov_names, paste0("a_", all_cov_names, "_pred")
  )
  setnames(
    simulated_another_a, 
    all_cov_names, 
    paste0("another_a_", all_cov_names, "_pred")
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
    
    model_used = if (cov_name %in% gam_covariates) {
      "Autoregressive (Smooth Treat.)"
    } else {
      "Mean Trajectory"
    }
    
    metrics = data.table(
      covariate = cov_name, model_used = model_used,
      RMSE_a = sqrt(
        mean((validation_data[[true_a_col]]-validation_data[[pred_a_col]])^2)
      ),
      Bias_a = mean(
        validation_data[[pred_a_col]] - validation_data[[true_a_col]]
      ),
      Cor_a = cor(
        plot_data[[true_a_col]], 
        plot_data[[pred_a_col]]
      ),
      RMSE_another_a = sqrt(
        mean((validation_data[[true_another_a_col]] - 
                validation_data[[pred_another_a_col]])^2)
      ),
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
      "Validação para", cov_name, "| Modelo:", model_used
    )
    subtitle_text = paste0(
      "RMSE: ", round(metrics$RMSE_a, 3), " (Sempre) / ", 
      round(metrics$RMSE_another_a, 3), " (Nunca) | ",
      "Corr: ", round(metrics$Cor_a, 3), " (Sempre) / ",
      round(metrics$Cor_another_a, 3), " (Nunca)"
    )
    
    plot = ggplot(plot_data, aes(x=t)) +
      geom_line(
        aes(y = .data[[true_a_col]], color = "Verdadeiro (Sempre Tratar)"),
        linetype = "dashed", size = 1
      ) +
      geom_line(
        aes(y = .data[[pred_a_col]], color = "Simulado (Sempre Tratar)"),
        size = 1
      ) +
      geom_line(
        aes(
          y = .data[[true_another_a_col]], 
          color = "Verdadeiro (Nunca Tratar)"
        ),
        linetype = "dashed", size = 1
      ) +
      geom_line(
        aes(
          y = .data[[pred_another_a_col]], 
          color = "Simulado (Nunca Tratar)"
        ),
        size = 1
      ) +
      scale_color_manual(
        values = c(
          "Verdadeiro (Sempre Tratar)" = "blue", 
          "Simulado (Sempre Tratar)" = "lightblue", 
          "Verdadeiro (Nunca Tratar)" = "red", 
          "Simulado (Nunca Tratar)" = "pink"
        )
      ) +
      labs(
        title = plot_title, 
        subtitle = subtitle_text,
        x = "Tempo", 
        y = "Valor", 
        color = "Trajetória"
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

