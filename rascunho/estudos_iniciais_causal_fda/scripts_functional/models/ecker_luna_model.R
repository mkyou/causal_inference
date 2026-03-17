#' @title (Corrigida) Prepara Dados para o Modelo Sequencial de Ecker & Luna
#' @description Extrai covariáveis do tempo imediatamente anterior ao primeiro
#' tratamento definido na configuração do cenário.
#' @param scenario_data Um cenário da lista all_scenarios_data.
#' @param scenario_config A configuração correspondente.
#' @return Uma data.table pronta para modelagem.

prepare_data_sequential_eckerluna = function(scenario_data, scenario_config) {
  dt_long = as.data.table(scenario_data$outcome$long)
  treatment_times = scenario_config$treatment_times
  
  treatment_history_dt = dt_long[t %in% treatment_times, .(
    treatment_history = paste(A, collapse = "_")
  ), by = id]
  
  # --- LÓGICA CORRIGIDA E ROBUSTA ---
  first_treatment_time = min(treatment_times)
  all_times = sort(unique(dt_long$t))
  # Pega o último ponto de tempo ESTRITAMENTE ANTERIOR ao primeiro tratamento.
  baseline_time = max(all_times[all_times < first_treatment_time])
  
  cov_names = grep("^X[0-9]+$", names(dt_long), value = TRUE)
  base_cov_names = paste0(cov_names, "_base")
  
  baseline_data = dt_long[t == baseline_time, .SD, .SDcols = c("id", cov_names)]
  setnames(baseline_data, cov_names, base_cov_names)
  
  model_data = merge(dt_long, baseline_data, by = "id")
  model_data = merge(model_data, treatment_history_dt, by = "id")
  
  return(model_data)
}


#' @title (NOVA VERSÃO) Ajusta Modelos GAM para Regimes Específicos
#' @description Ajusta modelos GAM apenas para os dois regimes de tratamento
#' de interesse (ex: "sempre tratar" e "nunca tratar").
#' @param model_data Dados da função `prepare_data_sequential_eckerluna`.
#' @param regime_a_label O rótulo do primeiro regime (ex: "1_1").
#' @param regime_b_label O rótulo do segundo regime (ex: "0_0").
#' @return Uma lista nomeada com os dois modelos GAM ajustados.

fit_regime_specific_gam_eckerluna = function(
    model_data, 
    regime_a_label, 
    regime_b_label
  ) {
  require(mgcv)
  
  # Constrói a fórmula (mesma de antes)
  base_cov_names = grep("_base$", names(model_data), value = TRUE)
  smooth_terms = paste0("s(t, by = ", base_cov_names, ", k = 10)")
  model_formula = as.formula(paste("y ~ s(t, k = 15) +", paste(
    smooth_terms, 
    collapse = " + "
  )))
  
  # Subconjunto dos dados para os regimes de interesse
  data_a = model_data[treatment_history == regime_a_label]
  data_b = model_data[treatment_history == regime_b_label]
  
  if(nrow(data_a) == 0 | nrow(data_b) == 0) {
    stop("Dados insuficientes para um ou ambos os regimes de interesse.")
  }
  
  # --- Ajuste dos dois modelos específicos ---
  
  # Modelo A (usa dados B como "dummy")
  dummy_from_b = data_b[sample(.N, min(3, .N))]
  fit_a = bam(
    model_formula, 
    data = rbind(data_a, dummy_from_b), 
    discrete = TRUE, 
    select = TRUE
  )
  
  # Modelo B (usa dados A como "dummy")
  dummy_from_a = data_a[sample(.N, min(3, .N))]
  fit_b = bam(
    model_formula, 
    data = rbind(data_b, dummy_from_a), 
    discrete = TRUE, 
    select = TRUE
  )
  
  # Retorna uma lista nomeada com os modelos
  models_list = list(fit_a, fit_b)
  names(models_list) = c(regime_a_label, regime_b_label)
  
  return(models_list)
}

#' @title (Renomeada) Executa a Suíte de Estimação Sequencial de Ecker & Luna
#' @description Orquestra a estimação para cenários sequenciais, comparando dois regimes.
#' @param scenario_data Um único cenário da lista `all_scenarios_data`.
#' @param scenario_config A configuração correspondente.
#' @param regime_a O regime de tratamento alvo (ex: c(1, 1)).
#' @param regime_b O regime de tratamento de comparação (ex: c(0, 0)).
#' @return Uma lista com o gráfico, métricas e modelos.

run_sequential_eckerluna_suite = function(
    scenario_data, 
    scenario_config, 
    regime_a, 
    regime_b
  ) {
  
  regime_a_label = paste(regime_a, collapse = "_")
  regime_b_label = paste(regime_b, collapse = "_")
  
  # --- 1. Preparação dos Dados ---
  model_data = prepare_data_sequential_eckerluna(scenario_data, scenario_config)
  
  # --- 2. Ajuste dos Modelos Específicos ---
  fitted_models = fit_regime_specific_gam_eckerluna(
    model_data, 
    regime_a_label, 
    regime_b_label
  )
  
  # --- 3. Predição Contrafactual ---
  pred_Y_a = predict(fitted_models[[regime_a_label]], newdata = model_data)
  pred_Y_b = predict(fitted_models[[regime_b_label]], newdata = model_data)
  
  model_data[, effect_estimate := pred_Y_a - pred_Y_b]
  causal_effect = model_data[, .(estimated = mean(effect_estimate)), by = t]
  
  # --- 4. Coleta de Efeitos Reais e Métricas (sem mudanças na lógica) ---
  wide_data = scenario_data$outcome$wide
  real_total_dt = data.table(
    t = wide_data$t_vec, 
    real_total = colMeans(wide_data$effect_a_vs_another)
  )
  
  plot_data = merge(causal_effect, real_total_dt, by = "t")
  metrics = data.table(
    RMSE_vs_Total = sqrt(mean((plot_data$estimated - plot_data$real_total)^2)),
    Bias_vs_Total = mean(plot_data$estimated - plot_data$real_total),
    Cor_vs_Total = cor(plot_data$estimated, plot_data$real_total)
  )
  
  # --- 5. Geração do Gráfico (sem mudanças na lógica) ---
  subtitle_text = paste(
    "RMSE (vs Total):", round(metrics$RMSE_vs_Total, 3),
    "| Bias (vs Total):", round(metrics$Bias_vs_Total, 3),
    "| Corr (vs Total):", round(metrics$Cor_vs_Total, 3)
  )
  
  effect_plot = ggplot(plot_data, aes(x = t)) +
    geom_line(
      aes(y = real_total, color = "Real Total"), 
      linetype = "dashed", 
      linewidth = 1.1
    ) +
    geom_line(aes(y = estimated, color = "Estimado (GAM)"), linewidth = 1.2) +
    labs(
      title = "Efeito Causal Estimado (Proxy E&L) vs. Real",
      subtitle = subtitle_text,
      x = "Tempo", y = "Efeito Causal", color = ""
    ) +
    scale_color_manual(
      values = c("Estimado (GAM)" = "#E69F00", "Real Total" = "black")
    ) +
    theme_minimal(base_size = 14)
  
  return(list(
    metrics = metrics,
    plot = effect_plot,
    models = fitted_models,
    plot_data = plot_data
  ))
}
