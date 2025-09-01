#' @title (VERSÃO FINAL) Ajusta Modelos para a G-Formula Funcional
#' @description Ajusta os modelos estruturais (GAMM para X, GAM para Y) usando
#' um modelo de intercepto e inclinação aleatórios para as covariáveis.
#' @param observed_data A data.table com os dados observados de um cenário.
#' @param scenario_config A configuração do cenário para obter os tempos de tratamento.
#' @return Uma lista nomeada com todos os modelos ajustados.

#' @title (VERSÃO FINAL E CORRETA) Ajusta Modelos para a G-Formula Funcional
#' @description Usa uma especificação de GAM de coeficientes variantes (s(t, by=...))
#' para estimar os efeitos funcionais dos tratamentos.

fit_gformula_models = function(observed_data, scenario_config) {
  require(mgcv)
  
  model_data = as.data.table(observed_data)
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  treatment_times = scenario_config$treatment_times
  
  # Prepara as colunas A_s, representando a história de tratamento
  if (!is.null(treatment_times)) {
    for (s in treatment_times) {
      status = model_data[t == s, .(id, A_s = A)]
      setnames(status, "A_s", paste0("A_", s))
      model_data = merge(model_data, status, by = "id", all.x = TRUE)
      model_data[is.na(get(paste0("A_", s))), (paste0("A_", s)) := 0]
    }
  }
  
  # --- 1. Ajustar um modelo GAM para cada covariável ---
  confounder_models = list()
  
  # Termos de efeito funcional do tratamento (feedback), REGULARIZADOS com k=6
  treatment_effect_terms = if (!is.null(treatment_times)) {
    paste0("s(t, by = A_", treatment_times, ", k = 6)")
  } else { "1" }
  
  for (cov in cov_names) {
    cat(paste("Ajustando modelo para", cov, "\n"))
    # Omitimos s(id) e lags de X em favor de um modelo puramente funcional e robusto
    formula_x = as.formula(
      paste(cov, "~ s(t) +", paste(treatment_effect_terms, collapse=" + "))
    )
    confounder_models[[cov]] = bam(formula_x, data = model_data,
                                   discrete = TRUE, select = TRUE)
  }
  
  # --- 2. Ajustar um modelo GAM para o desfecho Y ---
  cat("Ajustando modelo para o desfecho Y\n")
  
  # SUA SUGESTÃO: Usando produtos tensores para modelar beta(t)X(t)
  cov_effect_terms = paste0("te(t, ", cov_names, ", k = 10)")
  
  formula_y = as.formula(
    paste("y ~ s(t) +",
          paste(cov_effect_terms, collapse = " + "), "+",
          paste(treatment_effect_terms, collapse = " + "))
  )
  outcome_model = bam(formula_y, data = model_data, discrete = TRUE, select = TRUE)
  
  return(list(
    confounder_models = confounder_models,
    outcome_model = outcome_model,
    noise_sds = list(
      y = sd(outcome_model$residuals),
      X = sapply(confounder_models, function(m) sd(m$residuals))
    )
  ))
}


#' @title (VERSÃO FINAL E DEFINITIVA) Simula Intervenção Contrafactual
simulate_gformula_intervention = function(gformula_models, observed_data,
                                          scenario_config, intervention_regime) {
  all_times = sort(unique(observed_data$t))
  ids = unique(observed_data$id)
  n_ids = length(ids)
  cov_names = names(gformula_models$confounder_models)
  treatment_times = scenario_config$treatment_times
  
  sim_data = observed_data[t == min(all_times),
                           .SD, .SDcols = c("id", "t", "A", cov_names, "y")]
  
  # Guarda a história de tratamento COMPLETA do regime
  treatment_history = data.table(id = ids)
  if (!is.null(treatment_times)) {
    for (i in seq_along(treatment_times)) {
      s = treatment_times[i]
      treatment_history[, (paste0("A_", s)) := intervention_regime[i]]
    }
  }
  
  for (t_idx in 2:length(all_times)) {
    current_t = all_times[t_idx]
    newdata = data.table(id = ids, t = current_t)
    
    # Junta a história de tratamento COMPLETA, garantindo que A_1, A_2, etc. existam
    if (ncol(treatment_history) > 1) {
      newdata = merge(newdata, treatment_history, by = "id")
    }
    
    # Simula X(t)
    for(cov in cov_names) {
      pred = predict(gformula_models$confounder_models[[cov]],
                     newdata = newdata)
      newdata[, (cov) := pred + rnorm(n_ids, 0,
                                      gformula_models$noise_sds$X[cov])]
    }
    
    # Intervém em A(t) CONCORRENTE
    s_index = which(treatment_times == current_t)
    current_A = if (length(s_index) > 0) intervention_regime[s_index] else 0
    newdata[, A := current_A]
    
    # Simula Y(t)
    pred_y = predict(gformula_models$outcome_model, newdata = newdata)
    newdata[, y := pred_y + rnorm(n_ids, 0, gformula_models$noise_sds$y)]
    
    sim_data = rbind(sim_data, newdata[, .SD, .SDcols = names(sim_data)])
  }
  return(sim_data)
}

#' @title (FINAL) Executa a Suíte Completa de Estimação via G-Formula Funcional
#' @description Orquestra o ajuste dos modelos dinâmicos e a simulação de
#' intervenções para estimar o efeito causal.
#' @param scenario_data Um cenário da lista all_scenarios_data.
#' @param scenario_config A configuração correspondente.
#' @return Uma lista com o gráfico, métricas e modelos.

run_gformula_suite = function(scenario_data, scenario_config) {
  
  observed_data = as.data.table(scenario_data$outcome$long)
  regime_a = scenario_config$a_seq %||% c(1)
  regime_b = scenario_config$another_a_seq %||% c(0)
  
  # --- 1. ATO I: Aprender a "Física" do Mundo ---
  cat("--- Ajustando modelos da G-formula ---\n")
  gformula_models = fit_gformula_models(observed_data, scenario_config)
  
  # --- 2. ATO II: Simular Futuros Alternativos ---
  cat("\n--- Simulando intervenções ---\n")
  sim_data_a = simulate_gformula_intervention(gformula_models, observed_data, scenario_config, regime_a)
  sim_data_b = simulate_gformula_intervention(gformula_models, observed_data, scenario_config, regime_b)
  
  # --- 3. Calcular Efeito Causal e Métricas ---
  cat("--- Calculando efeito e métricas ---\n")
  mean_y_a = sim_data_a[, .(estimated_a = mean(y)), by = t]
  mean_y_b = sim_data_b[, .(estimated_b = mean(y)), by = t]
  causal_effect = merge(mean_y_a, mean_y_b, by = "t")
  causal_effect[, estimated := estimated_a - estimated_b]
  
  wide_data = scenario_data$outcome$wide
  real_total_dt = data.table(t = wide_data$t_vec, real_total = colMeans(wide_data$effect_a_vs_another))
  
  plot_data = merge(causal_effect, real_total_dt, by = "t")
  metrics = data.table(
    RMSE_vs_Total = sqrt(mean((plot_data$estimated - plot_data$real_total)^2)),
    Bias_vs_Total = mean(plot_data$estimated - plot_data$real_total),
    Cor_vs_Total = cor(plot_data$estimated, plot_data$real_total)
  )
  
  # --- 4. Geração do Gráfico ---
  subtitle_text = paste(
    "RMSE (vs Total):", round(metrics$RMSE_vs_Total, 3),
    "| Bias (vs Total):", round(metrics$Bias_vs_Total, 3),
    "| Corr (vs Total):", round(metrics$Cor_vs_Total, 3)
  )
  effect_plot = ggplot(plot_data, aes(x = t)) +
    geom_line(aes(y = real_total, color = "Real Total"), linetype = "dashed", linewidth = 1.1) +
    geom_line(aes(y = estimated, color = "Estimado (G-formula)"), linewidth = 1.2) +
    labs(
      title = "Efeito Causal Estimado (G-formula) vs. Real",
      subtitle = subtitle_text, x = "Tempo", y = "Efeito Causal", color = ""
    ) +
    scale_color_manual(values = c("Estimado (G-formula)" = "forestgreen", "Real Total" = "black")) +
    theme_minimal(base_size = 14)
  
  return(list(metrics = metrics, plot = effect_plot, models = gformula_models))
}
