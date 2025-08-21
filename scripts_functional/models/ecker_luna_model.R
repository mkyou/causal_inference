#' @title (Renomeada) Prepara Dados para o Modelo de Ecker & Luna
#' @description Extrai covariáveis do baseline de um cenário com tratamento 
#' pontual.
#' @param scenario_data Um único cenário da lista `all_scenarios_data`.
#' @return Uma data.table pronta para modelagem.

prepare_data_eckerluna = function(scenario_data) {
  dt_long = as.data.table(scenario_data$outcome$long)
  
  treatment_time = dt_long[A == 1, min(t, na.rm = TRUE)]
  if (is.infinite(treatment_time)) stop("Nenhum tratamento encontrado")
  
  all_times = sort(unique(dt_long$t))
  baseline_time = max(all_times[all_times <= treatment_time])
  
  cov_names = grep("^X[0-9]+$", names(dt_long), value = TRUE)
  base_cov_names = paste0(cov_names, "_base")
  
  baseline_data = dt_long[
    t == baseline_time, .SD, .SDcols = c("id", cov_names)
  ]
  setnames(baseline_data, cov_names, base_cov_names)
  
  model_data = merge(dt_long[t >= baseline_time], baseline_data, by = "id")
  
  treated_ids = dt_long[t == treatment_time & A == 1, unique(id)]
  model_data[, treated := (id %in% treated_ids)]
  
  return(model_data)
}


#' @title Ajusta Modelos GAM para o método de Ecker & Luna
#' @description Implementa um modelo de coeficientes variantes usando a base 
#' 'fs' do mgcv.
#' @param model_data Dados da função `prepare_data_eckerluna`.
#' @return Uma lista com os dois modelos GAM ajustados.

#' @title (VERSÃO 2) Ajusta Modelos GAM Separados de Forma Robusta
#' @description Usa a especificação padrão de coeficientes variantes (s(t, by=...))
#' que é mais estável para 'by' numérico. Mantém o truque de "dummy data".
#' @param model_data Dados da função `prepare_data_eckerluna`.
#' @return Uma lista com os dois modelos GAM ajustados.

fit_ecker_luna_gam = function(model_data) {
  require(mgcv)
  
  # --- Construção da Fórmula (Simplificada e Mais Robusta) ---
  base_cov_names = grep("_base$", names(model_data), value = TRUE)
  
  # Usando a sintaxe padrão e mais estável para coeficientes variantes
  # y ~ s(t) + s(t, by = X1_base) + s(t, by = X2_base) + ...
  smooth_terms = paste0("s(t, by = ", base_cov_names, ", k = 15)")
  model_formula = as.formula(paste("y ~ s(t, k = 20) +", paste(smooth_terms, collapse = " + ")))
  
  # --- Preparação dos Dados de Treinamento (Mantém o truque) ---
  data_treated = model_data[treated == TRUE]
  data_control = model_data[treated == FALSE]
  
  dummy_from_control = data_control[sample(.N, min(3, .N))]
  dummy_from_treated = data_treated[sample(.N, min(3, .N))]
  
  train_data_treated = rbind(data_treated, dummy_from_control)
  train_data_control = rbind(data_control, dummy_from_treated)
  
  # --- Ajuste dos Modelos ---
  fit_treated = bam(
    model_formula, 
    data = train_data_treated, 
    discrete = TRUE,
    select = TRUE
  )
  
  fit_control = bam(
    model_formula, 
    data = train_data_control, 
    discrete = TRUE,
    select = TRUE
  )
  
  return(list(treated = fit_treated, control = fit_control))
}


#' @title (Renomeada) Executa a Suíte de Estimação para o Modelo de 
#' Ecker & Luna
#' @description Orquestra todo o processo de ajuste, predição, cálculo de 
#' métricas e plotagem para o modelo proxy de Ecker & Luna.
#' @param scenario_data Um único cenário da lista `all_scenarios_data`.
#' @return Uma lista contendo o gráfico, as métricas e os modelos.

run_eckerluna_gam_suite = function(scenario_data) {
  
  # --- 1. Preparação dos Dados ---
  model_data = prepare_data_eckerluna(scenario_data)
  
  # --- 2. Ajuste dos Modelos ---
  fitted_models = fit_ecker_luna_gam(model_data)
  
  # --- 3. Predição e Cálculo do Efeito ---
  model_data[, pred_Y1 := predict(fitted_models$treated, newdata = model_data)]
  model_data[, pred_Y0 := predict(fitted_models$control, newdata = model_data)]
  model_data[, effect_estimate := pred_Y1 - pred_Y0]
  causal_effect = model_data[, .(estimated = mean(effect_estimate)), by = t]
  
  # --- 4. Coleta dos Efeitos Reais ---
  wide_data = scenario_data$outcome$wide
  real_total_dt = data.table(
    t = wide_data$t_vec, 
    real_total = colMeans(wide_data$effect_a_vs_another)
  )
  real_direct_dt = data.table(
    t = wide_data$t_vec, 
    real_direct = Reduce("+", wide_data$true_effect_list)
  )
  
  # --- 5. Cálculo das Métricas ---
  plot_data = merge(
    merge(causal_effect, real_total_dt, by = "t"), 
    real_direct_dt, by = "t"
  )
  
  metrics = data.table(
    RMSE_vs_Total = sqrt(mean((plot_data$estimated - plot_data$real_total)^2)),
    Bias_vs_Total = mean(plot_data$estimated - plot_data$real_total),
    Cor_vs_Total = cor(plot_data$estimated, plot_data$real_total)
  )
  
  # --- 6. Geração do Gráfico ---
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

# Supondo que `all_scenarios_data` foi gerado com sucesso.

# Forma correta de executar a análise para todos os cenários:
# Usamos lapply para aplicar a função a CADA elemento da lista de cenários.
all_results = lapply(all_scenarios_data, run_eckerluna_gam_suite)

# O resultado, `all_results`, será uma nova lista, com a mesma estrutura de nomes,
# onde cada elemento contém as métricas e o gráfico para aquele cenário.

# Para ver os resultados do Cenário 1:
print(all_results$`Cenario 1 (Ideal E&L)`$metrics)
print(all_results$`Cenario 1 (Ideal E&L)`$plot)

print(all_results$`Cenario 2 (Cov Dinamica)`$metrics)
print(all_results$`Cenario 2 (Cov Dinamica)`$plot)

# Para ver os resultados do Cenário 3 (o de feedback):
print(all_results$`Cenario 3 (Feedback)`$metrics)
print(all_results$`Cenario 3 (Feedback)`$plot)
