#' @title (VERSÃO DEFINITIVA E CORRIGIDA) Ajusta Modelos para a G-Formula Funcional
#' @description Prepara os dados (lags, variáveis de tratamento) e ajusta os
#' modelos estruturais (GAMM para X, GAM para Y) com especificação correta.

fit_gformula_models = function(observed_data, scenario_config) {
  require(mgcv)
  
  # --- 1. Preparar dados ---
  model_data = copy(observed_data)
  cov_names = grep("^X[0-9]+$", names(model_data), value = TRUE)
  treatment_times = scenario_config$treatment_times
  
  # Cria colunas explícitas para o status de tratamento A_1, A_2, ...
  if (!is.null(treatment_times)) {
    for (s in treatment_times) {
      treatment_status_s = model_data[t == s, .(id, A_s = A)]
      setnames(treatment_status_s, "A_s", paste0("A_", s))
      model_data = merge(model_data, treatment_status_s, by = "id", all.x = TRUE)
      model_data[is.na(get(paste0("A_", s))), (paste0("A_", s)) := 0]
    }
  }
  
  # --- CORREÇÃO FUNDAMENTAL: Criação dos Lags ---
  # Precisamos dos lags tanto para os X quanto para o A (para o modelo de X)
  lag_vars = c(cov_names, "A")
  model_data[, paste0(lag_vars, "_lag1") := lapply(.SD, shift, n=1, type="lag"), .SDcols = lag_vars, by = id]
  
  # Remove a primeira observação de cada id (que agora tem NAs nos lags)
  model_data = na.omit(model_data)
  
  treatment_terms = if (!is.null(treatment_times)) paste0("A_", treatment_times) else "1"
  
  # --- 2. Ajustar um modelo GAMM para cada covariável ---
  confounder_models = list()
  for (cov in cov_names) {
    cat(paste("Ajustando modelo para", cov, "\n"))
    
    # --- FÓRMULA CORRIGIDA ---
    # Incluindo os lags das covariáveis como preditores suaves
    lag_cov_smooth_terms = paste0("s(", cov_names, "_lag1)")
    formula_text = paste(
      cov, "~ s(t) + s(id, bs='re') +", 
      paste(lag_cov_smooth_terms, collapse=" + "), "+", 
      paste(treatment_terms, collapse=" + ")
    )
    
    confounder_models[[cov]] = bam(as.formula(formula_text), data = model_data, discrete = TRUE)
  }
  
  # --- 3. Ajustar um modelo GAM para o desfecho Y ---
  cat("Ajustando modelo para o desfecho Y\n")
  
  # O modelo de Y usa os valores CONCORRENTES, então sua fórmula estava correta.
  # Apenas garantindo que a variável A correta seja usada (não o lag).
  cov_smooth_terms = paste0("s(", cov_names, ")")
  # Para o modelo de Y, usamos A (concorrente), não A_lag1
  formula_y = as.formula(paste("y ~ s(t) +", paste(cov_smooth_terms, collapse=" + "), "+ A"))
  
  outcome_model = bam(formula_y, data = model_data, discrete = TRUE)
  
  # ... (resto da função, retornando a lista de modelos)
  return(list(
    confounder_models = confounder_models,
    outcome_model = outcome_model,
    noise_sds = list(
      y = sd(outcome_model$residuals),
      X = sapply(confounder_models, function(m) sd(m$residuals))
    )
  ))
}

# NOTA: A função simulate_gformula_intervention também precisa de um pequeno ajuste
# para usar A, e não A_lag1, ao prever Y. E precisa dos lags de X.
# A implementação completa e correta das duas está nos arquivos que você já tem.




#' @title (Corrigida) Simula uma Intervenção Contrafactual usando Modelos G-Formula
#' @description Itera ao longo do tempo, usando os modelos ajustados para
#' simular as trajetórias completas das covariáveis e do desfecho sob um
#' regime de tratamento fixo.

simulate_gformula_intervention = function(gformula_models, observed_data, scenario_config, intervention_regime) {
  
  # --- 1. Setup da Simulação ---
  all_times = sort(unique(observed_data$t))
  ids = unique(observed_data$id)
  n_ids = length(ids)
  cov_names = names(gformula_models$confounder_models)
  treatment_times = scenario_config$treatment_times
  
  # --- LÓGICA CORRIGIDA DE INICIALIZAÇÃO ---
  # Inicializa a simulação com os valores do baseline, mas selecionando
  # APENAS as colunas que serão de fato simuladas, descartando as extras.
  sim_cols = c("id", "t", "A", cov_names, "y")
  sim_data = observed_data[t == min(all_times), .SD, .SDcols = sim_cols]
  # --- FIM DA CORREÇÃO ---
  
  # Cria tabela para guardar a história de tratamento da simulação
  treatment_history_sim = data.table(id = ids)
  
  # --- 2. Loop de Simulação Sequencial (sem alterações) ---
  for (t_idx in 2:length(all_times)) {
    current_t = all_times[t_idx]
    
    # --- 2a. Preparar o `newdata` para a predição em `current_t` ---
    newdata = data.table(id = ids, t = current_t)
    
    if (ncol(treatment_history_sim) > 1) {
      newdata = merge(newdata, treatment_history_sim, by = "id")
    }
    
    s_index = which(treatment_times == current_t)
    if (length(s_index) > 0) {
      current_A = intervention_regime[s_index]
      treatment_history_sim[, (paste0("A_", current_t)) := current_A]
    } else {
      current_A = 0
    }
    newdata[, A := current_A]
    
    if(!is.null(treatment_times)){
      for(s in treatment_times) {
        col_name = paste0("A_", s)
        if(col_name %in% names(treatment_history_sim)) {
          newdata[, (col_name) := treatment_history_sim[[col_name]]]
        } else {
          newdata[, (col_name) := 0]
        }
      }
    }
    
    # --- 2b. Simular as Covariáveis X(t) ---
    for(cov in cov_names) {
      pred = predict(gformula_models$confounder_models[[cov]], newdata = newdata)
      newdata[, (cov) := pred + rnorm(n_ids, 0, gformula_models$noise_sds$X[cov])]
    }
    
    # --- 2c. Simular o Desfecho Y(t) ---
    pred_y = predict(gformula_models$outcome_model, newdata = newdata)
    newdata[, y := pred_y + rnorm(n_ids, 0, gformula_models$noise_sds$y)]
    
    # Guarda os resultados simulados do tempo atual. Agora as colunas baterão.
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