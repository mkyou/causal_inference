# ===================================================================
# Suíte de Análise de Cenário de Simulação (V2.1 - Corrigida)
#
# Contém a suíte completa de funções de análise com a correção
# na função `plot_confounding_structure` para lidar com a busca
# de tempo de forma mais robusta.
# ===================================================================

library(data.table)
library(ggplot2)

# -------------------------------------------------------------------
# Função 1: Checar Estrutura Causal (sem alterações)
# -------------------------------------------------------------------
check_causal_structure = function(scenario_data, scenario_config) {
  dt = scenario_data
  
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  feedback_active_on = c()
  for (cov in cov_names) {
    a_cov = paste0("a_", cov)
    b_cov = paste0("b_", cov)
    if (max(abs(dt[[a_cov]] - dt[[b_cov]]), na.rm = TRUE) > 1e-9) {
      feedback_active_on = c(feedback_active_on, cov)
    }
  }
  feedback_summary = if (length(feedback_active_on) > 0) 
    paste(feedback_active_on, collapse = ", ") else "Nenhum"
  
  confounding_summary = "Nenhum (por desenho)"
  if (!is.null(scenario_config$theta_matrix)) {
    if (any(scenario_config$theta_matrix[, -1] != 0) || 
        !is.null(scenario_config$theta_integral_matrix)) {
      confounding_summary = "Presente (por desenho)"
    }
  }
  
  return(data.table(
    Feedback_Ativo_em = feedback_summary,
    Confundimento_Desenhado = confounding_summary
  ))
}

# -------------------------------------------------------------------
# Função 2: Plotar Trajetórias (sem alterações)
# -------------------------------------------------------------------
plot_scenario_trajectories = function(scenario_data, n_plot = 10) {
  dt = scenario_data
  ids_to_plot = unique(dt$id)[1:min(n_plot, length(unique(dt$id)))]
  dt_plot = dt[id %in% ids_to_plot]
  
  measure_cols = names(dt_plot)[grep("^(X|a_X|b_X|y|a_y|b_y)", names(dt_plot))]
  
  dt_melt = melt(
    dt_plot, 
    id.vars = c("id", "t"), 
    measure.vars = measure_cols,
    variable.name = "Variavel", 
    value.name = "Valor"
  )
  
  ggplot(dt_melt, aes(x = t, y = Valor, group = id, color = factor(id))) +
    geom_line(alpha = 0.7) +
    facet_wrap(~ Variavel, scales = "free_y") +
    labs(
      title = "Trajetórias Individuais Geradas (Amostra)",
      subtitle = "Factual (X, y) e Contrafactual (a_*, b_*)",
      x = "Tempo", y = "Valor", color = "ID"
    ) +
    theme_bw() +
    theme(legend.position = "none")
}

# -------------------------------------------------------------------
# Função 3: Plotar Efeito Causal Total (sem alterações)
# -------------------------------------------------------------------
plot_scenario_effects = function(scenario_data) {
  dt = scenario_data
  
  total_effect = dt[, .(efeito_medio=mean(effect_a_vs_b, na.rm=TRUE)), by = t]
  
  ggplot(total_effect, aes(x = t, y = efeito_medio)) +
    geom_line(color = "#0072B2", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(
      title = "Efeito Causal Médio Total no Desfecho (Y)",
      subtitle = "Diferença Média: E[Y | Regime A] - E[Y | Regime B]",
      x = "Tempo",
      y = "Efeito Causal Médio (a_y - b_y)"
    ) +
    theme_bw()
}


# -------------------------------------------------------------------
# Função 4: Plotar Estrutura de Confundimento (CORRIGIDA)
# -------------------------------------------------------------------
#' @title (V2.2 - Corrigida) Visualiza a Força do Confundimento
plot_confounding_structure = function(scenario_data, scenario_config) {
  dt = scenario_data
  
  plot_data_list = list()
  
  # --- INÍCIO DA CORREÇÃO ---
  # Itera diretamente sobre os valores de treatment_times para evitar o
  # bug com seq_along() quando o comprimento é 1.
  for (s_time in scenario_config$treatment_times) {
    # --- FIM DA CORREÇÃO ---
    
    t_idx = which.min(abs(scenario_config$time_grid - s_time))
    prev_t_idx = t_idx - 1
    
    if (prev_t_idx < 1) next
    
    prev_time = scenario_config$time_grid[prev_t_idx]
    
    cov_name = "X1"
    
    temp_dt = dt[t %in% c(prev_time, s_time), .(id, t, A, X1)]
    # Verifica se há dados suficientes para o dcast
    if(nrow(temp_dt) == 0 || temp_dt[, .N, by=id][, all(N < 2)]){
      warning(paste("Dados insuficientes para dcast no tempo", s_time))
      next
    }
    
    temp_dt_wide = dcast(
      temp_dt, id ~ t, value.var = c("A", "X1")
    )
    
    old_names = c(
      paste0("X1_", prev_time), 
      paste0("A_", s_time)
    )
    # Checa se as colunas esperadas existem antes de renomear
    if (!all(old_names %in% names(temp_dt_wide))) {
      warning(paste("Colunas esperadas não encontradas para o tempo", s_time))
      next
    }
    
    new_names = c("Cov_pre_treatment", "Treatment_status")
    
    a_prev_name = paste0("A_", prev_time)
    if (a_prev_name %in% names(temp_dt_wide)) {
      temp_dt_wide[, (a_prev_name) := NULL]
    }
    
    setnames(temp_dt_wide, old_names, new_names)
    temp_dt_wide[, Treatment_Time := as.factor(s_time)]
    # Adiciona a s_idx para o caso de a lista ter nomes
    s_idx_char = as.character(s_time) 
    plot_data_list[[s_idx_char]] = temp_dt_wide
  }
  
  if (length(plot_data_list) == 0) {
    warning("Nenhum dado de confundimento para plotar.")
    return(NULL)
  }
  
  plot_data = rbindlist(plot_data_list, fill = TRUE)
  
  # Checa se a coluna Treatment_status existe e não é toda NA
  if (is.null(plot_data$Treatment_status) || all(is.na(plot_data$Treatment_status))) {
    warning("Coluna Treatment_status está vazia ou toda NA.")
    return(NULL)
  }
  
  ggplot(plot_data, aes(x = Cov_pre_treatment, y = Treatment_status)) +
    geom_jitter(
      aes(color = Treatment_Time), 
      width = 0, height = 0.05, alpha = 0.5
    ) +
    geom_smooth(method = "glm", method.args = list(family = "binomial")) +
    facet_wrap(~ Treatment_Time, labeller = label_both) +
    labs(
      title = "Verificação do Confundimento (X -> A)",
      subtitle = "Relação entre Covariável e Status do Tratamento",
      x = "Valor da Covariável (ex: X1) no Tempo s-1",
      y = "Probabilidade de Receber Tratamento em s",
      color = "Tempo de Trat."
    ) +
    theme_bw()
}


# -------------------------------------------------------------------
# Função 5: Ponto de Entrada (sem alterações)
# -------------------------------------------------------------------
analyze_scenario = function(
    scenario_data,
    scenario_config,
    scenario_name = "Cenário"
) {
  
  cat(paste("\n--- Analisando Cenário:", scenario_name, "---\n"))
  
  cat("  - Verificando estrutura causal...\n")
  summary_table = check_causal_structure(scenario_data, scenario_config)
  summary_table[, Scenario := scenario_name]
  
  cat("  - Gerando gráfico de trajetórias...\n")
  p_traj = plot_scenario_trajectories(scenario_data)
  
  cat("  - Gerando gráfico do efeito causal...\n")
  p_effects = plot_scenario_effects(scenario_data)
  
  cat("  - Gerando gráfico da estrutura de confundimento...\n")
  p_confounding = plot_confounding_structure(
    scenario_data,
    scenario_config
  )
  
  p_traj = p_traj + labs(title = paste(
    "Trajetórias Individuais | Cenário:", scenario_name
  ))
  p_effects = p_effects + labs(title = paste(
    "Efeito Causal Médio | Cenário:", scenario_name
  ))
  if(!is.null(p_confounding)) {
    p_confounding = p_confounding + labs(title = paste(
      "Estrutura de Confundimento | Cenário:", scenario_name
    ))
  }
  
  results = list(
    summary_table = summary_table,
    plots = list(
      trajectories = p_traj,
      causal_effect = p_effects,
      confounding_structure = p_confounding
    )
  )
  
  cat(paste("--- Análise para", scenario_name, "concluída. ---\n"))
  return(results)
}