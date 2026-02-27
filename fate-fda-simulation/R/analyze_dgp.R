library(ggplot2)
library(data.table)
library(patchwork)

#' Valida Distribuições Estáticas (Overlap e Balanceamento)
#' @param dgp_data O objeto retornado por generate_sequential_dgp
validate_static_distributions = function(dgp_data) {
  
  # Pega o t=0 para analisar as covariáveis basais e tratamento
  dt_static = dgp_data[t == min(t)]
  
  # 1. Checagem de Overlap em X1
  p1 = ggplot(dt_static, aes(x = X1, fill = factor(Trat_1))) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(
      values = c("firebrick", "royalblue"), 
      name="Tratamento"
    ) +
    labs(
      title = "Distribuição de X1 por Tratamento", 
      subtitle = "Verificação de Overlap/Suporte Comum",
      x = "Covariável X1"
    ) +
    theme_bw() + theme(legend.position = "bottom")
  
  # 2. Balanço do Tratamento
  prop = mean(dt_static$Trat_1)
  p2 = ggplot(dt_static, aes(x = factor(Trat_1))) +
    geom_bar(fill = "steelblue", alpha = 0.8) +
    labs(
      title = sprintf("Proporção de Tratados: %.1f%%", prop*100), 
      x = "Status (0/1)", y = "Contagem"
    ) +
    theme_bw()
  
  return(p1 + p2)
}

#' Valida Efeito Causal (True vs Naive) usando dados gerados
#' @param dgp_outcome O data.table retornado por generate_functional_outcome
#' @param t_treat O tempo do tratamento (para linha vertical)
validate_dgp_effects = function(dgp_outcome, t_treat) {
  
  # Verifica se as colunas necessárias existem
  cols_req = c("t", "y", "Trat_1", "effect_a_vs_b")
  if (!all(cols_req %in% names(dgp_outcome))) {
    stop("Colunas necessárias não encontradas no dataframe.")
  }
  
  # 1. ATE Verdadeiro (Ground Truth)
  # Média direta da coluna de efeito contrafactual
  ate_true = dgp_outcome[, .(Value = mean(effect_a_vs_b)), by = t]
  ate_true[, Type := "Verdadeiro (E[Ya - Yb])"]
  
  # 2. ATE Naive (Observacional / Biased)
  # Calcula médias por grupo
  naive_dt = dgp_outcome[, .(mean_y = mean(y)), by = .(t, Trat_1)]
  
  # Cria coluna de nome seguro para o dcast
  naive_dt[, Group := ifelse(Trat_1 == 1, "Treated", "Control")]
  
  # Pivotar (Wide)
  naive_wide = dcast(naive_dt, t ~ Group, value.var = "mean_y")
  
  # Garante que as colunas existam (para evitar erro se N for pequeno)
  if(!"Treated" %in% names(naive_wide)) naive_wide[, Treated := NA_real_]
  if(!"Control" %in% names(naive_wide)) naive_wide[, Control := NA_real_]
  
  # Calcula a diferença
  naive_wide[, Value := Treated - Control]
  naive_wide[, Type := "Naive (E[Y|1] - E[Y|0])"]
  
  # Seleciona colunas para empilhar
  ate_naive = naive_wide[, .(t, Value, Type)]
  
  # Combina para plot
  dt_plot = rbind(ate_true, ate_naive)
  
  # Plot
  p = ggplot(dt_plot, aes(x = t, y = Value, color = Type, linetype = Type)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = t_treat, linetype = "dotted", color = "black") +
    scale_color_manual(values = c("orange", "black")) +
    scale_linetype_manual(values = c("dashed", "solid")) +
    labs(
      title = "Validação Física: Efeito Verdadeiro vs. Naive",
      subtitle = "Discrepância = Viés de Confundimento (X -> Y)",
      y = "Efeito Médio do Tratamento (ATE)", 
      x = "Tempo"
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  return(p)
}

#' Valida Trajetórias Médias (Observado vs Contrafactuais)
#' Verifica se a "forma" global do desfecho é dominada pela 
#' tendência suave ou pelo efeito.
#' @param dgp_outcome Data.table com colunas y, a_y, b_y e t
validate_outcome_means = function(dgp_outcome, t_treat) {
  
  # Médias por tempo (agregando todos os indivíduos)
  dt_means = dgp_outcome[, .(
    Y_Factual = mean(y), # Média do que observamos
    Y0_Contrafactual = mean(b_y), # Se ninguém fosse tratado (Baseline puro)
    Y1_Contrafactual = mean(a_y) # Se todos fossem tratados (Base + Efeito)
  ), by = t]
  
  dt_long = melt(
    dt_means, id.vars = "t", 
    variable.name = "Curva", value.name = "Valor"
    )
  
  p = ggplot(
    dt_long, aes(x = t, y = Valor, color = Curva, linetype = Curva)
    ) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = t_treat, linetype = "dotted", color = "black") +
    scale_color_manual(values = c(
      "Y_Factual" = "black",
      "Y0_Contrafactual" = "royalblue", # Azul = Controle
      "Y1_Contrafactual" = "firebrick"  # Vermelho = Tratado
    )) +
    scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
    labs(
      title = "Trajetórias Médias Absolutas",
      subtitle = "Visualize a suavidade global e o impacto relativo do 
      tratamento",
      y = "Valor Médio do Desfecho (Y)", x = "Tempo"
    ) +
    theme_bw() +
    theme(legend.position = "top")
  
  return(p)
}