#' @title Decompõe e Plota os Efeitos Causais de um Cenário
#' @description Plota o efeito causal total (simulado) e o decompõe em suas
#' partes direta (teórica) e indireta (calculada).
#' @param scenario_data Um único elemento da lista `all_scenarios_data`.
#' @return Um objeto ggplot.

plot_scenario_effects = function(scenario_data) {
  require(ggplot2)
  require(data.table)
  
  dt_long = scenario_data$outcome$long
  time_grid = scenario_data$outcome$wide$t_vec
  
  # 1. Efeito Causal Total (Médio da Simulação)
  total_effect = dt_long[, .(value = mean(effect_a_vs_another)), by = t]
  total_effect[, type := "Total (Simulado)"]
  
  # 2. Efeito Direto (Teórico, a partir dos parâmetros)
  # Nota: Para S > 1, a soma pode ser mais complexa. Aqui, S=1.
  direct_effect_list = scenario_data$outcome$wide$true_effect_list
  direct_effect_vals = Reduce("+", direct_effect_list)
  direct_effect = data.table(t = time_grid, value = direct_effect_vals, type = "Direto (Teórico)")
  
  # 3. Efeito Indireto (Diferença)
  indirect_effect = data.table(
    t = time_grid,
    value = total_effect$value - direct_effect$value,
    type = "Indireto (via X)"
  )
  
  dt_plot = rbindlist(list(total_effect, direct_effect, indirect_effect))
  
  ggplot(dt_plot, aes(x = t, y = value, color = type, linetype = type)) +
    geom_line(linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    scale_color_manual(values = c("Total (Simulado)" = "black", "Direto (Teórico)" = "#0072B2", "Indireto (via X)" = "#D55E00")) +
    scale_linetype_manual(values = c("Total (Simulado)" = "solid", "Direto (Teórico)" = "dashed", "Indireto (via X)" = "dotted")) +
    labs(
      title = "Decomposição do Efeito Causal Médio",
      subtitle = "Total = Direto + Indireto (mediado pelas covariáveis)",
      x = "Tempo",
      y = "Efeito Causal",
      color = "Componente",
      linetype = "Componente"
    ) +
    theme_light()
}