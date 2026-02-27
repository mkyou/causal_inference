#' @title Plota Trajetórias Individuais de um Cenário (Versão Corrigida)
#' @description Gera gráficos das trajetórias de covariáveis e desfechos para
#' uma amostra de indivíduos, ajudando a validar visualmente o PGD.
#' @param scenario_data Um único elemento da lista `all_scenarios_data`.
#' @param n_plot Número de indivíduos a serem plotados.
#' @return Um objeto ggplot.

plot_scenario_trajectories = function(scenario_data, n_plot = 10) {
  require(ggplot2)
  require(data.table)
  
  # Usamos a saída completa (que contém a tabela 'long')
  dt = scenario_data$outcome$long
  ids_to_plot = unique(dt$id)[1:min(n_plot, length(unique(dt$id)))]
  dt_plot = dt[id %in% ids_to_plot]
  
  # --- BLOCO CORRIGIDO E MELHORADO ---
  # Em vez de 'patterns()', selecionamos explicitamente as colunas de interesse.
  # Isso é mais robusto e agora inclui as covariáveis contrafactuais.
  measure_cols = names(dt_plot)[grep("^(X|a_X|another_a_X|y|a_y|another_a_y)", names(dt_plot))]
  
  # Reformatar para formato longo para facetting
  dt_melt = melt(dt_plot, 
                  id.vars = c("id", "t"), 
                  measure.vars = measure_cols,
                  variable.name = "Variable", 
                  value.name = "Value" # Garantimos que a coluna se chame 'Value'
  )
  # --- FIM DO BLOCO CORRIGIDO ---
  
  ggplot(dt_melt, aes(x = t, y = Value, group = id, color = factor(id))) +
    geom_line(alpha = 0.8) +
    facet_wrap(~ Variable, scales = "free_y") +
    labs(
      title = "Trajetórias Individuais (Amostra)",
      subtitle = "Covariáveis e Desfechos (Observados e Contrafactuais)",
      x = "Tempo",
      y = "Valor",
      color = "ID do Indivíduo"
    ) +
    theme_light() +
    theme(legend.position = "none")
}
