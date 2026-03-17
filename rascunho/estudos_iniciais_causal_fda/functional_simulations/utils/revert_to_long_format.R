#' @title Converte Dados Processados de volta para Formato Longo 
#' @description Pega a lista de dados processados (matrizes) e a converte
#'   para um único data.table longo.
#' @param processed_data A lista de saída de `process_and_smooth_data`.
#' @return Um data.table em formato longo.

revert_to_long_format = function(processed_data) {
  
  ids = processed_data$id
  time_grid = processed_data$time_grid
  
  dt_long = CJ(id = ids, t = time_grid)
  
  # Itera e "derrete" cada variável funcional da lista para colunas
  for (var_name in names(processed_data$functional_vars)) {
    dt_long[, (var_name) := as.vector(
      t(processed_data$functional_vars[[var_name]])
    )]
  }
  
  # Converte e junta as variáveis de tratamento (que são escalares por id)
  scalar_vars_long = melt(
    processed_data$scalar_vars,
    id.vars = "id",
    variable.name = "t_char",
    value.name = "A"
  )
  scalar_vars_long[, t := as.numeric(gsub("A_s", "", t_char))]
  
  # Junta a informação de tratamento de volta ao data.table principal
  dt_long[scalar_vars_long, on = .(id, t), A := i.A]
  dt_long[is.na(A), A := 0]
  
  # Se os desfechos contrafactuais existirem, calcula o efeito verdadeiro
  if (all(c("a_y", "b_y") %in% names(dt_long))) {
    dt_long[, effect_a_vs_b := a_y - b_y]
  }
  
  return(dt_long)
}