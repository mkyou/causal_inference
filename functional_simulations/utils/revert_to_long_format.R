#' @title Converte Dados do Formato Refund de volta para Longo
#' @description Pega a lista de matrizes da `prepare_refund_data` e a
#'   "re-hidrata" para um único data.table longo, compatível com `mgcv`.
#' @param refund_data A lista de saída da `prepare_refund_data`.
#' @return Um data.table em formato longo.

revert_to_long_format = function(refund_data) {
  
  ids = refund_data$id
  time_grid = refund_data$time_grid
  
  dt_long = CJ(id = ids, t = time_grid)
  
  # Itera e "derrete" cada variável funcional
  for (var_name in names(refund_data$functional_vars)) {
    dt_long[, (var_name) := as.vector(
      t(refund_data$functional_vars[[var_name]])
    )]
  }
  
  # Converte e junta as variáveis de tratamento escalares
  scalar_vars_long = melt(
    refund_data$scalar_vars, 
    id.vars = "id", 
    variable.name = "t_char", 
    value.name = "A"
  )
  scalar_vars_long[, t := as.numeric(gsub("A_s", "", t_char))]
  
  dt_long[scalar_vars_long, on = .(id, t), A := i.A]
  dt_long[is.na(A), A := 0]
  
  # Recalcula colunas derivadas para compatibilidade com a suíte de análise
  if (all(c("a_y", "b_y") %in% names(dt_long))) {
    dt_long[, effect_a_vs_b := a_y - b_y]
  }
  
  return(dt_long)
}