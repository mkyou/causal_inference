library(data.table)

#' @title Prepara os Dados para Análise da G-Formula
#' @description Função centralizada que recebe os dados brutos em formato
#'   longo e adiciona todas as colunas necessárias para o ajuste dos
#'   modelos (lags e tratamentos por tempo).
#'
#' @param data_long O data.table no formato longo.
#' @param treatment_times Vetor com os tempos de tratamento.
#' @return Um data.table pronto para análise.

prepare_analysis_data = function(data_long, treatment_times) {
  dt = copy(data_long)
  setkey(dt, id, t)
  
  # Adiciona variáveis de lag para as covariáveis
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  lagged_names = paste0(cov_names, "_lag1")
  dt[
    , 
    (lagged_names) := shift(.SD, n = 1, type = "lag"), 
    by = id, 
    .SDcols = cov_names
  ]
  
  # Adiciona colunas de tratamento (A_s)
  for (s_time in treatment_times) {
    A_s_name = paste0("A_s", s_time)
    treatment_status_at_s = dt[t == s_time, .(id, A)]
    setnames(treatment_status_at_s, "A", A_s_name)
    dt = merge(dt, treatment_status_at_s, by = "id", all.x = TRUE)
  }
  
  # Adiciona colunas de "tempo desde o tratamento" para alinhar com o PGD
  for (s_time in treatment_times) {
    # Variável de tempo desde o tratamento
    time_since_col = paste0("time_since_s", s_time)
    dt[, (time_since_col) := fifelse(t >= s_time, t - s_time, 0)]
    
    # Variável 'by' para o spline
    by_col = paste0("by_s", s_time)
    A_s_col = paste0("A_s", s_time)
    dt[
      , (by_col) := .SD[[A_s_col]] * fifelse(t >= s_time, 1, 0), 
      .SDcols = A_s_col
    ]
  }
  
  return(dt)
}