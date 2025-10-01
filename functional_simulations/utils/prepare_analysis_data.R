library(data.table)

prepare_analysis_data = function(data_long, treatment_times, ...) {
  dt = copy(data_long)
  setkey(dt, id, t)
  
  cov_names = grep("^X[0-9]+$", names(dt), value = TRUE)
  
  # Revertendo para a criação de features de lag simples
  lagged_names = paste0(cov_names, "_lag1")
  dt[
    , (lagged_names) := shift(.SD, 1, type = "lag"), 
    by = id, .SDcols = cov_names
  ]
  
  treatment_status = dt[t %in% treatment_times, .(id, t, A)]
  wide_treatment = dcast(treatment_status, id ~ t, value.var = "A")
  new_names = paste0("A_s", treatment_times)
  setnames(wide_treatment, as.character(treatment_times), new_names)
  dt = wide_treatment[dt, on = "id"]
  
  for (s_time in treatment_times) {
    by_col = paste0("by_s", s_time)
    A_s_col = paste0("A_s", s_time)
    dt[
      , (by_col) := .SD[[A_s_col]] * fifelse(t >= s_time, 1, 0), 
      .SDcols = A_s_col
    ]
  }
  
  return(dt)
}