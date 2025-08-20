source(file = "scripts_functional/data_generation/no_direct_effect_scenario.R")

sanity_check_functional = function(sim_data, n_check = 5, covariates = NULL) {
  library(data.table)
  
  # Se for lista "both", pega o long
  if (is.list(sim_data) && !is.data.table(sim_data)) {
    if (!is.null(sim_data$long)) sim_data = sim_data$long
  }
  
  dt = copy(as.data.table(sim_data))
  
  # Detecta covariáveis automaticamente se não informado
  if (is.null(covariates)) {
    covariates = grep("^X", names(dt), value = TRUE)
    covariates = grep(
      "^(?!a_|another_a_)", 
      covariates, 
      perl = TRUE, 
      value = TRUE
    )
    
  }
  
  ids_to_check = unique(dt$id)[1:min(n_check, length(unique(dt$id)))]
  
  results = list()
  
  for (id_i in ids_to_check) {
    dt_i = dt[id==id_i]
    for (cov in covariates) {
      a_cov = paste0("a_", cov)
      another_cov = paste0("another_a_", cov)
      # 1) Diferença entre a_X e another_a_X
      diff_a = max(abs(dt_i[[a_cov]] - dt_i[[another_cov]]), na.rm = TRUE)
      # 2) Diferença entre tempos na covariável natural
      diff_t = max(abs(diff(dt_i[[cov]])))
      # 3) Diferença entre baseline (tempo inicial) e outros tempos (natural)
      diff_base = max(abs(dt_i[[cov]] - dt_i[[cov]][1]))
      
      results[[paste0("id",id_i,"_",cov)]] = list(
        max_diff_a_vs_another = diff_a,
        max_diff_over_time = diff_t,
        max_diff_from_baseline = diff_base
      )
    }
  }
  
  return(results)
}

check = sanity_check_functional(stress_test$`Stress 5`$confounders)

check$id1_X1
check$id1_X2
check$id1_X3
check$id1_X4
