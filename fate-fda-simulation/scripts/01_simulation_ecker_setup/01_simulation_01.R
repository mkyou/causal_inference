library(data.table)
library(future.apply)
library(progressr)
library(ggplot2)
library(dplyr)

# Carrega tudo
source("R/scenario_runner.R")
source("R/process_functional_data.R")
source("R/estimate_fate.R")

# CONFIG
N_SIM = 1000
BATCH_SIZE = 50
N_SUBJ = 2000
T_PTS = 64

scenarios = list(
  # Sigmoid: DD4 + Soft + Alpha baixo (0.1) para preservar suavidade
  list(
    name = "Sigmoid", mag = 2, 
    family = "cdf97", wt = "DaubLeAsymm",
    wt_filter = 8,
    thresh = "soft", alpha = 0.5, beta = 1
  ),
  
  # Doppler: CDF97 + Semisoft + Beta padrão (1.2)
  list(
    name = "Doppler", mag = 5, 
    family = "cdf97", wt = "DaubLeAsymm", wt_filter = 6,
    thresh = "semisoft", alpha = 0.3, beta = 1.2
  ),
  
  # Bumps: DB2 + Semisoft
  list(
    name = "Bumps", mag = 3, family = "db2", 
    wt = "DaubExPhase", wt_filter = 2, 
    thresh = "semisoft", alpha = 0.3, beta = 1.2
  ),
  
  # Heavisine: Haar + HARD Threshold (Melhor para degraus puros) + Beta baixo
  list(
    name = "Heavisine", mag = 4, 
    family = "haar", wt = "Haar", wt_filter = 1,
    thresh = "hard", alpha = 0.3, beta = 0.8
  )
)

methods = c("spline", "wavethresh", "rlifting_offline")

plan(multisession)
handlers(global = TRUE)

run_full_simulation = function() {
  n_batches = ceiling(N_SIM / BATCH_SIZE)
  
  # Create localized checkpoint directory
  ckpt_dir = file.path('output/checkpoints', sprintf("N%d", N_SUBJ))
  dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (scen in scenarios) {
    cat(sprintf("\n>>> Cenário: %s <<<\n", scen$name))
    
    for (b in 1:n_batches) {
      batch_file = file.path(
        ckpt_dir, 
        sprintf("%s_batch_%03d.rds", scen$name, b)
        )
      
      if (file.exists(batch_file)) {
        cat(sprintf("  Batch %d/%d já existe. Pulando.\n", b, n_batches))
        next
      }
      
      # Define intervalo
      start_idx = (b - 1) * BATCH_SIZE + 1
      end_idx = min(b * BATCH_SIZE, N_SIM)
      current_n = end_idx - start_idx + 1
      
      cat(
        sprintf(
          "  Processando Batch %d/%d (Sims %d-%d)... ", b, n_batches, 
          start_idx, end_idx
          )
        )
      
      # Executa apenas este lote
      p <- progressor(steps = current_n)
      
      sim_results_list = future_lapply(start_idx:end_idx, function(i) {
        dgp_data = generate_single_scenario_dataset(
          scen$name, scen$mag, 
          N_SUBJ, T_PTS
        )
        
        true_curve = attr(dgp_data, "true_effect_curve")
        true_range = diff(range(true_curve))
        if(true_range == 0) true_range = 1.0
        
        time_grid = attr(dgp_data, "time_grid")
        X_mat = attr(dgp_data, "X_mat")
        A_vec = attr(dgp_data, "A_vec")
        
        y_wide = as.matrix(dcast(dgp_data, id ~ t, value.var = "y")[, -1])
        
        # Limpa objeto pesado da memória do worker imediatamente
        rm(dgp_data) 
        
        res_sim_long = list()
        
        for (m in methods) {
          extra = list()
          
          if(m == "rlifting_offline") {
            extra = list(
              scheme_name = scen$family,
              alpha = scen$alpha, 
              beta = scen$beta,
              threshold_method = scen$thresh
            )
          }
          if(m == "wavethresh") {
            extra = list(
              family = scen$wt,
              filter.number = if(!is.null(scen$wt_filter)) 
                scen$wt_filter else 4
            )
          }
          
          y_smooth = do.call(
            smooth_data_matrix, c(list(time_grid, y_wide, m), extra)
          )
          
          if(is.null(y_smooth)) next
          
          exec_time = attr(y_smooth, "exec_time")
          
          fate = estimate_fate_pointwise(y_smooth, X_mat, A_vec, time_grid)
          
          error = fate$fate_est - true_curve
          se_sq = error^2
          abs_err = abs(error)
          covered = (
            true_curve >= fate$ci_lower & true_curve <= fate$ci_upper
            )
          
          res_sim_long[[m]] = data.table(
            SimID = i,
            Scenario = scen$name,
            Method = m,
            t = time_grid,
            Estimate = fate$fate_est,
            TrueVal = true_curve,
            SignalRange = true_range,
            SE_sq = se_sq,
            AbsErr = abs_err,
            ExecTime = exec_time
          )
        }
        
        p()
        return(rbindlist(res_sim_long))
        
      }, future.seed = TRUE)
      
      # Salva Lote
      saveRDS(rbindlist(sim_results_list), batch_file)
      
      # Limpeza de RAM
      rm(sim_results_list)
      gc() 
      cat("OK.\n")
    }
  }
  
  cat("\n--- Agregando Resultados Finais ---\n")
  
  # Lê e agrega tudo
  files = list.files(
    ckpt_dir, pattern = "\\.rds$", full.names = TRUE
    )
  dt_all = rbindlist(lapply(files, readRDS))
  
  # Agregação por Tempo
  stats_scen = dt_all[, .(
    RMSE_t = sqrt(mean(SE_sq, na.rm = TRUE)),
    NRMSE_t = sqrt(mean(SE_sq, na.rm = TRUE)) / mean(SignalRange),
    MAE_t = mean(AbsErr, na.rm = TRUE),
    
    RMSE_sd = sd(sqrt(SE_sq), na.rm = TRUE),
    NRMSE_sd = sd(sqrt(SE_sq)/mean(SignalRange), na.rm = TRUE),
    MAE_sd = sd(AbsErr, na.rm = TRUE),
    
    Mean_Estimate = mean(Estimate, na.rm=TRUE),
    True_Effect = mean(TrueVal, na.rm=TRUE)
  ), by = .(Scenario, Method, t)]
  
  # Agregação Global de Tempo
  time_stats = dt_all[, .(
    Time_Mean = mean(unique(ExecTime), na.rm=TRUE),
    Time_SD = sd(unique(ExecTime), na.rm=TRUE)
  ), by = .(Scenario, Method)]
  
  stats_final = merge(stats_scen, time_stats, by = c("Scenario", "Method"))
  
  # Salva com N_SUBJ e T_PTS no nome e pasta
  out_dir = sprintf("output/sim1_N%d_T%d", N_SUBJ, T_PTS)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(
    stats_final,
    file.path(
      out_dir, 
      sprintf("sim1_time_series_results_N%d_T%d.csv", N_SUBJ, T_PTS))
    )
  
  # Remove temporários apenas deste batch
  unlink(ckpt_dir, recursive = TRUE)
  
  return(stats_final)
}

# Executar
run_full_simulation()