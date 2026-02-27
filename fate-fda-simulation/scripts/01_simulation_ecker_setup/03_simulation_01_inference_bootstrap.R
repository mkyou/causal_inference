library(data.table)
library(future.apply)
library(progressr)
library(rLifting)
library(fda)
library(kableExtra)

source("R/scenario_runner.R")
source("R/process_functional_data.R")
source("R/estimate_fate.R")

# Configuração

QUICK_TEST = FALSE
RUN_PARALLEL = TRUE

if(QUICK_TEST) {
  CONFIG = list(
    N_SIM = 100, N_BOOT = 200, N_SUBJ = 500, 
    T_PTS = 128, BATCH_SIZE = 20
    )
} else {
  CONFIG = list(
    N_SIM = 1000, N_BOOT = 200, N_SUBJ = 200, 
    T_PTS = 128, BATCH_SIZE = 35
    )
}

SCENARIOS = list(
  list(name = "Sigmoid", mag = 2, family = "cdf97"),
  list(name = "Bumps", mag = 3, family = "db2"),
  list(name = "Doppler", mag = 5, family = "db2"),
  list(name = "Heavisine", mag = 4, family = "haar")
)

METHODS = c("hybrid", "ecker_spline")

dir.create("output/checkpoints_ic", recursive = TRUE, showWarnings = FALSE)
dir.create("output/sim1_artifacts", recursive = TRUE, showWarnings = FALSE)

n_cores = 4
if(RUN_PARALLEL) {
  plan(multisession, workers = n_cores)
} else {
    plan(sequential)
}

handlers(global = TRUE)


# Funções para Bootstrap
run_raw_bootstrap = function(y_wide, X_mat, A_vec, time_grid, n_boot) {
  
  n_subj = nrow(y_wide)
  boot_curves = matrix(NA_real_, nrow = n_boot, ncol = ncol(y_wide))
  
  for(b in 1:n_boot) {
    idx = sample(1:n_subj, n_subj, replace = TRUE)
    
    if(sum(A_vec[idx] == 1) < 15 || sum(A_vec[idx] == 0) < 15) next
    
    fate = estimate_fate_pointwise(
      y_wide[idx,,drop=FALSE], X_mat[idx,,drop=FALSE], 
      A_vec[idx], time_grid
    )
    if(!is.null(fate)) boot_curves[b, ] = fate$fate_est
  }
  return(boot_curves)
}

smooth_spline = function(y_wide, time_grid) {
  basis = create.bspline.basis(
    range(time_grid), 
    nbasis = min(50, length(time_grid)/4)
    )
  
  tryCatch(
    t(
      eval.fd(
        time_grid, 
        smooth.basis(
          time_grid, t(y_wide), 
          fdPar(basis, 2, 1)
          )$fd
        )
      ), 
    error = function(e) y_wide
    )
}

run_ecker_bootstrap = function(y_wide, X_mat, A_vec, time_grid, n_boot) {
  n_subj = nrow(y_wide)
  boot_curves = matrix(NA_real_, nrow = n_boot, ncol = ncol(y_wide))
  
  for(b in 1:n_boot) {
    idx = sample(1:n_subj, n_subj, replace = TRUE)
    if(sum(A_vec[idx] == 1) < 15 || sum(A_vec[idx] == 0) < 15) next
    y_smooth = smooth_spline(y_wide[idx,,drop=FALSE], time_grid)
    
    fate = estimate_fate_pointwise(
      y_smooth, X_mat[idx,,drop=FALSE], 
      A_vec[idx], time_grid
      )
    
    if(!is.null(fate)) boot_curves[b, ] = fate$fate_est
  }
  return(boot_curves)
}

build_band_hybrid = function(theta_hat, boot_curves, alpha = 0.05) {
  valid = complete.cases(boot_curves)
  
  if(sum(valid) < 50) 
    return(
      list(
        lower = rep(-Inf, length(theta_hat)), 
        upper = rep(Inf, length(theta_hat))
        )
      )
  
  boot_valid = boot_curves[valid, , drop = FALSE]
  deviations = apply(
    boot_valid, 1, 
    function(curve) 
      max(abs(curve - theta_hat), na.rm = TRUE)
    )
  q = quantile(deviations, 1 - alpha, na.rm = TRUE)
  list(lower = theta_hat - q, upper = theta_hat + q)
}

build_pointwise_band = function(boot_curves, alpha = 0.05) {
  valid = complete.cases(boot_curves)
  
  if(sum(valid) < 50) 
    return(
      list(
        lower = rep(-Inf, ncol(boot_curves)), 
        upper = rep(Inf, ncol(boot_curves))
        )
      )
  boot_valid = boot_curves[valid, , drop = FALSE]
  
  lower = apply(boot_valid, 2, quantile, probs = alpha/2, na.rm = TRUE)
  upper = apply(boot_valid, 2, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  list(lower = lower, upper = upper)
}

# Worker
run_single_sim = function(sim_id, scen, method, config) {
  dgp_data = generate_single_scenario_dataset(
    scen$name, scen$mag, 
    config$N_SUBJ, config$T_PTS
    )
  
  true_curve = attr(dgp_data, "true_effect_curve")
  true_range = diff(range(true_curve)); if(true_range == 0) true_range = 1.0
  time_grid = attr(dgp_data, "time_grid")
  X_mat = attr(dgp_data, "X_mat")
  A_vec = attr(dgp_data, "A_vec")
  true_scalar = mean(true_curve)
  
  y_wide = as.matrix(dcast(dgp_data, id ~ t, value.var = "y")[, -1])
  rm(dgp_data); gc(verbose = FALSE)
  
  if(method == "hybrid") {
    y_smooth = smooth_data_matrix(
      time_grid, y_wide, "rlifting_offline",
      scheme_name = scen$family, threshold_method = "semisoft", alpha = 0.3
      )
    
    if(is.null(y_smooth)) return(NULL)
    fate = estimate_fate_pointwise(y_smooth, X_mat, A_vec, time_grid)
    
    if(is.null(fate)) return(NULL)
    theta_hat = fate$fate_est
    
    boot_curves = run_raw_bootstrap(
      y_wide, X_mat, A_vec, time_grid, config$N_BOOT
      )
    
  } else if(method == "ecker_spline") {
    y_smooth = smooth_spline(y_wide, time_grid)
    fate = estimate_fate_pointwise(y_smooth, X_mat, A_vec, time_grid)
    if(is.null(fate)) return(NULL)
    theta_hat = fate$fate_est
    boot_curves = run_ecker_bootstrap(
      y_wide, X_mat, A_vec, time_grid, config$N_BOOT
      )
  }
  
  band_simul = build_band_hybrid(theta_hat, boot_curves, 0.05)
  band_point = build_pointwise_band(boot_curves, 0.05)
  
  cov_simul = all(
    true_curve >= band_simul$lower & 
      true_curve <= band_simul$upper
    )
  
  cov_point = mean(
    true_curve >= band_point$lower & 
      true_curve <= band_point$upper
    )
  
  scalar_boots = rowMeans(boot_curves, na.rm = TRUE)
  scalar_boots = scalar_boots[!is.na(scalar_boots)]
  
  cov_scalar = if(length(scalar_boots) >= 50) {
    ci = quantile(scalar_boots, c(0.025, 0.975))
    as.integer(true_scalar >= ci[1] & true_scalar <= ci[2])
  } else NA_integer_
  
  nrmse = sqrt(mean((theta_hat - true_curve)^2)) / true_range
  width = mean(band_simul$upper - band_simul$lower, na.rm = TRUE)
  width_point = mean(band_point$upper - band_point$lower, na.rm = TRUE)
  
  data.table(
    SimID = sim_id, Scenario = scen$name, Method = method,
    Cov_Simul = as.integer(cov_simul), Cov_Point = cov_point,
    Cov_Scalar = cov_scalar, Width = width, Width_Point = width_point, 
    NRMSE = nrmse
    )
}

# Executor com Batch
run_inference_simulation = function() {
  tasks = as.data.table(
    expand.grid(
      SimID = 1:CONFIG$N_SIM, ScenIdx = seq_along(SCENARIOS),
      Method = METHODS, stringsAsFactors = FALSE
      )
    )
  
  n_batches = ceiling(nrow(tasks) / CONFIG$BATCH_SIZE)
  
  cat(
    sprintf(
      "\n>>> IC Simulation (N=%d, Boot=%d, Sims=%d, Batches=%d) <<<\n",
      CONFIG$N_SUBJ, CONFIG$N_BOOT, 
      CONFIG$N_SIM, n_batches
      )
    )
  
  ckpt_dir = file.path("output/checkpoints_ic", sprintf("N%d", CONFIG$N_SUBJ))
  dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)
  
  for(b in 1:n_batches) {
    batch_file = file.path(ckpt_dir, sprintf("batch_%04d.rds", b))
    if(file.exists(batch_file)) {
      cat(sprintf("  Batch %d/%d existe. Pulando.\n", b, n_batches))
      next
      }
    
    start_idx = (b - 1) * CONFIG$BATCH_SIZE + 1
    end_idx = min(b * CONFIG$BATCH_SIZE, nrow(tasks))
    batch_tasks = tasks[start_idx:end_idx]
    
    cat(
      sprintf(
        "  Batch %d/%d (tasks %d-%d)... ", b, n_batches, start_idx, end_idx
        )
      )
    
    batch_results = future_lapply(1:nrow(batch_tasks), function(i) {
      row = batch_tasks[i]
      scen = SCENARIOS[[row$ScenIdx]]
      
      tryCatch(
        run_single_sim(row$SimID, scen, row$Method, CONFIG), 
        error = function(e) NULL
        )
    }, future.seed = TRUE)
    
    batch_dt = rbindlist(batch_results[!sapply(batch_results, is.null)])
    if(nrow(batch_dt) > 0) {
      saveRDS(batch_dt, batch_file); cat(sprintf("OK (%d)\n", nrow(batch_dt)))
    }
    else {
      cat("VAZIO\n")
      }
    
    rm(batch_results, batch_dt); gc(verbose = FALSE)
  }
  
  cat("\n--- Agregando Resultados ---\n")
  files = list.files(
    ckpt_dir, 
    pattern = "\\.rds$", full.names = TRUE
    )
  if(length(files) == 0) { cat("ERRO: Sem checkpoints!\n"); return(NULL) }
  
  all_results = rbindlist(lapply(files, readRDS))
  # Output directory
  base_dir = sprintf("output/sim1_N%d_T%d/inference", CONFIG$N_SUBJ, CONFIG$T_PTS)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  
  csv_file = file.path(base_dir, sprintf(
    "inference_%s_N%d_T%d.csv", 
    ifelse(QUICK_TEST, "quick", "full"),
    CONFIG$N_SUBJ, CONFIG$T_PTS
    ))
  
  fwrite(all_results, csv_file)
  summary_dt = all_results[, .(N = .N,
    `Cov Simul (%)` = round(100 * mean(Cov_Simul, na.rm = TRUE), 1),
    `Cov Point (%)` = round(100 * mean(Cov_Point, na.rm = TRUE), 1),
    `Cov Scalar (%)` = round(100 * mean(Cov_Scalar, na.rm = TRUE), 1),
    `Width Simul` = round(mean(Width, na.rm = TRUE), 2),
    `Width Point` = round(mean(Width_Point, na.rm = TRUE), 2),
    NRMSE = round(mean(NRMSE, na.rm = TRUE), 4)
  ), by = .(Scenario, Method)]
  print(summary_dt)
  
  # Pretty table output
  kbl_table = kable(
    summary_dt, format = "html", booktabs = TRUE,
    col.names = c(
      "Cenário", "Método", "N", "Cov Simul (%)", "Cov Point (%)",
      "Cov Scalar (%)", "Width Simul", "Width Point", "NRMSE"
      )
    ) |>
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
  
  html_file = file.path(base_dir, sprintf(
    "table_%s_ics_pretty_N%d_T%d.html", 
    ifelse(QUICK_TEST, "quick", "full"),
    CONFIG$N_SUBJ, CONFIG$T_PTS
    ))
  save_kable(kbl_table, html_file)
  
  cat(sprintf("\nSalvo: %s\n", csv_file))
  cat(sprintf("Tabela pretty: %s\n", html_file))
  
  # Limpar checkpoints temporários
  unlink(ckpt_dir, recursive = TRUE)
  cat("Checkpoints removidos.\n")
  
  return(all_results)
}

results = run_inference_simulation()