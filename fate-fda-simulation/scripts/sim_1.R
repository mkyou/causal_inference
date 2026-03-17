# =============================================================================
# scripts/sim_1.R
# Simulação 1: Comparação Offline (Wavelet Offline vs Wavelet Causal vs Spline)
# Grid: N ∈ {200, 500, 2000}, T ∈ {128, 256, 512}, 4 cenários, 1000 sims
# =============================================================================

library(data.table)
library(future.apply)
library(rLifting)
library(fda)
library(MASS)

# --- Setup portável ---
script_dir = tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) {
    args = commandArgs(trailingOnly = FALSE)
    file_arg = grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg)))
    else getwd()
  }
)
if (basename(script_dir) == "scripts") setwd(dirname(script_dir)) else setwd(script_dir)

source("R/config_ecker.R")
source("R/generate_sequential_dgp.R")
source("R/generate_ecker_dgp.R")
source("R/estimate_fate_optimized.R")
source("R/ecker_analytical_bands.R")
source("R/process_functional_data.R")
source("R/utils_dgp_functions.R")

# =============================================================================
# CONFIGURAÇÃO
# =============================================================================

QUICK_TEST = FALSE  # TRUE = 10 sims para depuração

N_SIM = ifelse(QUICK_TEST, 10, 1000)
BATCH_SIZE = ifelse(QUICK_TEST, 5, 25)
N_CORES = 4
ALPHA = 0.05

t_start = 20.0
t_end = 63.0
t_treat = 20.0

# Grid de cenários
GRID_N = c(200, 500, 2000)
GRID_T = c(128, 256, 512)
GRID_W = c(64,  64,  128)   # W para wavelet causal, indexado por T

# Cenários de efeito (ondaleta padronizada por cenário)
SCENARIOS = list(
  list(name = "Sigmoid",   func = effect_sigmoid,   mag = 2, family = "haar"),
  list(name = "Bumps",     func = effect_bumps,     mag = 3, family = "db2"),
  list(name = "Doppler",   func = effect_doppler,   mag = 5, family = "db2"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

# =============================================================================
# WORKER: executa uma simulação para um cenário + N + T
# =============================================================================

run_single_sim = function(
    sim_id, scen, N, T_pts, time_grid, L_chol,
    W, n_boot, scheme
) {
  # 1. Gera dataset
  dat = generate_ecker_dataset(scen, N, T_pts, time_grid, L_chol,
                                 t_start = t_start, t_end = t_end, t_treat = t_treat)
  true_curve = dat$true_curve
  true_range = diff(range(true_curve))
  if (true_range == 0) true_range = 1

  # 2. Projeções pré-computadas (compartilhadas entre os 3 métodos)
  proj = precompute_projection(dat$X_reg, dat$A_vec)

  # 3. Suavização
  # Wavelet offline
  y_off = tryCatch(
    smooth_data_matrix(time_grid, dat$Y_wide, "rlifting_offline",
                       scheme_name = scen$family, threshold_method = "semisoft", alpha = 0.3),
    error = function(e) NULL
  )

  # Wavelet causal (batch)
  max_levels = floor(log2(W))
  scheme_obj = lifting_scheme(scen$family)
  y_cau = tryCatch({
    t(apply(dat$Y_wide, 1, function(sig) {
      denoise_signal_causal(sig, scheme_obj, levels = max_levels, window_size = W,
                            method = "semisoft", alpha = 0.3)
    }))
  }, error = function(e) NULL)

  # 4. Estimação pontual
  fate_off = if (!is.null(y_off)) estimate_fate_fast(y_off, proj) else NULL
  fate_cau = if (!is.null(y_cau)) estimate_fate_fast(y_cau, proj) else NULL

  # 5. Bootstrap compartilhado (dados ruidosos — sem suavização)
  boot_curves = run_bootstrap_raw(dat$Y_wide, dat$X_reg, dat$A_vec, n_boot)

  # 6. Bandas bootstrap (usadas por wavelet offline E causal)
  compute_metrics = function(fate_est, method_name) {
    if (is.null(fate_est)) return(list(metrics = NULL, curves = NULL))
    band = build_simultaneous_band(fate_est$fate_est, boot_curves, ALPHA)
    cov_simul = as.integer(all(true_curve >= band$lower & true_curve <= band$upper))
    cov_point = mean(true_curve >= band$lower & true_curve <= band$upper)
    nrmse = sqrt(mean((fate_est$fate_est - true_curve)^2)) / true_range
    bias = mean(fate_est$fate_est - true_curve)
    width = mean(band$upper - band$lower, na.rm = TRUE)

    metrics_dt = data.table(
      SimID = sim_id, Scenario = scen$name, Method = method_name,
      N = N, T_pts = T_pts,
      Cov_Simul = cov_simul, Cov_Point = round(cov_point, 4),
      NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4)
    )

    curves_dt = data.table(
      SimID = sim_id, Scenario = scen$name, Method = method_name,
      N = N, T_pts = T_pts, t_val = time_grid,
      TrueCurve = true_curve, Est = fate_est$fate_est,
      Lower = band$lower, Upper = band$upper
    )

    list(metrics = metrics_dt, curves = curves_dt)
  }

  res_off = compute_metrics(fate_off, "wavelet_offline")
  res_cau = compute_metrics(fate_cau, "wavelet_causal")

  # 7. Splines Ecker (bandas analíticas)
  ecker_res = tryCatch(
    estimate_fate_ecker(dat$Y_wide, dat$X_reg, dat$A_vec, time_grid,
                         nbasis = 50, alpha = ALPHA, n_sim = 5000),
    error = function(e) NULL
  )

  res_spl = list(metrics = NULL, curves = NULL)
  if (!is.null(ecker_res)) {
    cov_simul = as.integer(all(true_curve >= ecker_res$lower & true_curve <= ecker_res$upper))
    cov_point = mean(true_curve >= ecker_res$lower & true_curve <= ecker_res$upper)
    nrmse = sqrt(mean((ecker_res$fate_est - true_curve)^2)) / true_range
    bias = mean(ecker_res$fate_est - true_curve)
    width = mean(ecker_res$upper - ecker_res$lower, na.rm = TRUE)

    res_spl$metrics = data.table(
      SimID = sim_id, Scenario = scen$name, Method = "spline_ecker",
      N = N, T_pts = T_pts,
      Cov_Simul = cov_simul, Cov_Point = round(cov_point, 4),
      NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4)
    )

    res_spl$curves = data.table(
      SimID = sim_id, Scenario = scen$name, Method = "spline_ecker",
      N = N, T_pts = T_pts, t_val = time_grid,
      TrueCurve = true_curve, Est = ecker_res$fate_est,
      Lower = ecker_res$lower, Upper = ecker_res$upper
    )
  }

  # Limpeza
  rm(dat, y_off, y_cau, boot_curves); gc(verbose = FALSE)

  list(
    metrics = rbindlist(list(res_off$metrics, res_cau$metrics, res_spl$metrics)[!sapply(list(res_off$metrics, res_cau$metrics, res_spl$metrics), is.null)]),
    curves = rbindlist(list(res_off$curves,  res_cau$curves,  res_spl$curves)[!sapply(list(res_off$curves, res_cau$curves, res_spl$curves), is.null)])
  )
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_1 = function() {

  for (T_idx in seq_along(GRID_T)) {
    T_pts = GRID_T[T_idx]
    W = GRID_W[T_idx]
    time_grid = seq(t_start, t_end, length.out = T_pts)

    # Cholesky pré-computado para este T
    cat(sprintf("Pré-computando Cholesky MGP para T=%d...\n", T_pts))
    L_chol = make_mgp_chol(time_grid, sigma = 0.5, rho = 0.8)

    for (N in GRID_N) {
      n_boot = N  # N_BOOT = N

      for (scen in SCENARIOS) {
        config_label = sprintf("N%d_T%d_%s", N, T_pts, scen$name)
        ckpt_dir = file.path("output", "sim1_checkpoints", config_label)
        dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

        # Tarefas
        sim_ids = 1:N_SIM
        n_batches = ceiling(N_SIM / BATCH_SIZE)

        cat(sprintf("\n>>> %s (N_BOOT=%d, Batches=%d) <<<\n", config_label, n_boot, n_batches))

        for (b in seq_len(n_batches)) {
          batch_file = file.path(ckpt_dir, sprintf("batch_%04d.rds", b))
          if (file.exists(batch_file)) {
            cat(sprintf("  Batch %d/%d existe. Pulando.\n", b, n_batches))
            next
          }

          start_idx = (b - 1) * BATCH_SIZE + 1
          end_idx = min(b * BATCH_SIZE, N_SIM)
          batch_ids = sim_ids[start_idx:end_idx]

          cat(sprintf("  Batch %d/%d (sims %d-%d)... ", b, n_batches, start_idx, end_idx))

          batch_results = lapply(batch_ids, function(sid) {
            tryCatch(
              run_single_sim(sid, scen, N, T_pts, time_grid, L_chol, W, n_boot, NULL),
              error = function(e) { message(sprintf("Sim %d erro: %s", sid, e$message)); NULL }
            )
          })

          batch_dt = rbindlist(lapply(batch_results, function(x) x$metrics))
          batch_curves = rbindlist(lapply(batch_results, function(x) x$curves))

          if (nrow(batch_dt) > 0) {
            saveRDS(list(metrics = batch_dt, curves = batch_curves), batch_file)
            cat(sprintf("OK (%d rows metrics, %d rows curves)\n", nrow(batch_dt), nrow(batch_curves)))
          } else {
            cat("VAZIO\n")
          }

          rm(batch_results, batch_dt, batch_curves); gc(verbose = FALSE)
        }
      }
    }
  }

  # --- Agregação final ---
  cat("\n\n========== AGREGANDO RESULTADOS ==========\n")

  all_files = list.files("output/sim1_checkpoints", pattern = "\\.rds$",
                           recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) { cat("ERRO: Sem checkpoints!\n"); return(NULL) }

  all_data = lapply(all_files, readRDS)
  all_results = rbindlist(lapply(all_data, function(x) x$metrics))
  all_curves = rbindlist(lapply(all_data, function(x) x$curves))

  # Summary
  summary_dt = all_results[, .(
    N_sims = .N,
    Cov_Simul = round(mean(Cov_Simul, na.rm = TRUE) * 100, 1),
    Cov_Point = round(mean(Cov_Point, na.rm = TRUE) * 100, 1),
    NRMSE = round(mean(NRMSE, na.rm = TRUE), 4),
    Bias = round(mean(Bias, na.rm = TRUE), 4),
    Width = round(mean(Width, na.rm = TRUE), 3)
  ), by = .(Scenario, Method, N, T_pts)]

  out_dir = "output/sim1_final"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(all_results, file.path(out_dir, "sim1_raw.csv"))
  fwrite(summary_dt,  file.path(out_dir, "sim1_summary.csv"))
  fwrite(all_curves,  file.path(out_dir, "sim1_curves.csv"))

  cat("\nSumário:\n")
  print(summary_dt)
  cat(sprintf("\nSalvo em: %s/\n", out_dir))

  return(summary_dt)
}

# =============================================================================
# RUN
# =============================================================================
results = run_simulation_1()
