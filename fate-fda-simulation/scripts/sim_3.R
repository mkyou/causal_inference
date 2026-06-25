# =============================================================================
# scripts/sim_3.R
# Simulação 3: Robustez a grade irregular (uniform vs. irregular)
# Métodos: wavelet_offline, wavelet_causal (ambos via rLifting com t=time_grid)
# Grid: N ∈ {200, 500}, T ∈ {128, 256}, 4 cenários, 2 ruídos, 1000 sims
# =============================================================================

library(data.table)
library(parallel)
library(rLifting)
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
source("R/process_functional_data.R")
source("R/utils_dgp_functions.R")

# =============================================================================
# CONFIGURAÇÃO
# =============================================================================

QUICK_TEST = TRUE

N_SIM = ifelse(QUICK_TEST, 1, 1000)
BATCH_SIZE = ifelse(QUICK_TEST, 1, 25)
N_CORES = max(1, detectCores() - 2)
ALPHA = 0.05

t_start = 20.0
t_end = 80.0
t_treat = 50.0

GRID_N = if (QUICK_TEST) c(500) else c(200, 500)
GRID_T = if (QUICK_TEST) c(256) else c(128, 256)
GRID_W = if (QUICK_TEST) c(64) else c(64, 64)

SCENARIOS = list(
  list(name = "Sigmoid",   func = effect_sigmoid,   mag = 2, family = "haar"),
  list(name = "Bumps",     func = effect_bumps,     mag = 3, family = "haar"),
  list(name = "Doppler",   func = effect_doppler,   mag = 5, family = "haar"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

NOISE_DISTS = c("gaussian", "t")
DF_T = 5
GRID_TYPES = c("uniform", "irregular")

safe_cores = function(N, T_pts, mem_safety = 0.75) {
  data_bytes = 12 * N * T_pts * 8
  fork_base = 400e6
  bytes_per_worker = data_bytes + fork_base
  avail_bytes = tryCatch(
    as.numeric(
      system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)
    ) * 1024,
    error = function(e) 8e9
  )
  max(1L, min(N_CORES, floor(avail_bytes * mem_safety / bytes_per_worker)))
}

# =============================================================================
# WORKER
# =============================================================================

run_single_sim_3 = function(
    sim_id, scen, N, T_pts, time_grid, L_chol,
    W, n_boot, error_dist = "gaussian", grid_type = "uniform"
) {
  dat = generate_ecker_dataset(
    scen, N, T_pts, time_grid, L_chol,
    t_start = t_start, t_end = t_end, t_treat = t_treat,
    error_dist = error_dist, df_t = DF_T
  )
  true_curve = dat$true_curve
  idx_post = time_grid >= t_treat
  true_range = diff(range(true_curve[idx_post]))
  if (true_range == 0) true_range = 1
  raw_est = colMeans(dat$Y_wide[dat$A_vec == 1, ]) -
    colMeans(dat$Y_wide[dat$A_vec == 0, ])
  cfate_true = cfate_cumulative(true_curve, time_grid, t_treat)

  proj = precompute_projection(dat$X_reg, dat$A_vec)

  # Wavelet offline — passa time_grid (suporta grade irregular)
  y_off = tryCatch(
    smooth_data_matrix(
      time_grid, dat$Y_wide, "rlifting_offline",
      scheme_name = scen$family, threshold_method = "semisoft", alpha = 0.3
    ),
    error = function(e) NULL
  )

  # Wavelet causal — passa t = time_grid (suporta grade irregular)
  max_levels = floor(log2(W))
  scheme_obj = lifting_scheme(scen$family)
  y_cau = tryCatch({
    t(apply(dat$Y_wide, 1, function(sig) {
      denoise_signal_causal(
        sig, scheme_obj, levels = max_levels, window_size = W,
        method = "semisoft", alpha = 0.3, t = time_grid
      )
    }))
  }, error = function(e) NULL)

  fate_off = if (!is.null(y_off)) estimate_fate_fast(y_off, proj) else NULL
  fate_cau = if (!is.null(y_cau)) estimate_fate_fast(y_cau, proj) else NULL

  boot_curves = run_bootstrap_raw(dat$Y_wide, dat$X_reg, dat$A_vec, n_boot)

  compute_metrics = function(fate_est, method_name) {
    if (is.null(fate_est)) return(list(metrics = NULL, curves = NULL))
    band = build_simultaneous_band(fate_est$fate_est, boot_curves, ALPHA)
    cov_simul = as.integer(all(true_curve >= band$lower & true_curve <= band$upper))
    cov_point = mean(true_curve >= band$lower & true_curve <= band$upper)
    post_diff = fate_est$fate_est[idx_post] - true_curve[idx_post]
    nrmse = sqrt(mean(post_diff^2)) / true_range
    bias = mean(post_diff)
    width = mean(band$upper - band$lower, na.rm = TRUE)

    cfate_est_vec = cfate_cumulative(fate_est$fate_est, time_grid, t_treat)
    cfate_lower = numeric(length(time_grid))
    cfate_upper = numeric(length(time_grid))
    valid_boot = which(complete.cases(boot_curves))
    if (length(valid_boot) >= 50) {
      cfate_boot_mat = t(apply(
        boot_curves[valid_boot, , drop = FALSE], 1,
        function(r) cfate_cumulative(r, time_grid, t_treat)
      ))
      cb = build_simultaneous_band(
        cfate_est_vec[idx_post],
        cfate_boot_mat[, idx_post, drop = FALSE],
        ALPHA
      )
      cfate_lower[idx_post] = cb$lower
      cfate_upper[idx_post] = cb$upper
    }

    metrics_dt = data.table(
      SimID = sim_id, GridType = grid_type, ErrorDist = error_dist,
      Scenario = scen$name, Method = method_name,
      N = N, T_pts = T_pts,
      Cov_Simul = cov_simul, Cov_Point = round(cov_point, 4),
      NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4)
    )
    curves_dt = data.table(
      SimID = sim_id, GridType = grid_type, ErrorDist = error_dist,
      Scenario = scen$name, Method = method_name,
      N = N, T_pts = T_pts, t_val = time_grid,
      TrueCurve = true_curve, RawEst = raw_est,
      Est = fate_est$fate_est, Lower = band$lower, Upper = band$upper,
      CFATE_True = cfate_true, CFATE_Est = cfate_est_vec,
      CFATE_Lower = cfate_lower, CFATE_Upper = cfate_upper
    )
    list(metrics = metrics_dt, curves = curves_dt)
  }

  res_off = compute_metrics(fate_off, "wavelet_offline")
  res_cau = compute_metrics(fate_cau, "wavelet_causal")

  rm(dat, y_off, y_cau, boot_curves); gc(verbose = FALSE)

  all_metrics = list(res_off$metrics, res_cau$metrics)
  all_curves = list(res_off$curves, res_cau$curves)
  list(
    metrics = rbindlist(all_metrics[!sapply(all_metrics, is.null)]),
    curves  = rbindlist(all_curves[!sapply(all_curves,  is.null)])
  )
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_3 = function() {
  sim_tag = if (QUICK_TEST) "sim3_test" else "sim3"
  ckpt_tag = paste0(sim_tag, "_checkpoints")
  out_tag = paste0(sim_tag, "_final")

  for (grid_type in GRID_TYPES) {
  for (error_dist in NOISE_DISTS) {
  for (T_idx in seq_along(GRID_T)) {
    T_pts = GRID_T[T_idx]
    W = GRID_W[T_idx]

    set.seed(42 + T_idx)
    time_grid = if (grid_type == "uniform") {
      seq(t_start, t_end, length.out = T_pts)
    } else {
      sort(runif(T_pts, t_start, t_end))
    }

    cat(sprintf("Pré-computando Cholesky MGP para T=%d (%s)...\n", T_pts, grid_type))
    L_chol = make_mgp_chol(time_grid, sigma = 0.5, rho = 0.8)

    for (N in GRID_N) {
      n_boot = N

      for (scen in SCENARIOS) {
        config_label = sprintf(
          "%s_%s_N%d_T%d_%s", grid_type, error_dist, N, T_pts, scen$name
        )
        ckpt_dir = file.path("output", ckpt_tag, config_label)
        dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

        sim_ids = 1:N_SIM
        n_batches = ceiling(N_SIM / BATCH_SIZE)
        n_cores_this = safe_cores(N, T_pts)

        cat(sprintf(
          "\n>>> %s (N_BOOT=%d, Batches=%d, cores=%d) <<<\n",
          config_label, n_boot, n_batches, n_cores_this
        ))

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

          batch_results = mclapply(batch_ids, function(sid) {
            tryCatch(
              run_single_sim_3(
                sid, scen, N, T_pts, time_grid, L_chol, W, n_boot, error_dist, grid_type
              ),
              error = function(e) {
                message(sprintf("Sim %d erro: %s", sid, e$message)); NULL
              }
            )
          }, mc.cores = n_cores_this)

          batch_dt = rbindlist(lapply(batch_results, function(x) x$metrics))
          batch_curves = rbindlist(lapply(batch_results, function(x) x$curves))

          if (nrow(batch_dt) > 0) {
            saveRDS(list(metrics = batch_dt, curves = batch_curves), batch_file)
            cat(sprintf(
              "OK (%d rows metrics, %d rows curves)\n",
              nrow(batch_dt), nrow(batch_curves)
            ))
          } else {
            cat("VAZIO\n")
          }

          rm(batch_results, batch_dt, batch_curves); gc(verbose = FALSE)
        }
      }
    }
  }
  } # error_dist
  } # grid_type

  # --- Agregação final ---
  cat("\n\n========== AGREGANDO RESULTADOS ==========\n")

  all_files = list.files(
    file.path("output", ckpt_tag), pattern = "\\.rds$",
    recursive = TRUE, full.names = TRUE
  )
  if (length(all_files) == 0) { cat("ERRO: Sem checkpoints!\n"); return(NULL) }

  all_data = lapply(all_files, readRDS)
  all_results = rbindlist(lapply(all_data, function(x) x$metrics))
  all_curves  = rbindlist(lapply(all_data, function(x) x$curves))

  summary_dt = all_results[, .(
    N_sims = .N,
    Cov_Simul = round(mean(Cov_Simul, na.rm = TRUE) * 100, 1),
    Cov_Point = round(mean(Cov_Point, na.rm = TRUE) * 100, 1),
    NRMSE = round(mean(NRMSE, na.rm = TRUE), 4),
    Bias = round(mean(Bias, na.rm = TRUE), 4),
    Width = round(mean(Width, na.rm = TRUE), 3)
  ), by = .(GridType, ErrorDist, Scenario, Method, N, T_pts)]


  out_dir = file.path("output", out_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(all_results, file.path(out_dir, paste0(sim_tag, "_raw.csv")))
  fwrite(summary_dt,  file.path(out_dir, paste0(sim_tag, "_summary.csv")))
  fwrite(all_curves,  file.path(out_dir, paste0(sim_tag, "_curves.csv")))

  cat("\nSumário:\n")
  print(summary_dt)
  cat(sprintf("\nSalvo em: %s/\n", out_dir))

  return(summary_dt)
}

# =============================================================================
# RUN
# =============================================================================
results = run_simulation_3()
