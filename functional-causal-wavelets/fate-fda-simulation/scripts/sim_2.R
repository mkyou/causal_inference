# =============================================================================
# scripts/sim_2.R
# Simulação 2: Latência do estimador online (streaming)
# Grid: idêntico à sim_1 (2 dists × 4 cenários × 3N × 3T = 72 configs)
# K = 100 réplicas, execução sequencial
# Mede: latência de denoising por passo, latência da etapa FATE, NRMSE final
# =============================================================================

library(data.table)
library(rLifting)
library(MASS)

script_dir = tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) {
    args = commandArgs(trailingOnly = FALSE)
    file_arg = grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    } else {
      getwd()
    }
  }
)
if (basename(script_dir) == "scripts") {
  setwd(dirname(script_dir))
} else {
  setwd(script_dir)
}

source("R/config_ecker.R")
source("R/generate_sequential_dgp.R")
source("R/generate_ecker_dgp.R")
source("R/estimate_fate_optimized.R")
source("R/utils_dgp_functions.R")

# =============================================================================
# CONFIGURAÇÃO
# =============================================================================

QUICK_TEST = FALSE
N_SIM = ifelse(QUICK_TEST, 5, 100)
BATCH_SIZE = 25

t_start = 20.0
t_end = 80.0
t_treat = 50.0

GRID_N = if (QUICK_TEST) c(200) else c(200, 500, 2000)
GRID_T = if (QUICK_TEST) c(128) else c(128, 256, 512)
GRID_W = if (QUICK_TEST) c(64) else c(64, 64, 128)
GRID_DIST = if (QUICK_TEST) "gaussian" else c("gaussian", "t")

SCENARIOS = list(
  list(name = "Sigmoid", func = effect_sigmoid, mag = 2, family = "haar"),
  list(name = "Bumps", func = effect_bumps, mag = 3, family = "db2"),
  list(name = "Doppler", func = effect_doppler, mag = 5, family = "db2"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

# =============================================================================
# WORKER: Uma réplica de streaming com medição de latência
# =============================================================================

run_single_sim_online = function(
  sim_id, scen, N, T_pts, time_grid, L_chol, W, error_dist
) {
  dat = generate_ecker_dataset(
    scen, N, T_pts, time_grid, L_chol,
    t_start = t_start, t_end = t_end, t_treat = t_treat,
    error_dist = error_dist
  )

  true_curve = dat$true_curve
  idx_post = which(time_grid >= t_treat)
  true_range = diff(range(true_curve[idx_post]))
  if (true_range == 0) true_range = 1

  proj = precompute_projection(dat$X_reg, dat$A_vec)

  # Setup streaming: N processadores independentes
  scheme_obj = lifting_scheme(scen$family)
  max_levels = floor(log2(W))
  processors = lapply(seq_len(N), function(i) {
    new_wavelet_stream(
      scheme_obj,
      window_size = W, levels = max_levels,
      shrinkage = "semisoft", alpha = 0.3
    )
  })

  Y_online = matrix(0, N, T_pts)
  latency_denoise_us = numeric(T_pts)

  # Processa cada instante t em sequência estrita
  for (k in seq_len(T_pts)) {
    t0 = proc.time()[["elapsed"]]
    for (i in seq_len(N)) {
      Y_online[i, k] = processors[[i]](dat$Y_wide[i, k])
    }
    # Latência média por indivíduo neste instante (µs)
    latency_denoise_us[k] = (proc.time()[["elapsed"]] - t0) / N * 1e6
  }

  # Etapa FATE: só após todos os T instantes (idêntica ao modo causal)
  t0 = proc.time()[["elapsed"]]
  fate_final = estimate_fate_fast(Y_online, proj)
  latency_fate_us = (proc.time()[["elapsed"]] - t0) * 1e6

  post_diff = fate_final$fate_est[idx_post] - true_curve[idx_post]
  nrmse_final = sqrt(mean(post_diff^2)) / true_range

  metrics_dt = data.table(
    SimID = sim_id, Scenario = scen$name, ErrorDist = error_dist,
    N = N, T_pts = T_pts,
    NRMSE_final = round(nrmse_final, 6),
    Latency_denoise_mean_us = round(mean(latency_denoise_us), 4),
    Latency_fate_us = round(latency_fate_us, 4)
  )

  latency_dt = data.table(
    SimID = sim_id, Scenario = scen$name, ErrorDist = error_dist,
    N = N, T_pts = T_pts,
    t_val = time_grid,
    Latency_denoise_us = latency_denoise_us
  )

  rm(dat, Y_online, processors, fate_final)
  gc(verbose = FALSE)

  list(metrics = metrics_dt, latency = latency_dt)
}

# =============================================================================
# AGREGAÇÃO (dois passes — OOM-safe)
# =============================================================================

aggregate_simulation_2 = function() {
  cat("\n\n========== AGREGANDO SIM 2 ==========\n")
  config_dirs = list.dirs(
    "output/sim2_checkpoints",
    recursive = FALSE, full.names = TRUE
  )
  if (length(config_dirs) == 0) {
    cat("ERRO: Nenhum checkpoint encontrado.\n")
    return(invisible(NULL))
  }

  out_dir = "output/sim2_final"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Passe 1: métricas (leve) ---
  metrics_list = vector("list", length(config_dirs))
  for (i in seq_along(config_dirs)) {
    batch_files = list.files(
      config_dirs[i],
      pattern = "\\.rds$",
      full.names = TRUE
    )
    metrics_list[[i]] = rbindlist(
      lapply(batch_files, function(f) readRDS(f)$metrics)
    )
  }
  all_metrics = rbindlist(metrics_list)

  summary_dt = all_metrics[, .(
    N_sims = .N,
    NRMSE_final = round(mean(NRMSE_final), 4),
    Latency_denoise_mean_us = round(mean(Latency_denoise_mean_us), 2),
    Latency_fate_us = round(mean(Latency_fate_us), 2)
  ), by = .(Scenario, ErrorDist, N, T_pts)]

  fwrite(all_metrics, file.path(out_dir, "sim2_raw.csv"))
  fwrite(summary_dt, file.path(out_dir, "sim2_summary.csv"))
  rm(all_metrics, metrics_list)
  gc(verbose = FALSE)

  # --- Passe 2: curvas de latência, config por config ---
  latency_summary_list = vector("list", length(config_dirs))
  for (i in seq_along(config_dirs)) {
    batch_files = list.files(
      config_dirs[i], pattern = "\\.rds$",
      full.names = TRUE
    )
    lat_i = rbindlist(lapply(batch_files, function(f) readRDS(f)$latency))
    latency_summary_list[[i]] = lat_i[, .(
      Latency_mean_us = round(mean(Latency_denoise_us), 4),
      Latency_p95_us = round(quantile(Latency_denoise_us, 0.95), 4)
    ), by = .(Scenario, ErrorDist, N, T_pts, t_val)]
    rm(lat_i)
    gc(verbose = FALSE)
  }
  fwrite(
    rbindlist(latency_summary_list),
    file.path(out_dir, "sim2_latency_curves.csv")
  )

  cat(sprintf("Salvo em: %s/\n", out_dir))
  cat("\nSumário:\n")
  print(summary_dt)
  invisible(summary_dt)
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_2 = function() {
  for (error_dist in GRID_DIST) {
    for (T_idx in seq_along(GRID_T)) {
      T_pts = GRID_T[T_idx]
      W = GRID_W[T_idx]
      time_grid = seq(t_start, t_end, length.out = T_pts)

      cat(
        sprintf(
          "\nPré-computando Cholesky MGP para T=%d, dist=%s...\n", 
          T_pts, error_dist
        )
      )
      L_chol = make_mgp_chol(time_grid, sigma = 0.5, rho = 0.8)

      for (N in GRID_N) {
        for (scen in SCENARIOS) {
          config_label = sprintf(
            "%s_N%d_T%d_%s", error_dist,
            N, T_pts, scen$name
          )
          ckpt_dir = file.path("output", "sim2_checkpoints", config_label)
          dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

          n_batches = ceiling(N_SIM / BATCH_SIZE)
          cat(sprintf("\n>>> %s (Batches=%d) <<<\n", config_label, n_batches))

          for (b in seq_len(n_batches)) {
            batch_file = file.path(ckpt_dir, sprintf("batch_%04d.rds", b))
            if (file.exists(batch_file)) {
              cat(sprintf("  Batch %d/%d existe. Pulando.\n", b, n_batches))
              next
            }

            start_idx = (b - 1) * BATCH_SIZE + 1
            end_idx = min(b * BATCH_SIZE, N_SIM)
            batch_ids = start_idx:end_idx

            cat(sprintf(
              "  Batch %d/%d (sims %d-%d)... ", b, n_batches, start_idx, end_idx
            ))

            batch_results = lapply(batch_ids, function(sid) {
              tryCatch(
                run_single_sim_online(
                  sid, scen, N, T_pts, time_grid, L_chol, W, error_dist
                ),
                error = function(e) {
                  message(sprintf("Sim %d erro: %s", sid, e$message))
                  NULL
                }
              )
            })
            batch_results = Filter(Negate(is.null), batch_results)

            if (length(batch_results) > 0) {
              saveRDS(
                list(
                  metrics = rbindlist(lapply(batch_results, `[[`, "metrics")),
                  latency = rbindlist(lapply(batch_results, `[[`, "latency"))
                ),
                batch_file
              )
              cat(sprintf("OK (%d sims)\n", length(batch_results)))
            } else {
              cat("VAZIO\n")
            }

            rm(batch_results)
            gc(verbose = FALSE)
          }
        }
      }
    }
  }

  aggregate_simulation_2()
}

# =============================================================================
# RUN
# =============================================================================

results = run_simulation_2()
