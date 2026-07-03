# =============================================================================
# scripts/sim_3.R
# Simulação 3: rLifting causal em grade irregular compartilhada
# Grid: 2 dists × 3 irr × 2N × 2T × 4 cenários = 96 configs
# K = 1000, mclapply por batch (mesmo padrão da sim_1)
# Grade irregular fixa por (irr_level, T_pts), compartilhada pelos N sujeitos
# =============================================================================

library(data.table)
library(parallel)
library(rLifting)
library(MASS)

script_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
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

QUICK_TEST <- FALSE
N_SIM <- ifelse(QUICK_TEST, 5, 1000)
BATCH_SIZE <- ifelse(QUICK_TEST, 5, 25)
ALPHA <- 0.05

t_start <- 20.0
t_end <- 80.0
t_treat <- 50.0

GRID_N <- if (QUICK_TEST) c(200) else c(200, 2000)
GRID_T <- if (QUICK_TEST) c(128) else c(128, 512)
GRID_W <- if (QUICK_TEST) c(64) else c(64, 128)
GRID_DIST <- if (QUICK_TEST) "gaussian" else c("gaussian", "t")
IRR_LEVELS <- if (QUICK_TEST) "moderado" else c("leve", "moderado", "severo")

SCENARIOS <- list(
  list(name = "Sigmoid", func = effect_sigmoid, mag = 2, family = "haar"),
  list(name = "Bumps", func = effect_bumps, mag = 3, family = "haar"),
  list(name = "Doppler", func = effect_doppler, mag = 5, family = "haar"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

N_CORES <- max(1L, detectCores() - 2L)

safe_cores <- function(N, T_pts, mem_safety = 0.75) {
  data_bytes <- 10 * N * T_pts * 8
  fork_base <- 400e6
  avail_bytes <- tryCatch(
    as.numeric(
      system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)
    ) * 1024,
    error = function(e) 8e9
  )
  max(
    1L,
    min(N_CORES, floor(avail_bytes * mem_safety / (data_bytes + fork_base)))
  )
}

# =============================================================================
# GRADES IRREGULARES (fixas por irr_level × T_pts, seed determinístico)
# =============================================================================

make_irregular_grid <- function(irr_level, T_pts, t_start, t_end, t_treat) {
  seed <- switch(irr_level,
    leve = 1001L,
    moderado = 1002L,
    severo = 1003L
  ) + T_pts
  set.seed(seed)
  switch(irr_level,
    leve = {
      reg <- seq(t_start, t_end, length.out = T_pts)
      delta <- (t_end - t_start) / (T_pts - 1) * 0.25
      g <- reg + runif(T_pts, -delta, delta)
      g[1] <- t_start
      g[T_pts] <- t_end
      sort(g)
    },
    moderado = sort(runif(T_pts, t_start, t_end)),
    severo = {
      n_pre <- round(2 / 3 * T_pts)
      n_post <- T_pts - n_pre
      sort(c(runif(n_pre, t_start, t_treat), runif(n_post, t_treat, t_end)))
    }
  )
}

GRIDS <- list()
for (irr in IRR_LEVELS) {
  for (T_pts in GRID_T) {
    GRIDS[[paste0(irr, "_T", T_pts)]] <- make_irregular_grid(
      irr, T_pts, t_start, t_end, t_treat
    )
  }
}

# =============================================================================
# WORKER: Uma réplica — rLifting causal + bootstrap
# =============================================================================

run_single_sim_s3 <- function(
  sim_id, scen, N, T_pts, time_grid, L_chol, W,
  error_dist, irr_level, n_boot
) {
  dat <- generate_ecker_dataset(
    scen, N, T_pts, time_grid, L_chol,
    t_start = t_start, t_end = t_end, t_treat = t_treat,
    error_dist = error_dist
  )
  true_curve <- dat$true_curve
  idx_post <- which(time_grid >= t_treat)
  true_range <- diff(range(true_curve[idx_post]))
  if (true_range == 0) true_range <- 1

  proj <- precompute_projection(dat$X_reg, dat$A_vec)

  scheme_obj <- lifting_scheme("haar")
  max_levels <- floor(log2(W))

  y_cau <- tryCatch(
    t(apply(dat$Y_wide, 1, function(sig) {
      denoise_signal_causal(
        sig, scheme_obj,
        levels = max_levels, window_size = W,
        shrinkage = "semisoft", alpha = 0.3
      )
    })),
    error = function(e) NULL
  )
  if (is.null(y_cau)) {
    return(NULL)
  }

  fate_cau <- estimate_fate_fast(y_cau, proj)
  boot_crv <- run_bootstrap_raw(dat$Y_wide, dat$X_reg, dat$A_vec, n_boot)
  band <- build_simultaneous_band(fate_cau$fate_est, boot_crv, ALPHA)

  cov_simul <- as.integer(all(true_curve >= band$lower & true_curve <= band$upper))
  cov_point <- mean(true_curve >= band$lower & true_curve <= band$upper)
  post_diff <- fate_cau$fate_est[idx_post] - true_curve[idx_post]
  nrmse <- sqrt(mean(post_diff^2)) / true_range
  bias <- mean(post_diff)
  width <- mean(band$upper - band$lower, na.rm = TRUE)

  cfate_hat <- cfate_cumulative(fate_cau$fate_est, time_grid, t_treat)
  cfate_boot <- t(
    apply(
      boot_crv, 1,
      function(r) cfate_cumulative(r, time_grid, t_treat)
    )
  )
  cfate_band <- build_simultaneous_band(cfate_hat, cfate_boot, ALPHA)

  true_cfate <- cfate_cumulative(true_curve, time_grid, t_treat)
  cfate_range <- max(abs(true_cfate[idx_post]))
  if (cfate_range == 0) cfate_range <- 1
  cfate_diff <- cfate_hat[idx_post] - true_cfate[idx_post]
  cfate_nrmse <- sqrt(mean(cfate_diff^2)) / cfate_range
  cfate_cov_simul <- as.integer(
    all(true_cfate >= cfate_band$lower & true_cfate <= cfate_band$upper)
  )
  cfate_cov_point <- mean(
    true_cfate >= cfate_band$lower & true_cfate <= cfate_band$upper
  )

  metrics_dt <- data.table(
    SimID = sim_id, Scenario = scen$name, ErrorDist = error_dist,
    Irregularity = irr_level, N = N, T_pts = T_pts,
    Cov_Simul = cov_simul * 100, Cov_Point = cov_point * 100,
    NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4),
    CFATE_NRMSE = round(cfate_nrmse, 6),
    CFATE_Cov_Simul = cfate_cov_simul * 100,
    CFATE_Cov_Point = round(cfate_cov_point * 100, 2)
  )

  curves_dt <- data.table(
    SimID = sim_id, Scenario = scen$name, ErrorDist = error_dist,
    Irregularity = irr_level, N = N, T_pts = T_pts,
    t_val = time_grid,
    TrueCurve = true_curve, Est = fate_cau$fate_est,
    Lower = band$lower, Upper = band$upper,
    TrueCFATE = true_cfate, CFATE_Est = cfate_hat,
    CFATE_Lower = cfate_band$lower, CFATE_Upper = cfate_band$upper
  )

  rm(dat, y_cau, boot_crv, fate_cau)
  gc(verbose = FALSE)
  list(metrics = metrics_dt, curves = curves_dt)
}

# =============================================================================
# AGREGAÇÃO (dois passes — OOM-safe)
# =============================================================================

aggregate_simulation_3 <- function(n_exemplars = 10) {
  cat("\n\n========== AGREGANDO SIM 3 ==========\n")
  config_dirs <- list.dirs(
    "output/sim3_checkpoints",
    recursive = FALSE, full.names = TRUE
  )
  if (length(config_dirs) == 0) {
    cat("ERRO: Nenhum checkpoint encontrado.\n")
    return(invisible(NULL))
  }
  out_dir <- "output/sim3_final"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Passe 1: métricas
  metrics_list <- vector("list", length(config_dirs))
  for (i in seq_along(config_dirs)) {
    bf <- list.files(config_dirs[i], pattern = "\\.rds$", full.names = TRUE)
    metrics_list[[i]] <- rbindlist(lapply(bf, function(f) readRDS(f)$metrics))
  }
  all_metrics <- rbindlist(metrics_list)

  summary_dt <- all_metrics[, .(
    N_sims = .N,
    Cov_Simul = round(mean(Cov_Simul), 2),
    Cov_Point = round(mean(Cov_Point), 2),
    NRMSE = round(mean(NRMSE), 4),
    NRMSE_SD = round(sd(NRMSE), 4),
    Bias = round(mean(Bias), 4),
    Width = round(mean(Width), 4),
    CFATE_NRMSE = round(mean(CFATE_NRMSE), 4),
    CFATE_Cov_Simul = round(mean(CFATE_Cov_Simul), 2),
    CFATE_Cov_Point = round(mean(CFATE_Cov_Point), 2)
  ), by = .(ErrorDist, Scenario, Irregularity, N, T_pts)]

  fwrite(all_metrics, "output/sim3_final/sim3_raw.csv")
  fwrite(summary_dt, "output/sim3_final/sim3_summary.csv")
  rm(all_metrics, metrics_list)
  gc(verbose = FALSE)

  cat("Sumário:\n")
  print(summary_dt)

  # Passe 2: curvas, config por config
  csumm_list <- vector("list", length(config_dirs))
  exemplar_list <- vector("list", length(config_dirs))

  for (i in seq_along(config_dirs)) {
    bf <- list.files(config_dirs[i], pattern = "\\.rds$", full.names = TRUE)
    curves <- rbindlist(lapply(bf, function(f) readRDS(f)$curves))

    csumm_list[[i]] <- curves[, .(
      RMSE_t = sqrt(mean((Est - TrueCurve)^2)),
      Bias_t = mean(Est - TrueCurve),
      Cov_t = mean(Lower <= TrueCurve & TrueCurve <= Upper) * 100,
      Width_t = mean(Upper - Lower, na.rm = TRUE),
      CFATE_RMSE_t = sqrt(mean((CFATE_Est - TrueCFATE)^2)),
      CFATE_Bias_t = mean(CFATE_Est - TrueCFATE),
      CFATE_Cov_t = mean(
        CFATE_Lower <= TrueCFATE &
          TrueCFATE <= CFATE_Upper
      ) * 100,
      CFATE_Width_t = mean(CFATE_Upper - CFATE_Lower, na.rm = TRUE)
    ), by = .(ErrorDist, Scenario, Irregularity, N, T_pts, t_val)]

    eids <- sort(
      sample(
        unique(curves$SimID),
        min(n_exemplars, length(unique(curves$SimID)))
      )
    )
    exemplar_list[[i]] <- curves[SimID %in% eids]
    rm(curves)
    gc(verbose = FALSE)
  }

  fwrite(rbindlist(csumm_list), "output/sim3_final/sim3_curve_summary.csv")
  fwrite(rbindlist(exemplar_list), "output/sim3_final/sim3_exemplars.csv")

  cat(sprintf("\nSalvo em: output/sim3_final/\n"))
  invisible(summary_dt)
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_3 <- function() {
  for (error_dist in GRID_DIST) {
    for (irr_level in IRR_LEVELS) {
      for (T_idx in seq_along(GRID_T)) {
        T_pts <- GRID_T[T_idx]
        W <- GRID_W[T_idx]
        time_grid <- GRIDS[[paste0(irr_level, "_T", T_pts)]]

        cat(
          sprintf(
            "\nCholesky MGP: dist=%s, irr=%s, T=%d...\n",
            error_dist, irr_level, T_pts
          )
        )
        L_chol <- make_mgp_chol(time_grid, sigma = 0.5, rho = 0.8)

        for (N in GRID_N) {
          n_boot <- N
          for (scen in SCENARIOS) {
            config_label <- sprintf(
              "%s_%s_N%d_T%d_%s", error_dist, irr_level, N, T_pts, scen$name
            )
            ckpt_dir <- file.path("output", "sim3_checkpoints", config_label)
            dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

            n_batches <- ceiling(N_SIM / BATCH_SIZE)
            n_cores_this <- safe_cores(N, T_pts)
            cat(
              sprintf(
                "\n>>> %s (Batches=%d, Cores=%d) <<<\n",
                config_label, n_batches, n_cores_this
              )
            )

            for (b in seq_len(n_batches)) {
              batch_file <- file.path(ckpt_dir, sprintf("batch_%04d.rds", b))
              if (file.exists(batch_file)) {
                cat(sprintf("  Batch %d/%d existe. Pulando.\n", b, n_batches))
                next
              }

              start_idx <- (b - 1) * BATCH_SIZE + 1
              end_idx <- min(b * BATCH_SIZE, N_SIM)
              batch_ids <- start_idx:end_idx
              cat(
                sprintf(
                  "  Batch %d/%d (sims %d-%d)... ", b, n_batches,
                  start_idx, end_idx
                )
              )

              batch_results <- mclapply(batch_ids, function(sid) {
                tryCatch(
                  run_single_sim_s3(
                    sid, scen, N, T_pts, time_grid, L_chol, W,
                    error_dist, irr_level, n_boot
                  ),
                  error = function(e) {
                    message(sprintf("Sim %d erro: %s", sid, e$message))
                    NULL
                  }
                )
              }, mc.cores = n_cores_this)

              batch_results <- Filter(Negate(is.null), batch_results)
              batch_dt <- rbindlist(lapply(batch_results, `[[`, "metrics"))
              batch_curves <- rbindlist(lapply(batch_results, `[[`, "curves"))

              if (nrow(batch_dt) > 0) {
                saveRDS(
                  list(metrics = batch_dt, curves = batch_curves),
                  batch_file
                )
                cat(sprintf("OK (%d sims)\n", nrow(batch_dt)))
              } else {
                cat("VAZIO\n")
              }
              rm(batch_results, batch_dt, batch_curves)
              gc(verbose = FALSE)
            }
          }
        }
      }
    }
  }

  aggregate_simulation_3()
}

# =============================================================================
# RUN
# =============================================================================

results <- run_simulation_3()
