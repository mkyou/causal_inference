# =============================================================================
# scripts/sim_1.R
# Simulação 1: Comparação Offline (Wavelet Offline vs Wavelet Causal vs Spline)
# Grid: N ∈ {200, 500, 2000}, T ∈ {128, 256, 512}, 4 cenários, 1000 sims
# =============================================================================

library(data.table)
library(parallel)
library(rLifting)
library(fda)
library(MASS)

# --- Setup portável ---
script_dir = tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) {
    args = commandArgs(trailingOnly = FALSE)
    file_arg = grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0)
      dirname(normalizePath(sub("^--file=", "", file_arg)))
    else getwd()
  }
)
if (basename(script_dir) == "scripts")
  setwd(dirname(script_dir)) else setwd(script_dir)

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

QUICK_TEST = FALSE

N_SIM = ifelse(QUICK_TEST, 1, 1000)
BATCH_SIZE = ifelse(QUICK_TEST, 1, 25)
N_CORES = max(1, detectCores() - 2)

safe_cores = function(N, T_pts, mem_safety = 0.75) {
  data_bytes = 12 * N * T_pts * 8 # Y_wide + suavizados + bootstraps + cfate
  fork_base  = 400e6 # overhead do processo pai ao fazer fork
  bytes_per_worker = data_bytes + fork_base
  avail_bytes = tryCatch(
    as.numeric(
      system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern = TRUE)
    ) * 1024,
    error = function(e) 8e9
  )
  max(1L, min(N_CORES, floor(avail_bytes * mem_safety / bytes_per_worker)))
}
ALPHA = 0.05

t_start = 20.0
t_end = 80.0
t_treat = 50.0

# Grid de cenários
GRID_N = if (QUICK_TEST) c(500) else c(200, 500, 2000)
GRID_T = if (QUICK_TEST) c(256) else c(128, 256, 512)
GRID_W = if (QUICK_TEST) c(64) else c(64,  64,  128)

# Cenários de efeito (ondaleta padronizada por cenário)
SCENARIOS = list(
  list(name = "Sigmoid", func = effect_sigmoid, mag = 2, family = "haar"),
  list(name = "Bumps", func = effect_bumps, mag = 3, family = "haar"),
  list(name = "Doppler", func = effect_doppler, mag = 5, family = "haar"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

CONF_LEVELS = list(
  list(name = "conf_low",  cols = 1:8),
  list(name = "conf_med",  cols = 1:6),
  list(name = "conf_high", cols = 1:4)
)

NOISE_DISTS = c("gaussian", "t")
DF_T = 5

# =============================================================================
# WORKER: executa uma simulação para um cenário + N + T
# =============================================================================

run_single_sim = function(
  sim_id, scen, N, T_pts, time_grid, L_chol,
  W, n_boot, scheme, error_dist = "gaussian"
) {
  # 1. Gera dataset
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

  # 2. Projeção pré-computada
  proj = precompute_projection(dat$X_reg, dat$A_vec)

  # 3. Suavização
  # Wavelet offline
  y_off = tryCatch(
    smooth_data_matrix(
      time_grid, dat$Y_wide, "rlifting_offline",
      scheme_name = scen$family, threshold_method = "semisoft", alpha = 0.3
    ),
    error = function(e) NULL
  )

  # Wavelet causal (batch)
  max_levels = floor(log2(W))
  scheme_obj = lifting_scheme(scen$family)
  y_cau = tryCatch({
    t(
      apply(
        dat$Y_wide, 1, 
        function(sig) {
          denoise_signal_causal(
            sig, scheme_obj, levels = max_levels, window_size = W,
            method = "semisoft", alpha = 0.3
          )
        }
      )
    )
  }, error = function(e) NULL)

  # 4. Estimação pontual
  fate_off = if (!is.null(y_off)) estimate_fate_fast(y_off, proj) else NULL
  fate_cau = if (!is.null(y_cau)) estimate_fate_fast(y_cau, proj) else NULL

  # 5. Bootstrap (dados ruidosos — sem suavização)
  boot_curves = run_bootstrap_raw(dat$Y_wide, dat$X_reg, dat$A_vec, n_boot)

  # 6. Bandas bootstrap (usadas por wavelet offline E causal)
  compute_metrics = function(fate_est, method_name, boot = boot_curves) {
    if (is.null(fate_est)) return(list(metrics = NULL, curves = NULL))
    band = build_simultaneous_band(fate_est$fate_est, boot, ALPHA)
    cov_simul = as.integer(
      all(
        true_curve >= band$lower & true_curve <= band$upper
      )
    )
    cov_point = mean(true_curve >= band$lower & true_curve <= band$upper)
    post_diff = fate_est$fate_est[idx_post] - true_curve[idx_post]
    nrmse = sqrt(mean(post_diff^2)) / true_range
    bias  = mean(post_diff)
    width = mean(band$upper - band$lower, na.rm = TRUE)

    cfate_est_vec = cfate_cumulative(fate_est$fate_est, time_grid, t_treat)
    cfate_lower = numeric(length(time_grid))
    cfate_upper = numeric(length(time_grid))
    valid_boot = which(complete.cases(boot))
    if (length(valid_boot) >= 50) {
      cfate_boot_mat = t(
        apply(
          boot[valid_boot, , drop = FALSE], 1,
          function(r) cfate_cumulative(r, time_grid, t_treat)
        )
      )
      cb = build_simultaneous_band(
        cfate_est_vec[idx_post],
        cfate_boot_mat[, idx_post, drop = FALSE],
        ALPHA
      )
      cfate_lower[idx_post] = cb$lower
      cfate_upper[idx_post] = cb$upper
    }

    metrics_dt = data.table(
      SimID = sim_id, ErrorDist = error_dist,
      Scenario = scen$name, Method = method_name,
      N = N, T_pts = T_pts,
      Cov_Simul = cov_simul, Cov_Point = round(cov_point, 4),
      NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4)
    )

    curves_dt = data.table(
      SimID = sim_id, ErrorDist = error_dist,
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

  # 6b. Estimadores confundidos por nível
  res_conf_list = lapply(
    CONF_LEVELS,
    function(lev) {
      X_rc = cbind(1, dat$X_mat[, lev$cols])
      pc = precompute_projection(X_rc, dat$A_vec)
      fate_off_c = if (!is.null(y_off)) estimate_fate_fast(y_off, pc) else NULL
      fate_cau_c = if (!is.null(y_cau)) estimate_fate_fast(y_cau, pc) else NULL
      boot_c = run_bootstrap_raw(dat$Y_wide, X_rc, dat$A_vec, n_boot)
      list(
        off = compute_metrics(
          fate_off_c, paste0("wavelet_offline_", lev$name),
          boot = boot_c
        ),
        cau = compute_metrics(
          fate_cau_c, paste0("wavelet_causal_", lev$name),
          boot = boot_c
        )
      )
    }
  )

  # 7. Splines Ecker (bandas analíticas)
  ecker_res = tryCatch(
    estimate_fate_ecker(dat$Y_wide, dat$X_reg, dat$A_vec, time_grid,
                         alpha = ALPHA),
    error = function(e) NULL
  )

  res_spl = list(metrics = NULL, curves = NULL)
  if (!is.null(ecker_res)) {
    cov_simul = as.integer(
      all(
        true_curve >= ecker_res$lower &
          true_curve <= ecker_res$upper
      )
    )
    cov_point = mean(
      true_curve >= ecker_res$lower & 
        true_curve <= ecker_res$upper
    )
    post_diff = ecker_res$fate_est[idx_post] - true_curve[idx_post]
    nrmse = sqrt(mean(post_diff^2)) / true_range
    bias  = mean(post_diff)
    width = mean(ecker_res$upper - ecker_res$lower, na.rm = TRUE)

    spl_cfate_est = cfate_cumulative(ecker_res$fate_est, time_grid, t_treat)
    spl_cfate_lower = cfate_cumulative(ecker_res$lower, time_grid, t_treat)
    spl_cfate_upper = cfate_cumulative(ecker_res$upper, time_grid, t_treat)

    res_spl$metrics = data.table(
      SimID = sim_id, ErrorDist = error_dist,
      Scenario = scen$name, Method = "spline_ecker",
      N = N, T_pts = T_pts,
      Cov_Simul = cov_simul, Cov_Point = round(cov_point, 4),
      NRMSE = round(nrmse, 6), Bias = round(bias, 6), Width = round(width, 4)
    )

    res_spl$curves = data.table(
      SimID = sim_id, ErrorDist = error_dist,
      Scenario = scen$name, Method = "spline_ecker",
      N = N, T_pts = T_pts, t_val = time_grid,
      TrueCurve = true_curve, RawEst = raw_est,
      Est = ecker_res$fate_est, Lower = ecker_res$lower, 
      Upper = ecker_res$upper,
      CFATE_True = cfate_true, CFATE_Est = spl_cfate_est,
      CFATE_Lower = spl_cfate_lower, CFATE_Upper = spl_cfate_upper
    )
  }

  # Limpeza
  rm(dat, y_off, y_cau, boot_curves); gc(verbose = FALSE)

  conf_metrics = unlist(
    lapply(res_conf_list, function(x) list(x$off$metrics, x$cau$metrics)),
    recursive = FALSE
  )
  conf_curves = unlist(
    lapply(res_conf_list, function(x) list(x$off$curves, x$cau$curves)),
    recursive = FALSE
  )
  all_metrics = c(
    list(
      res_off$metrics, res_cau$metrics, res_spl$metrics
    ), 
    conf_metrics
  )
  
  all_curves = c(
    list(
      res_off$curves, res_cau$curves, res_spl$curves
    ), 
    conf_curves
  )
  list(
    metrics = rbindlist(all_metrics[!sapply(all_metrics, is.null)]),
    curves  = rbindlist(all_curves[!sapply(all_curves,  is.null)])
  )
}

# =============================================================================
# AGREGAÇÃO FINAL
# =============================================================================

aggregate_simulation_1 = function(sim_tag, ckpt_tag, out_tag, n_exemplars = 10) {
  out_dir   = file.path("output", out_tag)
  ckpt_base = file.path("output", ckpt_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  config_dirs = sort(list.dirs(ckpt_base, recursive = FALSE))
  n_configs   = length(config_dirs)
  if (n_configs == 0) { cat("ERRO: Sem configs!\n"); return(NULL) }

  set.seed(42)

  # --- Passo 1/2: métricas brutas (sem curvas) ---
  cat(sprintf("\n--- Passo 1/2: métricas brutas (%d configs) ---\n", n_configs))
  metrics_list = vector("list", n_configs)
  for (i in seq_along(config_dirs)) {
    batch_files = list.files(config_dirs[i], pattern = "\\.rds$", full.names = TRUE)
    cat(sprintf("  [%d/%d] %s\r", i, n_configs, basename(config_dirs[i])))
    metrics_list[[i]] = rbindlist(lapply(batch_files, function(f) readRDS(f)$metrics))
  }
  cat("\n")

  all_metrics = rbindlist(metrics_list)
  rm(metrics_list); gc(verbose = FALSE)

  summary_dt = all_metrics[, .(
    N_sims    = .N,
    Cov_Simul = round(mean(Cov_Simul, na.rm = TRUE) * 100, 1),
    Cov_Point = round(mean(Cov_Point, na.rm = TRUE) * 100, 1),
    NRMSE     = round(mean(NRMSE, na.rm = TRUE), 4),
    NRMSE_SD  = round(sd(NRMSE, na.rm = TRUE), 4),
    Bias      = round(mean(Bias, na.rm = TRUE), 4),
    Width     = round(mean(Width, na.rm = TRUE), 3)
  ), by = .(ErrorDist, Scenario, Method, N, T_pts)]

  fwrite(all_metrics, file.path(out_dir, paste0(sim_tag, "_raw.csv")))
  fwrite(summary_dt,  file.path(out_dir, paste0(sim_tag, "_summary.csv")))
  cat(sprintf(
    "  raw: %d linhas | summary: %d linhas\n",
    nrow(all_metrics), nrow(summary_dt)
  ))
  rm(all_metrics); gc(verbose = FALSE)

  # --- Passo 2/2: curvas-resumo + exemplares (config por config) ---
  cat(sprintf("\n--- Passo 2/2: curvas-resumo + exemplares ---\n"))
  curve_summ_list = vector("list", n_configs)
  exemplar_list   = vector("list", n_configs)

  for (i in seq_along(config_dirs)) {
    batch_files = list.files(config_dirs[i], pattern = "\\.rds$", full.names = TRUE)
    cat(sprintf("  [%d/%d] %s...\n", i, n_configs, basename(config_dirs[i])))

    curves = rbindlist(lapply(batch_files, function(f) readRDS(f)$curves))

    curve_summ_list[[i]] = curves[, .(
      RMSE_t  = sqrt(mean((Est - TrueCurve)^2, na.rm = TRUE)),
      Bias_t  = mean(Est - TrueCurve, na.rm = TRUE),
      Cov_t   = mean(TrueCurve >= Lower & TrueCurve <= Upper, na.rm = TRUE),
      Width_t = mean(Upper - Lower, na.rm = TRUE)
    ), by = .(ErrorDist, Scenario, Method, N, T_pts, t_val)]

    eids = sort(sample(unique(curves$SimID), min(n_exemplars, length(unique(curves$SimID)))))
    exemplar_list[[i]] = curves[SimID %in% eids]

    rm(curves); gc(verbose = FALSE)
  }

  all_curve_summ = rbindlist(curve_summ_list)
  all_exemplars  = rbindlist(exemplar_list)

  fwrite(all_curve_summ, file.path(out_dir, paste0(sim_tag, "_curve_summary.csv")))
  fwrite(all_exemplars,  file.path(out_dir, paste0(sim_tag, "_exemplars.csv")))
  cat(sprintf(
    "  curve_summary: %d linhas | exemplars: %d linhas\n",
    nrow(all_curve_summ), nrow(all_exemplars)
  ))

  cat("\nSumário:\n")
  print(summary_dt)
  cat(sprintf("\nSalvo em: %s/\n", out_dir))
  return(summary_dt)
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_1 = function() {
  sim_tag  = if (QUICK_TEST) "sim1_test" else "sim1"
  ckpt_tag = paste0(sim_tag, "_checkpoints")
  out_tag  = paste0(sim_tag, "_final")

  for (error_dist in NOISE_DISTS) {
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
          config_label = sprintf(
            "%s_N%d_T%d_%s", error_dist,
            N, T_pts, scen$name
          )
          ckpt_dir = file.path("output", ckpt_tag, config_label)
          dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

          # Tarefas
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

            cat(
              sprintf(
                "  Batch %d/%d (sims %d-%d)... ", b, n_batches, 
                start_idx, end_idx
              )
            )

            batch_results = mclapply(batch_ids, function(sid) {
              tryCatch(
                run_single_sim(
                  sid, scen, N, T_pts, time_grid, 
                  L_chol, W, n_boot, NULL, error_dist
                ),
                error = function(e) { 
                  message(sprintf("Sim %d erro: %s", sid, e$message)); 
                  NULL
                }
              )
            }, mc.cores = n_cores_this)

            batch_dt = rbindlist(lapply(batch_results, function(x) x$metrics))
            batch_curves = rbindlist(
              lapply(batch_results, function(x) x$curves)
            )

            if (nrow(batch_dt) > 0) {
              saveRDS(
                list(metrics = batch_dt, curves = batch_curves), 
                batch_file
              )
              cat(
                sprintf(
                  "OK (%d rows metrics, %d rows curves)\n", nrow(batch_dt), 
                  nrow(batch_curves)
                )
              )
            } else {
              cat("VAZIO\n")
            }

            rm(batch_results, batch_dt, batch_curves); gc(verbose = FALSE)
          }
        }
      }
    }
  } # error_dist

  # --- Agregação final ---
  cat("\n\n========== AGREGANDO RESULTADOS ==========\n")
  return(aggregate_simulation_1(sim_tag, ckpt_tag, out_tag, n_exemplars = 10))
}

# =============================================================================
# RUN
# =============================================================================
results = run_simulation_1()

if (QUICK_TEST) {
  library(ggplot2)
  dt_c = data.table::fread("output/sim1_test_final/sim1_test_curves.csv")

  p = ggplot(dt_c, aes(x = t_val)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.2) +
    geom_line(aes(y = Est, color = Method), linewidth = 0.8) +
    geom_line(aes(y = TrueCurve),
              color = "black", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = t_treat, linetype = "dotted", color = "gray40") +
    facet_wrap(~ Method, ncol = 1) +
    labs(x = "t", y = "FATE(t)", title = "Validação: janela pré-tratamento") +
    theme_bw() + theme(legend.position = "none")

  out_path = "output/sim1_test_final/test_curves.png"
  ggsave(out_path, p, width = 8, height = 9, dpi = 150)
  cat(sprintf("\nGráfico salvo em: %s\n", out_path))
}
