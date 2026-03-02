# =============================================================================
# scripts/sim_2.R
# Simulação 2: Estimador Online (Streaming via new_wavelet_stream)
# Grid: N ∈ {200, 500, 2000}, T ∈ {128, 256, 512}, 4 cenários, 1000 sims
# =============================================================================

library(data.table)
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
source("R/utils_dgp_functions.R")

# =============================================================================
# CONFIGURAÇÃO
# =============================================================================

QUICK_TEST = FALSE  # TRUE = 10 sims para depuração

N_SIM = ifelse(QUICK_TEST, 10, 1000)
BATCH_SIZE = ifelse(QUICK_TEST, 5, 25)

t_start = 20.0
t_end = 63.0
t_treat = 20.0

# Grid
GRID_N = c(200, 500, 2000)
GRID_T = c(128, 256, 512)
GRID_W = c(64,  64,  128)

# Cenários (mesmas ondaletas da sim 1)
SCENARIOS = list(
  list(name = "Sigmoid",   func = effect_sigmoid,   mag = 2, family = "haar"),
  list(name = "Bumps",     func = effect_bumps,     mag = 3, family = "db2"),
  list(name = "Doppler",   func = effect_doppler,   mag = 5, family = "db2"),
  list(name = "Heavisine", func = effect_heavisine, mag = 4, family = "haar")
)

# =============================================================================
# WORKER: Uma simulação streaming com medição de latência ponto a ponto
# =============================================================================

run_single_sim_online = function(sim_id, scen, N, T_pts, time_grid, L_chol, W) {
  # 1. Gera dataset
  dat = generate_ecker_dataset(scen, N, T_pts, time_grid, L_chol,
                                 t_start = t_start, t_end = t_end, t_treat = t_treat)
  true_curve = dat$true_curve
  true_range = diff(range(true_curve))
  if (true_range == 0) true_range = 1

  proj = precompute_projection(dat$X_reg, dat$A_vec)

  # 2. Suavização STREAMING via new_wavelet_stream (gravando tempo)
  scheme_obj = lifting_scheme(scen$family)
  max_levels = floor(log2(W))

  Y_online = matrix(0, N, T_pts)

  # Vamos gravar o tempo médio que levou para processar todos os indivíduos em cada instante de tempo
  # latencia_t[k] = (tempo total gasto pelos N processadores no ponto k) / N
  # Retorna ms (milissegundos) ou us (microssegundos) por ponto
  latencia_t_s = numeric(T_pts)

  # Instancia N processadores independendes
  processors = lapply(1:N, function(i) {
    new_wavelet_stream(scheme_obj, window_size = W, levels = max_levels,
                       method = "semisoft", alpha = 0.3)
  })

  # Processamento estritamente causal ao longo do tempo (simulando a chegada conjunta de dados)
  for (k in seq_len(T_pts)) {
    # Em R, Sys.time tem resolução baixa. Como calcularemos para N indivíduos no loop,
    # a resolução dilui. Usamos microbenchmark.
    t0 = Sys.time()
    for (i in seq_len(N)) {
      # Alimenta 1 novo ponto no processador do indivíduo i
      Y_online[i, k] = processors[[i]](dat$Y_wide[i, k])
    }
    t1 = Sys.time()

    # Tempo em microssegundos por indivíduo no instante k
    latencia_t_s[k] = as.numeric(difftime(t1, t0, units = "secs")) / N
  }

  # 3. Estimação FATE OLS
  fate_online = estimate_fate_fast(Y_online, proj)
  nrmse_online = sqrt(mean((fate_online$fate_est - true_curve)^2)) / true_range

  # 4. Saída
  metrics_dt = data.table(
    SimID = sim_id, Scenario = scen$name, N = N, T_pts = T_pts,
    NRMSE_online = round(nrmse_online, 6),
    Latency_mean_us = mean(latencia_t_s) * 1e6
  )

  curves_dt = data.table(
    SimID = sim_id, Scenario = scen$name, N = N, T_pts = T_pts, t_val = time_grid,
    TrueCurve = true_curve, Est = fate_online$fate_est,
    Latency_us = latencia_t_s * 1e6  # Vetor T de latências (microssegundos)
  )

  # Limpeza
  rm(dat, Y_online, processors); gc(verbose = FALSE)

  list(metrics = metrics_dt, curves = curves_dt)
}

# =============================================================================
# EXECUTOR PRINCIPAL
# =============================================================================

run_simulation_2 = function() {

  for (T_idx in seq_along(GRID_T)) {
    T_pts = GRID_T[T_idx]
    W = GRID_W[T_idx]
    time_grid = seq(t_start, t_end, length.out = T_pts)

    cat(sprintf("Pré-computando Cholesky MGP para T=%d...\n", T_pts))
    L_chol = make_mgp_chol(time_grid, sigma = 0.5, rho = 0.8)

    for (N in GRID_N) {
      for (scen in SCENARIOS) {
        config_label = sprintf("N%d_T%d_%s", N, T_pts, scen$name)
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

          cat(sprintf("  Batch %d/%d (sims %d-%d)... ", b, n_batches, start_idx, end_idx))

          batch_results = lapply(batch_ids, function(sid) {
            tryCatch(
              run_single_sim_online(sid, scen, N, T_pts, time_grid, L_chol, W),
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

  # --- Agregação ---
  cat("\n\n========== AGREGANDO RESULTADOS SIM 2 ==========\n")
  all_files = list.files("output/sim2_checkpoints", pattern = "\\.rds$",
                           recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) { cat("ERRO: Sem checkpoints!\n"); return(NULL) }

  all_data = lapply(all_files, readRDS)
  all_results = rbindlist(lapply(all_data, function(x) x$metrics))
  all_curves = rbindlist(lapply(all_data, function(x) x$curves))

  summary_dt = all_results[, .(
    N_sims = .N,
    NRMSE_online = round(mean(NRMSE_online,  na.rm = TRUE), 4),
    Latency_mean_us = round(mean(Latency_mean_us, na.rm = TRUE), 2)
  ), by = .(Scenario, N, T_pts)]

  out_dir = "output/sim2_final"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(all_results, file.path(out_dir, "sim2_raw.csv"))
  fwrite(summary_dt,  file.path(out_dir, "sim2_summary.csv"))
  fwrite(all_curves,  file.path(out_dir, "sim2_curves.csv"))

  cat("\nSumário:\n")
  print(summary_dt)
  cat(sprintf("\nSalvo em: %s/\n", out_dir))

  return(summary_dt)
}

# =============================================================================
# RUN
# =============================================================================
results = run_simulation_2()
