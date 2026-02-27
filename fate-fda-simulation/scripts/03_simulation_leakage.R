library(data.table)
library(ggplot2)
library(patchwork)
library(dplyr)
library(kableExtra)
library(rLifting)

source("R/utils_dgp_functions.R")
source("R/generate_sequential_dgp.R")
source("R/generate_functional_outcomes.R")
source("R/process_functional_data.R")
source("R/config_ecker.R")

# Para não ter que editar a original; copiei e coloquei o drop = FALSE
smooth_spline_robust = function(time_grid, wide_matrix, nbasis = 50) {
  # Cria base B-Spline
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), 
    nbasis = nbasis
    )
  
  log_lambdas = seq(-10, 1, length.out = 20)
  gcvs = numeric(length(log_lambdas))
  
  sample_indices = if(nrow(wide_matrix) > 50) 
    sample(nrow(wide_matrix), 50) else 1:nrow(wide_matrix)
  test_data = t(wide_matrix[sample_indices, , drop = FALSE]) 
  
  for(i in seq_along(log_lambdas)) {
    fdPar_obj = fdPar(bspline_basis, Lfdobj = 2, lambda = 10^log_lambdas[i])
    fd_try = smooth.basis(time_grid, test_data, fdPar_obj)
    gcvs[i] = sum(fd_try$gcv)
  }
  
  best_lambda = 10^log_lambdas[which.min(gcvs)]
  
  final_fdPar = fdPar(bspline_basis, Lfdobj = 2, lambda = best_lambda)
  final_fd = smooth.basis(time_grid, t(wide_matrix), final_fdPar)$fd
  
  return(t(eval.fd(time_grid, final_fd)))
}

LEAK_CONFIG = list(
  N_SIM = 1000,
  BATCH_SIZE = 50,
  N_SUBJ = 500,          
  T_PTS = 64,          
  T_START = 20,
  T_END = 80,
  TIME_T1 = 35.0,        
  TIME_T2 = 60.0,        
  LEAK_WINDOW = 4.0      
)

# FUNÇÕES AUXILIARES
smooth_causal_wrapper = function(y_vec, scheme="db2", window=64) {
  scheme_obj = rLifting::lifting_scheme(scheme)
  lvl = floor(log2(window))
  y_smooth = rLifting::denoise_signal_causal(
    y_vec,
    scheme = scheme_obj,
    levels = lvl,
    window_size = window,
    method = "semisoft",
    alpha = 0.3, 
    beta = 1.2
  )
  
  return(y_smooth)
}

theme_paper = theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 12),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# Simulação
run_leakage_experiment = function() {
  
  # Cria diretório de checkpoints
  dir_checkpoints = file.path(
    "output/checkpoints_sim3", sprintf("N%d", LEAK_CONFIG$N_SUBJ)
    )
  dir.create(dir_checkpoints, recursive = TRUE, showWarnings = FALSE)
  
  # Preparação
  time_grid = seq(
    LEAK_CONFIG$T_START, LEAK_CONFIG$T_END, 
    length.out = LEAK_CONFIG$T_PTS
    )
  
  # Funções de efeito (definição global)
  func_sigmoid = function(t, s) 
    effect_sigmoid(t, s, LEAK_CONFIG$T_START, LEAK_CONFIG$T_END, magnitude = 2)
  func_step = function(t, s) 
    effect_heavisine(
      t, s, LEAK_CONFIG$T_START, 
      LEAK_CONFIG$T_END, magnitude = -3
      )
  
  true_curve_a = func_sigmoid(time_grid, LEAK_CONFIG$TIME_T1)
  true_curve_b = func_step(time_grid, LEAK_CONFIG$TIME_T2)
  true_effect = true_curve_a - true_curve_b
  
  # Índices de leakage
  get_leak_indices = function(t_event) {
    which(
      time_grid >= (t_event - LEAK_CONFIG$LEAK_WINDOW) & 
        time_grid < t_event
      )
  }
  idx_leak_t1 = get_leak_indices(LEAK_CONFIG$TIME_T1)
  idx_leak_t2 = get_leak_indices(LEAK_CONFIG$TIME_T2)
  
  # Loop do batch
  n_batches = ceiling(LEAK_CONFIG$N_SIM / LEAK_CONFIG$BATCH_SIZE)
  
  cat(
    sprintf(
      ">>> Iniciando simulação de Leakage (Total: %d, Lotes: %d) <<<\n",
      LEAK_CONFIG$N_SIM, n_batches
      )
    )
  
  # Barra de progresso global
  pb = txtProgressBar(min = 0, max = LEAK_CONFIG$N_SIM, style = 3)
  
  for(b in 1:n_batches) {
    # Arquivo do lote atual
    batch_file = file.path(dir_checkpoints, sprintf("batch_%03d.rds", b))
    
    # Define intervalo de simulações deste lote
    idx_start = (b - 1) * LEAK_CONFIG$BATCH_SIZE + 1
    idx_end = min(b * LEAK_CONFIG$BATCH_SIZE, LEAK_CONFIG$N_SIM)
    
    # Verifica se já existe
    if(file.exists(batch_file)) {
      setTxtProgressBar(pb, idx_end)
      next
    }
    
    # Lista temporária para armazenar resultados deste lote
    batch_results_list = vector("list", length = idx_end - idx_start + 1)
    counter = 1
    
    # Loop interno (simulações do lote)
    for(sim in idx_start:idx_end) {
      cfg = get_ecker_configuration(LEAK_CONFIG$N_SUBJ, LEAK_CONFIG$T_PTS)
      cfg$theta_matrix = rbind(cfg$theta_matrix, cfg$theta_matrix)
      cfg$gamma_matrix = cbind(cfg$gamma_matrix, cfg$gamma_matrix)
      
      dgp_covs = generate_sequential_dgp(
        LEAK_CONFIG$N_SUBJ, time_grid, 
        c(LEAK_CONFIG$TIME_T1, LEAK_CONFIG$TIME_T2),
        cfg$baseline_func_templates, cfg$baseline_params_list,
        cfg$beta_funcs_matrix, cfg$gamma_matrix, cfg$theta_matrix,
        rep(0, 9), cfg$static_indices,
        regime_a = c(1, 0), regime_b = c(0, 1)
      )
      
      dgp_out = generate_functional_outcome(
        dgp_covs, cfg$baseline_y_func, 
        list(func_sigmoid, func_step),
        cfg$concurrent_effect_funcs, 
        c(LEAK_CONFIG$TIME_T1, LEAK_CONFIG$TIME_T2), 
        noise_sd = 0.8
      )
      
      mean_curves = dgp_out[
        , .(Mean_Ya = mean(a_y), Mean_Yb = mean(b_y)), by = t
        ]
      y_a_vec = mean_curves$Mean_Ya
      y_b_vec = mean_curves$Mean_Yb
      
      # Suavização
      # 1. Spline
      spline_a = smooth_spline_robust(
        time_grid, matrix(y_a_vec, nrow = 1), nbasis = 60
        )[1,]
      spline_b = smooth_spline_robust(
        time_grid, matrix(y_b_vec, nrow = 1), nbasis = 60
        )[1,]
      eff_spline = spline_a - spline_b
      
      # 2. rLifting Offline
      rlift_off_a = smooth_data_matrix(
        time_grid, matrix(y_a_vec, nrow = 1), 
        "rlifting_offline", scheme_name = "db2"
        )[1,]
      rlift_off_b = smooth_data_matrix(
        time_grid, matrix(y_b_vec, nrow = 1), 
        "rlifting_offline", scheme_name = "db2"
        )[1,]
      eff_offline = rlift_off_a - rlift_off_b
      
      # 3. rLifting Causal
      rlift_causal_a = smooth_causal_wrapper(y_a_vec, scheme = "db2")
      rlift_causal_b = smooth_causal_wrapper(y_b_vec, scheme = "db2")
      eff_causal = rlift_causal_a - rlift_causal_b
      
      # Métricas
      calc_leak = function(est, truth, idx) mean(est[idx] - truth[idx])
      
      res_dt = data.table(
        SimID = sim,
        Method = c("Spline", "Offline Wavelet", "Causal Wavelet"),
        NRMSE = c(
          sqrt(mean((eff_spline - true_effect)^2)),
          sqrt(mean((eff_offline - true_effect)^2)),
          sqrt(mean((eff_causal - true_effect)^2))
        ) / diff(range(true_effect)),
        MaxError = c(
          max(abs(eff_spline - true_effect)),
          max(abs(eff_offline - true_effect)),
          max(abs(eff_causal - true_effect))
        ),
        Leakage_T1 = c(
          calc_leak(eff_spline, true_effect, idx_leak_t1),
          calc_leak(eff_offline, true_effect, idx_leak_t1),
          calc_leak(eff_causal, true_effect, idx_leak_t1)
        ),
        Leakage_T2 = c(
          calc_leak(eff_spline, true_effect, idx_leak_t2),
          calc_leak(eff_offline, true_effect, idx_leak_t2),
          calc_leak(eff_causal, true_effect, idx_leak_t2)
        )
      )
      
      curve_dt = data.table(
        SimID = sim, t = time_grid, True_Effect = true_effect,
        Spline = eff_spline, Offline = eff_offline, Causal = eff_causal
      )
      
      batch_results_list[[counter]] = list(
        metrics = res_dt, curves = curve_dt
        )
      counter = counter + 1
      
      setTxtProgressBar(pb, sim)
    } 
    
    # Salva batch e limpa memória
    saveRDS(batch_results_list, batch_file)
    rm(batch_results_list)
    gc() # Garbage Collection
    
  } # Fim Loop Batch
  
  close(pb)
  
  # Agregação
  cat("\n>>> Agregando resultados dos batches... <<<\n")
  files = list.files(dir_checkpoints, full.names = TRUE, pattern = "\\.rds$")
  
  # Lê todos os batches e 'desenrola' a lista (unlist recursive=FALSE)
  # Isso cria uma lista plana onde cada elemento é uma simulação
  all_sims_flat = unlist(lapply(files, readRDS), recursive = FALSE)
  
  # Limpa checkpoints temporários
  unlink(dir_checkpoints, recursive = TRUE)
  
  return(all_sims_flat)
}

# Artefatos
generate_leakage_artifacts = function(res_list) {
  artifacts_dir = sprintf(
    "output/sim3_N%d_T%d/artifacts", 
    LEAK_CONFIG$N_SUBJ, 
    LEAK_CONFIG$T_PTS
    )
  dir.create(artifacts_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_metrics = rbindlist(lapply(res_list, `[[`, "metrics"))
  all_curves = rbindlist(lapply(res_list, `[[`, "curves"))
  
  # 1. Tabela comparativa completa
  tbl_leakage = all_metrics[, .(
    NRMSE_M = mean(NRMSE), 
    NRMSE_sd = sd(NRMSE),
    MaxErr_M = mean(MaxError),
    MaxErr_sd = sd(MaxError),
    Leak_T1_M = mean(Leakage_T1), Leak_T1_SD = sd(Leakage_T1),
    Leak_T2_M = mean(Leakage_T2), Leak_T2_SD = sd(Leakage_T2)
  ), by = Method]
  
  tbl_pretty = tbl_leakage |>
    dplyr::mutate(
      `NRMSE (SD)` = sprintf("%.4f (%.4f)", NRMSE_M, NRMSE_sd),
      `Max Error (SD)` = sprintf("%.4f (%.4f)", MaxErr_M, MaxErr_sd),
      `Leakage T1 (Suave)` = sprintf("%.4f (%.4f)", Leak_T1_M, Leak_T1_SD),
      `Leakage T2 (Abrupto)` = sprintf("%.4f (%.4f)", Leak_T2_M, Leak_T2_SD)
    ) |>
    dplyr::select(
      `Metodo` = Method, `NRMSE (SD)`, 
      `Max Error (SD)`, `Leakage T1 (Suave)`, 
      `Leakage T2 (Abrupto)`
      )
  
  kbl = kable(
    tbl_pretty, format = "html", booktabs = TRUE,
    caption = "Viés de antecipação (leakage) em eventos suaves vs. abruptos"
    ) |>
    kable_styling(bootstrap_options = c("striped", "hover")) |>
    row_spec(which(tbl_pretty$Method == "Spline")) |>
    row_spec(which(tbl_pretty$Method == "Causal wavelet"))
  
  save_kable(
    kbl, file.path(artifacts_dir, sprintf(
      "table_leakage_comparison_N%d_T%d.html", 
      LEAK_CONFIG$N_SUBJ, LEAK_CONFIG$T_PTS
      ))
    )
  
  # 2. Gráfico multi-painel (global + zoom T1 + zoom T2)
  avg_curves = all_curves[, .(
    Spline = mean(Spline), 
    Offline = mean(Offline), 
    Causal = mean(Causal), 
    Truth = mean(True_Effect)
  ), by = t]
  
  long_curves = melt(
    avg_curves, id.vars = c("t", "Truth"),
    variable.name = "Method", value.name = "Estimate"
    )
  
  cols = c(
    "Spline" = "firebrick", 
    "Offline" = "forestgreen", 
    "Causal" = "royalblue"
    )
  
  # A. Global
  p_main = ggplot(long_curves, aes(x = t)) +
    geom_line(
      aes(y = Truth), 
      color = "black", 
      linetype = "dashed", 
      alpha = 0.5
      ) +
    geom_line(aes(y = Estimate, color = Method), linewidth = 0.7) +
    scale_color_manual(values = cols) +
    labs(x = NULL, y = "Efeito") +
    theme_paper + theme(legend.position = "top")
  
  # Função para gerar Zooms
  plot_zoom = function(center_time, title, y_lims) {
    ggplot(long_curves, aes(x = t)) +
      # Zona de leakage
      geom_rect(
        aes(
          xmin = center_time - LEAK_CONFIG$LEAK_WINDOW, 
          xmax = center_time, ymin = -Inf, 
          ymax = Inf
          ),
        fill = "grey90", alpha = 0.01, inherit.aes = FALSE
      ) +
      geom_vline(xintercept = center_time, linetype = "solid") +
      geom_line(aes(y = Truth), color = "black", linetype = "dashed") +
      geom_line(aes(y = Estimate, color = Method), linewidth = 1.0) +
      scale_color_manual(values = cols) +
      coord_cartesian(
        xlim = c(center_time - 5, center_time + 3), 
        ylim = y_lims
        ) +
      labs(title = NULL, x = "Tempo", y = NULL) +
      theme_paper + theme(legend.position = "none")
  }
  
  # B. Zoom T1 (Sigmoid) - Espera-se pouco erro
  p_zoom_t1 = plot_zoom(
    LEAK_CONFIG$TIME_T1, "Zoom T1: sigmoide (suave)", c(-0.2, 0.5))
  
  # C. Zoom T2 (Heavisine) - Espera-se mais erro no Spline
  p_zoom_t2 = plot_zoom(
    LEAK_CONFIG$TIME_T2, "Zoom T2: degrau (abrupto)", c(-0.5, 2.2)
    )
  
  # Montagem: main em cima, zooms em baixo
  layout = "
  AAAA
  BCCD
  "
  # Ajuste simples: main em cima, dois zooms lado a lado em baixo
  p_final = p_main / (p_zoom_t1 | p_zoom_t2) + 
    plot_layout(heights = c(1, 1))
  
  ggsave(
    file.path(artifacts_dir, sprintf(
      "plot_leakage_panels_N%d_T%d.png", 
      LEAK_CONFIG$N_SUBJ, LEAK_CONFIG$T_PTS
      )), 
    p_final, width = 12, height = 8)
  
  cat("\nArtefatos gerados: Tabela comparativa e Painel de Zooms (T1 e T2).\n")
}

# EXECUÇÃO
results = run_leakage_experiment()

# Cria diretorio artifacts com T no nome da pasta
artifacts_dir = sprintf(
  "output/sim3_N%d_T%d/artifacts", 
  LEAK_CONFIG$N_SUBJ, LEAK_CONFIG$T_PTS
  )
dir.create(artifacts_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  results, 
  file.path(
    artifacts_dir, 
    sprintf(
      "sim3_raw_results_N%d_T%d.rds", LEAK_CONFIG$N_SUBJ, 
      LEAK_CONFIG$T_PTS
      )
    )
  )

generate_leakage_artifacts(results)

