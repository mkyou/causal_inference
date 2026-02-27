library(data.table)
library(future.apply)
library(progressr)
library(ggplot2)
library(rLifting)
library(patchwork)
library(kableExtra)

source("R/scenario_runner.R")
source("R/process_functional_data.R")
source("R/estimate_fate.R")
source("R/utils_dgp_functions.R")
source("R/generate_sequential_dgp.R")
source("R/generate_functional_outcomes.R")
source("R/config_ecker.R")

# CONFIGURAÇÕES
RT_CONFIG = list(
  N_SIM = 1000,       
  BATCH_SIZE = 50,
  N_SUBJ = 200,
  T_TOTAL = 1024,
  T_TREAT = 256,
  WINDOW_SZ = 64,
  
  # Configuração rLifting
  SCHEME = "db2",
  THRESH = "semisoft",
  ALPHA = 0.3,
  BETA = 1.2
)

SCENARIOS = list(
  list(name = "Bumps", mag = 3),
  list(name = "Doppler", mag = 5)
)

# Diretórios
dir.create("output/checkpoints_sim2", recursive = TRUE, showWarnings = FALSE)
dir.create("output/sim2_artifacts", recursive = TRUE, showWarnings = FALSE)

# Paralelismo
plan(multisession)
handlers(global = TRUE)

# Estilo visual
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

# Prepara pesos estáticos (estimador linear)
prepare_static_weights = function(dgp_data, t_idx) {
  X_cols = grep("^X[0-9]+", names(dgp_data), value=TRUE)
  X_raw = as.matrix(dgp_data[t == t_idx, ..X_cols])
  X_des = cbind(1, X_raw)
  A_vec = dgp_data[t == t_idx, Trat_1]
  
  idx1 = which(A_vec == 1); idx0 = which(A_vec == 0)
  x_bar = colMeans(X_des)
  
  get_w = function(idx) {
    if(length(idx) == 0) return(numeric(0))
    X_sub = X_des[idx, , drop=FALSE]
    XtX = crossprod(X_sub); diag(XtX) = diag(XtX) + 1e-8
    inv = tryCatch(solve(XtX), error = function(e) return(NULL))
    if(is.null(inv)) return(rep(0, length(idx)))
    as.numeric(t(x_bar) %*% inv %*% t(X_sub))
  }
  
  list(w1 = get_w(idx1), w0 = get_w(idx0), idx1 = idx1, idx0 = idx0)
}

# Worker Unitário
run_single_simulation = function(sim_id, scen, conf) {
  time_grid = seq(20, 100, length.out = conf$T_TOTAL)
  t_treat_val = time_grid[conf$T_TREAT]
  
  cfg = get_ecker_configuration(conf$N_SUBJ, conf$T_TOTAL)
  
  dgp_covs = generate_sequential_dgp(
    conf$N_SUBJ, time_grid, c(t_treat_val),
    cfg$baseline_func_templates, cfg$baseline_params_list,
    cfg$beta_funcs_matrix, cfg$gamma_matrix, cfg$theta_matrix,
    rep(0, 9), cfg$static_indices
  )
  
  eff_func = switch(
    scen$name, 
    "Bumps" = effect_bumps, 
    "Doppler" = effect_doppler
    )
  target_f = function(t, s) eff_func(t, s, 20, 100, scen$mag)
  
  dgp_out = generate_functional_outcome(
    dgp_covs, cfg$baseline_y_func, list(target_f), 
    cfg$concurrent_effect_funcs, c(t_treat_val), noise_sd = 0.5
  )
  
  # Matriz de entrada (N x T)
  Y_full = as.matrix(dcast(dgp_out, id ~ t, value.var="y")[, -1])
  
  # Verdade
  True_Curve = target_f(time_grid, t_treat_val)
  True_Range = diff(range(True_Curve)); if(True_Range == 0) True_Range = 1.0
  
  # Inicialização streaming (rLifting)
  weights = prepare_static_weights(dgp_covs, time_grid[1])
  scheme = rLifting::lifting_scheme(conf$SCHEME)
  max_levels = floor(log2(conf$WINDOW_SZ))
  
  # Cria os processadores de fluxo
  streams = lapply(1:conf$N_SUBJ, function(i) {
    rLifting::new_wavelet_stream(
      scheme, window_size = conf$WINDOW_SZ, levels = max_levels,
      method = conf$THRESH, alpha = conf$ALPHA, beta = conf$BETA
    )
  })
  
  # Pré-alocação de métricas
  metrics_mat = matrix(NA, nrow = conf$T_TOTAL, ncol = 3) 
  colnames(metrics_mat) = c("Estimate", "Latency_us", "Time_Idx")
  
  # Loop de tempo real
  for(t in 1:conf$T_TOTAL) {
    y_input_vec = Y_full[, t]
    
    start_clock = Sys.time()
    
    # Processamento paralelo simulado (vetorizado no R)
    # Na prática, cada stream processaria seu dado
    y_smooth_t = numeric(conf$N_SUBJ)
    
    # O gargalo aqui é o R loop, 
    # mas serve para medir a latência do ALGORITMO
    for(i in 1:conf$N_SUBJ) {
      y_smooth_t[i] = streams[[i]](as.numeric(y_input_vec[i]))
    }
    
    # Estimativa instantânea (regressão ponderada)
    est_val = sum(weights$w1 * y_smooth_t[weights$idx1]) - 
      sum(weights$w0 * y_smooth_t[weights$idx0])
    
    latency = as.numeric(Sys.time() - start_clock) * 1e6 # microsegundos
    
    metrics_mat[t, ] = c(est_val, latency, t)
  }
  
  # Formatação de saída
  dt_metrics = as.data.table(metrics_mat)
  dt_metrics[, `:=`(
    SimID = sim_id,
    Scenario = scen$name,
    Truth = True_Curve,
    Error = Estimate - True_Curve
  )]
  
  # Reduzindo tamanho: não precisamos salvar tudo para a métrica global
  global_stats = data.table(
    SimID = sim_id,
    Scenario = scen$name,
    NMRSE = sqrt(mean(dt_metrics$Error^2)) / True_Range,
    MaxError = max(abs(dt_metrics$Error)),
    Mean_Latency_us = mean(dt_metrics$Latency_us),
    Max_Latency_us = max(dt_metrics$Latency_us)
  )
  
  return(list(ts = dt_metrics, global = global_stats))
}

# EXECUTOR COM BATCH
run_realtime_simulation = function() {
  grid_tasks = expand.grid(
    SimID = 1:RT_CONFIG$N_SIM, 
    ScenIdx = 1:length(SCENARIOS)
  )
  setDT(grid_tasks)
  
  n_total = nrow(grid_tasks)
  n_batches = ceiling(n_total / RT_CONFIG$BATCH_SIZE)
  
  cat(
    sprintf(
      "\n>>> Iniciando Simulação Real-Time (Total: %d, Lotes: %d) <<<\n", 
      n_total, 
      n_batches
      )
    )
  
  with_progress({
    p <- progressor(steps = n_total)
    
    ckpt_dir = file.path("output/checkpoints_sim2", sprintf("N%d", RT_CONFIG$N_SUBJ))
    dir.create(ckpt_dir, recursive=TRUE, showWarnings=FALSE)
    
    for(b in 1:n_batches) {
      ckpt_file = file.path(ckpt_dir, sprintf("batch_%03d.rds", b))
      
      if(file.exists(ckpt_file)) {
        p(amount = RT_CONFIG$BATCH_SIZE)
        next
      }
      
      # Índices do lote atual
      idx_start = (b - 1) * RT_CONFIG$BATCH_SIZE + 1
      idx_end = min(b * RT_CONFIG$BATCH_SIZE, n_total)
      current_indices = idx_start:idx_end
      
      # Execução paralela do lote
      batch_results = future_lapply(current_indices, function(i) {
        row = grid_tasks[i]
        scen = SCENARIOS[[row$ScenIdx]]
        res = run_single_simulation(row$SimID, scen, RT_CONFIG)
        p()
        return(res)
      }, future.seed = TRUE)
      
      # Separa e comprime
      batch_ts = rbindlist(lapply(batch_results, `[[`, "ts"))
      batch_global = rbindlist(lapply(batch_results, `[[`, "global"))
      
      saveRDS(list(ts = batch_ts, global = batch_global), ckpt_file)
      
      # Limpeza de memória
      rm(batch_results, batch_ts, batch_global)
      gc()
    }
  })
  
  cat("\n>>> Agregando Resultados Finais... <<<\n")
  
  files = list.files(
    ckpt_dir, 
    full.names = TRUE, 
    pattern = "\\.rds$"
    )
  all_data = lapply(files, readRDS)
  
  # Junta tudo
  final_ts = rbindlist(lapply(all_data, `[[`, "ts"))
  final_global = rbindlist(lapply(all_data, `[[`, "global"))
  
  # Salva CSVs brutos
  # Cria diretorio artifacts com base em N e T
  artifacts_dir = sprintf(
    "output/sim2_N%d_T%d/artifacts", 
    RT_CONFIG$N_SUBJ, 
    RT_CONFIG$WINDOW_SZ
    )
  dir.create(artifacts_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(
    final_ts, 
    file.path(
      artifacts_dir, 
      sprintf(
        "sim2_time_series_N%d_T%d.csv", 
        RT_CONFIG$N_SUBJ, 
        RT_CONFIG$WINDOW_SZ
        )
      )
    )
  
  fwrite(
    final_global, 
    file.path(
      artifacts_dir, 
      sprintf(
        "sim2_global_stats_N%d_T%d.csv", 
        RT_CONFIG$N_SUBJ, RT_CONFIG$WINDOW_SZ
        )
      )
    )
  
  # Limpa checkpoints temporários
  unlink(ckpt_dir, recursive = TRUE)
  
  return(list(ts = final_ts, global = final_global))
}

# Gerador de artefatos da simulação
generate_artifacts_paper = function(res_list) {
  dt_ts = res_list$ts
  dt_global = res_list$global
  colors_scen = c("Bumps" = "firebrick", "Doppler" = "royalblue")
  
  # output dir já deve existir ou criamos de novo para garantir
  artifacts_dir = sprintf(
    "output/sim2_N%d_T%d/artifacts", 
    RT_CONFIG$N_SUBJ, 
    RT_CONFIG$WINDOW_SZ
    )
  
  dir.create(artifacts_dir, showWarnings = FALSE, recursive = TRUE)

  # Após processar todos os dados, qual a performance média?
  tbl_summary = dt_global[, .(
    # Média e SD do NRMSE Global de cada simulação
    Mean_NRMSE = mean(NMRSE),
    SD_NRMSE = sd(NMRSE),
    
    # Média e SD do Erro Máximo Global
    Mean_MaxErr = mean(MaxError),
    SD_MaxErr = sd(MaxError),
    
    # Latência
    Mean_Lat = mean(Mean_Latency_us),
    SD_Lat = sd(Mean_Latency_us)
  ), by = Scenario]
  
  tbl_formatted = tbl_summary |>
    dplyr::mutate(
      `NRMSE (SD)` = sprintf("%.4f (%.4f)", Mean_NRMSE, SD_NRMSE),
      `Max Error (SD)` = sprintf("%.2f (%.2f)", Mean_MaxErr, SD_MaxErr),
      `Latência µs (SD)` = sprintf("%.0f (%.0f)", Mean_Lat, SD_Lat)
    ) |>
    dplyr::select(
      Scenario, `NRMSE (SD)`, `Max Error (SD)`, `Latência µs (SD)`
      )
  
  kbl_obj = kable(
    tbl_formatted, format = "html", booktabs = TRUE,
    caption = "Performance global da estimação em tempo real"
  ) |>
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed"), 
      full_width = F
      )
  
  save_kable(
    kbl_obj, 
    file.path(
      artifacts_dir, 
      sprintf(
        "table_realtime_metrics_N%d_T%d.html", 
        RT_CONFIG$N_SUBJ, 
        RT_CONFIG$WINDOW_SZ
        )
      )
    )
  
  # Snapshots da trajetória
  time_map = seq(20, 100, length.out = max(dt_ts$Time_Idx))
  dt_traj_mean = dt_ts[, .(
    Mean_Est = mean(Estimate),
    True_Val = mean(Truth)
  ), by = .(Scenario, Time_Idx)]
  
  dt_traj_mean[, Time_Real := time_map[Time_Idx]]
  
  cutoffs_idx = c(250, 500, 750, 1024)
  plot_list = list()
  
  for(cut in cutoffs_idx) {
    data_view = dt_traj_mean[Time_Idx <= cut]
    t_real_cut = time_map[cut]
    
    p = ggplot(data_view, aes(x = Time_Real)) +
      geom_line(
        aes(y = True_Val), color = "black", 
        linetype = "dashed", linewidth = 0.6
        ) +
      geom_line(
        aes(y = Mean_Est, color = Scenario), 
        linewidth = 0.9, alpha = 0.7
        ) +
      scale_color_manual(values = colors_scen) +
      coord_cartesian(xlim = c(20, 100), ylim = c(-2, 6)) +
      facet_wrap(~Scenario, ncol = 1, scales = "free_y") +
      labs(
        subtitle = sprintf(
          "Tempo: %.1f (nº pontos observados=%d)", 
          t_real_cut, cut
          ),
        x = NULL, y = "FATE"
      ) +
      theme_paper +
      theme(
        legend.position = "none",
        plot.subtitle = element_text(size = 9, color = "grey30")
      )
    
    plot_list[[as.character(cut)]] = p
  }
  
  p_combined = (
    plot_list[[1]] | plot_list[[2]] | plot_list[[3]] | plot_list[[4]]
    ) +
    plot_annotation(
      theme = theme(text = element_text(family = "serif", size = 14))
    )
  
  ggsave(
    file.path(
      artifacts_dir, 
      sprintf(
        "plot_snapshots_trajectory_N%d_T%d.png", 
        RT_CONFIG$N_SUBJ, RT_CONFIG$WINDOW_SZ
        )
      ), 
    p_combined, width = 14, height = 7
    )
  
  # NRMSE acumulado; evolução da qualidade média da simulação
  dt_ts[, SqErr := Error^2]
  dt_mse_t = dt_ts[, .(MSE_t = mean(SqErr)), by = .(Scenario, Time_Idx)]
  dt_mse_t[, CumMSE := cumsum(MSE_t) / (1:.N), by = Scenario]
  
  ranges = dt_ts[, .(R = diff(range(Truth))), by = Scenario]
  dt_mse_t = merge(dt_mse_t, ranges, by="Scenario")
  dt_mse_t[, CumNRMSE := sqrt(CumMSE) / R]
  dt_mse_t[, Time_Real := time_map[Time_Idx]]
  
  p_cum_nrmse = ggplot(
    dt_mse_t, aes(x = Time_Real, y = CumNRMSE, color = Scenario)
    ) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = colors_scen) +
    labs(
      y = "NRMSE acumulado", x = "Tempo"
    ) +
    facet_wrap(~Scenario) +
    theme_paper +
    theme(legend.position = "none")
  
  ggsave(
    file.path(
      artifacts_dir, 
      sprintf(
        "plot_nrmse_cumulative_N%d_T%d.png", 
        RT_CONFIG$N_SUBJ, RT_CONFIG$WINDOW_SZ
        )
      ), 
    p_cum_nrmse, width = 10, height = 5
    )
  
  # Latência
  limit = quantile(dt_ts$Latency_us, 0.995)
  dt_lat = dt_ts[Latency_us < limit]
  dt_lat_agg = dt_lat[
    , 
    .(Mean_Lat = mean(Latency_us)), 
    by = .(Scenario, Time_Idx)
    ]
  dt_lat_agg[, Time_Real := time_map[Time_Idx]]
  
  p_lat = ggplot(dt_lat_agg, aes(x = Time_Real, y = Mean_Lat)) +
    geom_point(aes(color = Scenario), size = 0.5, alpha = 0.3) +
    geom_smooth(
      aes(color = Scenario), 
      method = "loess", span = 0.1, se = FALSE, 
      color="black", linewidth=0.6
      ) +
    scale_color_manual(values = colors_scen) +
    labs(
      y = "Latência (µs)", x = "Tempo"
    ) +
    facet_wrap(~Scenario, ncol=1, scales="free_y") +
    theme_paper +
    theme(legend.position = "none")
  
  ggsave(
    file.path(
      artifacts_dir, 
      sprintf(
        "plot_latency_step_N%d_T%d.png", 
        RT_CONFIG$N_SUBJ, RT_CONFIG$WINDOW_SZ
        )
      ), 
    p_lat, width = 8, height = 6
    )
}

results = run_realtime_simulation()
generate_artifacts_paper(results)
