# scripts/sim1_generate_artifacts.R

library(data.table)
library(ggplot2)

# =============================================================================
# 1. SETUP
# =============================================================================

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

dir.create(
  "output/sim1_final/plots/pt",
  recursive = TRUE, showWarnings = FALSE
)
dir.create(
  "output/sim1_final/plots/eng",
  recursive = TRUE, showWarnings = FALSE
)

T_TREAT = 50.0

# =============================================================================
# 2. SISTEMA VISUAL
# =============================================================================

METHOD_COLORS = c(
  wavelet_causal = "royalblue",
  wavelet_causal_conf_low = "royalblue",
  wavelet_causal_conf_med = "royalblue",
  wavelet_causal_conf_high = "royalblue",
  wavelet_offline = "forestgreen",
  wavelet_offline_conf_low = "forestgreen",
  wavelet_offline_conf_med = "forestgreen",
  wavelet_offline_conf_high = "forestgreen",
  spline_ecker = "firebrick"
)

LT_MAIN = c(
  wavelet_causal = "solid",
  wavelet_offline = "dashed",
  spline_ecker = "dotted"
)

LT_CONF = c(
  wavelet_causal = "solid",
  wavelet_causal_conf_low = "dashed",
  wavelet_causal_conf_med = "dotdash",
  wavelet_causal_conf_high = "dotted",
  wavelet_offline = "solid",
  wavelet_offline_conf_low = "dashed",
  wavelet_offline_conf_med = "dotdash",
  wavelet_offline_conf_high = "dotted"
)

SH_MAIN = c(
  wavelet_causal = 17,
  wavelet_offline = 4,
  spline_ecker = 19
)

SH_CONF = c(
  wavelet_causal = 17,
  wavelet_causal_conf_low = 19,
  wavelet_causal_conf_med = 4,
  wavelet_causal_conf_high = 15,
  wavelet_offline = 17,
  wavelet_offline_conf_low = 19,
  wavelet_offline_conf_med = 4,
  wavelet_offline_conf_high = 15
)

LABELS = list(
  pt = list(
    methods = c(
      wavelet_causal = "Ond. causal",
      wavelet_causal_conf_low = "Ond. causal (conf. baixo)",
      wavelet_causal_conf_med = "Ond. causal (conf. médio)",
      wavelet_causal_conf_high = "Ond. causal (conf. alto)",
      wavelet_offline = "Ond. offline",
      wavelet_offline_conf_low = "Ond. offline (conf. baixo)",
      wavelet_offline_conf_med = "Ond. offline (conf. médio)",
      wavelet_offline_conf_high = "Ond. offline (conf. alto)",
      spline_ecker = "Spline (Ecker)"
    ),
    scenarios = c(
      Bumps = "Bumps",
      Doppler = "Doppler",
      Heavisine = "Escada",
      Sigmoid = "Sigmoide"
    ),
    dists = c(gaussian = "Gaussiano", t = "t (5 g.l.)"),
    x_T = "Pontos na grade (T)",
    x_t = "Tempo (t)",
    y_nrmse = "NRMSE (mediana ± quartis)",
    y_cov = "Cobertura simultânea empírica",
    y_ptwise_nrmse = "NRMSE pontual",
    y_ptwise_bias = "Viés pontual",
    y_ptwise_width = "Largura pontual da banda",
    y_ptwise_cov = "Cobertura pontual (%)",
    y_cfate_nrmse = "NRMSE pontual (CFATE)",
    y_cfate_cov = "Cobertura pontual CFATE (%)",
    y_fate = "FATE(t)",
    y_cfate = "CFATE(t)"
  ),
  eng = list(
    methods = c(
      wavelet_causal = "Wav. causal",
      wavelet_causal_conf_low = "Wav. causal (conf. low)",
      wavelet_causal_conf_med = "Wav. causal (conf. med.)",
      wavelet_causal_conf_high = "Wav. causal (conf. high)",
      wavelet_offline = "Wav. offline",
      wavelet_offline_conf_low = "Wav. offline (conf. low)",
      wavelet_offline_conf_med = "Wav. offline (conf. med.)",
      wavelet_offline_conf_high = "Wav. offline (conf. high)",
      spline_ecker = "Spline (Ecker)"
    ),
    scenarios = c(
      Bumps = "Bumps",
      Doppler = "Doppler",
      Heavisine = "Staircase",
      Sigmoid = "Sigmoid"
    ),
    dists = c(gaussian = "Gaussian", t = "t (5 d.f.)"),
    x_T = "Grid points (T)",
    x_t = "Time (t)",
    y_nrmse = "NRMSE (median ± quartiles)",
    y_cov = "Empirical simultaneous coverage",
    y_ptwise_nrmse = "Pointwise NRMSE",
    y_ptwise_bias = "Pointwise bias",
    y_ptwise_width = "Pointwise band width",
    y_ptwise_cov = "Pointwise coverage (%)",
    y_cfate_nrmse = "Pointwise NRMSE (CFATE)",
    y_cfate_cov  = "Pointwise CFATE coverage (%)",
    y_fate = "FATE(t)",
    y_cfate = "CFATE(t)"
  )
)

MAIN_METHODS = c("wavelet_causal", "wavelet_offline", "spline_ecker")

OFFLINE_CONF = c(
  "wavelet_offline", "wavelet_offline_conf_low",
  "wavelet_offline_conf_med", "wavelet_offline_conf_high"
)

CAUSAL_CONF = c(
  "wavelet_causal", "wavelet_causal_conf_low",
  "wavelet_causal_conf_med", "wavelet_causal_conf_high"
)

CASES = list(
  worst = list(N = 200,  T = 128),
  med   = list(N = 500,  T = 256),
  best  = list(N = 2000, T = 512)
)

base_theme = theme_bw(base_size = 11) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 9, color = "gray40"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(1.5, "cm"),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

# =============================================================================
# 3. DADOS
# =============================================================================

cat("Carregando dados...\n")
dt_raw = fread("output/sim1_final/sim1_raw.csv")
dt_csumm = fread("output/sim1_final/sim1_curve_summary.csv")
dt_ex = fread("output/sim1_final/sim1_exemplars.csv")

# NRMSE_t = RMSE_t / range(TrueCurve)
true_ranges = dt_ex[
  Method == "wavelet_offline",
  .(true_range = diff(range(TrueCurve))),
  by = .(Scenario, N, T_pts, ErrorDist)
]
dt_csumm = merge(
  dt_csumm, true_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist"), all.x = TRUE
)
dt_csumm[, NRMSE_t := RMSE_t / true_range]
dt_csumm[, true_range := NULL]

# CFATE_NRMSE_t (se disponível)
if ("CFATE_RMSE_t" %in% names(dt_csumm)) {
  cfate_ranges = dt_ex[
    t_val > T_TREAT,
    .(cfate_range = diff(range(CFATE_True))),
    by = .(Scenario, N, T_pts, ErrorDist)
  ]
  dt_csumm = merge(
    dt_csumm, cfate_ranges,
    by = c("Scenario", "N", "T_pts", "ErrorDist"), all.x = TRUE
  )
  dt_csumm[, CFATE_NRMSE_t := CFATE_RMSE_t / cfate_range]
  dt_csumm[, c("CFATE_RMSE_t", "cfate_range") := NULL]
}

has_cfate = all(c("CFATE_NRMSE", "CFATE_Cov_Simul") %in% names(dt_raw))
if (!has_cfate) {
  cat(
    "
    AVISO: métricas CFATE ausentes em sim1_raw.csv. 
    Execute sim1_add_cfate.R primeiro.\n
    "
  )
}

prep = function(dt, lbl) {
  dt[, Method_label := factor(
    lbl$methods[Method],
    levels = unname(lbl$methods)
  )]
  dt[, Scenario_label := factor(
    lbl$scenarios[Scenario],
    levels = unname(lbl$scenarios)
  )]
  dt[, ErrorDist_label := factor(
    lbl$dists[ErrorDist],
    levels = unname(lbl$dists)
  )]
  dt[, N_label := factor(
    paste0("N = ", N),
    levels = paste0("N = ", c(200, 500, 2000))
  )]
  dt[, T_fac := factor(T_pts, levels = c(128, 256, 512))]
  dt
}

# =============================================================================
# 4. FUNÇÕES
# =============================================================================

save_plot = function(p, fname, lang, width = 9, height = 6) {
  path = file.path("output/sim1_final/plots", lang, fname)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  cat(sprintf("  -> %s/%s\n", lang, fname))
}

plot_nrmse = function(dt, methods, lt, sh, lbl) {
  dt_s = dt[Method %in% methods]
  dt_s[, Method := factor(Method, levels = methods)]
  ggplot(
    dt_s, 
    aes(
      x = T_fac, y = NRMSE,
      color = Method, linetype = Method, shape = Method, group = Method
    )
  ) +
  stat_summary(
    fun = median, geom = "line",
    position = position_dodge(0.4), linewidth = 0.7
  ) +
  stat_summary(
    fun = median, geom = "point",
    position = position_dodge(0.4), size = 2.5
  ) +
  stat_summary(
    fun.min = \(x) quantile(x, 0.25),
    fun.max = \(x) quantile(x, 0.75),
    geom = "errorbar",
    position = position_dodge(0.4), width = 0.25
  ) +
  facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
  scale_color_manual(
    values = METHOD_COLORS[methods], labels = lbl$methods[methods]
  ) +
  scale_linetype_manual(
    values = lt[methods], labels = lbl$methods[methods]
  ) +
  scale_shape_manual(
    values = sh[methods], labels = lbl$methods[methods]
  ) +
  labs(x = lbl$x_T, y = lbl$y_nrmse) +
  base_theme
}

plot_cfate_nrmse = function(dt, methods, lt, sh, lbl) {
  dt_s = dt[Method %in% methods]
  dt_s[, Method := factor(Method, levels = methods)]
  ggplot(dt_s, aes(
    x = T_fac, y = CFATE_NRMSE,
    color = Method, linetype = Method, shape = Method, group = Method
  )) +
    stat_summary(
      fun = median, geom = "line",
      position = position_dodge(0.4), linewidth = 0.7
    ) +
    stat_summary(
      fun = median, geom = "point",
      position = position_dodge(0.4), size = 2.5
    ) +
    stat_summary(
      fun.min = \(x) quantile(x, 0.25),
      fun.max = \(x) quantile(x, 0.75),
      geom = "errorbar",
      position = position_dodge(0.4), width = 0.25
    ) +
    facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
    scale_color_manual(
      values = METHOD_COLORS[methods], labels = lbl$methods[methods]
    ) +
    scale_linetype_manual(
      values = lt[methods], labels = lbl$methods[methods]
    ) +
    scale_shape_manual(
      values = sh[methods], labels = lbl$methods[methods]
    ) +
    labs(x = lbl$x_T, y = lbl$y_nrmse) +
    base_theme
}

plot_coverage = function(dt, methods, lt, sh, lbl, col = "Cov_Simul") {
  dt_cov = dt[
    Method %in% methods, 
    .(
      Cov = mean(get(col), na.rm = TRUE)
    ),
    by = .(ErrorDist_label, N_label, T_fac, Method)
  ]
  dt_cov[, Method := factor(Method, levels = methods)]
  
  ggplot(
    dt_cov,
    aes(
      x = T_fac, y = Cov,
      color = Method, linetype = Method, shape = Method, group = Method
    )
  ) +
  geom_hline(
    yintercept = 0.95, color = "gray40",
    linewidth = 0.4, linetype = "dashed"
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
  scale_color_manual(
    values = METHOD_COLORS[methods], labels = lbl$methods[methods]
  ) +
  scale_linetype_manual(
    values = lt[methods], labels = lbl$methods[methods]
  ) +
  scale_shape_manual(
    values = sh[methods], labels = lbl$methods[methods]
  ) +
  scale_y_continuous(
    limits = c(0, 1.01),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(x = lbl$x_T, y = lbl$y_cov) +
  base_theme
}

plot_pointwise = function(
  dt_cs, metric, ylab, methods, lt, sh, lbl,
  N_sel, T_sel, add_hline = NULL
) {
  dt_s = dt_cs[Method %in% methods & N == N_sel & T_pts == T_sel]
  dt_s[, Method := factor(Method, levels = methods)]

  n_pts = dt_s[, .N, by = .(Scenario_label, ErrorDist_label, Method)][, max(N)]
  every = max(1L, floor(n_pts / 12))
  dt_pts = dt_s[, .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, ErrorDist_label, Method)
  ]

  p = ggplot(
    dt_s, 
    aes(
      x = t_val, y = .data[[metric]],
      color = Method, linetype = Method
    )
  ) +
  geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
  geom_line(linewidth = 0.7) +
  geom_point(data = dt_pts, aes(shape = Method), size = 2.2) +
  facet_grid(Scenario_label ~ ErrorDist_label, scales = "free_y") +
  scale_color_manual(
    values = METHOD_COLORS[methods], labels = lbl$methods[methods]
  ) +
  scale_linetype_manual(
    values = lt[methods], labels = lbl$methods[methods]
  ) +
  scale_shape_manual(
    values = sh[methods], labels = lbl$methods[methods]
  ) +
  labs(
    x = lbl$x_t, y = ylab,
    subtitle = sprintf("N = %d,  T = %d", N_sel, T_sel)
  ) +
  base_theme
  
  if (!is.null(add_hline)) {
    p = p + geom_hline(
      yintercept = add_hline, color = "gray40",
      linewidth = 0.4, linetype = "dashed"
    )
  }
  p
}

plot_exemplar = function(
  dt_e, methods, lt, sh, lbl, N_sel, T_sel,
  dist_sel = "gaussian",
  curve_col = "TrueCurve", est_col = "Est",
  lower_col = "Lower", upper_col = "Upper",
  ylab = NULL) {
  if (is.null(ylab)) ylab = lbl$y_fate
  dt_s = dt_e[
    Method %in% methods & N == N_sel &
      T_pts == T_sel & ErrorDist == dist_sel
  ]
  dt_s = dt_s[
    dt_s[, .I[SimID == min(SimID)], by = .(Scenario, Method)]$V1
  ]
  dt_s[, Method := factor(Method, levels = methods)]

  n_pts = dt_s[, .N, by = .(Scenario_label, Method)][, max(N)]
  every = max(1L, floor(n_pts / 12))
  dt_pts = dt_s[, .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, Method)
  ]

  ggplot(dt_s, aes(x = t_val)) +
    geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
    geom_ribbon(
      aes(ymin = .data[[lower_col]], ymax = .data[[upper_col]], fill = Method),
      alpha = 0.2
    ) +
    geom_line(
      aes(y = .data[[curve_col]]),
      color = "black", linewidth = 0.55, linetype = "longdash"
    ) +
    geom_line(
      aes(y = .data[[est_col]], color = Method, linetype = Method),
      linewidth = 0.8
    ) +
    geom_point(
      data = dt_pts,
      aes(y = .data[[est_col]], color = Method, shape = Method), size = 2.2
    ) +
    facet_grid(Scenario_label ~ Method_label, scales = "free_y") +
    scale_color_manual(values = METHOD_COLORS[methods]) +
    scale_fill_manual(values = METHOD_COLORS[methods]) +
    scale_linetype_manual(values = lt[methods]) +
    scale_shape_manual(values = sh[methods]) +
    labs(
      x = lbl$x_t, y = ylab,
      subtitle = sprintf(
        "%s  |  N = %d,  T = %d",
        lbl$dists[dist_sel], N_sel, T_sel
      )
    ) +
    base_theme +
    theme(legend.position = "none")
}

# =============================================================================
# 5. GERAÇÃO
# =============================================================================

for (lang in c("pt", "eng")) {
  lbl = LABELS[[lang]]

  dt_r = prep(copy(dt_raw), lbl)
  dt_cs = prep(copy(dt_csumm), lbl)
  dt_e = prep(copy(dt_ex), lbl)

  # --- 5.1 NRMSE e cobertura (3 métodos principais) ---
  cat(sprintf("\n[%s | 5.1] NRMSE e cobertura...\n", lang))
  for (scen in names(lbl$scenarios)) {
    dt_s = dt_r[Scenario == scen]
    save_plot(
      plot_nrmse(dt_s, MAIN_METHODS, LT_MAIN, SH_MAIN, lbl),
      sprintf("nrmse_%s.png", scen), lang
    )
    save_plot(
      plot_coverage(dt_s, MAIN_METHODS, LT_MAIN, SH_MAIN, lbl),
      sprintf("coverage_%s.png", scen), lang
    )
    if (has_cfate) {
      save_plot(
        plot_cfate_nrmse(dt_s, MAIN_METHODS, LT_MAIN, SH_MAIN, lbl),
        sprintf("cfate_nrmse_%s.png", scen), lang
      )
      save_plot(
        plot_coverage(
          dt_s, MAIN_METHODS, LT_MAIN,
          SH_MAIN, lbl, col = "CFATE_Cov_Simul"
        ),
        sprintf("cfate_coverage_%s.png", scen), lang
      )
    }
  }

  # --- 5.2 Curvas ponto-a-ponto ---
  cat(sprintf("\n[%s | 5.2] Ptwise NRMSE, viés, largura...\n", lang))
  for (case_name in names(CASES)) {
    N_sel = CASES[[case_name]]$N
    T_sel = CASES[[case_name]]$T
    save_plot(
      plot_pointwise(
        dt_cs, "NRMSE_t", lbl$y_ptwise_nrmse,
        MAIN_METHODS, LT_MAIN, SH_MAIN, lbl, N_sel, T_sel
      ),
      sprintf("ptwise_nrmse_%s.png", case_name), lang, width = 9, height = 7
    )
    save_plot(
      plot_pointwise(
        dt_cs, "Bias_t", lbl$y_ptwise_bias,
        MAIN_METHODS, LT_MAIN, SH_MAIN, lbl, N_sel, T_sel
      ),
      sprintf("ptwise_bias_%s.png", case_name), lang, width = 9, height = 7
    )
    save_plot(
      plot_pointwise(
        dt_cs, "Width_t", lbl$y_ptwise_width,
        MAIN_METHODS, LT_MAIN, SH_MAIN, lbl, N_sel, T_sel
      ),
      sprintf("ptwise_width_%s.png", case_name), lang, width = 9, height = 7
    )
    if (has_cfate && "CFATE_NRMSE_t" %in% names(dt_cs)) {
      save_plot(
        plot_pointwise(
          dt_cs, "CFATE_NRMSE_t", lbl$y_cfate_nrmse,
          MAIN_METHODS, LT_MAIN, SH_MAIN, lbl, N_sel, T_sel
        ),
        sprintf("ptwise_cfate_nrmse_%s.png", case_name),
        lang, width = 9, height = 7
      )
      save_plot(
        plot_pointwise(
          dt_cs, "CFATE_Cov_t", lbl$y_cfate_cov,
          MAIN_METHODS, LT_MAIN, SH_MAIN, lbl, N_sel, T_sel,
          add_hline = 0.95
        ),
        sprintf("ptwise_cfate_cov_%s.png", case_name), 
        lang, width = 9, height = 7
      )
    }
  }

  # --- 5.3 Exemplares: 3 métodos principais (FATE) ---
  cat(sprintf("\n[%s | 5.3] Exemplares (FATE)...\n", lang))
  for (case_name in names(CASES)) {
    for (dist_sel in c("gaussian", "t")) {
      save_plot(
        plot_exemplar(
          dt_e, MAIN_METHODS, LT_MAIN, SH_MAIN, lbl,
          CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel,
          ylab = lbl$y_fate
        ),
        sprintf("exemplar_%s_%s.png", case_name, dist_sel), lang,
        width = 10, height = 9
      )
    }
  }

  # --- 5.4 Exemplares: CFATE ---
  if (has_cfate) {
    cat(sprintf("\n[%s | 5.4] Exemplares (CFATE)...\n", lang))
    for (case_name in names(CASES)) {
      for (dist_sel in c("gaussian", "t")) {
        save_plot(
          plot_exemplar(
            dt_e, MAIN_METHODS, LT_MAIN, SH_MAIN, lbl,
            CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel,
            curve_col = "CFATE_True", est_col = "CFATE_Est",
            lower_col = "CFATE_Lower", upper_col = "CFATE_Upper",
            ylab = lbl$y_cfate
          ),
          sprintf("cfate_exemplar_%s_%s.png", case_name, dist_sel), lang,
          width = 10, height = 9
        )
      }
    }
  }

  # --- 5.5 Exemplares: confundimento ---
  cat(sprintf("\n[%s | 5.5] Exemplares (confundimento)...\n", lang))
  for (case_name in names(CASES)) {
    for (dist_sel in c("gaussian", "t")) {
      save_plot(
        plot_exemplar(
          dt_e, OFFLINE_CONF, LT_CONF, SH_CONF, lbl,
          CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel,
          ylab = lbl$y_fate
        ),
        sprintf("conf_off_exemplar_%s_%s.png", case_name, dist_sel), lang,
        width = 12, height = 9
      )
      save_plot(
        plot_exemplar(
          dt_e, CAUSAL_CONF, LT_CONF, SH_CONF, lbl,
          CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel,
          ylab = lbl$y_fate
        ),
        sprintf("conf_cau_exemplar_%s_%s.png", case_name, dist_sel), lang,
        width = 12, height = 9
      )
    }
  }
}

cat("\nGeração de artefatos concluída.\n")
