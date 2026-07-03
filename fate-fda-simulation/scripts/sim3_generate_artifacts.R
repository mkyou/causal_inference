# scripts/sim3_generate_artifacts.R

library(data.table)
library(ggplot2)

# =============================================================================
# 1. SETUP
# =============================================================================

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

dir.create(
  "output/sim3_final/plots/pt",
  recursive = TRUE, showWarnings = FALSE
)
dir.create(
  "output/sim3_final/plots/eng",
  recursive = TRUE, showWarnings = FALSE
)

T_TREAT <- 50.0

# =============================================================================
# 2. SISTEMA VISUAL
# =============================================================================

IRR_COLORS <- c(
  leve = "forestgreen",
  moderado = "royalblue",
  severo = "firebrick"
)

IRR_LINETYPES <- c(
  leve = "solid",
  moderado = "dashed",
  severo = "dotted"
)

IRR_SHAPES <- c(
  leve = 17,
  moderado = 19,
  severo = 15
)

IRR_LEVELS <- c("leve", "moderado", "severo")

CASES <- list(
  worst = list(N = 200, T = 128),
  best = list(N = 2000, T = 512)
)

LABELS <- list(
  pt = list(
    scenarios = c(
      Bumps = "Bumps",
      Doppler = "Doppler",
      Heavisine = "Escada",
      Sigmoid = "Sigmoide"
    ),
    dists = c(gaussian = "Gaussiano", t = "t (5 g.l.)"),
    irr = c(
      leve = "Irregularidade leve",
      moderado = "Irregularidade moderada",
      severo = "Irregularidade severa"
    ),
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
    scenarios = c(
      Bumps = "Bumps",
      Doppler = "Doppler",
      Heavisine = "Staircase",
      Sigmoid = "Sigmoid"
    ),
    dists = c(gaussian = "Gaussian", t = "t (5 d.f.)"),
    irr = c(
      leve = "Mild irregularity",
      moderado = "Moderate irregularity",
      severo = "Severe irregularity"
    ),
    x_T = "Grid points (T)",
    x_t = "Time (t)",
    y_nrmse = "NRMSE (median ± quartiles)",
    y_cov = "Empirical simultaneous coverage",
    y_ptwise_nrmse = "Pointwise NRMSE",
    y_ptwise_bias = "Pointwise bias",
    y_ptwise_width = "Pointwise band width",
    y_ptwise_cov = "Pointwise coverage (%)",
    y_cfate_nrmse = "Pointwise NRMSE (CFATE)",
    y_cfate_cov = "Pointwise CFATE coverage (%)",
    y_fate = "FATE(t)",
    y_cfate = "CFATE(t)"
  )
)

base_theme <- theme_bw(base_size = 11) +
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
dt_raw <- fread("output/sim3_final/sim3_raw.csv")
dt_csumm <- fread("output/sim3_final/sim3_curve_summary.csv")
dt_ex <- fread("output/sim3_final/sim3_exemplars.csv")

# NRMSE_t = RMSE_t / range(TrueCurve)
true_ranges <- dt_ex[
  ,
  .(true_range = diff(range(TrueCurve))),
  by = .(Scenario, N, T_pts, ErrorDist, Irregularity)
]
dt_csumm <- merge(
  dt_csumm, true_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist", "Irregularity"), all.x = TRUE
)
dt_csumm[, NRMSE_t := RMSE_t / true_range]
dt_csumm[, true_range := NULL]

# CFATE_NRMSE_t = CFATE_RMSE_t / range(TrueCFATE)
cfate_ranges <- dt_ex[
  t_val > T_TREAT,
  .(cfate_range = diff(range(TrueCFATE))),
  by = .(Scenario, N, T_pts, ErrorDist, Irregularity)
]
dt_csumm <- merge(
  dt_csumm, cfate_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist", "Irregularity"), all.x = TRUE
)
dt_csumm[, CFATE_NRMSE_t := CFATE_RMSE_t / cfate_range]
dt_csumm[, c("CFATE_RMSE_t", "cfate_range") := NULL]

prep <- function(dt, lbl) {
  dt[, Irregularity := factor(Irregularity, levels = IRR_LEVELS)]
  dt[, Irr_label := factor(
    lbl$irr[as.character(Irregularity)],
    levels = unname(lbl$irr)
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
    levels = paste0("N = ", c(200, 2000))
  )]
  dt[, T_fac := factor(T_pts, levels = c(128, 512))]
  dt
}

# =============================================================================
# 4. FUNÇÕES
# =============================================================================

save_plot <- function(p, fname, lang, width = 9, height = 6) {
  path <- file.path("output/sim3_final/plots", lang, fname)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  cat(sprintf("  -> %s/%s\n", lang, fname))
}

plot_nrmse <- function(dt, lbl, metric = "NRMSE") {
  ggplot(dt, aes(
    x = T_fac, y = .data[[metric]],
    color = Irregularity, linetype = Irregularity,
    shape = Irregularity, group = Irregularity
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
    scale_color_manual(values = IRR_COLORS, labels = lbl$irr) +
    scale_linetype_manual(values = IRR_LINETYPES, labels = lbl$irr) +
    scale_shape_manual(values = IRR_SHAPES, labels = lbl$irr) +
    labs(x = lbl$x_T, y = lbl$y_nrmse) +
    base_theme
}

plot_coverage <- function(dt, lbl, metric = "Cov_Simul") {
  dt_cov <- dt[, .(Cov = mean(get(metric), na.rm = TRUE)),
    by = .(ErrorDist_label, N_label, T_fac, Irregularity)
  ]
  ggplot(
    dt_cov,
    aes(
      x = T_fac, y = Cov / 100,
      color = Irregularity, linetype = Irregularity,
      shape = Irregularity, group = Irregularity
    )
  ) +
    geom_hline(
      yintercept = 0.95, color = "gray40",
      linewidth = 0.4, linetype = "dashed"
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
    scale_color_manual(values = IRR_COLORS, labels = lbl$irr) +
    scale_linetype_manual(values = IRR_LINETYPES, labels = lbl$irr) +
    scale_shape_manual(values = IRR_SHAPES, labels = lbl$irr) +
    scale_y_continuous(
      limits = c(0, 1.01),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(x = lbl$x_T, y = lbl$y_cov) +
    base_theme
}

plot_pointwise <- function(
  dt_cs, metric, ylab, lbl,
  N_sel, T_sel, add_hline = NULL
) {
  dt_s <- dt_cs[N == N_sel & T_pts == T_sel]

  n_pts <- dt_s[
    , .N,
    by = .(Scenario_label, ErrorDist_label, Irregularity)
  ][, max(N)]
  every <- max(1L, floor(n_pts / 12))
  dt_pts <- dt_s[
    , .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, ErrorDist_label, Irregularity)
  ]

  p <- ggplot(
    dt_s,
    aes(
      x = t_val, y = .data[[metric]],
      color = Irregularity, linetype = Irregularity
    )
  ) +
    geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
    geom_line(linewidth = 0.7) +
    geom_point(data = dt_pts, aes(shape = Irregularity), size = 2.2) +
    facet_grid(Scenario_label ~ ErrorDist_label, scales = "free_y") +
    scale_color_manual(values = IRR_COLORS, labels = lbl$irr) +
    scale_linetype_manual(values = IRR_LINETYPES, labels = lbl$irr) +
    scale_shape_manual(values = IRR_SHAPES, labels = lbl$irr) +
    labs(
      x = lbl$x_t, y = ylab,
      subtitle = sprintf("N = %d,  T = %d", N_sel, T_sel)
    ) +
    base_theme
  if (!is.null(add_hline)) {
    p <- p + geom_hline(
      yintercept = add_hline, color = "gray40",
      linewidth = 0.4, linetype = "dashed"
    )
  }
  p
}

plot_exemplar <- function(
  dt_e, lbl, N_sel, T_sel,
  dist_sel = "gaussian",
  curve_col = "TrueCurve", est_col = "Est",
  lower_col = "Lower", upper_col = "Upper",
  ylab = NULL
) {
  if (is.null(ylab)) ylab <- lbl$y_fate
  dt_s <- dt_e[N == N_sel & T_pts == T_sel & ErrorDist == dist_sel]
  dt_s <- dt_s[
    dt_s[, .I[SimID == min(SimID)], by = .(Scenario, Irregularity)]$V1
  ]

  n_pts <- dt_s[, .N, by = .(Scenario_label, Irr_label)][, max(N)]
  every <- max(1L, floor(n_pts / 12))
  dt_pts <- dt_s[
    , .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, Irr_label)
  ]

  ggplot(dt_s, aes(x = t_val)) +
    geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
    geom_ribbon(
      aes(
        ymin = .data[[lower_col]], ymax = .data[[upper_col]],
        fill = Irregularity
      ),
      alpha = 0.2
    ) +
    geom_line(
      aes(y = .data[[curve_col]]),
      color = "black", linewidth = 0.55, linetype = "longdash"
    ) +
    geom_line(
      aes(y = .data[[est_col]], color = Irregularity, linetype = Irregularity),
      linewidth = 0.8
    ) +
    geom_point(
      data = dt_pts,
      aes(y = .data[[est_col]], color = Irregularity, shape = Irregularity),
      size = 2.2
    ) +
    facet_grid(Scenario_label ~ Irr_label, scales = "free_y") +
    scale_color_manual(values = IRR_COLORS, labels = lbl$irr) +
    scale_fill_manual(values = IRR_COLORS, labels = lbl$irr) +
    scale_linetype_manual(values = IRR_LINETYPES, labels = lbl$irr) +
    scale_shape_manual(values = IRR_SHAPES, labels = lbl$irr) +
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
  lbl <- LABELS[[lang]]

  dt_r <- prep(copy(dt_raw), lbl)
  dt_cs <- prep(copy(dt_csumm), lbl)
  dt_e <- prep(copy(dt_ex), lbl)

  # --- 5.1 NRMSE e cobertura por cenário ---
  cat(sprintf("\n[%s | 5.1] NRMSE e cobertura...\n", lang))
  for (scen in names(lbl$scenarios)) {
    dt_s <- dt_r[Scenario == scen]
    save_plot(plot_nrmse(dt_s, lbl), sprintf("nrmse_%s.png", scen), lang)
    save_plot(plot_coverage(dt_s, lbl), sprintf("coverage_%s.png", scen), lang)
    save_plot(
      plot_nrmse(dt_s, lbl, "CFATE_NRMSE"),
      sprintf("cfate_nrmse_%s.png", scen), lang
    )
    save_plot(
      plot_coverage(dt_s, lbl, "CFATE_Cov_Simul"),
      sprintf("cfate_coverage_%s.png", scen), lang
    )
  }

  # --- 5.2 Curvas ponto-a-ponto ---
  cat(
    sprintf(
      "\n[%s | 5.2] Ptwise NRMSE, viés, largura e cobertura...\n", lang
    )
  )
  for (case_name in names(CASES)) {
    N_sel <- CASES[[case_name]]$N
    T_sel <- CASES[[case_name]]$T
    save_plot(
      plot_pointwise(dt_cs, "NRMSE_t", lbl$y_ptwise_nrmse, lbl, N_sel, T_sel),
      sprintf("ptwise_nrmse_%s.png", case_name), lang,
      width = 9, height = 7
    )
    save_plot(
      plot_pointwise(dt_cs, "Bias_t", lbl$y_ptwise_bias, lbl, N_sel, T_sel),
      sprintf("ptwise_bias_%s.png", case_name), lang,
      width = 9, height = 7
    )
    save_plot(
      plot_pointwise(dt_cs, "Width_t", lbl$y_ptwise_width, lbl, N_sel, T_sel),
      sprintf("ptwise_width_%s.png", case_name), lang,
      width = 9, height = 7
    )
    save_plot(
      plot_pointwise(
        dt_cs, "Cov_t", lbl$y_ptwise_cov, lbl,
        N_sel, T_sel,
        add_hline = 0.95
      ),
      sprintf("ptwise_cov_%s.png", case_name), lang,
      width = 9, height = 7
    )
    save_plot(
      plot_pointwise(
        dt_cs, "CFATE_NRMSE_t", lbl$y_cfate_nrmse,
        lbl, N_sel, T_sel
      ),
      sprintf("ptwise_cfate_nrmse_%s.png", case_name),
      lang,
      width = 9, height = 7
    )
    save_plot(
      plot_pointwise(
        dt_cs, "CFATE_Cov_t", lbl$y_cfate_cov,
        lbl, N_sel, T_sel,
        add_hline = 0.95
      ),
      sprintf("ptwise_cfate_cov_%s.png", case_name), lang,
      width = 9, height = 7
    )
  }

  # --- 5.3 Exemplares: FATE ---
  cat(sprintf("\n[%s | 5.3] Exemplares (FATE)...\n", lang))
  for (case_name in names(CASES)) {
    for (dist_sel in c("gaussian", "t")) {
      save_plot(
        plot_exemplar(
          dt_e, lbl, CASES[[case_name]]$N,
          CASES[[case_name]]$T, dist_sel,
          ylab = lbl$y_fate
        ),
        sprintf("exemplar_%s_%s.png", case_name, dist_sel), lang,
        width = 10, height = 9
      )
    }
  }

  # --- 5.4 Exemplares: CFATE ---
  cat(sprintf("\n[%s | 5.4] Exemplares (CFATE)...\n", lang))
  for (case_name in names(CASES)) {
    for (dist_sel in c("gaussian", "t")) {
      save_plot(
        plot_exemplar(
          dt_e, lbl, CASES[[case_name]]$N,
          CASES[[case_name]]$T, dist_sel,
          curve_col = "TrueCFATE", est_col = "CFATE_Est",
          lower_col = "CFATE_Lower", upper_col = "CFATE_Upper",
          ylab = lbl$y_cfate
        ),
        sprintf("cfate_exemplar_%s_%s.png", case_name, dist_sel), lang,
        width = 10, height = 9
      )
    }
  }
}

cat("\nGeração de artefatos concluída.\n")
