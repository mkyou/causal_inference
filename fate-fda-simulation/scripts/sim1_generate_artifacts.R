# scripts/sim1_generate_artifacts.R

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

out_dir <- "output/sim1_final/plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. SISTEMA VISUAL
# =============================================================================

METHOD_COLORS <- c(
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

# Linha: diferencia métodos no gráfico de comparação
LT_MAIN <- c(
  wavelet_causal = "solid",
  wavelet_offline = "dashed",
  spline_ecker = "dotted"
)

# Linha: diferencia níveis de confundimento dentro de uma família
LT_CONF <- c(
  wavelet_causal = "solid",
  wavelet_causal_conf_low = "dashed",
  wavelet_causal_conf_med = "dotdash",
  wavelet_causal_conf_high = "dotted",
  wavelet_offline = "solid",
  wavelet_offline_conf_low = "dashed",
  wavelet_offline_conf_med = "dotdash",
  wavelet_offline_conf_high = "dotted"
)

# Forma: causal=triângulo, offline=X, spline=círculo
# Confundimento: sem conf=triângulo, baixo=círculo, médio=X, alto=quadrado
SH_MAIN <- c(
  wavelet_causal = 17,
  wavelet_offline = 4,
  spline_ecker = 19
)

SH_CONF <- c(
  wavelet_causal = 17,
  wavelet_causal_conf_low = 19,
  wavelet_causal_conf_med = 4,
  wavelet_causal_conf_high = 15,
  wavelet_offline = 17,
  wavelet_offline_conf_low = 19,
  wavelet_offline_conf_med = 4,
  wavelet_offline_conf_high = 15
)

METHOD_LABELS <- c(
  wavelet_causal = "Ond. causal",
  wavelet_causal_conf_low = "Ond. causal (conf. baixo)",
  wavelet_causal_conf_med = "Ond. causal (conf. médio)",
  wavelet_causal_conf_high = "Ond. causal (conf. alto)",
  wavelet_offline = "Ond. offline",
  wavelet_offline_conf_low = "Ond. offline (conf. baixo)",
  wavelet_offline_conf_med = "Ond. offline (conf. médio)",
  wavelet_offline_conf_high = "Ond. offline (conf. alto)",
  spline_ecker = "Spline (Ecker)"
)

SCENARIO_LABELS <- c(
  Bumps = "Bumps", Doppler = "Doppler",
  Heavisine = "Escada", Sigmoid = "Sigmoide"
)

DIST_LABELS <- c(gaussian = "Gaussiano", t = "t (5 g.l.)")

MAIN_METHODS <- c("wavelet_causal", "wavelet_offline", "spline_ecker")

OFFLINE_CONF <- c(
  "wavelet_offline", "wavelet_offline_conf_low",
  "wavelet_offline_conf_med", "wavelet_offline_conf_high"
)

CAUSAL_CONF <- c(
  "wavelet_causal", "wavelet_causal_conf_low",
  "wavelet_causal_conf_med", "wavelet_causal_conf_high"
)

CASES <- list(
  worst = list(N = 200, T = 128),
  med = list(N = 500, T = 256),
  best = list(N = 2000, T = 512)
)

T_TREAT <- 50.0

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
dt_raw <- fread("output/sim1_final/sim1_raw.csv")
dt_csumm <- fread("output/sim1_final/sim1_curve_summary.csv")
dt_ex <- fread("output/sim1_final/sim1_exemplars.csv")

# NRMSE_t = RMSE_t / range(TrueCurve): range é determinístico,
# extrai dos exemplares
true_ranges <- dt_ex[
  Method == "wavelet_offline",
  .(
    true_range = diff(range(TrueCurve))
  ),
  by = .(Scenario, N, T_pts, ErrorDist)
]

dt_csumm <- merge(
  dt_csumm, true_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist"), all.x = TRUE
)
dt_csumm[, NRMSE_t := RMSE_t / true_range]

prep <- function(dt) {
  dt[, Method_label := factor(
    METHOD_LABELS[Method],
    levels = unname(METHOD_LABELS)
  )]
  dt[, Scenario_label := factor(
    SCENARIO_LABELS[Scenario],
    levels = unname(SCENARIO_LABELS)
  )]
  dt[, ErrorDist_label := factor(
    DIST_LABELS[ErrorDist],
    levels = unname(DIST_LABELS)
  )]
  dt[, N_label := factor(
    paste0("N = ", N),
    levels = paste0("N = ", c(200, 500, 2000))
  )]
  dt[, T_fac := factor(T_pts, levels = c(128, 256, 512))]
  dt
}
dt_raw <- prep(dt_raw)
dt_csumm <- prep(dt_csumm)
dt_ex <- prep(dt_ex)

# =============================================================================
# 4. FUNÇÕES
# =============================================================================

save_plot <- function(p, fname, width = 9, height = 6) {
  path <- file.path(out_dir, fname)
  ggsave(path, p, width = width, height = height, dpi = 300, bg = "white")
  cat(sprintf("  -> %s\n", fname))
}

# NRMSE: mediana + quartis, faceta ErrorDist × N
plot_nrmse <- function(dt, methods, lt, sh) {
  dt_s <- dt[Method %in% methods]
  dt_s[, Method := factor(Method, levels = methods)]
  ggplot(dt_s, aes(
    x = T_fac, y = NRMSE,
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
      values = METHOD_COLORS[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_linetype_manual(
      values = lt[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_shape_manual(
      values = sh[methods],
      labels = METHOD_LABELS[methods]
    ) +
    labs(x = "Pontos na grade (T)", y = "NRMSE (mediana ± quartis)") +
    base_theme
}

# Cobertura simultânea: linhas, faceta ErrorDist × N
plot_coverage <- function(dt, methods, lt, sh) {
  dt_cov <- dt[Method %in% methods, .(
    Cov = mean(Cov_Simul, na.rm = TRUE)
  ), by = .(ErrorDist_label, N_label, T_fac, Method)]
  dt_cov[, Method := factor(Method, levels = methods)]
  ggplot(dt_cov, aes(
    x = T_fac, y = Cov,
    color = Method, linetype = Method, shape = Method, group = Method
  )) +
    geom_hline(
      yintercept = 0.95, color = "gray40",
      linewidth = 0.4, linetype = "dashed"
    ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
    scale_color_manual(
      values = METHOD_COLORS[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_linetype_manual(
      values = lt[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_shape_manual(
      values = sh[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_y_continuous(
      limits = c(0, 1.01),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(
      x = "Pontos na grade (T)",
      y = "Cobertura simultânea empírica"
    ) +
    base_theme
}

# Curvas ponto-a-ponto (Bias_t ou Width_t por t)
# Faceta: cenário (linhas) × distribuição (colunas), escala livre por linha
plot_pointwise <- function(dt_cs, metric, ylab, methods, lt, sh,
                           N_sel, T_sel, add_hline = NULL) {
  dt_s <- dt_cs[Method %in% methods & N == N_sel & T_pts == T_sel]
  dt_s[, Method := factor(Method, levels = methods)]

  n_pts <- dt_s[, .N, by = .(Scenario_label, ErrorDist_label, Method)][, max(N)]
  every <- max(1L, floor(n_pts / 12))
  dt_pts <- dt_s[, .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, ErrorDist_label, Method)
  ]

  p <- ggplot(dt_s, aes(
    x = t_val, y = .data[[metric]],
    color = Method, linetype = Method
  )) +
    geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
    geom_line(linewidth = 0.7) +
    geom_point(data = dt_pts, aes(shape = Method), size = 2.2) +
    facet_grid(Scenario_label ~ ErrorDist_label, scales = "free_y") +
    scale_color_manual(
      values = METHOD_COLORS[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_linetype_manual(
      values = lt[methods],
      labels = METHOD_LABELS[methods]
    ) +
    scale_shape_manual(
      values = sh[methods],
      labels = METHOD_LABELS[methods]
    ) +
    labs(
      x = "Tempo (t)", y = ylab,
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

# Exemplar: curva verdadeira + estimativa + banda (uma run por grupo)
# Faceta: cenário (linhas) × método (colunas)
plot_exemplar <- function(dt_e, methods, lt, sh, N_sel, T_sel,
                          dist_sel = "gaussian") {
  dt_s <- dt_e[
    Method %in% methods & N == N_sel &
      T_pts == T_sel & ErrorDist == dist_sel
  ]
  dt_s <- dt_s[
    dt_s[, .I[SimID == min(SimID)], by = .(Scenario, Method)]$V1
  ]
  dt_s[, Method := factor(Method, levels = methods)]

  n_pts <- dt_s[, .N, by = .(Scenario_label, Method)][, max(N)]
  every <- max(1L, floor(n_pts / 12))
  dt_pts <- dt_s[, .SD[seq(1, .N, by = every)],
    by = .(Scenario_label, Method)
  ]

  ggplot(dt_s, aes(x = t_val)) +
    geom_vline(xintercept = T_TREAT, linewidth = 0.4, color = "black") +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper, fill = Method),
      alpha = 0.2
    ) +
    geom_line(
      aes(y = TrueCurve),
      color = "black", linewidth = 0.55, linetype = "longdash"
    ) +
    geom_line(
      aes(y = Est, color = Method, linetype = Method),
      linewidth = 0.8
    ) +
    geom_point(
      data = dt_pts,
      aes(y = Est, color = Method, shape = Method), size = 2.2
    ) +
    facet_grid(Scenario_label ~ Method_label, scales = "free_y") +
    scale_color_manual(values = METHOD_COLORS[methods]) +
    scale_fill_manual(values = METHOD_COLORS[methods]) +
    scale_linetype_manual(values = lt[methods]) +
    scale_shape_manual(values = sh[methods]) +
    labs(
      x = "Tempo (t)", y = "FATE(t)",
      subtitle = sprintf(
        "%s  |  N = %d,  T = %d",
        DIST_LABELS[dist_sel], N_sel, T_sel
      )
    ) +
    base_theme +
    theme(legend.position = "none")
}

# =============================================================================
# 5. GERAÇÃO
# =============================================================================

# --- 5.1 Comparação principal: NRMSE e cobertura ---
cat("\n[5.1] NRMSE e cobertura (3 métodos)...\n")
for (scen in names(SCENARIO_LABELS)) {
  dt_s <- dt_raw[Scenario == scen]
  save_plot(
    plot_nrmse(dt_s, MAIN_METHODS, LT_MAIN, SH_MAIN),
    sprintf("nrmse_%s.png", scen)
  )
  save_plot(
    plot_coverage(dt_s, MAIN_METHODS, LT_MAIN, SH_MAIN),
    sprintf("coverage_%s.png", scen)
  )
}

# --- 5.2 Curvas ponto-a-ponto: NRMSE, viés e largura (3 métodos) ---
cat("\n[5.2] Ptwise NRMSE, viés e largura (3 métodos)...\n")
for (case_name in names(CASES)) {
  N_sel <- CASES[[case_name]]$N
  T_sel <- CASES[[case_name]]$T
  save_plot(
    plot_pointwise(
      dt_csumm, "NRMSE_t", "NRMSE pontual",
      MAIN_METHODS, LT_MAIN, SH_MAIN, N_sel, T_sel
    ),
    sprintf("ptwise_nrmse_%s.png", case_name),
    width = 9, height = 7
  )
  save_plot(
    plot_pointwise(
      dt_csumm, "Bias_t", "Viés pontual",
      MAIN_METHODS, LT_MAIN, SH_MAIN, N_sel, T_sel
    ),
    sprintf("ptwise_bias_%s.png", case_name),
    width = 9, height = 7
  )
  save_plot(
    plot_pointwise(
      dt_csumm, "Width_t", "Largura pontual da banda",
      MAIN_METHODS, LT_MAIN, SH_MAIN, N_sel, T_sel
    ),
    sprintf("ptwise_width_%s.png", case_name),
    width = 9, height = 7
  )
}

# --- 5.3 Exemplares: 3 métodos principais ---
cat("\n[5.3] Exemplares (3 métodos)...\n")
for (case_name in names(CASES)) {
  for (dist_sel in c("gaussian", "t")) {
    save_plot(
      plot_exemplar(
        dt_ex, MAIN_METHODS, LT_MAIN, SH_MAIN,
        CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel
      ),
      sprintf("exemplar_%s_%s.png", case_name, dist_sel),
      width = 10, height = 9
    )
  }
}

# --- 5.4 Exemplares: confundimento offline e causal ---
cat("\n[5.4] Exemplares (confundimento)...\n")
for (case_name in names(CASES)) {
  for (dist_sel in c("gaussian", "t")) {
    save_plot(
      plot_exemplar(
        dt_ex, OFFLINE_CONF, LT_CONF, SH_CONF,
        CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel
      ),
      sprintf("conf_off_exemplar_%s_%s.png", case_name, dist_sel),
      width = 12, height = 9
    )
    save_plot(
      plot_exemplar(
        dt_ex, CAUSAL_CONF, LT_CONF, SH_CONF,
        CASES[[case_name]]$N, CASES[[case_name]]$T, dist_sel
      ),
      sprintf("conf_cau_exemplar_%s_%s.png", case_name, dist_sel),
      width = 12, height = 9
    )
  }
}

cat("\nGeração de artefatos concluída.\n")
