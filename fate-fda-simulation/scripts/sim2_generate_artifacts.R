# scripts/sim2_generate_artifacts.R

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

input_dir <- "output/sim2_final"
dir.create(
  file.path(input_dir, "plots/pt"),
  recursive = TRUE, showWarnings = FALSE
)
dir.create(
  file.path(input_dir, "plots/eng"),
  recursive = TRUE, showWarnings = FALSE
)

SCENARIOS <- c("Bumps", "Doppler", "Heavisine", "Sigmoid")

# =============================================================================
# 2. SISTEMA VISUAL
# =============================================================================

LABELS <- list(
  pt = list(
    scenarios = c(
      Bumps = "Bumps", Doppler = "Doppler",
      Heavisine = "Escada", Sigmoid = "Sigmoide"
    ),
    dists = c(gaussian = "Gaussiano", t = "t (5 g.l.)"),
    x_T = "Pontos na grade (T)",
    x_t = "Tempo (t)",
    y_nrmse_main = "NRMSE (Mediana e Quartis)",
    y_nrmse_box = "NRMSE (Boxplot)",
    y_latency = "Latência média (µs)",
    method_on = "Ondaletas (Online)"
  ),
  eng = list(
    scenarios = c(
      Bumps = "Bumps", Doppler = "Doppler",
      Heavisine = "Staircase", Sigmoid = "Sigmoid"
    ),
    dists = c(gaussian = "Gaussian", t = "t (5 d.f.)"),
    x_T = "Grid points (T)",
    x_t = "Time (t)",
    y_nrmse_main = "NRMSE (Median and Quartiles)",
    y_nrmse_box = "NRMSE (Boxplot)",
    y_latency = "Mean latency (µs)",
    method_on = "Wavelets (Online)"
  )
)

base_theme <- theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(margin = margin(5, 5, 5, 5), face = "bold")
  )

# =============================================================================
# 3. DADOS
# =============================================================================

cat("Carregando dados...\n")
df_raw <- fread(file.path(input_dir, "sim2_raw.csv"))
df_latency <- fread(file.path(input_dir, "sim2_latency_curves.csv"))

df_raw[, N_label := factor(
  paste("N =", N),
  levels = paste("N =", sort(unique(N)))
)]
df_raw[, T_factor := factor(T_pts, levels = sort(unique(T_pts)))]

df_latency[
  ,
  N_label := factor(
    paste("N =", N),
    levels = paste("N =", sort(unique(N)))
  )
]
df_latency[, T_factor := factor(T_pts, levels = sort(unique(T_pts)))]
df_latency[
  ,
  ErrorDist_label := factor(
    ifelse(ErrorDist == "gaussian", "Gaussiano", "t (5 g.l.)"),
    levels = c("Gaussiano", "t (5 g.l.)")
  )
]

# =============================================================================
# 4. GERAÇÃO
# =============================================================================

for (lang in c("pt", "eng")) {
  lbl <- LABELS[[lang]]
  out_dir <- file.path(input_dir, "plots", lang)

  df_raw[
    ,
    ErrorDist_label := factor(
      lbl$dists[ErrorDist],
      levels = unname(lbl$dists)
    )
  ]
  df_latency[
    , ErrorDist_label := factor(
      lbl$dists[ErrorDist],
      levels = unname(lbl$dists)
    )
  ]
  df_raw[, Method := lbl$method_on]

  # ── 4.1 NRMSE ──────────────────────────────────────────────────────────────
  cat(sprintf("\n[%s | 4.1] NRMSE...\n", lang))
  for (scen in SCENARIOS) {
    p_data <- df_raw[Scenario == scen]
    if (nrow(p_data) == 0) next

    p_main <- ggplot(
      p_data,
      aes(x = T_factor, y = NRMSE_final, color = Method)
    ) +
      stat_summary(fun.data = "median_hilow", geom = "pointrange", size = 0.8) +
      facet_wrap(~N_label, ncol = 3) +
      scale_color_manual(values = setNames("royalblue", lbl$method_on)) +
      labs(x = lbl$x_T, y = lbl$y_nrmse_main) +
      base_theme
    ggsave(
      file.path(out_dir, sprintf("nrmse_%s.png", scen)),
      p_main,
      width = 9, height = 5, dpi = 300, bg = "white"
    )
    cat(sprintf("  -> %s/nrmse_%s.png\n", lang, scen))

    p_box <- ggplot(p_data, aes(x = T_factor, y = NRMSE_final, fill = Method)) +
      geom_boxplot(color = "black", alpha = 0.8) +
      facet_wrap(~N_label, ncol = 3) +
      scale_fill_manual(values = setNames("royalblue", lbl$method_on)) +
      labs(x = lbl$x_T, y = lbl$y_nrmse_box) +
      base_theme
    ggsave(
      file.path(out_dir, sprintf("boxplot_nrmse_%s.png", scen)),
      p_box,
      width = 9, height = 5, dpi = 300, bg = "white"
    )
    cat(sprintf("  -> %s/boxplot_nrmse_%s.png\n", lang, scen))
  }

  # ── 4.2 Curvas de latência ─────────────────────────────────────────────────
  cat(sprintf("\n[%s | 4.2] Curvas de latência...\n", lang))
  for (scen in SCENARIOS) {
    lat_data <- df_latency[Scenario == scen]
    if (nrow(lat_data) == 0) next

    p <- ggplot(
      lat_data,
      aes(x = t_val, y = Latency_mean_us, color = T_factor, group = T_factor)
    ) +
      geom_line(linewidth = 0.7) +
      facet_grid(ErrorDist_label ~ N_label, scales = "fixed") +
      labs(x = lbl$x_t, y = lbl$y_latency, color = "T") +
      base_theme
    ggsave(
      file.path(out_dir, sprintf("latency_%s.png", scen)),
      p,
      width = 9, height = 6, dpi = 300, bg = "white"
    )
    cat(sprintf("  -> %s/latency_%s.png\n", lang, scen))
  }
}

cat("\nGeração de artefatos da Simulação 2 concluída.\n")