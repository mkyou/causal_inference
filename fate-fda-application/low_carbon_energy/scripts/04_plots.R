# 04_plots.R
# Visualização dos resultados FATE/CFATE e análise de sensibilidade causal.
#
# Saídas:
#   output/lcl_final/plots/texto_principal/ — figuras principais
#   output/lcl_final/plots/apendice/ — gráficos complementares/diagnóstico
#
# Convenções visuais idênticas a 02_preprocess.R:
#   theme_bw(14), serif, plot.title = element_blank()
#   Modo offline: firebrick  |  Modo causal: royalblue
#   ACORN: paleta padrão ggplot2 (mesma ordem de 02)
#   dpi = 300, bg = "white"

library(data.table)
library(ggplot2)
library(here)

source(
  normalizePath(
    file.path(
      dirname(here()), "fate-fda-simulation", 
      "R", "estimate_fate_optimized.R"
    )
  )
)

Sys.setlocale("LC_TIME", "pt_BR.UTF-8")

DATA_DIR = here("low_carbon_energy", "data", "processed")
OUT_MAIN = here(
  "low_carbon_energy", "output", "lcl_final", 
  "plots", "texto_principal"
)
OUT_APP = here("low_carbon_energy", "output", "lcl_final", "plots", "apendice")

for (d in c(OUT_MAIN, OUT_APP))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# §1  Carregar estimativas e metadados
# =============================================================================

message("Carregando estimativas...")
fate_res = readRDS(file.path(DATA_DIR, "fate_estimates.rds"))
meta = fread(file.path(DATA_DIR, "X_meta.csv"))
setkey(meta, LCLid)
A_vec = meta$A
X_reg = model.matrix(~ Acorn_grouped, data = meta)
x_bar = colMeans(X_reg)
N = nrow(meta)
T_pts = fate_res[mode == "offline", .N]
ref_grid = fate_res[mode == "offline", tstp]
fate_offline = fate_res[mode == "offline", fate]
fate_causal = fate_res[mode == "causal",  fate]

# Tema base — idêntico ao 02_preprocess.R
base_theme = theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(margin = margin(5, 5, 5, 5), face = "bold")
  )

mode_colors = c("Offline (Haar)" = "firebrick", "Causal (Haar)" = "royalblue")

acorn_colors = c(
  "ACORN-U"     = "firebrick",
  "Adversity"   = "salmon",
  "Affluent"    = "forestgreen",
  "Comfortable" = "royalblue"
)
trim_labels = c(
  "1" = "Q1 (jan–mar)", "2" = "Q2 (abr–jun)",
  "3" = "Q3 (jul–set)", "4" = "Q4 (out–dez)"
)

fate_res[
  ,
  mode_label := factor(
    mode,
    levels = c("offline", "causal"),
    labels = c("Offline (Haar)", "Causal (Haar)")
  )
]

# =============================================================================
# §2  FATE por grupo ACORN — lido do arquivo gerado em script 03
# =============================================================================

message("Carregando fate_by_acorn...")
acorn_res      = readRDS(file.path(DATA_DIR, "fate_by_acorn.rds"))
fate_acorn_off = acorn_res$offline
fate_acorn_cau = acorn_res$causal
acorn_levels   = names(fate_acorn_off)

cumtrapz = function(y) c(0, cumsum((head(y, -1L) + y[-1L]) / 2))

# =============================================================================
# §3  FATE bruto vs estimado
#     intraday → texto_principal   |   anual → apendice
# =============================================================================

message("Carregando raw_wide para FATE bruto...")
raw_wide = readRDS(file.path(DATA_DIR, "raw_wide.rds"))
fate_raw = colMeans(raw_wide[A_vec == 1L, ], na.rm = TRUE) - colMeans(raw_wide[A_vec == 0L, ], na.rm = TRUE)
rm(raw_wide); gc()

dt_compare = rbindlist(
  list(
    data.table(tstp = ref_grid, fate = fate_raw,     serie = "Diferença bruta"),
    data.table(tstp = ref_grid, fate = fate_offline, serie = "Offline (Haar)"),
    data.table(tstp = ref_grid, fate = fate_causal,  serie = "Causal (Haar)")
  )
)
dt_compare[, hora := as.numeric(format(tstp, "%H")) + as.numeric(format(tstp, "%M")) / 60]
dt_compare[, date := as.Date(tstp)]
dt_compare[
  , 
  serie := factor(
    serie,
    levels = c("Diferença bruta", "Offline (Haar)", "Causal (Haar)")
  )
]

compare_colors = c(
  "Diferença bruta" = "black",
  "Offline (Haar)"  = "firebrick",
  "Causal (Haar)"   = "royalblue"
)

# 3a. Intraday → texto_principal
intraday_comp = dt_compare[, .(fate = mean(fate)), by = .(hora, serie)]

p3a = ggplot(intraday_comp, aes(x = hora, y = fate, colour = serie)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = compare_colors) +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  labs(
    x = "Hora do dia",
    y = expression(hat(theta)(t) ~ "(kWh/hh)")
  ) +
  base_theme

ggsave(
  file.path(OUT_MAIN, "fate_intraday_comparison.png"), p3a,
  width = 10, height = 5, dpi = 300, bg = "white"
)

# 3b. Anual (LOESS) → apendice
annual_comp = dt_compare[, .(fate = mean(fate)), by = .(date, serie)]

p3b = ggplot(annual_comp, aes(x = date, y = fate, colour = serie)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_line(
    data = annual_comp[serie != "Diferença bruta"],
    linewidth = 0.8
  ) +
  geom_line(
    data = annual_comp[serie == "Diferença bruta"],
    linewidth = 0.3, alpha = 0.3
  ) +
  geom_smooth(
    data = annual_comp[serie == "Diferença bruta"],
    method = "loess", formula = y ~ x, span = 0.15,
    se = FALSE, linewidth = 0.9
  ) +
  scale_colour_manual(values = compare_colors) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Data (2013)",
    y = expression(hat(theta)(t) ~ "(kWh/hh)")
  ) +
  base_theme

ggsave(
  file.path(OUT_APP, "fate_annual_comparison.png"), p3b,
  width = 10, height = 5, dpi = 300, bg = "white"
)

message("§3 feito.")

# =============================================================================
# §4  FATE por timestamp com bandas simultâneas bootstrap 95% → texto_principal
# =============================================================================

p4 = ggplot(fate_res, aes(x = tstp)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = mode_label), alpha = 0.2) +
  geom_line(aes(y = fate, colour = mode_label), linewidth = 0.45) +
  facet_wrap(~ mode_label, ncol = 1L) +
  scale_colour_manual(values = mode_colors) +
  scale_fill_manual(values = mode_colors) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Data (2013)",
    y = expression(hat(theta)(t) ~ "(kWh/hh)")
  ) +
  base_theme + theme(legend.position = "none")

ggsave(
  file.path(OUT_MAIN, "fate_bands.png"), p4,
  width = 10, height = 7, dpi = 300, bg = "white"
)

message("§4 feito.")

# =============================================================================
# §5  FATE intradiário por grupo ACORN → texto_principal
# =============================================================================

dt_acorn = rbindlist(
  lapply(
    acorn_levels,
    function(grp) {
      rbindlist(
        list(
          data.table(
            tstp = ref_grid, fate = fate_acorn_off[[grp]],
            acorn = grp, mode_label = "Offline (Haar)"
          ),
          data.table(
            tstp = ref_grid, fate = fate_acorn_cau[[grp]],
            acorn = grp, mode_label = "Causal (Haar)"
          )
        )
      )
    }
  )
)

dt_acorn[
  ,
  hora := as.numeric(format(tstp, "%H")) + as.numeric(format(tstp, "%M")) / 60
]

dt_acorn[
  ,
  mode_label := factor(
    mode_label,
    levels = c("Offline (Haar)", "Causal (Haar)")
  )
]

intraday_acorn = dt_acorn[
  ,
  .(fate = mean(fate)),
  by = .(hora, acorn, mode_label)
]

p5 = ggplot(
  intraday_acorn, 
  aes(x = hora, y = fate, colour = acorn)
) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ mode_label, ncol = 2L) +
  scale_colour_manual(values = acorn_colors) +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  labs(
    x = "Hora do dia",
    y = expression(hat(theta)(t) ~ "(kWh/hh)"),
    colour = "ACORN"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, family = "serif")
  )

ggsave(
  file.path(OUT_MAIN, "fate_by_acorn_intraday.png"), p5,
  width = 10, height = 5, dpi = 300, bg = "white"
)

message("§5 feito.")

# =============================================================================
# §6  CFATE ao longo do ano com bandas → texto_principal
# =============================================================================

p6 = ggplot(fate_res, aes(x = tstp)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_ribbon(
    aes(ymin = cfate_low, ymax = cfate_up, fill = mode_label),
    alpha = 0.2
  ) +
  geom_line(aes(y = cfate, colour = mode_label), linewidth = 0.7) +
  facet_wrap(~ mode_label, ncol = 1L) +
  scale_colour_manual(values = mode_colors) +
  scale_fill_manual(values = mode_colors) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Data (2013)",
    y = expression(hat(Theta)(t) ~ "(kWh acumulado/hh)")
  ) +
  base_theme + theme(legend.position = "none")

ggsave(
  file.path(OUT_MAIN, "cfate_annual.png"), p6,
  width = 10, height = 7, dpi = 300, bg = "white"
)

message("§6 feito.")

# =============================================================================
# §7  CFATE acumulado por grupo ACORN → texto_principal
# =============================================================================

dt_acorn_cfate = rbindlist(
  lapply(
    acorn_levels, 
    function(grp) {
      rbindlist(
        list(
          data.table(
            tstp  = ref_grid,
            cfate = cumtrapz(fate_acorn_off[[grp]]),
            acorn = grp, mode_label = "Offline (Haar)"
          ),
          data.table(
            tstp  = ref_grid,
            cfate = cumtrapz(fate_acorn_cau[[grp]]),
            acorn = grp, mode_label = "Causal (Haar)"
          )
        )
      )
    }
  )
)

dt_acorn_cfate[
  ,
  mode_label := factor(
    mode_label,
    levels = c("Offline (Haar)", "Causal (Haar)")
  )
]

p7 = ggplot(dt_acorn_cfate, aes(x = tstp, y = cfate, colour = acorn)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ mode_label, ncol = 2L) +
  scale_colour_manual(values = acorn_colors) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Data (2013)",
    y = expression(hat(Theta)(t) ~ "(kWh acumulado/hh)"),
    colour = "ACORN"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, family = "serif")
  )

ggsave(
  file.path(OUT_MAIN, "cfate_by_acorn.png"), p7,
  width = 10, height = 5, dpi = 300, bg = "white"
)

message("§7 feito.")

# =============================================================================
# §8  Tabela CFATE(T) por ACORN
# =============================================================================

cfate_table = rbindlist(
  lapply(
    acorn_levels, 
    function(grp) {
      data.table(
        acorn = grp,
        n = meta[Acorn_grouped == grp, .N],
        cfate_off = round(cumtrapz(fate_acorn_off[[grp]])[T_pts], 2),
        cfate_cau = round(cumtrapz(fate_acorn_cau[[grp]])[T_pts], 2)
      )
    }
  )
)

cfate_table = rbindlist(
  list(
    cfate_table, 
    data.table(
      acorn = "Global",
      n = N,
      cfate_off = round(fate_res[mode == "offline", cfate[T_pts]], 2),
      cfate_cau = round(fate_res[mode == "causal",  cfate[T_pts]], 2)
    )
  )
)

setnames(
  cfate_table,
  c("Grupo ACORN", "N", "CFATE offline (kWh/hh)", "CFATE causal (kWh/hh)")
)
fwrite(cfate_table, file.path(OUT_MAIN, "cfate_table_by_acorn.csv"))
print(cfate_table)
message("§8 feito.")

# =============================================================================
# §9  Sensibilidade causal
#
# 9a. Estabilidade sazonal: FATE intradiário por trimestre → apendice
# 9b. Teste de placebo por permutação (B = 200) → apendice
# =============================================================================

message("\n--- §9 Sensibilidade causal ---")

# 9a. Perfil intradiário por trimestre (Q1–Q4)
dt_seas = data.table(
  tstp = ref_grid,
  fate_off = fate_offline,
  fate_cau = fate_causal
)
dt_seas[, trim  := as.integer(ceiling(as.integer(format(tstp, "%m")) / 3))]
dt_seas[
  ,
  hora  := as.numeric(format(tstp, "%H")) + as.numeric(format(tstp, "%M")) / 60
]

seas_long = melt(
  dt_seas[
    , 
    .(fate_off = mean(fate_off), fate_cau = mean(fate_cau)),
    by = .(hora, trim)
  ],
  id.vars = c("hora", "trim"),
  measure.vars = c("fate_off", "fate_cau"),
  variable.name = "mode", value.name = "fate"
)
seas_long[
  , 
  mode := factor(
    mode,
    levels = c("fate_off", "fate_cau"),
    labels = c("Offline (Haar)", "Causal (Haar)")
  )
]
seas_long[
  ,
  trim_label := factor(
    trim_labels[as.character(trim)],
    levels = trim_labels
  )
]

p9a = ggplot(seas_long, aes(x = hora, y = fate, colour = mode)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", colour = "gray50",
    linewidth = 0.4
  ) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ trim_label, ncol = 2L) +
  scale_colour_manual(values = mode_colors) +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  labs(
    x = "Hora do dia",
    y = expression(hat(theta)(t) ~ "(kWh/hh)")
  ) +
  base_theme

ggsave(
  file.path(OUT_APP, "sensitivity_seasonal.png"), p9a,
  width = 10, height = 7, dpi = 300, bg = "white"
)
message("  §9a feito.")

# 9b. Teste de placebo por permutação
perm_cache = file.path(DATA_DIR, "perm_stats.rds")

if (file.exists(perm_cache)) {
  message("  §9b: carregando cache.")
  perm_stats = readRDS(perm_cache)
} else {
  N_PERM = 200L
  message(sprintf("  §9b: %d permutações...", N_PERM))
  smooth_cau_mat = readRDS(file.path(DATA_DIR, "smooth_causal.rds"))
  set.seed(2024L)
  perm_stats = vapply(seq_len(N_PERM), function(b) {
    proj_p = precompute_projection(X_reg, sample(A_vec))
    max(abs(estimate_fate_fast(smooth_cau_mat, proj_p)$fate_est))
  }, numeric(1L))
  rm(smooth_cau_mat); gc()
  saveRDS(perm_stats, perm_cache, compress = FALSE)
}

obs_stat = max(abs(fate_causal))
p_val = mean(perm_stats >= obs_stat)
message(sprintf("  Observado: %.4f | p = %.4f", obs_stat, p_val))

p9b = ggplot(data.table(stat = perm_stats), aes(x = stat)) +
  geom_histogram(bins = 30L, fill = "gray70", colour = "white") +
  geom_vline(xintercept = obs_stat, colour = "royalblue", linewidth = 1.2) +
  annotate(
    "text", x = obs_stat * 0.985, y = Inf,
    vjust = 1.5, hjust = 1, family = "serif", size = 4,
    label = sprintf("Observado = %.3f\np = %.3f", obs_stat, p_val),
    colour = "royalblue"
  ) +
  labs(
    x = expression(max[t] ~ "|" ~ hat(theta)(t) ~ "|"),
    y = "Frequência"
  ) +
  base_theme

ggsave(
  file.path(OUT_APP, "sensitivity_placebo.png"), p9b,
  width = 8, height = 5, dpi = 300, bg = "white"
)
message("  §9b feito.")

message("\nScript 04 concluído.")
message("  texto_principal: ", OUT_MAIN)
message("  apendice:        ", OUT_APP)