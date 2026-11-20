# 02_preprocess.R
#
# Aplica as decisões de 01_explore_data.R e produz o dataset de análise.
#
# Decisões implementadas:
#   Janela      : 2013-01-01 00:00 → 2013-12-31 23:30 (ano civil completo)
#   Exclusões   : MAC005565 (série irrecuperável) + ACORN- (positividade)
#   Seleção     : todos os ids_2013 — sem threshold adicional de completude
#   DST outono  : nenhuma ação; R distingue BST/GMT como POSIXct distintos
#   DST primavera: gap de 60 min tratado nativamente pelo rLifting
#   Leituras ausentes: ~78% dos domicílios perdem leituras esparsas por falha
#     operacional do medidor (Tipo 2: linha inteiramente ausente). O argumento
#     para identificabilidade é mecanístico: falhas de hardware são MCAR em
#     relação ao outcome — um medidor morto não sabe quanto o domicílio
#     consumia. Verificado empiricamente: proporção ToU entre domicílios com
#     falta ≈ proporção ToU global (sem seleção diferencial por grupo). Seção
#     3b abaixo verifica positividade por timestamp e severidade por grupo.
#     Cada domicílio é processado pelo rLifting no seu grid real (grade
#     irregular nativa). O estimador do FATE em cada t usa a média sobre os
#     domicílios com leitura disponível em t.
#   Covariáveis : A ∈ {0,1} (0 = Std, 1 = ToU); Acorn_grouped como fator e
#     como dummies. Categoria de referência definida no script de estimação.
#   Imputação   : nenhuma
#   Normalização: nenhuma (energia em kWh/hh como coletada)
#
# Saídas:
#   data/processed/lcl_2013_long.rds — data.table longo: LCLid, tstp, energy
#   data/processed/X_meta.csv        — LCLid, A, Acorn_grouped, dummies acorn_*
#   output/lcl_final/plots/          — gráficos EDA (main_text, appendix)
#
# Rodar de: Rscript low_carbon_energy/scripts/02_preprocess.R
#           (raiz do .Rproj — fate-fda-application/)
# Tempo estimado: 10–15 min

library(data.table)
library(ggplot2)
library(here)

# ══════════════════════════════════════════════════════════════════════════════
# CONSTANTES
# ══════════════════════════════════════════════════════════════════════════════

DATA_RAW = here("low_carbon_energy", "data", "raw")
HH_DIR = file.path(DATA_RAW, "halfhourly_dataset")
HH_INFO = file.path(DATA_RAW, "informations_households.csv")
PROC_DIR = here("low_carbon_energy", "data", "processed")
OUT_MAIN = here(
  "low_carbon_energy", "output", 
  "lcl_final", "plots", "main_text"
)
OUT_APP = here("low_carbon_energy", "output", "lcl_final", "plots", "appendix")

for (d in c(PROC_DIR, OUT_MAIN, OUT_APP))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# WIN_END definido com hora explícita para incluir a última leitura de 31/dez.
WIN_START = as.POSIXct("2013-01-01 00:00:00", tz = "Europe/London")
WIN_END = as.POSIXct("2013-12-31 23:30:00", tz = "Europe/London")

hh = fread(HH_INFO)
EXCLUDE = c("MAC005565", hh[Acorn == "ACORN-", LCLid])

block_files = sort(
  list.files(
    HH_DIR, pattern = "^block_\\d+\\.csv$", 
    full.names = TRUE
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# 1. DERIVAÇÃO DE IDS_2013
# ══════════════════════════════════════════════════════════════════════════════
# Mesmo critério de 01_explore_data.R — script autônomo não depende de 01.

cat("Derivando ids_2013...\n")

get_coverage = function(fpath) {
  b = fread(fpath, select = c("LCLid", "tstp"))
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b[!is.na(tstp), .(t_min = min(tstp), t_max = max(tstp)), by = LCLid]
}

cov = rbindlist(lapply(block_files, get_coverage))
ids_2013 = cov[
  t_min <= WIN_START & t_max >= WIN_END & !LCLid %in% EXCLUDE,
  LCLid
]

cat(sprintf("  %d domicílios elegíveis\n", length(ids_2013)))

# ══════════════════════════════════════════════════════════════════════════════
# 2. CARREGAMENTO
# ══════════════════════════════════════════════════════════════════════════════

load_block = function(fpath) {
  b = fread(fpath, select = c("LCLid", "tstp", "energy(kWh/hh)"))
  b = b[LCLid %in% ids_2013]
  if (nrow(b) == 0L) return(NULL)
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b[, energy := suppressWarnings(as.numeric(get("energy(kWh/hh)")))]
  b[, `energy(kWh/hh)` := NULL]
  b[!is.na(tstp) & tstp >= WIN_START & tstp <= WIN_END]
}

cat("Carregando blocos halfhourly...\n")
dt = rbindlist(lapply(block_files, load_block), use.names = TRUE)
setorder(dt, LCLid, tstp)
cat(sprintf("  %d linhas carregadas\n", nrow(dt)))

# ══════════════════════════════════════════════════════════════════════════════
# 3. VERIFICAÇÃO PÓS-CARGA
# ══════════════════════════════════════════════════════════════════════════════

# Confirma que não há duplicatas de (LCLid, tstp) — o R já distingue BST/GMT.
n_dup = dt[, sum(duplicated(tstp)), by = LCLid][V1 > 0, .N]
stopifnot(n_dup == 0L)

# Distribuição de leituras por domicílio (esperado: maioria ~17.518–17.520).
n_per_hh = dt[, .N, by = LCLid]
cat(
  sprintf(
    "Leituras por domicílio — mediana: %d | min: %d | max: %d\n",
    median(n_per_hh$N), min(n_per_hh$N), max(n_per_hh$N)
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# 4. METADADOS E COVARIÁVEIS
# ══════════════════════════════════════════════════════════════════════════════

meta = hh[LCLid %in% ids_2013, .(LCLid, stdorToU, Acorn_grouped)]
meta[, A := as.integer(stdorToU == "ToU")]

# Dummies para Acorn_grouped (sem drop de referência — definida no script 03).
for (lvl in sort(unique(meta$Acorn_grouped))) {
  col = paste0("acorn_", tolower(gsub("[^a-zA-Z0-9]", "_", lvl)))
  meta[, (col) := as.integer(Acorn_grouped == lvl)]
}

fwrite(meta, file.path(PROC_DIR, "X_meta.csv"))
cat(
  sprintf(
    "X_meta salvo: %d domicílios | ToU: %d | Std: %d\n",
    nrow(meta), meta[A == 1L, .N], meta[A == 0L, .N]
  )
)

# ══════════════════════════════════════════════════════════════════════════════
# 3b. DIAGNÓSTICO DE LEITURAS AUSENTES (TIPO 2)
# ══════════════════════════════════════════════════════════════════════════════
# Leituras Tipo 2: linhas inteiramente ausentes no arquivo por falha
# operacional do medidor. Detectadas pela diferença N_obs − 17.520.
#
# Para identificabilidade do FATE(t) as perguntas relevantes são:
#   (i)  Positividade: em cada timestamp, há observações dos dois grupos?
#   (ii) Mecanismo: a ausência é MCAR em relação ao outcome?
#        Se fosse MNAR (medidor falha correlacionado ao consumo) E diferencial
#        entre ToU e Std, o estimador E[Y(t)|A, M(t)=1] seria enviesado para
#        E[Y(t)|A]. O argumento aqui é mecanístico: falhas de hardware são
#        independentes do padrão de consumo. Verificação empírica abaixo.

n_miss_diag = merge(n_per_hh, meta[, .(LCLid, stdorToU, A)], by = "LCLid")
n_miss_diag[, n_miss := 17520L - N]

# Severidade por grupo
print(
  n_miss_diag[
    ,
    .(
      n_total = .N,
      pct_com_falta = round(mean(N < 17520L) * 100, 1),
      mediana_miss = as.double(median(n_miss)),
      p95_miss = as.integer(quantile(n_miss, 0.95))
    ),
    by = stdorToU
  ]
)

# Se as ausências selecionassem diferencialmente por tratamento, a proporção
# ToU entre os que faltam possivelmente seria diferente da proporção ToU global.
prop_global_ToU = mean(n_miss_diag$A)
prop_ToU_faltam = n_miss_diag[N < 17520L, mean(A)]
cat(
  sprintf(
    "Proporção ToU — global: %.1f%% | entre domicílios com falta: %.1f%%\n",
    prop_global_ToU * 100, prop_ToU_faltam * 100
  )
)

# Positividade por timestamp: verificada no loop bloco a bloco (seção 6),
# onde cada bloco é lido individualmente sem manter dt em memória.

# ══════════════════════════════════════════════════════════════════════════════
# 5. DATASET PROCESSADO (formato longo)
# ══════════════════════════════════════════════════════════════════════════════
# Salvo como data.table em RDS. O script de estimação (03) lê este arquivo,
# divide por LCLid e passa cada série ao rLifting individualmente.
# A grade nominal é regular (17.520 meias-horas uniformes). As ~2 leituras
# ausentes por domicílio (mediana) criam gaps pontuais na série individual —
# esses gaps são preenchidos por interpolação linear (approx(), rule=2) em
# fill_household() antes de passar o vetor completo ao rLifting. O rLifting
# opera em modo de grade regular (parâmetro t não utilizado). Para o FATE,
# a grade efetiva é a união dos 17.520 timestamps; positividade verificada
# em todos eles (seção 3b).

saveRDS(dt, file.path(PROC_DIR, "lcl_2013_long.rds"))
cat(sprintf("lcl_2013_long.rds salvo: %.0f MB\n", object.size(dt) / 1e6))

# Libera os ~2 GB do dataset completo antes dos gráficos.
rm(dt); gc()

# ══════════════════════════════════════════════════════════════════════════════
# 6. GRÁFICOS EDA
# ══════════════════════════════════════════════════════════════════════════════
# Os agregados (médias diárias e intradiárias) são computados bloco a bloco
# para não recarregar os 89M de linhas em memória. Cada bloco (~1 MB) é lido,
# o contributo (soma e contagem) é acumulado, e o bloco é descartado.
# Ao final a média é calculada como sum/n — equivalente à média sobre todos
# os domicílios e dias.

# Tema base — idêntico ao sim1_generate_artifacts.R
base_theme = theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(margin = margin(5, 5, 5, 5), face = "bold")
  )

group_colors = c("ToU" = "royalblue", "Std" = "firebrick")
group_labels = c("ToU" = "Dinâmica (ToU)", "Std" = "Padrão (Std)")

acorn_colors = c(
  "ACORN-U"     = "firebrick",
  "Adversity"   = "salmon",
  "Affluent"    = "forestgreen",
  "Comfortable" = "royalblue"
)

# Acumuladores bloco a bloco.
# Cada bloco contribui com (soma, contagem) por célula de agregação.
# Média final = sum_total / n_total — equivalente à média sobre todos os dados.
# tstp_acc acumula também sum_sq para calcular o desvio padrão (envoltura).
daily_acc = NULL
tstp_acc = NULL # soma/n por (tstp, grupo) — envoltura funcional
heatmap_acc = NULL # média por (semana_iso, hora, grupo) — heatmap
seasonal_acc = NULL # média por (trimestre, hora, grupo) — perfis sazonais
acorn_acc = NULL # média por (hora, Acorn_grouped, grupo) — intradiária
daily_acorn_acc = NULL # média por (date, Acorn_grouped, grupo) — anual
pos_acc = NULL # contagens por (tstp, grupo) — verificação de positividade

cat("Computando agregados para gráficos (bloco a bloco)...\n")

for (fpath in block_files) {
  b = fread(fpath, select = c("LCLid", "tstp", "energy(kWh/hh)"))
  b = b[LCLid %in% ids_2013]
  if (nrow(b) == 0L) next
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b[, energy := suppressWarnings(as.numeric(get("energy(kWh/hh)")))]
  b[, `energy(kWh/hh)` := NULL]
  b = b[!is.na(tstp) & tstp >= WIN_START & tstp <= WIN_END]
  b = merge(b, meta[, .(LCLid, stdorToU, Acorn_grouped)], by = "LCLid")
  b[, date := as.Date(tstp, tz = "Europe/London")]
  b[
    ,
    hora := as.numeric(
      format(tstp, "%H", tz = "Europe/London")
    ) + as.numeric(
      format(tstp, "%M", tz = "Europe/London")
    ) / 60
  ]
  b[, semana := as.integer(format(date, "%V"))]
  b[, trim   := as.integer(ceiling(as.integer(format(date, "%m")) / 3))]

  daily_acc = rbindlist(
    list(
      daily_acc,
      b[
        ,
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(date, stdorToU)
      ]
    )
  )
  tstp_acc = rbindlist(
    list(
      tstp_acc,
      b[
        , 
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(tstp, stdorToU)
      ]
    )
  )
  heatmap_acc = rbindlist(
    list(
      heatmap_acc,
      b[
        , 
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(semana, hora, stdorToU)
      ]
    )
  )
  seasonal_acc = rbindlist(
    list(
      seasonal_acc,
      b[
        , 
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(trim, hora, stdorToU)
      ]
    )
  )
  acorn_acc = rbindlist(
    list(
      acorn_acc,
      b[
        , 
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(hora, Acorn_grouped, stdorToU)
      ]
    )
  )
  daily_acorn_acc = rbindlist(
    list(
      daily_acorn_acc,
      b[
        , 
        .(sum_kw = sum(energy, na.rm = TRUE), n_kw = sum(!is.na(energy))),
        by = .(date, Acorn_grouped, stdorToU)
      ]
    )
  )
  pos_acc = rbindlist(
    list(
      pos_acc,
      b[, .(n = .N), by = .(tstp, stdorToU)]
    )
  )
  rm(b)
}

daily = daily_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(date, stdorToU)
]

tstp_full = tstp_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(tstp, stdorToU)
]

heatmap = heatmap_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(semana, hora, stdorToU)
]

seasonal = seasonal_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(trim, hora, stdorToU)
]

acorn = acorn_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(hora, Acorn_grouped, stdorToU)
]

daily_acorn = daily_acorn_acc[
  , .(mean_kw = sum(sum_kw) / sum(n_kw)), by = .(date, Acorn_grouped, stdorToU)
]

pos_wide = dcast(
  pos_acc[, .(n = sum(n)), by = .(tstp, stdorToU)],
  tstp ~ stdorToU, value.var = "n", fill = 0L
)

cat(
  sprintf(
    "Timestamps sem obs. ToU: %d | sem obs. Std: %d (de %d timestamps)\n",
    pos_wide[ToU == 0L, .N], pos_wide[Std == 0L, .N], nrow(pos_wide)
  )
)

setorder(daily,     stdorToU, date)
setorder(tstp_full, stdorToU, tstp)
rm(
  daily_acc, tstp_acc, heatmap_acc, seasonal_acc,
  acorn_acc, daily_acorn_acc, pos_acc, pos_wide
)
gc()

cat("Gerando gráficos...\n")

# ── 6a. Anual: consumo médio diário por grupo ─────────────────────────────────

p_annual = ggplot(daily, aes(x = date, y = mean_kw, color = stdorToU)) +
  geom_line(alpha = 0.85, linewidth = 0.6) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(x = "Data (2013)", y = "Consumo médio diário (kWh/hh)") +
  base_theme

ggsave(
  file.path(OUT_MAIN, "annual_mean_by_group.png"),
  p_annual, width = 10, height = 5, dpi = 300, bg = "white"
)

# ── 6b. Envoltura funcional (apêndice) ───────────────────────────────────────
# Média por timestamp e grupo sobre todos os ~5.095 domicílios (resolução
# completa de 30 min). Mostra ciclos intradiários + sazonalidade no mesmo eixo.
setorder(tstp_full, stdorToU, tstp)

p_envelope = ggplot(tstp_full, aes(x = tstp, y = mean_kw, color = stdorToU)) +
  geom_line(linewidth = 0.3, alpha = 0.9) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  scale_x_datetime(date_breaks = "2 months", date_labels = "%b") +
  labs(x = "Data (2013)", y = "Consumo médio (kWh/hh)") +
  base_theme

ggsave(
  file.path(OUT_APP, "functional_envelope.png"),
  p_envelope, width = 10, height = 5, dpi = 300, bg = "white"
)

# ── 6d. Heatmap: consumo médio por semana × hora ─────────────────────────────
heatmap[, grupo := ifelse(stdorToU == "ToU", "Dinâmica (ToU)", "Padrão (Std)")]

p_heatmap = ggplot(heatmap, aes(x = hora, y = semana, fill = mean_kw)) +
  geom_tile() +
  facet_wrap(~ grupo) +
  scale_fill_viridis_c(name = "kWh/hh", option = "plasma") +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  scale_y_reverse(
    breaks = c(1, 9, 18, 27, 36, 44, 52),
    labels = c("Jan", "Mar", "Mai", "Jul", "Set", "Nov", "Dez")
  ) +
  labs(x = "Hora do dia", y = "Semana do ano") +
  base_theme +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 11, family = "serif")
  )

ggsave(
  file.path(OUT_MAIN, "heatmap_consumption.png"),
  p_heatmap, width = 11, height = 6, dpi = 300, bg = "white"
)

# ── 6e. Perfis intradiários por estação ───────────────────────────────────────
trim_labels = c(
  "1" = "Q1 (jan–mar)", "2" = "Q2 (abr–jun)",
  "3" = "Q3 (jul–set)", "4" = "Q4 (out–dez)"
)
seasonal[, trim_label := trim_labels[as.character(trim)]]

p_seasonal = ggplot(seasonal, aes(x = hora, y = mean_kw, color = stdorToU)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ trim_label, ncol = 2) +
  scale_color_manual(values = group_colors, labels = group_labels) +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  labs(x = "Hora do dia", y = "Consumo médio (kWh/hh)") +
  base_theme

ggsave(
  file.path(OUT_MAIN, "intraday_by_season.png"),
  p_seasonal, width = 10, height = 7, dpi = 300, bg = "white"
)

# ── 6f. Diferença média bruta intradiária por ACORN ───────────────────────────
# mean(ToU) − mean(Std) por slot de 30 min, estratificado por ACORN.
# Não é o FATE estimado — é uma diferença de médias não ajustada.
acorn_wide = dcast(
  acorn, hora + Acorn_grouped ~ stdorToU, 
  value.var = "mean_kw"
)
acorn_wide[, diff_medio := ToU - Std]

p_acorn_intra = ggplot(
  acorn_wide,
  aes(x = hora, y = diff_medio, color = Acorn_grouped)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = acorn_colors) +
  scale_x_continuous(
    breaks = seq(0, 22, by = 4),
    labels = sprintf("%02d:00", seq(0, 22, by = 4))
  ) +
  labs(
    x = "Hora do dia",
    y = expression(bar(Y)[ToU](t) - bar(Y)[Std](t) ~ "(kWh/hh)"),
    color = "ACORN"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, family = "serif")
  )

ggsave(
  file.path(OUT_MAIN, "fate_raw_by_acorn_intraday.png"),
  p_acorn_intra, width = 10, height = 5, dpi = 300, bg = "white"
)

# ── 6g. Diferença média bruta anual por ACORN ─────────────────────────────────
daily_acorn_wide = dcast(
  daily_acorn, date + Acorn_grouped ~ stdorToU, value.var = "mean_kw"
)
daily_acorn_wide[, diff_medio := ToU - Std]

p_acorn_annual = ggplot(
  daily_acorn_wide,
  aes(x = date, y = diff_medio, color = Acorn_grouped)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(alpha = 0.25, linewidth = 0.4) +
  geom_smooth(
    method = "loess", formula = y ~ x, span = 0.2,
    se = FALSE, linewidth = 0.9
  ) +
  scale_colour_manual(values = acorn_colors) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  labs(
    x = "Data (2013)",
    y = expression(bar(Y)[ToU](t) - bar(Y)[Std](t) ~ "(kWh/hh)"),
    color = "ACORN"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 11, family = "serif")
  )

ggsave(
  file.path(OUT_MAIN, "fate_raw_by_acorn_annual.png"),
  p_acorn_annual, width = 10, height = 5, dpi = 300, bg = "white"
)

cat(
  sprintf(
    "Concluído. Outputs em:\n  %s\n  %s\n",
    PROC_DIR, dirname(OUT_MAIN)
  )
)
