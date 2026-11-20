# scripts/sim1_add_cfate.R
# Agrega métricas de CFATE dos checkpoints e atualiza os CSVs da sim 1.
# Idempotente: pode ser re-executado sem problema.

library(data.table)

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

ckpt_base <- "output/sim1_checkpoints"
out_dir <- "output/sim1_final"
T_TREAT <- 50.0

config_dirs <- sort(list.dirs(ckpt_base, recursive = FALSE))
n_configs <- length(config_dirs)
cat(sprintf("Configs encontradas: %d\n", n_configs))

# =============================================================================
# 1. Lê checkpoints → métricas CFATE por run e por ponto
# =============================================================================

cfate_raw_list <- vector("list", n_configs)
cfate_curve_list <- vector("list", n_configs)

for (i in seq_along(config_dirs)) {
  batch_files <- list.files(
    config_dirs[i],
    pattern = "\\.rds$", full.names = TRUE
  )
  cat(sprintf("  [%d/%d] %s\r", i, n_configs, basename(config_dirs[i])))

  curves <- rbindlist(lapply(batch_files, function(f) readRDS(f)$curves))

  # Métricas por run (t > T_TREAT, onde CFATE é não-trivial)
  cfate_raw_list[[i]] <- curves[t_val > T_TREAT, .(
    CFATE_RMSE = sqrt(mean((CFATE_Est - CFATE_True)^2, na.rm = TRUE)),
    CFATE_Cov_Simul = as.numeric(
      all(
        CFATE_True >= CFATE_Lower & CFATE_True <= CFATE_Upper
      )
    ),
    CFATE_Cov_Point = mean(
      CFATE_True >= CFATE_Lower & CFATE_True <= CFATE_Upper,
      na.rm = TRUE
    )
  ), by = .(SimID, ErrorDist, Scenario, Method, N, T_pts)]

  # Resumo ponto-a-ponto (todos os t para manter alinhamento com curve_summary)
  cfate_curve_list[[i]] <- curves[, .(
    CFATE_RMSE_t = sqrt(mean((CFATE_Est - CFATE_True)^2, na.rm = TRUE)),
    CFATE_Bias_t = mean(CFATE_Est - CFATE_True, na.rm = TRUE),
    CFATE_Cov_t = mean(
      CFATE_True >= CFATE_Lower & CFATE_True <= CFATE_Upper,
      na.rm = TRUE
    ),
    CFATE_Width_t = mean(CFATE_Upper - CFATE_Lower, na.rm = TRUE)
  ), by = .(ErrorDist, Scenario, Method, N, T_pts, t_val)]

  rm(curves)
  gc(verbose = FALSE)
}
cat("\n")

all_cfate_raw <- rbindlist(cfate_raw_list)
all_cfate_curve <- rbindlist(cfate_curve_list)
rm(cfate_raw_list, cfate_curve_list)
gc(verbose = FALSE)

# =============================================================================
# 2. Normaliza CFATE_RMSE pelo range da curva verdadeira (via exemplares)
# =============================================================================

dt_ex <- fread(file.path(out_dir, "sim1_exemplars.csv"))

cfate_ranges <- dt_ex[
  t_val > T_TREAT,
  .(cfate_range = diff(range(CFATE_True))),
  by = .(Scenario, N, T_pts, ErrorDist)
]

all_cfate_raw <- merge(
  all_cfate_raw, cfate_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist"), all.x = TRUE
)
all_cfate_raw[, CFATE_NRMSE := round(CFATE_RMSE / cfate_range, 4)]
all_cfate_raw[, c("CFATE_RMSE", "cfate_range") := NULL]

all_cfate_curve <- merge(
  all_cfate_curve, cfate_ranges,
  by = c("Scenario", "N", "T_pts", "ErrorDist"), all.x = TRUE
)
all_cfate_curve[, CFATE_NRMSE_t := CFATE_RMSE_t / cfate_range]
all_cfate_curve[, c("CFATE_RMSE_t", "cfate_range") := NULL]

# =============================================================================
# 3. Atualiza sim1_raw.csv
# =============================================================================

cat("Atualizando sim1_raw.csv...\n")
dt_raw <- fread(file.path(out_dir, "sim1_raw.csv"))

# Remove colunas CFATE se já existirem (idempotência)
cfate_cols <- c("CFATE_NRMSE", "CFATE_Cov_Simul", "CFATE_Cov_Point")
dt_raw[, (intersect(cfate_cols, names(dt_raw))) := NULL]

dt_raw <- merge(
  dt_raw, all_cfate_raw,
  by = c("SimID", "ErrorDist", "Scenario", "Method", "N", "T_pts"),
  all.x = TRUE
)
fwrite(dt_raw, file.path(out_dir, "sim1_raw.csv"))

# =============================================================================
# 4. Recomputa sim1_summary.csv
# =============================================================================

cat("Recomputando sim1_summary.csv...\n")
summary_dt <- dt_raw[, .(
  N_sims = .N,
  Cov_Simul = round(mean(Cov_Simul, na.rm = TRUE) * 100, 1),
  Cov_Point = round(mean(Cov_Point, na.rm = TRUE) * 100, 1),
  NRMSE = round(mean(NRMSE, na.rm = TRUE), 4),
  NRMSE_SD = round(sd(NRMSE, na.rm = TRUE), 4),
  Bias = round(mean(Bias, na.rm = TRUE), 4),
  Width = round(mean(Width, na.rm = TRUE), 3),
  CFATE_NRMSE = round(mean(CFATE_NRMSE, na.rm = TRUE), 4),
  CFATE_Cov_Simul = round(mean(CFATE_Cov_Simul, na.rm = TRUE) * 100, 1),
  CFATE_Cov_Point = round(mean(CFATE_Cov_Point, na.rm = TRUE) * 100, 1)
), by = .(ErrorDist, Scenario, Method, N, T_pts)]
fwrite(summary_dt, file.path(out_dir, "sim1_summary.csv"))

# =============================================================================
# 5. Atualiza sim1_curve_summary.csv
# =============================================================================

cat("Atualizando sim1_curve_summary.csv...\n")
dt_csumm <- fread(file.path(out_dir, "sim1_curve_summary.csv"))

cfate_curve_cols <- c(
  "CFATE_NRMSE_t", "CFATE_Bias_t",
  "CFATE_Cov_t", "CFATE_Width_t"
)
dt_csumm[, (intersect(cfate_curve_cols, names(dt_csumm))) := NULL]

dt_csumm <- merge(
  dt_csumm, all_cfate_curve,
  by = c("ErrorDist", "Scenario", "Method", "N", "T_pts", "t_val"),
  all.x = TRUE
)
fwrite(dt_csumm, file.path(out_dir, "sim1_curve_summary.csv"))

cat("\nConcluído.\n")
print(
  summary_dt[
    ,
    .(ErrorDist, Scenario, Method, N, T_pts, CFATE_NRMSE, CFATE_Cov_Simul)
  ]
)
