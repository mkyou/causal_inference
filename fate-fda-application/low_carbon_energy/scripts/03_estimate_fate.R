# =============================================================================
# 03_estimate_fate.R
# Estimação do FATE e do CFATE com rLifting 1.0.0 (Haar, offline e causal)
#
# Estrutura:
#   §1  Setup + instalação rLifting 1.0.0
#   §2  Dados processados + grade de referência
#   §3  Validação de configuração (extension × levels) e tunagem de α/β
#   §4  Fase 1a – suavização offline → smooth_offline.rds
#   §5  Fase 1b – suavização causal → smooth_causal.rds
#   §6  Fase 1c – matriz bruta → raw_wide.rds  (libera dt)
#   §7  Fase 2  – estimação FATE + delta_beta + FATE/ACORN + bootstrap + bandas
#   §8  Cálculo do CFATE + bandas
#   §9  Guardar resultados
#
# Memória máxima em cada fase: dt (2.1 GB) + uma matriz N×T (~714 MB) = 2.8 GB
# Bootstrap: raw_wide (714 MB) + boot_mat (5095 × 17520 × 8 ≈ 714 MB) = 1.4 GB
# =============================================================================

library(data.table)
library(here)
library(parallel)

# =============================================================================
# §1  Instalar / actualizar rLifting 1.0.0 a partir do código-fonte
# =============================================================================

pkg_src = normalizePath(file.path(dirname(dirname(here())), "rLifting"))
installed_ver = tryCatch(
  packageVersion("rLifting"),
  error = function(e) package_version("0.0.0")
)

if (installed_ver < package_version("1.0.0")) {
  message("Instalando rLifting 1.0.0 a partir do código-fonte...")
  install.packages(pkg_src, repos = NULL, type = "source")
}

library(rLifting)
message("rLifting: ", packageVersion("rLifting"))

# =============================================================================
# Funções de estimação
# =============================================================================

# Imputa média da coluna em matrizes N×T com NAs esparsos (MCAR).
impute_col_mean = function(mat) {
  if (!anyNA(mat)) return(mat)
  na_cols = which(colSums(is.na(mat)) > 0L)
  col_means = colMeans(mat[, na_cols, drop = FALSE], na.rm = TRUE)
  for (k in seq_along(na_cols)) {
    j = na_cols[k]
    mat[is.na(mat[, j]), j] = col_means[k]
  }
  mat
}

precompute_projection = function(X_reg, A_vec, ridge_lambda = 1e-7) {
  x_bar = colMeans(X_reg)
  idx_1 = which(A_vec == 1)
  idx_0 = which(A_vec == 0)

  calc = function(X) {
    XtX = crossprod(X)
    diag(XtX) = diag(XtX) + ridge_lambda
    XtX_inv = solve(XtX)
    list(
      H = XtX_inv %*% t(X), 
      var_factor = as.numeric(t(x_bar) %*% XtX_inv %*% x_bar)
    )
  }

  m1 = calc(X_reg[idx_1, , drop = FALSE])
  m0 = calc(X_reg[idx_0, , drop = FALSE])

  list(
    H1 = m1$H, H0 = m0$H,
    var_factor_1 = m1$var_factor, var_factor_0 = m0$var_factor,
    x_bar = x_bar,
    idx_1 = idx_1, idx_0 = idx_0,
    N1 = length(idx_1), N0 = length(idx_0),
    p = ncol(X_reg)
  )
}

# Dados chegam já imputados — OLS vetorizado direto.
# Retorna fate_est (T-vetor) e delta_beta (p×T) para uso por ACORN.
estimate_fate_fast = function(Y_wide, proj) {
  beta1 = proj$H1 %*% Y_wide[proj$idx_1, , drop = FALSE]
  beta0 = proj$H0 %*% Y_wide[proj$idx_0, , drop = FALSE]
  delta_beta = beta1 - beta0
  list(
    fate_est   = as.numeric(proj$x_bar %*% delta_beta),
    delta_beta = delta_beta
  )
}

# Bootstrap paralelo com hat-weights (evita submatrizes N×T).
# n_cores > 1 usa mclapply (fork-based, Linux/macOS).
run_bootstrap_raw = function(
  Y_wide, X_reg, A_vec, n_boot,
  ridge_lambda = 1e-7, n_cores = 1L
) {
  N = nrow(Y_wide)
  T_pts = ncol(Y_wide)
  x_bar = colMeans(X_reg)
  p = ncol(X_reg)

  one_boot = function(b) {
    idx = sample(N, N, replace = TRUE)
    A_b = A_vec[idx]
    i1  = which(A_b == 1)
    i0  = which(A_b == 0)

    if (length(i1) < p + 5L || length(i0) < p + 5L) return(rep(NA_real_, T_pts))

    m1 = tabulate(idx[i1], nbins = N)
    m0 = tabulate(idx[i0], nbins = N)
    act1 = which(m1 > 0L)
    act0 = which(m0 > 0L)

    X1 = X_reg[act1, , drop = FALSE]
    X0 = X_reg[act0, , drop = FALSE]

    XtX1 = crossprod(X1 * sqrt(m1[act1]))
    diag(XtX1) = diag(XtX1) + ridge_lambda
    
    XtX0 = crossprod(X0 * sqrt(m0[act0]))
    diag(XtX0) = diag(XtX0) + ridge_lambda

    Omega1 = tryCatch(solve(XtX1), error = function(e) NULL)
    Omega0 = tryCatch(solve(XtX0), error = function(e) NULL)
    if (is.null(Omega1) || is.null(Omega0)) return(rep(NA_real_, T_pts))

    h1 = m1[act1] * as.numeric(x_bar %*% Omega1 %*% t(X1))
    h0 = m0[act0] * as.numeric(x_bar %*% Omega0 %*% t(X0))

    w_full = numeric(N)
    w_full[act1] = h1
    w_full[act0] = -h0

    as.numeric(w_full %*% Y_wide)
  }

  n_chunks  = 10L
  chunk_ids = split(seq_len(n_boot),
                    ceiling(seq_len(n_boot) * n_chunks / n_boot))
  res = vector("list", n_chunks)
  for (ci in seq_len(n_chunks)) {
    res[[ci]] = if (n_cores > 1L)
      parallel::mclapply(chunk_ids[[ci]], one_boot, mc.cores = n_cores)
    else
      lapply(chunk_ids[[ci]], one_boot)
    message(sprintf("  Bootstrap: %3d%%  (%d/%d)",
                    round(max(chunk_ids[[ci]]) / n_boot * 100),
                    max(chunk_ids[[ci]]), n_boot))
  }

  do.call(rbind, unlist(res, recursive = FALSE))
}

build_simultaneous_band = function(theta_hat, boot_curves, alpha = 0.05) {
  valid = complete.cases(boot_curves)
  if (sum(valid) < 50L) {
    return(
      list(
        lower = rep(-Inf, length(theta_hat)),
        upper = rep(Inf,  length(theta_hat))
      )
    )
  }

  devs = apply(
    boot_curves[valid, , drop = FALSE], 1,
    function(r) max(abs(r - theta_hat), na.rm = TRUE)
  )
  q = quantile(devs, 1 - alpha, na.rm = TRUE)

  list(lower = theta_hat - q, upper = theta_hat + q)
}

# =============================================================================
# §2  Dados processados + grade de referência
# =============================================================================
DATA_DIR = here("low_carbon_energy", "data", "processed")

message("Carregando lcl_2013_long.rds (~2 GB)...")
dt = readRDS(file.path(DATA_DIR, "lcl_2013_long.rds"))
meta = fread(file.path(DATA_DIR, "X_meta.csv"))

setkey(dt, LCLid)
setkey(meta, LCLid)

ref_grid = sort(unique(dt$tstp)) # 17 520 timestamps
T_pts = length(ref_grid) # 17 520
ref_dt = data.table(tstp = ref_grid)

lclids = meta$LCLid # ordem canónica (5 095)
N = length(lclids)
A_vec = meta$A # 1 = ToU, 0 = Std

# X_reg: intercepto + dummies ACORN (referência = Adversity)
X_reg = model.matrix(~ Acorn_grouped, data = meta)
x_bar = colMeans(X_reg)

message(
  sprintf(
    "N = %d  |  T = %d  |  ACORN cols = %d",
    N, T_pts, ncol(X_reg) - 1L
  )
)

# =============================================================================
# §3  Validação de configuração e tunagem de α/β
#
# (a) Compara (extension × levels) pelo SURE mínimo em 30 domicílios amostrais.
# (b) Para a config escolhida, compara SURE(defaults α=0.3, β=1.2) com
#     SURE(tunado). Adopta valores tunados só se melhoria ≥ 2%.
# =============================================================================

set.seed(42L)
tune_ids = c(
  lclids[sample(which(A_vec == 1L), 15L)],
  lclids[sample(which(A_vec == 0L), 15L)]
)

scheme_haar = lifting_scheme("haar")

# --- (a) Escolha de extension × levels ---------------------------------------

configs = CJ(extension = c("symmetric", "zero"), levels = c(3L, 4L, 5L))

message("\n--- Validação: extension × levels (n=30) ---")

extract_and_interp = function(id) {
  hh = dt[.(id), .(tstp, energy)]
  hh_full = hh[ref_dt, on = "tstp"]
  if (anyNA(hh_full$energy)) {
    ok = which(!is.na(hh_full$energy))
    hh_full[
      ,
      energy := approx(
        ok, energy[ok],
        xout = seq_len(T_pts),
        rule = 2L
      )$y
    ]
  }
  hh_full$energy
}

configs[
  ,
  sure_med := {
    ext = extension; lev = levels
    sures = vapply(
      tune_ids,
      function(id) {
        sig = extract_and_interp(id)
        tryCatch(
          tune_alpha_beta(sig, scheme_haar, levels = lev, extension = ext)$sure,
          error = function(e) NA_real_
        )
      },
      FUN.VALUE = numeric(1L)
    )
    median(sures, na.rm = TRUE)
  },
  by = .(extension, levels)
]

print(configs[order(sure_med)])

best = configs[which.min(sure_med)]
EXT_OFF = best$extension
EXT_CAU = "zero" # causal: mantém o path Lagrange (vignette v05)
LEVELS = best$levels

message(
  sprintf(
    "Config: levels=%d | offline ext='%s' | causal ext='%s'",
    LEVELS, EXT_OFF, EXT_CAU
  )
)

# --- (b) Comparar tunado vs defaults -----------------------------------------
message("\n--- Validação: tunado vs defaults ---")

alpha_default = 0.3
beta_default = 1.2

compare_tbl = vapply(
  tune_ids,
  function(id) {
    sig = extract_and_interp(id)
    tuned = tune_alpha_beta(
      sig, scheme_haar,
      levels = LEVELS, extension = EXT_OFF
    )
    sure_def = rLifting:::.sure_alpha_beta(
      sig, scheme_haar, levels = LEVELS, extension = EXT_OFF,
      alpha = alpha_default, beta = beta_default
    )
    c(
      sure_tuned = tuned$sure,
      sure_def = sure_def,
      alpha = tuned$alpha,
      beta = tuned$beta
    )
  },
  FUN.VALUE = numeric(4L)
)

med_sure_tuned = median(compare_tbl["sure_tuned", ])
med_sure_def = median(compare_tbl["sure_def",   ])
improvement_pct = (med_sure_def - med_sure_tuned) / abs(med_sure_def) * 100

message(
  sprintf(
    "  SURE mediana: defaults=%.2f  tunado=%.2f  melhoria=%.1f%%",
    med_sure_def, med_sure_tuned, improvement_pct
  )
)

if (med_sure_tuned < med_sure_def) {
  ALPHA_HAT = median(compare_tbl["alpha", ])
  BETA_HAT = median(compare_tbl["beta",  ])
  message(
    sprintf(
      "  Usando α/β tunados: α=%.3f  β=%.3f",
      ALPHA_HAT, BETA_HAT
    )
  )
} else {
  ALPHA_HAT = alpha_default
  BETA_HAT = beta_default
  message("  Tunado não melhora: usando defaults α=0.3, β=1.2")
}

# Window size para modo causal (W odd, W >= 16 * 2^levels)
W = max(337L, 16L * 2L^LEVELS)
if (W %% 2L == 0L) W = W + 1L
message(sprintf("Window causal: W=%d", W))

N_CORES = max(1L, parallel::detectCores() - 1L)
message(
  sprintf(
    "Paralelismo: %d / %d núcleos.", 
    N_CORES, 
    parallel::detectCores()
  )
)

# =============================================================================
# Auxiliares: alinhamento e suavização com grade irregular nativa
# =============================================================================

fill_household = function(id) {
  hh = dt[.(id), .(tstp, energy)]
  hh_full = hh[ref_dt, on = "tstp"]
  obs_idx = which(!is.na(hh_full$energy))
  list(
    signal = hh_full$energy[obs_idx],
    t = if (length(obs_idx) < T_pts) as.double(obs_idx) else NULL,
    obs_idx = obs_idx
  )
}

to_full_grid = function(sm, obs_idx) {
  if (length(obs_idx) == T_pts) return(sm)
  full = rep(NA_real_, T_pts)
  full[obs_idx] = sm
  full
}

fill_household_raw = function(id) {
  hh = dt[.(id), .(tstp, energy)]
  hh_full = hh[ref_dt, on = "tstp"]
  hh_full$energy
}

# =============================================================================
# §4  Fase 1a – suavização offline → smooth_offline.rds
# Pico: dt (2.1 GB) + smooth_offline (714 MB) = 2.8 GB
# =============================================================================

OUT_OFF = file.path(DATA_DIR, "smooth_offline.rds")

if (!file.exists(OUT_OFF)) {
  message(sprintf("\n--- Fase 1a: offline (N=%d, cores=%d) ---", N, N_CORES))
  setDTthreads(1L)
  rows = parallel::mclapply(seq_len(N), function(i) {
    hh = fill_household(lclids[i])
    sm = denoise_signal_offline(
      hh$signal, scheme_haar,
      levels = LEVELS,
      extension = EXT_OFF,
      threshold_method = "universal",
      shrinkage = "scad",
      alpha = ALPHA_HAT,
      beta = BETA_HAT,
      t = hh$t
    )
    to_full_grid(sm, hh$obs_idx)
  }, mc.cores = N_CORES)
  setDTthreads(0L)
  smooth_offline = do.call(rbind, rows)
  saveRDS(smooth_offline, OUT_OFF, compress = FALSE)
  message("  smooth_offline.rds salvo.")
  rm(smooth_offline, rows); gc()
} else {
  message("smooth_offline.rds já existe — pulando fase 1a.")
}

# =============================================================================
# §5  Fase 1b – suavização causal → smooth_causal.rds
# Pico: dt (2.1 GB) + smooth_causal (714 MB) = 2.8 GB
# =============================================================================

OUT_CAU = file.path(DATA_DIR, "smooth_causal.rds")

if (!file.exists(OUT_CAU)) {
  message(sprintf("\n--- Fase 1b: causal (N=%d, W=%d, L=%d, cores=%d) ---", N, W, LEVELS, N_CORES))
  setDTthreads(1L)
  rows = parallel::mclapply(seq_len(N), function(i) {
    hh = fill_household(lclids[i])
    sm = denoise_signal_causal(
      hh$signal, scheme_haar,
      levels = LEVELS,
      window_size = W,
      extension = EXT_CAU,
      threshold_method = "universal",
      shrinkage = "scad",
      alpha = ALPHA_HAT,
      beta = BETA_HAT,
      update_freq = 1L,
      t = hh$t
    )
    to_full_grid(sm, hh$obs_idx)
  }, mc.cores = N_CORES)
  setDTthreads(0L)
  smooth_causal = do.call(rbind, rows)
  saveRDS(smooth_causal, OUT_CAU, compress = FALSE)
  message("  smooth_causal.rds salvo.")
  rm(smooth_causal, rows); gc()
} else {
  message("smooth_causal.rds já existe — pulando fase 1b.")
}

# =============================================================================
# §6  Fase 1c – matriz bruta (sem suavização) → raw_wide.rds
# Necessária porque o bootstrap é sobre dados ruidosos (§sec-bootstrap-ic).
# Pico: dt (2.1 GB) + raw_wide (714 MB) = 2.8 GB
# =============================================================================

OUT_RAW = file.path(DATA_DIR, "raw_wide.rds")

if (!file.exists(OUT_RAW)) {
  message(sprintf("\n--- Fase 1c: raw_wide (N=%d, cores=%d) ---", N, N_CORES))
  setDTthreads(1L)
  rows = parallel::mclapply(seq_len(N), function(i) {
    fill_household_raw(lclids[i])
  }, mc.cores = N_CORES)
  setDTthreads(0L)
  raw_wide = do.call(rbind, rows)
  saveRDS(raw_wide, OUT_RAW, compress = FALSE)
  message("  raw_wide.rds salvo.")
  rm(raw_wide, rows); gc()
} else {
  message("raw_wide.rds já existe — pulando fase 1c.")
}

rm(dt)
gc() # dt não é mais necessário

# =============================================================================
# §7  Fase 2 – estimação FATE + bootstrap + bandas simultâneas
#
# - Estimativa pontual: Y_smooth (offline ou causal) via precompute_projection
#    e estimate_fate_fast.
# - Bootstrap: run_bootstrap_raw(raw_wide, ...) com B = N = 5095
# (§sec-bootstrap-ic da dissertação: reamostrar os dados ruidosos.
# banda: build_simultaneous_band(theta_hat, boot_mat, alpha = 0.05))
#
# Nota: boot_mat (B×T) é partilhado pelas duas estimativas (offline e causal)
# — as curvas bootstrap não dependem do modo de suavização. Apenas o D_b
# muda (depende de theta_hat que difere por modo).
# =============================================================================

N_BOOT = N # B = N = 5095
ALPHA_BAND = 0.05

proj = precompute_projection(X_reg, A_vec)

# --- 7.1 Estimativas pontuais + delta_beta + FATE por ACORN ------------------

message("\n--- Fase 2: estimativas pontuais ---")

smooth_off_mat = impute_col_mean(readRDS(OUT_OFF))
res_off = estimate_fate_fast(smooth_off_mat, proj)
fate_offline = res_off$fate_est
delta_beta_off = res_off$delta_beta
rm(smooth_off_mat, res_off)
gc()

saveRDS(
  delta_beta_off, 
  file.path(DATA_DIR, "delta_beta_offline.rds"), 
  compress = FALSE
)

smooth_cau_mat = impute_col_mean(readRDS(OUT_CAU))
res_cau = estimate_fate_fast(smooth_cau_mat, proj)
fate_causal = res_cau$fate_est
delta_beta_cau = res_cau$delta_beta
rm(smooth_cau_mat, res_cau) 
gc()

saveRDS(
  delta_beta_cau, 
  file.path(DATA_DIR, "delta_beta_causal.rds"),
  compress = FALSE
)

message(sprintf("  FATE offline: média=%.4f kWh/hh", mean(fate_offline)))
message(sprintf("  FATE causal:  média=%.4f kWh/hh", mean(fate_causal)))

# FATE por grupo ACORN: x_k' · delta_beta (estimativas pontuais)
acorn_levels = levels(factor(meta$Acorn_grouped))
x_acorn = lapply(acorn_levels, function(grp) {
  x = setNames(numeric(ncol(X_reg)), colnames(X_reg))
  x["(Intercept)"] = 1
  col_nm = paste0("Acorn_grouped", grp)
  if (col_nm %in% names(x)) x[col_nm] = 1
  x
})
names(x_acorn) = acorn_levels

fate_acorn_off = lapply(x_acorn, function(xk) as.numeric(xk %*% delta_beta_off))
fate_acorn_cau = lapply(x_acorn, function(xk) as.numeric(xk %*% delta_beta_cau))
saveRDS(
  list(offline = fate_acorn_off, causal = fate_acorn_cau),
  file.path(DATA_DIR, "fate_by_acorn.rds"), compress = FALSE
)
rm(delta_beta_off, delta_beta_cau); gc()
message("  delta_beta e fate_by_acorn salvos.")

# --- 7.2 Bootstrap sobre dados ruidosos (B = N) ------------------------------

message(
  sprintf(
    "\n--- Bootstrap: B=%d sobre raw_wide (cores=%d) ---", 
    N_BOOT, N_CORES
  )
)

raw_wide = impute_col_mean(readRDS(OUT_RAW))
message(
  sprintf(
    "  raw_wide imputado: %d NAs remanescentes.", 
    sum(is.na(raw_wide))
  )
)

set.seed(42L)
boot_mat = run_bootstrap_raw(raw_wide, X_reg, A_vec, N_BOOT, n_cores = N_CORES)
rm(raw_wide)
gc()

message(
  sprintf(
    "  boot_mat: %d x %d  (válidos: %d)",
    nrow(boot_mat), ncol(boot_mat),
    sum(complete.cases(boot_mat))
  )
)

# --- 7.3 Bandas simultâneas --------------------------------------------------

band_offline = build_simultaneous_band(fate_offline, boot_mat, ALPHA_BAND)
band_causal  = build_simultaneous_band(fate_causal,  boot_mat, ALPHA_BAND)

message(
  sprintf(
    "  Banda offline: largura média=%.4f kWh/hh",
    mean(band_offline$upper - band_offline$lower, na.rm = TRUE)
  )
)
message(
  sprintf(
    "  Banda causal:  largura média=%.4f kWh/hh",
    mean(band_causal$upper  - band_causal$lower,  na.rm = TRUE)
  )
)

# =============================================================================
# §8  CFATE – Cumulative Functional Average Treatment Effect
#
# CFATE(t) = ∫₀ᵗ θ(s) ds  ≈  regra trapezoidal uniforme
#          = c(0, cumsum((θ[1:(T-1)] + θ[2:T]) / 2))   [kWh/hh acumulados]
#
# Cada slot é 30 min = 0.5 h, e energy já está em kWh/slot (não kW),
# portanto o passo da integral é 1 slot (sem factor ×0.5).
# A regra trapezoidal trata θ como função linear entre slots consecutivos.
#
# Bandas: CFATE_b via mesma regra sobre cada linha de boot_mat.
# =============================================================================

cumtrapz = function(y) c(0, cumsum((head(y, -1L) + y[-1L]) / 2))

cfate_offline = cumtrapz(fate_offline)
cfate_causal  = cumtrapz(fate_causal)

# CFATE para cada curva bootstrap (B × T)
boot_cfate = t(apply(boot_mat, 1L, cumtrapz))   # B × T

band_cfate_off = build_simultaneous_band(cfate_offline, boot_cfate, ALPHA_BAND)
band_cfate_cau = build_simultaneous_band(cfate_causal,  boot_cfate, ALPHA_BAND)

message(
  sprintf(
    "\n  CFATE(T) offline = %.2f kWh/hh/ano  [%.2f, %.2f]",
    cfate_offline[T_pts],
    band_cfate_off$lower[T_pts],
    band_cfate_off$upper[T_pts]
  )
)
message(
  sprintf(
    "  CFATE(T) causal  = %.2f kWh/hh/ano  [%.2f, %.2f]",
    cfate_causal[T_pts],
    band_cfate_cau$lower[T_pts],
    band_cfate_cau$upper[T_pts]
  )
)

rm(boot_mat, boot_cfate); gc()

# =============================================================================
# §9  Guardar resultados
# =============================================================================

OUT_DIR = here("low_carbon_energy", "output")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

fate_results = rbindlist(
  list(
    data.table(
      tstp = ref_grid,
      fate = fate_offline,
      lower = band_offline$lower,
      upper = band_offline$upper,
      cfate = cfate_offline,
      cfate_low = band_cfate_off$lower,
      cfate_up = band_cfate_off$upper,
      mode = "offline"
    ),
    data.table(
      tstp = ref_grid,
      fate = fate_causal,
      lower = band_causal$lower,
      upper = band_causal$upper,
      cfate = cfate_causal,
      cfate_low = band_cfate_cau$lower,
      cfate_up = band_cfate_cau$upper,
      mode = "causal"
    )
  )
)

saveRDS(
  fate_results,
  file.path(DATA_DIR, "fate_estimates.rds"),
  compress = FALSE
)
fwrite(fate_results, file.path(OUT_DIR, "fate_estimates.csv"))

message("\nfate_estimates.rds + .csv salvos.")
message("Script 03 concluído.")
