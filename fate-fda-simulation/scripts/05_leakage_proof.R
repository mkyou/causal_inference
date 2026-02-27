library(data.table)
library(dplyr)
library(fda)
library(rLifting)

# Wrapper para Splines
smooth_spline_robust = function(time_grid, wide_matrix, nbasis = 50) {
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), 
    nbasis = nbasis
  )
  
  # Busca simples por GCV para lambda
  log_lambdas = seq(-10, 1, length.out = 20)
  gcvs = numeric(length(log_lambdas))
  test_data = t(wide_matrix)
  
  for(i in seq_along(log_lambdas)) {
    fdPar_obj = fdPar(
      bspline_basis, Lfdobj = 2, lambda = 10^log_lambdas[i]
    )
    fd_try = smooth.basis(time_grid, test_data, fdPar_obj)
    gcvs[i] = sum(fd_try$gcv)
  }
  
  best_lambda = 10^log_lambdas[which.min(gcvs)]
  final_fdPar = fdPar(bspline_basis, Lfdobj = 2, lambda = best_lambda)
  final_fd = smooth.basis(time_grid, t(wide_matrix), final_fdPar)$fd
  
  return(t(eval.fd(time_grid, final_fd)))
}

# Wrapper Genérico para Ondaletas (Offline/Causal + Esquema variável)
smooth_wavelet_wrapper = function(
    y_vec, scheme_name, type = c("causal", "offline"), window = 64
) {
  scheme_obj = rLifting::lifting_scheme(scheme_name)
  lvl = floor(log2(window))
  
  if (type == "causal") {
    y_smooth = rLifting::denoise_signal_causal(
      y_vec,
      scheme = scheme_obj,
      levels = lvl,
      window_size = window,
      method = "semisoft",
      alpha = 0.3, 
      beta = 1.2
    )
  } else {
    # Para offline, ajustamos levels com base no tamanho total n
    n = length(y_vec)
    lvl_off = floor(log2(n))
    y_smooth = rLifting::denoise_signal_offline(
      y_vec,
      scheme = scheme_obj,
      levels = lvl_off,
      method = "semisoft",
      alpha = 0.3, 
      beta = 1.2
    )
  }
  return(y_smooth)
}

# Experimento
run_leakage_metric = function() {
  
  # Parâmetros
  T_total = 128
  t_impulse = 64
  magnitude = 10
  
  # 1. Criação do Sinal Delta
  time_grid = 1:T_total
  y_impulse = rep(0, T_total)
  y_impulse[t_impulse] = magnitude
  
  # 2. Aplicação dos Métodos
  
  # A. Splines
  y_spline = smooth_spline_robust(
    time_grid, matrix(y_impulse, nrow=1), nbasis = 80
  )[1,]
  
  # B. Ondaletas (Curtas vs Longas | Offline vs Causal)
  # Curta = db2 (4 coeficientes efetivos aprox)
  # Longa = dd4 (Biorthogonal 9/7 - usada no JPEG2000, suporte longo)
  
  y_off_short = smooth_wavelet_wrapper(y_impulse, "db2", "offline")
  y_off_long  = smooth_wavelet_wrapper(y_impulse, "dd4", "offline") 
  
  y_causal_short = smooth_wavelet_wrapper(
    y_impulse, "db2", "causal", window = 64
  )
  y_causal_long  = smooth_wavelet_wrapper(
    y_impulse, "dd4", "causal", window = 64
  )
  
  # 3. Cálculo da Métrica: Vazamento de Energia (%)
  
  calc_leak_energy = function(y_est, t_evt) {
    # Energia antes do evento / Energia total
    # Se causal, energy_pre deve ser 0 (ou muito próximo de zero numérico)
    energy_pre = sum(y_est[1:(t_evt-1)]^2)
    energy_total = sum(y_est^2)
    
    if (energy_total == 0) return(0)
    return((energy_pre / energy_total) * 100)
  }
  
  metrics = data.table(
    Method = c(
      "Spline", 
      "Offline (Short - db2)", "Offline (Medium - dd4)", 
      "Causal (Short - db2)", "Causal (Medium - dd4)"
    ),
    Leakage_Energy_Percent = c(
      calc_leak_energy(y_spline, t_impulse),
      calc_leak_energy(y_off_short, t_impulse),
      calc_leak_energy(y_off_long, t_impulse),
      calc_leak_energy(y_causal_short, t_impulse),
      calc_leak_energy(y_causal_long, t_impulse)
    )
  )
  
  return(metrics)
}

# Execução e salvamento
results_table = run_leakage_metric()

# Exibe no console para conferência rápida
cat(">>> Métricas de Vazamento de Energia <<<\n")
print(results_table)

# Salva em CSV
dir.create("output", showWarnings = FALSE)
write.csv(
  results_table, 
  "output/table_leakage_metric.csv", 
  row.names = FALSE
)

cat("\nTabela salva em 'output/table_leakage_metric.csv'.\n")