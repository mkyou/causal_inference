library(data.table)
library(future.apply)
library(rLifting)
library(fda)
library(ggplot2)
library(patchwork)

source("R/scenario_runner.R")
source("R/process_functional_data.R")
source("R/estimate_fate.R")

# Configuração
N_SIM = 1
N_BOOT = 200
N_SUBJ = 200
T_PTS = 128

SCENARIOS = list(
  list(name = "Sigmoid", mag = 2, family = "cdf97"),
  list(name = "Bumps", mag = 3, family = "db2"),
  list(name = "Doppler", mag = 5, family = "db2"),
  list(name = "Heavisine", mag = 4, family = "haar")
)

# Configurar Paralelização
n_cores = 4
plan(multisession, workers = n_cores)

# Funções Auxiliares (Copiadas de 03_simulation_01_inference_bootstrap.R)

run_raw_bootstrap = function(y_wide, X_mat, A_vec, time_grid, n_boot) {
  n_subj = nrow(y_wide)
  boot_curves = matrix(NA_real_, nrow = n_boot, ncol = ncol(y_wide))
  
  for(b in 1:n_boot) {
    idx = sample(1:n_subj, n_subj, replace = TRUE)
    if(sum(A_vec[idx] == 1) < 15 || sum(A_vec[idx] == 0) < 15) next
    
    fate = estimate_fate_pointwise(
      y_wide[idx,,drop=FALSE], X_mat[idx,,drop=FALSE], 
      A_vec[idx], time_grid
    )
    if(!is.null(fate)) boot_curves[b, ] = fate$fate_est
  }
  return(boot_curves)
}

build_band_hybrid = function(theta_hat, boot_curves, alpha = 0.05) {
  valid = complete.cases(boot_curves)
  if(sum(valid) < 50) 
    return(
      list(
        lower = rep(-Inf, length(theta_hat)), 
        upper = rep(Inf, length(theta_hat))
        )
      )
  
  boot_valid = boot_curves[valid, , drop = FALSE]
  deviations = apply(
    boot_valid, 1, 
    function(curve) max(abs(curve - theta_hat), na.rm = TRUE)
    )
  q = quantile(deviations, 1 - alpha, na.rm = TRUE)
  list(lower = theta_hat - q, upper = theta_hat + q)
}

# Worker Function
run_viz_scenario = function(scen) {
  # Gerar dados
  dgp_data = generate_single_scenario_dataset(
    scen$name, scen$mag, 
    N_SUBJ, T_PTS
  )
  
  true_curve = attr(dgp_data, "true_effect_curve")
  time_grid = attr(dgp_data, "time_grid")
  X_mat = attr(dgp_data, "X_mat")
  A_vec = attr(dgp_data, "A_vec")
  
  y_wide = as.matrix(dcast(dgp_data, id ~ t, value.var = "y")[, -1])
  
  # Estimativa (rLifting Hybrid)
  y_smooth = smooth_data_matrix(
    time_grid, y_wide, "rlifting_offline",
    scheme_name = scen$family, threshold_method = "semisoft", alpha = 0.3
  )
  
  if(is.null(y_smooth)) return(NULL)
  
  fate = estimate_fate_pointwise(y_smooth, X_mat, A_vec, time_grid)
  if(is.null(fate)) return(NULL)
  
  theta_hat = fate$fate_est
  
  # Bootstrap
  boot_curves = run_raw_bootstrap(y_wide, X_mat, A_vec, time_grid, N_BOOT)
  
  # Intervalo de Confiança Simultâneo
  band = build_band_hybrid(theta_hat, boot_curves, 0.05)
  
  # Retornar dados para plot
  data.table(
    Scenario = scen$name,
    t = time_grid,
    True = true_curve,
    Estimate = theta_hat,
    Lower = band$lower,
    Upper = band$upper
  )
}

# Execução
cat("Rodando simulações de visualização para 4 cenários (Parallel)...\n")
results_list = future_lapply(SCENARIOS, function(scen) {
  run_viz_scenario(scen)
}, future.seed = TRUE)

dt_plot = rbindlist(results_list)

cat("Gerando gráficos...\n")

# Configuração Visual
theme_paper = theme_bw() +
  theme(
    text = element_text(family = "serif", size = 12),
    axis.title = element_blank(), # Sem títulos de eixos para limpeza no LaTeX
    plot.title = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_list = list()
scenario_names = sapply(SCENARIOS, function(x) x$name)

for(scen in scenario_names) {
  d = dt_plot[Scenario == scen]
  
  p = ggplot(d, aes(x = t)) +
    # Faixa de Confiança
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "grey50", alpha = 0.4) +
    # Verdade
    geom_line(
      aes(y = True), linetype = "dashed", 
      color = "black", linewidth = 0.8
      ) +
    # Estimativa
    geom_line(aes(y = Estimate), color = "black", linewidth = 0.8) +
    theme_paper +
    theme(
      plot.margin = margin(5, 5, 5, 5)
    )
    # Nota: Sem labels ou títulos, controlado pelo LaTeX
  
  plot_list[[scen]] = p
}

# Combinar com Patchwork
final_plot = (plot_list[["Sigmoid"]] | plot_list[["Bumps"]]) /
             (plot_list[["Doppler"]] | plot_list[["Heavisine"]])

# Create output dir
base_dir = sprintf("output/sim1_N%d_T%d/artifacts", N_SUBJ, T_PTS)
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(base_dir, sprintf("plot_inference_examples_N%d_T%d.png", N_SUBJ, T_PTS)), 
  final_plot, 
  width = 8, height = 6, dpi = 300
  )

cat(sprintf("Gráfico salvo em %s\n", base_dir))
