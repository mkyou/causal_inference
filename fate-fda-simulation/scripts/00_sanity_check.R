library(ggplot2)
library(data.table)
library(patchwork)

source("R/utils_dgp_functions.R")
source("R/generate_sequential_dgp.R")
source("R/generate_functional_outcomes.R")
source("R/process_functional_data.R")

# 1. Configuração do Cenário
N = 4 
time_grid = seq(0, 1, by = 0.001) 
t_treat = 0.5

func_doppler = function(t, s) effect_doppler(t, s, magnitude = 5)

# 2. Gerar Covariáveis
dgp_covs = generate_sequential_dgp(
  n_subjects = N,
  time_grid = time_grid,
  treatment_times = c(t_treat),
  baseline_func_templates = list(function(t, p) 0),
  baseline_params_list = list(matrix(0, N, 1)),
  beta_funcs_matrix = matrix(list(), 1, 1),
  gamma_matrix = matrix(list(function(t,s) 0), 1, 1),
  theta_matrix = matrix(100, 1, 2), # Força A=1
  noise_sd_vec = c(0.1)
)

# 3. Gerar Desfecho (Agora deve funcionar pois t=0.5 existe)
dgp_out = generate_functional_outcome(
  dgp_data = dgp_covs,
  baseline_y_func = function(t) 0,
  direct_effect_funcs = list(func_doppler),
  concurrent_effect_funcs = list(function(t) 0),
  treatment_times = c(t_treat),
  noise_sd = 0.5 
)

# 4. Suavização
y_wide = dcast(dgp_out, id ~ t, value.var = "y")
y_matrix = as.matrix(y_wide[, -1]) # Remove coluna id

# A. Splines (FDA)
res_spline = smooth_data_matrix(time_grid, y_matrix, method = "spline")

# B. rLifting
res_rlift = smooth_data_matrix(
  time_grid, 
  y_matrix, 
  method = "rlifting_offline",
  scheme_name = 'cdf97'
  )

# C. Wavethresh
res_wave = smooth_data_matrix(time_grid, y_matrix, method = "wavethresh")

# 5. Visualização
id_plot = 1
# Usamos a_y como referência (é o contrafactual sob tratamento, com ruído)
# Na simulação forçamos A=1, então y e a_y devem ser idênticos neste caso.
true_curve_noisy = dgp_out[id == id_plot, a_y] 

df_plot = data.table(
  Tempo = time_grid,
  Ruidoso = y_matrix[id_plot, ],
  Spline = res_spline[id_plot, ],
  rLifting = res_rlift[id_plot, ],
  Wavethresh = res_wave[id_plot, ]
)

df_long = melt(
  df_plot, 
  id.vars = "Tempo", 
  variable.name = "Tipo", 
  value.name = "Valor"
  )

cores = c(
  "Ruidoso" = "grey80", 
  "Spline" = "firebrick", 
  "rLifting" = "royalblue", 
  "Wavethresh" = "forestgreen"
)

df_long |> ggplot(aes(x = Tempo, y = Valor, color = Tipo)) +
  geom_line(aes(alpha = Tipo, linewidth = Tipo)) +
  scale_color_manual(values = cores) +
  scale_alpha_manual(
    values = c(
      "Ruidoso" = 0.5, "Spline" = 1,
      "rLifting" = 1, "Wavethresh" = 1
      )
    ) +
  scale_linewidth_manual(
    values = c(
      "Ruidoso" = 0.5, "Spline" = 0.8, 
      "rLifting" = 0.8, "Wavethresh" = 0.8
      )
    ) +
  geom_vline(xintercept = t_treat, linetype = "dotted", color = "black") +
  annotate(
    "text", x = t_treat + 0.02, 
    y = max(df_plot$Ruidoso), label = "Tratamento (Doppler)", 
    hjust = 0
    ) +
  labs(
    title = "Sanity Check: Recuperação do Efeito Doppler",
    subtitle = "Comparação de métodos de suavização sobre Y(t)",
    y = "Desfecho Y(t)"
  ) +
  theme_bw() +
  theme(legend.position = "top")