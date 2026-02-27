library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(kableExtra)
library(patchwork)

# Configurações visuais
theme_paper = theme_bw() +
  theme(
    legend.position = "bottom",
    text = element_text(family = "serif", size = 12),
    axis.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", colour = "black"),
    strip.text = element_text(face = "bold")
  )

method_colors = c(
  "spline" = "firebrick", 
  "wavethresh" = "forestgreen", 
  "rlifting_offline" = "royalblue"
)
method_labels = c(
  "spline" = "Spline (FDA)", 
  "wavethresh" = "Wavethresh (SURE)", 
  "rlifting_offline" = "rLifting"
)

# Configuração do tamanho da simulação
N_SUBJ = 2000
T_PTS = 64

# Configuração de diretórios de artefatos
base_dir = sprintf("output/sim1_N%d_T%d", N_SUBJ, T_PTS)
res_file = file.path(
  base_dir, sprintf("sim1_time_series_results_N%d_T%d.csv", N_SUBJ, T_PTS)
  )

if(!file.exists(res_file)) {
  rds_file = file.path(
    base_dir, sprintf("sim1_results_N%d_T%d.rds", N_SUBJ, T_PTS)
    )
  if(file.exists(rds_file)) {
     stop(sprintf("Arquivo CSV não encontrado em: %s", res_file))
  } else {
     stop(sprintf("Nenhum resultado encontrado em: %s", base_dir))
  }
}
res = fread(res_file)

artifacts_dir = file.path(base_dir, "artifacts")
dir.create(artifacts_dir, recursive=TRUE, showWarnings=FALSE)

# Tabela global
# Agregação: média das métricas ao longo do tempo t
table_global = res |>
  group_by(Scenario, Method) |>
  summarise(
    # Métricas de erro (média no tempo)
    Global_Max_Err = max(MAE_t, na.rm=TRUE), 
    
    Global_NRMSE = mean(NRMSE_t, na.rm=TRUE),
    Global_NRMSE_sd = mean(NRMSE_sd, na.rm=TRUE),
    
    Global_MAE = mean(MAE_t, na.rm=TRUE),
    Global_MAE_sd = mean(MAE_sd, na.rm=TRUE),
    
    # Viés global (média da diferença absoluta entre Estimate e True)
    Global_Bias = mean(Mean_Estimate - True_Effect, na.rm=TRUE),
    Global_Bias_sd = mean(RMSE_sd, na.rm=TRUE),
    
    # Tempo
    Exec_Time = mean(Time_Mean, na.rm=TRUE),
    Exec_Time_sd = mean(Time_SD, na.rm=TRUE),
    
    .groups = "drop"
  ) |>
  mutate(
    # Formatação "Média (SD)"
    NRMSE_Fmt = sprintf("%.3f (%.3f)", Global_NRMSE, Global_NRMSE_sd),
    MAE_Fmt = sprintf("%.3f (%.3f)", Global_MAE, Global_MAE_sd),
    Bias_Fmt = sprintf("%.3f (%.3f)", Global_Bias, Global_Bias_sd),
    MaxErr_Fmt = sprintf("%.3f", Global_Max_Err),
    Time_Fmt = sprintf("%.4f s (%.4f s)", Exec_Time, Exec_Time_sd)
  ) |>
  dplyr::select(
    Scenario, Method, MaxErr_Fmt, NRMSE_Fmt, 
    MAE_Fmt, Bias_Fmt, Time_Fmt
  ) |>
  arrange(Scenario, Method)

# Salvar CSV bruto da tabela
fwrite(
  table_global,
  file.path(
    artifacts_dir, 
    sprintf("table_global_metrics_N%d_T%d.csv", N_SUBJ, T_PTS))
  )

# Gerar HTML/LaTeX da tabela para visualização
kbl_table = kable(
  table_global, format = "html", booktabs = TRUE,
  col.names = c(
    "Cenário", "Método", "Max Err", "NRMSE (SD)", 
    "MAE (SD)", "Bias (SD)", "Tempo"
    )
  ) |>
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))

# Salvar tabela renderizada
save_kable(
  kbl_table, 
  file.path(
    artifacts_dir, 
    sprintf("table_global_pretty_N%d_T%d.html", N_SUBJ, T_PTS)
    )
  )

print(table_global)


# Gráfico 1: NRMSE ao longo do tempo
p_nrmse = ggplot(res, aes(x = t, y = NRMSE_t, color = Method)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  facet_wrap(~Scenario, scales = "free_y") +
  scale_color_manual(values = method_colors, labels = method_labels) +
  labs(x = "Time", y = "NRMSE", color = NULL) +
  theme_paper

ggsave(
  file.path(
    artifacts_dir, sprintf("plot_nrmse_time_N%d_T%d.png", N_SUBJ, T_PTS)
    ), 
  p_nrmse, width = 10, height = 6, dpi = 300
  )


# Gráfico 2: viés visual (estimativa vs. verdade)
plot_list = list()
scenarios_vec = unique(res$Scenario)

for(scen in scenarios_vec) {
  data_scen = res[Scenario == scen]
  
  # Cria camada da "verdade" (apenas uma linha preta tracejada)
  # Pegamos a verdade do primeiro método (é igual para todos)
  truth_line = data_scen[
    Method == method_labels[1] | Method == unique(Method)[1]
    ]
  
  p = ggplot(data_scen, aes(x = t)) +
    # Linha da verdade (fundo)
    geom_line(
      aes(y = True_Effect), 
      color = "black", 
      linetype = "dashed", 
      linewidth = 0.8
      ) +
    # Linhas dos métodos
    geom_line(
      aes(y = Mean_Estimate, color = Method), 
      linewidth = 0.7, 
      alpha = 0.9
      ) +
    scale_color_manual(values = method_colors, labels = method_labels) +
    labs(subtitle = scen, x = "Time", y = "Effect", color = NULL) +
    theme_paper +
    theme(plot.subtitle = element_text(hjust = 0.5, face = "bold"))
  
  plot_list[[scen]] = p
}

# Combina os plots (2x2)
p_bias_combined = wrap_plots(plot_list, ncol = 2) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

ggsave(
  file.path(
    artifacts_dir, sprintf("plot_bias_visual_N%d_T%d.png", N_SUBJ, T_PTS)
    ), 
  p_bias_combined, 
  width = 10, 
  height = 7, 
  dpi = 300
  )
