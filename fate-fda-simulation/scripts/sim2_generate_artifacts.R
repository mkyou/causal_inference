library(data.table)
library(ggplot2)

# Configurações de caminhos
input_dir = "output/sim2_final"
output_dir = "output/sim2_final/plots"

out_dir_main = file.path(output_dir, "texto_principal")
out_dir_app  = file.path(output_dir, "apendice")

if (!dir.exists(out_dir_main)) dir.create(out_dir_main, recursive = TRUE)
if (!dir.exists(out_dir_app))  dir.create(out_dir_app, recursive = TRUE)

# Tema base para todos os graficos harmonizado com LaTeX
base_theme = theme_bw(base_size = 14) +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(margin = margin(5, 5, 5, 5), face = "bold")
  )

# -----------------------------------------------------------------------------
# 1. Gráficos de Boxplot do NRMSE Online por Cenário
# -----------------------------------------------------------------------------
cat("Gerando boxplots do NRMSE...\n")
df_raw = fread(file.path(input_dir, "sim2_raw.csv"))

for(scen_name in c("Bumps", "Doppler", "Heavisine", "Sigmoid")) {
    p_data = df_raw[Scenario == scen_name]
    if(nrow(p_data) == 0) next
    
    # Label para o título do arquivo (Escada se Heavisine)
    scen_label = ifelse(scen_name == "Heavisine", "Escada", scen_name)

    # Preparação dos fatores
    p_data[, N_label := factor(paste("N =", N), levels = paste("N =", sort(unique(N))))]
    p_data[, T_factor := factor(T_pts, levels = sort(unique(T_pts)))]
    p_data[, Method := "Ondaletas (Online)"]

    # a. NRMSE
    # Versão Texto Principal: Mediana e Quartis
    p_main = ggplot(p_data, aes(x = T_factor, y = NRMSE_online, color = Method)) +
        stat_summary(fun.data = "median_hilow", geom = "pointrange", size = 0.8) +
        facet_wrap(~ N_label, ncol = 3) +
        scale_color_manual(values = c("Ondaletas (Online)" = "royalblue")) +
        labs(x = "Pontos na grade (T)", y = "NRMSE online (Mediana e Quartis)") +
        base_theme
    
    ggsave(
        filename = file.path(out_dir_main, sprintf("nrmse_%s.png", scen_label)),
        plot = p_main, width = 9, height = 5, dpi = 300
    )

    # Versão Apêndice: Boxplot padrão
    p_app = ggplot(p_data, aes(x = T_factor, y = NRMSE_online, fill = Method)) +
        geom_boxplot(color = "black", alpha = 0.8) +
        facet_wrap(~ N_label, ncol = 3) +
        scale_fill_manual(values = c("Ondaletas (Online)" = "royalblue")) +
        labs(x = "Pontos na grade (T)", y = "NRMSE online (Boxplot)") +
        base_theme

    ggsave(
        filename = file.path(out_dir_app, sprintf("boxplot_nrmse_%s.png", scen_label)),
        plot = p_app, width = 9, height = 5, dpi = 300
    )
}

# -----------------------------------------------------------------------------
# 2. Gráficos do FATE Ponto a Ponto (Verdadeiro vs. Estimado Online)
# -----------------------------------------------------------------------------
cat("Gerando curvas do FATE...\n")
df_curves = fread(file.path(input_dir, "sim2_curves.csv"))

# Filtrar para um cenário ilustrativo de alta resolução
p_data_curves = df_curves[N == 2000 & T_pts == 256]

# Obter estimativa média empírica
p_mean = p_data_curves[, .(
  mean_online = mean(Est),
  theta_true = mean(TrueCurve)
), by = .(Scenario, t_val)]

for (scen_name in c("Bumps", "Doppler", "Heavisine", "Sigmoid")) {
  scen_data = p_mean[Scenario == scen_name]
  if(nrow(scen_data) == 0) next
  
  scen_label = ifelse(scen_name == "Heavisine", "Escada", scen_name)

  p = ggplot(scen_data, aes(x = t_val)) +
    geom_line(aes(y = theta_true, color = "Verdadeiro"), linewidth = 0.8) +
    geom_line(aes(y = mean_online, color = "Ondaletas (Online)"), linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(
        values = c("Verdadeiro" = "black", "Ondaletas (Online)" = "royalblue")
    ) +
    labs(
         x = "Tempo (t)",
         y = expression(theta(t)),
         color = NULL
    ) +
    base_theme

  ggsave(
      filename = file.path(out_dir_main, sprintf("fate_curves_%s.png", scen_label)),
      plot = p,
      width = 9,
      height = 5,
      dpi = 300
  )
}

# =============================================================================
# 6. TABELA GERENCIAL SUMMARY
# =============================================================================

cat("Gerando tabela gerencial summary para Simulação 2...\\n")

dt_summary = dt_raw[, .(
  NRMSE_Mean = mean(NRMSE_online, na.rm = TRUE),
  NRMSE_Median = median(NRMSE_online, na.rm = TRUE),
  Latency_Mean_us = mean(Latency_mean_us, na.rm = TRUE),
  Latency_Median_us = median(Latency_mean_us, na.rm = TRUE)
), by = .(Scenario, N_label, T_factor)]

setorder(dt_summary, Scenario, N_label, T_factor)
fwrite(dt_summary, "output/sim2_final/sim2_summary.csv")

cat("Artifact generation para a SIMULAÇÃO 2 concluída com sucesso.\\n")

# -----------------------------------------------------------------------------
# 3. Tabela Gerencial Summary
# -----------------------------------------------------------------------------
cat("Gerando tabela gerencial summary...\n")

# Calcular RMSE Global e Bias Global por Simulação a partir das curvas
agg_sim = df_curves[, .(
  RMSE = sqrt(mean((Est - TrueCurve)^2)),
  Bias = mean(Est - TrueCurve)
), by = .(SimID, Scenario, N, T_pts)]

# Adicionar RMSE e Bias ao df_raw para consolidação
df_raw[agg_sim, c("RMSE", "Bias") := .(i.RMSE, i.Bias), on = .(SimID, Scenario, N, T_pts)]

dt_summary = df_raw[, .(
  NRMSE_Mean = mean(NRMSE_online, na.rm = TRUE),
  NRMSE_Median = median(NRMSE_online, na.rm = TRUE),
  NRMSE_SD = sd(NRMSE_online, na.rm = TRUE),
  RMSE_Mean = mean(RMSE, na.rm = TRUE),
  Bias_Mean = mean(Bias, na.rm = TRUE),
  Latency_Mean_us = mean(Latency_mean_us, na.rm = TRUE)
), by = .(Scenario, N, T_pts)]

setorder(dt_summary, Scenario, N, T_pts)
fwrite(dt_summary, file.path(input_dir, "sim2_summary.csv"))

cat("Pronto! Artefatos salvos em", output_dir, "e", input_dir, "\n")
