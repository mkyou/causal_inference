library(data.table)
library(ggplot2)

# Configurações de caminhos
input_dir = "output/sim2_final"
output_dir = "output/sim2_final/plots"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# 1. Gráficos de Boxplot do NRMSE Online por Cenário
# -----------------------------------------------------------------------------
cat("Gerando boxplots do NRMSE...\n")
df_raw = fread(file.path(input_dir, "sim2_raw.csv"))

scenarios = unique(df_raw$Scenario)

for(scen in scenarios) {
    p_data = df_raw[Scenario == scen]

    # Preparação dos fatores e labels
    p_data[, N_factor := as.factor(N)]
    p_data[, T_label := factor(paste("T =", T_pts), levels = c("T = 64", "T = 128", "T = 256", "T = 512"))]

    p = ggplot(p_data, aes(x = N_factor, y = NRMSE_online)) +
        geom_boxplot(fill = "royalblue", color = "black", alpha = 0.8) +
        facet_wrap(~ T_label, ncol = 3) +
        theme_bw(base_family = "serif", base_size = 12) +
        labs(
            x = "Tamanho da amostra",
            y = "NRMSE online"
        ) +
        theme(
            legend.position = "none",
            strip.background = element_rect(fill = "grey95")
        )

    ggsave(
        filename = file.path(output_dir, sprintf("boxplot_nrmse_%s.png", scen)),
        plot = p,
        width = 6.5,
        height = 3.5,
        dpi = 300
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

for (scen in scenarios) {
  scen_data = p_mean[Scenario == scen]

  p = ggplot(scen_data, aes(x = t_val)) +
    geom_line(aes(y = theta_true, color = "FATE Verdadeiro"), linewidth = 0.8) +
    geom_line(aes(y = mean_online, color = "FATE Estimado (Online)"), linewidth = 0.8, linetype = "dashed") +
    scale_color_manual(
        values = c("FATE Verdadeiro" = "black", "FATE Estimado (Online)" = "royalblue")
    ) +
    theme_bw(base_family = "serif", base_size = 12) +
    labs(
         x = "Tempo",
         y = expression(theta(t)),
         color = NULL
    ) +
    theme(
        legend.position = "bottom"
    )

  ggsave(
      filename = file.path(output_dir, sprintf("fate_curves_%s.png", scen)),
      plot = p,
      width = 5.5,
      height = 3.5,
      dpi = 300
  )
}

cat("Pronto! Artefatos salvos em", output_dir, "\n")
