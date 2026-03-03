# scripts/sim1_generate_artifacts.R

library(data.table)
library(ggplot2)

# =============================================================================
# 1. SETUP DE DIRETÓRIOS E CORES
# =============================================================================

out_dir = "output/sim1_final/plots"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

method_colors = c(
  "Ondaletas (Online)"  = "royalblue",
  "Ondaletas (Offline)" = "forestgreen",
  "Splines de Ecker"    = "firebrick"
)

method_labels_mapping = c(
  "wavelet_causal"  = "Ondaletas (Online)",
  "wavelet_offline" = "Ondaletas (Offline)",
  "spline_ecker"    = "Splines de Ecker"
)

scenarios = c("Bumps", "Doppler", "Heavisine", "Sigmoid")

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

# =============================================================================
# 2. CARREGAMENTO DOS DADOS (SIM 1 FINAL)
# =============================================================================

cat("Carregando resultados da Simulação 1...\\n")

f_raw = "output/sim1_final/sim1_raw.csv"
f_curves = "output/sim1_final/sim1_curves.csv"

if (!file.exists(f_raw) || !file.exists(f_curves)) {
  stop("Os arquivos finais não foram encontrados. Garanta que sim_1.R terminou sua execução.")
}

dt_raw = fread(f_raw)
dt_curves = fread(f_curves)

# Pre-processamento (Fatores e Labels)
dt_raw[, Method := factor(method_labels_mapping[Method], levels = names(method_colors))]
dt_raw[, Scenario := factor(Scenario, levels = scenarios, labels = c("Bumps", "Doppler", "Escada", "Sigmoid"))]
dt_raw[, N_fac := factor(N, levels = sort(unique(N)))]
dt_raw[, T_fac := factor(T_pts, levels = sort(unique(T_pts)))]

dt_curves[, Method := factor(method_labels_mapping[Method], levels = names(method_colors))]
dt_curves[, Scenario := factor(Scenario, levels = scenarios, labels = c("Bumps", "Doppler", "Escada", "Sigmoid"))]

# =============================================================================
# 3. TABELA GERENCIAL
# =============================================================================

cat("Gerando tabela gerencial summary...\\n")

# Calcular RMSE Global por Simulação a partir das curvas
rmse_dt = dt_curves[, .(
  RMSE = sqrt(mean((Est - TrueCurve)^2))
), by = .(SimID, Scenario, Method, N, T_pts)]

# Adicionar RMSE ao dt_raw para consolidação
dt_raw[rmse_dt, RMSE := i.RMSE, on = .(SimID, Scenario, Method, N, T_pts)]

dt_summary = dt_raw[, .(
  NRMSE_Mean = mean(NRMSE, na.rm = TRUE),
  NRMSE_Median = median(NRMSE, na.rm = TRUE),
  NRMSE_SD = sd(NRMSE, na.rm = TRUE),
  RMSE_Mean = mean(RMSE, na.rm = TRUE),
  Bias_Mean = mean(Bias, na.rm = TRUE),
  Width_Mean = mean(Width, na.rm = TRUE),
  Width_SD = sd(Width, na.rm = TRUE),
  Cov_Simul_Mean = mean(Cov_Simul, na.rm = TRUE)
), by = .(Scenario, N_fac, T_fac, Method)]

setorder(dt_summary, Scenario, N_fac, T_fac, Method)
fwrite(dt_summary, "output/sim1_final/sim1_summary.csv")

# =============================================================================
# 4. BOXPLOTS AGREGADOS (NRMSE, LARGURA e COBERTURA)
# =============================================================================

cat("Gerando gráficos de caixas agregados...\\n")

for (scen in levels(dt_raw$Scenario)) {
  scen_raw = dt_raw[Scenario == scen]
  
  # a. Boxplot: NRMSE
  p_nrmse = ggplot(scen_raw, aes(x = T_fac, y = NRMSE, fill = Method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    facet_wrap(~ N_fac, labeller = label_both) +
    scale_fill_manual(values = method_colors) +
    labs(x = "Pontos na grade (T)", y = "NRMSE global") +
    base_theme
  
  ggsave(sprintf("%s/boxplot_nrmse_%s.png", out_dir, scen), p_nrmse, width = 9, height = 5, bg = "white", dpi = 300)
  
  # b. Boxplot: Largura das Bandas (Width)
  p_width = ggplot(scen_raw, aes(x = T_fac, y = Width, fill = Method)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    facet_wrap(~ N_fac, labeller = label_both) +
    scale_fill_manual(values = method_colors) +
    labs(x = "Pontos na grade (T)", y = "Largura média da banda") +
    base_theme
  
  ggsave(sprintf("%s/boxplot_width_%s.png", out_dir, scen), p_width, width = 9, height = 5, bg = "white", dpi = 300)
  

  # c. Ponto/Barra: Cobertura Simultânea Média (Cov_Simul) - teórica 95%
  p_cov = ggplot(dt_summary[Scenario == scen], aes(x = T_fac, y = Cov_Simul_Mean, color = Method, group = Method)) +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    facet_wrap(~ N_fac, labeller = label_both) +
    scale_color_manual(values = method_colors) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = "Pontos na grade (T)", y = "Taxa empírica de cobertura simultânea") +
    base_theme
  
  ggsave(sprintf("%s/lines_covsimul_%s.png", out_dir, scen), p_cov, width = 9, height = 5, bg = "white", dpi = 300)
}

# =============================================================================
# 5. GRÁFICOS PONTO-A-PONTO (MÉDIA)
# =============================================================================

cat("Gerando gráficos ponto-a-ponto para Pior e Melhor caso...\\n")

best_case = list(N = 2000, T_pts = 512)
worst_case = list(N = 200, T_pts = 128)

generate_pointwise_aggregates = function(case_label, condition) {
  dt_sub = dt_curves[condition]
  
  if (nrow(dt_sub) == 0) return()
  
  dt_p = dt_sub[, .(
    Mean_NRMSE = sqrt(mean((Est - TrueCurve)^2)) / mean(diff(range(TrueCurve))),
    Mean_Width = mean(Upper - Lower),
    Mean_Cov = mean(TrueCurve >= Lower & TrueCurve <= Upper)
  ), by = .(Scenario, Method, t_val)]
  
  # a. Pontual RMSE over time
  p1 = ggplot(dt_p, aes(x = t_val, y = Mean_NRMSE, color = Method)) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ Scenario, scales = "free_y") +
    scale_color_manual(values = method_colors) +
    labs(x = "Tempo (t)", y = "NRMSE pontual") +
    base_theme
  ggsave(sprintf("%s/pointwise_nrmse_%s.png", out_dir, case_label), p1, width = 10, height = 6, bg = "white", dpi = 300)
}

generate_pointwise_aggregates("worst_case", dt_curves$N == worst_case$N & dt_curves$T_pts == worst_case$T_pts)
generate_pointwise_aggregates("best_case", dt_curves$N == best_case$N & dt_curves$T_pts == best_case$T_pts)

# =============================================================================
# 6. GRÁFICOS DE EXEMPLO (True vs. Est + Bands)
# =============================================================================

cat("Gerando gráficos de uma simulação exemplar (DGP fixo) com bandas de confiança...\\n")

generate_exemplary_plot = function(case_label, condition) {
  # Fixando o primeiro dataset gerado SimID == 1 para visualização
  dt_sub = dt_curves[condition & SimID == 1]
  
  if (nrow(dt_sub) == 0) return()
  
  p = ggplot(dt_sub, aes(x = t_val)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Method), alpha = 0.25) +
    geom_line(aes(y = TrueCurve, color = "Verdadeiro"), linewidth = 0.7, linetype = "solid") +
    geom_line(aes(y = Est, color = Method), linewidth = 0.8, linetype = "dashed") +
    facet_grid(Scenario ~ Method, scales = "free_y") +
    scale_color_manual(
      values = c(method_colors, "Verdadeiro" = "black")
    ) +
    scale_fill_manual(values = method_colors) +
    labs(x = "Tempo (t)", y = "Efeito causal FATE(t)") +
    base_theme +
    theme(legend.position = "none") # Já dividido por facetas, n precisa legenda por painel
  
  ggsave(sprintf("%s/exemplary_run_%s.png", out_dir, case_label), p, width = 12, height = 9, bg = "white", dpi = 300)
}

generate_exemplary_plot("worst_case", dt_curves$N == worst_case$N & dt_curves$T_pts == worst_case$T_pts)
generate_exemplary_plot("best_case", dt_curves$N == best_case$N & dt_curves$T_pts == best_case$T_pts)

cat("Artifact generation para a SIMULAÇÃO 1 concluída com sucesso.\\n")
