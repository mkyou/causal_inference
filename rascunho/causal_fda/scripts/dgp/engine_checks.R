library(data.table)
library(ggplot2)
library(patchwork)

COR_PRIMARIA = "navy"
COR_SECUNDARIA = "palevioletred"

# Funções de verificação

#' @title Verifica Parâmetros de Sanidade do DGP
#' @description Compara os parâmetros esperados da simulação 
#' com os dados realmente gerados.
#'
#' @param dgp_data (list) A saída de `generate_fdg_causal`.
#' @param params_esperados (list) Uma lista nomeada com os 
#'   parâmetros esperados (N, T, k, t_A_time, noise_sd, etc.)
#'
#' @return (invisible) NULL. Apenas imprime no console.
#' 
check_sanity = function(dgp_data, params_esperados) {
  
  # N (Número de indivíduos)
  N_real = nrow(dgp_data$static_data)
  message(sprintf(
    "-> N (Amostra): Esperado=%d, Real=%d",
    params_esperados$N, N_real
  ))
  if (N_real != params_esperados$N) {
    warning("N (Amostra) diverge do esperado!")
  }
  
  # T (Pontos de Tempo)
  T_real = length(unique(dgp_data$curve_data$t))
  message(sprintf(
    "-> T (Pontos): Esperado=%d, Real=%d",
    params_esperados$T, T_real
  ))
  if (T_real != params_esperados$T) {
    warning("T (Pontos) diverge do esperado!")
  }
  
  # k (Número de covariáveis)
  X_cols = grep("^X[0-9]+", names(dgp_data$static_data), 
                value = TRUE)
  k_real = length(X_cols)
  message(sprintf(
    "-> k (Covariáveis): Esperado=%d, Real=%d",
    params_esperados$k, k_real
  ))
  if (k_real != params_esperados$k) {
    warning("k (Covariáveis) diverge do esperado!")
  }
  
  # Outros parâmetros (apenas informativo)
  message(sprintf(
    "-> t_A (Tempo Trat.): %.2f", params_esperados$t_A_time
  ))
  message(sprintf(
    "-> Noise SD (Ruído): %.2f", params_esperados$noise_sd
  ))
  
  invisible(NULL)
}


#' @title Verifica a Distribuição dos Escores de Propensão
#' @description Estima P(A|X) e plota sua distribuição 
#' para avaliar o 'overlap' (sobreposição).
#'
#' @param dgp_data (list) A saída de `generate_fdg_causal`.
#'
#' @return (ggplot) Um objeto ggplot.
#' 
check_propensity_scores = function(dgp_data) {
  
  dt_analysis = copy(dgp_data$static_data)
  X_cols = grep("^X[0-9]+", names(dt_analysis), value = TRUE)
  
  if (length(X_cols) > 0) {
    # Se há covariáveis, estima o escore de propensão
    model_formula = as.formula(paste(
      "A ~", paste(X_cols, collapse = " + ")
    ))
    
    # Ajusta o modelo logístico
    prop_model = glm(
      model_formula,
      data = dt_analysis,
      family = "binomial"
    )
    
    dt_analysis[, p_score := predict(prop_model, type = "response")]
    plot_title = "Distribuição do Escore de Propensão (P(A|X))"
    
  } else {
    # Se não há covariáveis, a propensão é apenas a média
    dt_analysis[, p_score := mean(A)]
    plot_title = "Propensão Média (Sem Covariáveis)"
  }
  
  # Plotar
  p = ggplot(dt_analysis, aes(x = p_score, fill = factor(A))) +
    geom_density(linewidth = .5, alpha = .5) +
    scale_fill_manual(
      values = c("0" = COR_SECUNDARIA, "1" = COR_PRIMARIA),
      name = "Tratamento (A)"
    ) +
    labs(
      title = plot_title,
      subtitle = "Usado para verificar sobreposição (overlap)",
      x = "Escore de Propensão P(A=1 | X)",
      y = "Densidade"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(p)
}


#' @title Verifica o ATE (Efeito Médio do Tratamento)
#' @description Plota o ATE verdadeiro com ATE ingênuo 
#' (gerado a partid dos dados observados, que contém viés e ruído) 
#' E imprime estatísticas descritivas.
#'
#' @param dgp_data (list) A saída de `generate_fdg_causal`.
#' @param t_A_time (Numeric) O tempo de aplicação do tratamento.
#'
#' @return (ggplot) Um objeto ggplot combinado (patchwork).
#' 
check_ate = function(dgp_data, t_A_time) {
  
  # Cálculo de ATE verdadeiro (sem ruído) e ATE com ruído
  #    Usamos as curvas médias populacionais
  dt_true = dgp_data$curve_data[, .(
    ATE_true = mean(true_Y_cf_A1) - mean(true_Y_cf_A0)
  ), by = t]
  
  # Cálculo do ATE ingênuo
  dt_obs = dgp_data$curve_data[, .(id, t, Y_obs)] |>
    merge(dgp_data$static_data[, .(id, A)], by = "id")
  
  # Média separada por grupo de tratamento
  dt_naive_agg = dt_obs[, .(Y_mean = mean(Y_obs)), by = .(t, A)]
  
  # Pivotar para calcular a diferença
  dt_naive = dcast(dt_naive_agg, t ~ A,
                   value.var = "Y_mean", fill = NA)
  
  # Garantir que ambas as colunas '0' e '1' existam
  if (!"0" %in% names(dt_naive)) dt_naive[, `0` := NA_real_]
  if (!"1" %in% names(dt_naive)) dt_naive[, `1` := NA_real_]
  
  # ATE Ingênuo = E[Y|A=1] - E[Y|A=0]
  dt_naive[, ATE_naive := `1` - `0`]
  
  # Combinar e calcular a diferença (Viés + Ruído)
  dt_ate = dt_true[dt_naive[, .(t, ATE_naive)], on = "t"]
  dt_ate[, ATE_diff := ATE_naive - ATE_true]
  
  
  # Estatísticas Descritivas (Antes vs. Depois de t_A)
  message("--- Descritivas do ATE (Efeito Médio) ---")
  
  # Antes de t_A
  message(sprintf("Período Pré-Tratamento (t < %.2f):", t_A_time))
  desc_pre = dt_ate[t < t_A_time,
                    .(ATE_True_Pre = ATE_true,
                      ATE_Naive_Pre = ATE_naive)] |>
    summary()
  print(desc_pre)
  
  # Depois de t_A
  message(sprintf("Período Pós-Tratamento (t >= %.2f):", t_A_time))
  desc_post = dt_ate[t >= t_A_time,
                     .(ATE_True_Post = ATE_true,
                       ATE_Naive_Post = ATE_naive)] |>
    summary()
  print(desc_post)
  
  # Preparando dados para plot
  dt_plot_long = melt(
    dt_ate,
    id.vars = "t",
    measure.vars = c("ATE_true", "ATE_naive"),
    variable.name = "Tipo",
    value.name = "ATE"
  )
  
  # Plot 1: ATE Verdadeiro vs. Ingênuo
  p_ate = ggplot(
    dt_plot_long, aes(x = t, y = ATE, linetype = Tipo)
    ) +
    geom_line(linewidth = .8) +
    geom_vline(
      xintercept = t_A_time,
      linetype = "dashed",
      color = "black") +
    scale_linetype_manual(
      values = c(
        "ATE_true" = "solid",
        "ATE_naive" = "dotted"
      ),
      labels = c(
        "ATE Verdadeiro (Ground Truth)",
        "ATE Ingênuo (Observado E[Y|1]-E[Y|0])"
      )
    ) +
    labs(
      title = "Verificação do ATE: Verdadeiro vs. Ingênuo",
      subtitle = "ATE Ingênuo = Viés de Confundimento + Ruído",
      y = "ATE(t)",
      x = "Tempo (t)"
    ) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  # Plot 2: Diferença (Viés + Ruído)
  p_diff = ggplot(dt_ate, aes(x = t, y = ATE_diff)) +
    geom_line(color = "gray40", linewidth = 0.8) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    labs(
      title = "Diferença (ATE Ingênuo - ATE Verdadeiro)",
      subtitle = "Isto representa o (Viés de Confundimento + Ruído Amostral)",
      y = "Viés + Ruído",
      x = "Tempo (t)"
    ) +
    theme_bw()
  
  # Combinar gráficos
  return(p_ate / p_diff)
}