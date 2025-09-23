# ===================================================================
# SCRIPT MESTRE PARA EXECUÇÃO DE EXPERIMENTOS DE SIMULAÇÃO
# ===================================================================

# --- SETUP ---
# Limpa o ambiente e os gráficos para uma execução limpa
rm(list = ls()); graphics.off()

# Carrega todas as bibliotecas e scripts necessários
source("functional_simulations/data_examples/confouding_feedback.R")

source("functional_simulations/models/ecker_luna.R")
source("functional_simulations/models/g_formula.R")
source("functional_simulations/models/marginal_structural_model.R")

source("functional_simulations/utils/assess_estimation.R")
source("functional_simulations/utils/run_model_evaluation.R")

# ----- EXECUÇÃO DE MODELOS ------------------------------------------------


# ----- G-FORMULA ESPARSA -------------------------------------
results_msm = run_model_evaluation(
  all_scenarios_data = all_combined_scenarios,
  scenario_configs = scenario_configs_combined,
  model_function = run_fmsm,
  assessment_function = assess_full_msm_process,
  model_name = "MSM Funcional"
)

# --- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS (MSM FUNCIONAL) ---
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS (MSM FUNCIONAL) ---\n")
purrr::iwalk(
  results_msm, # Certifique-se de que a variável é a correta
  function(assessment, scenario_name) {
    # --- Cabeçalho para o Cenário ---
    cat("\n=====================================================\n")
    cat(">>> Análise do Cenário:", scenario_name, "\n")
    cat("=====================================================\n\n")
    
    # --- 1. Diagnósticos dos Modelos de Propensão ---
    cat("--- 1. Avaliação dos Modelos de Propensão ---\n")
    if (!is.null(assessment$propensity_model_assessment)) {
      # Itera sobre o diagnóstico de cada modelo de propensão (um por tempo de tratamento)
      purrr::iwalk(assessment$propensity_model_assessment, function(prop_model, time) {
        cat(paste("\n-- Diagnóstico de Propensões em t =", time, "--\n"))
        cat("Resumo das Probabilidades Previstas (Positividade):\n")
        print(prop_model$summary)
        print(prop_model$plot) # Gráfico de calibração
      })
    } else {
      cat("Avaliação do modelo de propensão não encontrada.\n")
    }
    cat("\n-----------------------------------------------------\n")
    
    # --- 2. Diagnósticos dos Pesos ---
    cat("\n--- 2. Avaliação dos Pesos (SW) ---\n")
    if (!is.null(assessment$weight_assessment)) {
      cat("Resumo da Distribuição dos Pesos:\n")
      print(assessment$weight_assessment$summary)
      print(assessment$weight_assessment$plot) # Histograma dos pesos
    } else {
      cat("Avaliação dos pesos não encontrada.\n")
    }
    cat("\n-----------------------------------------------------\n")
    
    # --- 3. Ajuste do Modelo de Desfecho ---
    cat("\n--- 3. Avaliação do Ajuste Factual do Modelo de Desfecho ---\n")
    if (!is.null(assessment$outcome_model_factual_fit)) {
      print(assessment$outcome_model_factual_fit)
    } else {
      cat("Avaliação do modelo de desfecho não encontrada.\n")
    }
    cat("\n-----------------------------------------------------\n")
    
    # --- 4. Estimação do Efeito Causal ---
    cat("\n--- 4. Avaliação da Estimação do Efeito Causal ---\n")
    if (!is.null(assessment$causal_effect_assessment)) {
      cat("Métricas Resumidas (Geral):\n")
      print(assessment$causal_effect_assessment$summary)
      
      cat("\nMétricas Funcionais (Viés e RMSE ao longo do tempo):\n")
      print(assessment$causal_effect_assessment$functional)
      
      cat("\nGráfico do Efeito Causal Médio:\n")
      print(assessment$causal_effect_assessment$plot)
    } else {
      cat("Avaliação do efeito causal não encontrada.\n")
    }
    
    # Pausa para permitir a visualização antes de passar para o próximo
    readline(prompt = "\nPressione [enter] para ver o próximo cenário...")
  }
)

# ----- G-FORMULA -------------------------------------
cat("\n\n--- AVALIANDO DADOS (G-Formula) ---\n")
results_gformula = run_model_evaluation(
  all_scenarios_data = all_combined_scenarios,
  scenario_configs = scenario_configs_combined,
  model_function = run_gformula_full,
  assessment_function = assess_full_gformula_process,
  model_name = "G-Formula"
)

# --- Tabela de Resumo para o EFEITO CAUSAL ---
causal_summary_table = purrr::map_dfr(
  results_gformula,
  # O caminho correto para as métricas de resumo do efeito causal
  ~ .x$causal_effect_assessment$summary,
  .id = "Scenario"
)

cat("\n\n--- TABELA DE RESUMO: EFEITO CAUSAL ---\n")
print(causal_summary_table)

outcome_model_summary_table = purrr::map_dfr(
  results_gformula,
  ~ .x$outcome_model_factual_fit,
  .id = "Scenario"
)

cat("\n\n--- TABELA DE RESUMO: MODELO DESFECHO ---\n")
print(outcome_model_summary_table)

potential_outcome_summary_table = purrr::map_dfr(
  results_gformula,
  ~ .x$outcome_assessment_metric,
  .id = "Scenario"
)

cat("\n\n--- TABELA DE RESUMO: RESULTADOS POTENCIAIS ---\n")
print(potential_outcome_summary_table)


# --- Tabela de Resumo para as COVARIÁVEIS ---
covariate_summary_table = purrr::map_dfr(
  results_gformula,
  function(scenario_results) {
    purrr::map_dfr(
      scenario_results$covariate_assessment$metrics,
      function(metrics_per_cov) {
        dplyr::bind_rows(
          "Regime A" = metrics_per_cov$regime_a,
          "Regime B" = metrics_per_cov$regime_b,
          .id = "Regime"
        )
      }, .id = "Covariate"
    )
  },
  .id = "Scenario"
)

cat("\n\n--- TABELA DE RESUMO: SIMULAÇÃO DE COVARIÁVEIS ---\n")
print(covariate_summary_table)

results_gformula$CF2_LowCross$covariate_assessment$plots

# Apresentação interativa dos gráficos
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS (G-FORMULA) ---\n")
purrr::iwalk(
  results_gformula,
  function(assessment, scenario_name) {
    cat("\n=====================================================\n")
    cat(">>> Análise do Cenário:", scenario_name, "\n")
    cat("=====================================================\n\n")
    
    print(assessment$outcome_assessment_metrics)
    print(assessment$causal_effect_assessment$summary)
    print(assessment$causal_effect_assessment$plot)
    
    
    # Pausa para permitir a visualização antes de passar para o próximo
    readline(prompt = "Pressione [enter] para ver o próximo cenário...")
  }
)

# ----- G-FORMULA PARA COVARIÁVEIS-------------------------------------
cat("\n\n--- AVALIANDO COVARIAVEIS (G-Formula) ---\n")
results_gformula_covariates = run_model_evaluation(
  all_scenarios_data = all_feedback_scenarios,
  scenario_configs = scenario_configs_feedback,
  model_function = run_gformula_covariates_for_assessment,
  assessment_function = assess_covariate_simulation,
  model_name = "G-Formula (Etapa de Covariáveis)"
)

# Apresentação interativa dos gráficos
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS (COVS G-FORMULA) ---\n")
purrr::iwalk(
  results_gformula_covariates,
  function(assessment, scenario_name) {
    cat("\n=====================================================\n")
    cat(">>> Análise do Cenário:", scenario_name, "\n")
    cat("=====================================================\n\n")
    
    print(assessment$plot)
    
    # Pausa para permitir a visualização antes de passar para o próximo
    readline(prompt = "Pressione [enter] para ver o próximo cenário...")
  }
)



# ----- ECKER E LUNA -----------------------------------------
cat("\n\n--- AVALIANDO O EFEITO FINAL (ECKER & LUNA) ---\n")
results_ecker_luna = run_model_evaluation(
  all_scenarios_data = all_scenarios,
  scenario_configs = scenario_configs,
  model_function = fit_ecker_luna_mgcv,
  assessment_function = assess_effect_estimation,
  model_name = "Ecker & Luna"
)

# --- APRESENTAÇÃO DOS RESULTADOS ---

# Tabela de Resumo Comparativa
summary_table = purrr::map_dfr(
  results_ecker_luna, 
  ~ .x$summary, 
  .id = "Scenario"
)
cat("\n\n--- TABELA DE RESUMO COMPLETA (ECKER & LUNA) ---\n")
print(summary_table)


# Apresentação interativa dos gráficos
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS (ECKER & LUNA) ---\n")
purrr::iwalk(
  results_ecker_luna,
  function(assessment, scenario_name) {
    cat("\n=====================================================\n")
    cat(">>> Análise do Cenário:", scenario_name, "\n")
    cat("=====================================================\n\n")
    
    print(assessment$plot)
    
    # Pausa para permitir a visualização antes de passar para o próximo
    readline(prompt = "Pressione [enter] para ver o próximo cenário...")
  }
)
