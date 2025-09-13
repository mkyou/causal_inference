# ===================================================================
# SCRIPT MESTRE PARA EXECUÇÃO DE EXPERIMENTOS DE SIMULAÇÃO
# ===================================================================

# --- 1. SETUP ---
# Limpa o ambiente e os gráficos para uma execução limpa
rm(list = ls()); graphics.off()

# Carrega todas as bibliotecas e scripts necessários
source("functional_simulations/data_examples/confouding.R")

source("functional_simulations/models/ecker_luna.R")
source("functional_simulations/models/parametric_g_formula.R")
source("functional_simulations/models/parametric_g_formula_covariate_stage.R")

source("functional_simulations/utils/assess_estimation.R")
source("functional_simulations/utils/run_model_evaluation.R")

# --- 2. EXECUÇÃO DOS MODELOS E AVALIAÇÃO ---
# Graças à função `run_model_evaluation`, esta parte é muito limpa.
# Podemos testar múltiplos modelos aqui.


cat("\n\n--- AVALIANDO O EFEITO FINAL NO DESFECHO (ECKER & LUNA) ---\n")
results_eckerluna_effect = run_model_evaluation(
  all_scenarios_data = all_scenarios,
  scenario_configs = scenario_configs,
  model_function = fit_ecker_luna, # <-- Modelo que estima Y
  assessment_function = assess_effect_estimation, # <-- Avalia Y
  model_name = "Ecker & Luna (Etapa de Desfecho)"
)

cat("\n--- AVALIANDO A QUALIDADE DA SIMULAÇÃO DAS COVARIÁVEIS ---\n")
results_gformula_covariates = run_model_evaluation(
  all_scenarios_data = all_scenarios[c(4, 7)],
  scenario_configs = scenario_configs,
  model_function = run_gformula_covariates,
  assessment_function = assess_covariate_simulation, # <-- Avalia X
  model_name = "G-Formula (Etapa de Covariáveis)"
)


results_g_formula = run_model_evaluation(
 all_scenarios_data = all_scenarios,
 scenario_configs = scenario_configs,
 model_function = fit_g_formula,
 model_name = "G-Formula Funcional"
)


# --- 4. APRESENTAÇÃO DOS RESULTADOS ---

# Melhoria Sugerida: Tabela de Resumo Comparativa
# Une todas as métricas de resumo em uma única tabela para fácil comparação
summary_table = purrr::map_dfr(
  results_ecker_luna, 
  ~ .x$summary, 
  .id = "Scenario"
)
cat("\n\n--- TABELA DE RESUMO COMPLETA (ECKER & LUNA) ---\n")
print(summary_table)

# Avaliação dos modelos das covariáveis
summary_table = purrr::map_dfr(
  results_gformula_covariates, 
  ~ .x$summary, 
  .id = "Scenario"
)
cat("\n\n--- TABELA DE RESUMO COMPLETA (G-FORMULA COVARIÁVEIS) ---\n")
print(summary_table)


# Apresentação interativa dos gráficos, como no seu script
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS ---\n")
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

# Apresentação interativa dos gráficos, como no seu script
cat("\n\n--- APRESENTAÇÃO INTERATIVA DOS GRÁFICOS ---\n")
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
