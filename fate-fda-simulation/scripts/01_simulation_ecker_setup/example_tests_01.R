# Carregar Infraestrutura
source("R/utils_dgp_functions.R")
source("R/generate_sequential_dgp.R")     
source("R/generate_functional_outcomes.R") 
source("R/config_ecker.R")                
source("R/analyze_dgp.R")  
source("R/scenario_runner.R") 

# Parâmetros
N = 2000
T_pts = 512

# Wrapper Rápido de Teste
test_scenario = function(name, mag) {
  cat(sprintf("\n>>> Testando Cenário: %s <<<\n", name))
  
  # Gera usando a mesma função da simulação
  dgp_data = generate_single_scenario_dataset(name, mag, N, T_pts)
  
  # Extrai parâmetros para validação
  t_treat = 20.0
  time_grid = attr(dgp_data, "time_grid")
  true_curve = attr(dgp_data, "true_effect_curve")
  
  # Validação Visual
  print(validate_static_distributions(dgp_data))
  print(validate_dgp_effects(dgp_data, t_treat))
  print(validate_outcome_means(dgp_data, t_treat))
}

# Execução
test_scenario("Sigmoid", 2)
test_scenario("Doppler", 5)
test_scenario("Bumps", 3)
test_scenario("Heavisine", 4)