library(data.table); library(fda); library(ggplot2)

# Funções úteis para geração de dados
generate_simulation_scenarios = function(
    n = 1000,
    time_grid = seq(0, 10, by = 0.1),
    base_y_fun = function(t) 1.5 + 0.5 * sin(0.3 * t),
    treatment_effect = list(function(t, s) 1.2 * exp(-0.05 * (t - s))),
    baseline_covariates = list(
      function(t) 0.5 + 0.15 * t^1.2,
      function(t) 0.3 + 0.1 * sqrt(t),
      function(t) 0.6 + 0.1 * log1p(t),
      function(t) 0.4 + 0.05 * t
    ),
    concurrent_effects_X = list(
      function(t) 0.8,  # Efeito constante para Cenário 1
      function(t) 0.6,  # Efeito constante para Cenário 1
      function(t) 0.75, # Efeito constante para Cenário 1
      function(t) 0.5   # Efeito constante para Cenário 1
    ),
    effects_covariates_from_A = lapply(1:4, function(k) {
      list(function(t, s) 0)  # Sem efeito de A em X para Cenário 1
    }),
    prop_func = NULL,  # Randomização para Cenário 1
    noise_sd = 0.5,
    p_treat = 0.5,
    lambda = 1e-6,
    scenario_labels = c("Cenário 1", "Cenário 2", "Cenário 3", "Cenário 4")
) {
  # Verificações iniciais
  if (length(scenario_labels) != 4) {
    stop("Deve haver exatamente 4 rótulos de cenário")
  }
  
  scenarios = list()
  
  # Configurações específicas por cenário
  scenario_configs = list(
    list(
      concurrent_effects_X = list(
        function(t) 0.8, function(t) 0.6, function(t) 0.75, function(t) 0.5
      ),
      effects_covariates_from_A = lapply(1:4, function(k) list(function(t, s) 0)),
      prop_func = NULL
    ),
    list(
      concurrent_effects_X = list(
        function(t) 0.8 * log1p(t),
        function(t) 0.6 * log1p(t),
        function(t) 0.75 * log1p(t),
        function(t) 0.5 * log1p(t)
      ),
      effects_covariates_from_A = lapply(1:4, function(k) list(function(t, s) 0)),
      prop_func = NULL
    ),
    list(
      concurrent_effects_X = list(
        function(t) 0.8 * log1p(t),
        function(t) 0.6 * log1p(t),
        function(t) 0.75 * log1p(t),
        function(t) 0.5 * log1p(t)
      ),
      effects_covariates_from_A = lapply(1:4, function(k) list(function(t, s) 0)),
      prop_func = list(function(X_list, t_index) {
        plogis(-2 + 2 * X_list[[1]][, t_index] + 1.5 * X_list[[2]][, t_index])
      })
    ),
    list(
      concurrent_effects_X = list(
        function(t) 0.8 * log1p(t),
        function(t) 0.6 * log1p(t),
        function(t) 0.75 * log1p(t),
        function(t) 0.5 * log1p(t)
      ),
      effects_covariates_from_A = lapply(1:4, function(k) {
        list(function(t, s) {
          if (t >= s) {
            amplitude = 0.1 + 0.05 * k
            decay_rate = 3.5
            shape_param = 0.7 
            time_diff = t - s
            effect = amplitude * (1 - exp(-decay_rate * time_diff)) * 
              exp(-shape_param * time_diff)
            return(effect)
          } else {
            return(0)
          }
        })
      }),
      prop_func = list(function(X_list, t_index) {
        plogis(-2 + 2 * X_list[[1]][, t_index] + 1.5 * X_list[[2]][, t_index])
      })
    )
  )
  
  for (i in 1:4) {
    # Gerar dados
    dt_conf = generate_confounder_pointwise(
      n = n,
      effects_function_list = scenario_configs[[i]]$effects_covariates_from_A,
      noise_sd = noise_sd,
      p_treat = p_treat,
      time_grid = time_grid,
      baseline_function_list = baseline_covariates,
      propensity_function_list = scenario_configs[[i]]$prop_func,
      return_format = "long"
    )
    
    dt_y = generate_y_concurrent(
      dt = dt_conf,
      base_function = base_y_fun,
      effects_function_list = treatment_effect,
      covariates_concurrent_effect_list = scenario_configs[[i]]$concurrent_effects_X,
      noise_sd = noise_sd,
      return_format = "long"
    )
    
    # Suavizar dados
    s_dt = smooth_simulated_data(dt_y, lambda = lambda)
    
    # Armazenar todos os dados
    scenarios[[i]] = list(
      smooth = s_dt,
      confounders = dt_conf,
      outcome = dt_y
    )
    
    names(scenarios)[i] = scenario_labels[i]
  }
  
  return(scenarios)
}

plot_effects_raw_from_dt = function(
    dt_y, 
    time_grid = sort(unique(dt_y$t)), 
    s_val = 1
) {
  s_index = which.min(abs(time_grid - s_val))
  
  # Identificar indivíduos tratados em s_val
  treated_ids = dt_y[A == 1 & t == s_val, unique(id)]
  dt_treated = dt_y[id %in% treated_ids]
  
  # Número de tempos
  T = length(time_grid)
  
  # Reformatar para matriz: uma linha por indivíduo, uma coluna por tempo
  y_mat = matrix(dt_treated$y, nrow = length(treated_ids), byrow = TRUE)
  control_mat = matrix(dt_treated$control_y, nrow = length(treated_ids), byrow = TRUE)
  
  # Efeito total: y - control_y
  total_effect = y_mat - control_mat
  mean_total_effect = colMeans(total_effect)
  
  # Efeito direto verdadeiro: usar coluna true_A_direct_effect_s1
  direct_mat = matrix(dt_treated$true_A_direct_effect_s1, nrow = length(treated_ids), byrow = TRUE)
  mean_direct_effect = colMeans(direct_mat)
  
  # Construir data.tables
  dt_total = data.table(
    t = time_grid,
    Mean_Effect = mean_total_effect,
    Type = "Total"
  )
  
  dt_direct = data.table(
    t = time_grid,
    Mean_Effect = mean_direct_effect,
    Type = "Direto"
  )
  
  plot_dt = rbind(dt_total, dt_direct)
  
  # Gráfico
  ggplot(plot_dt, aes(x = t, y = Mean_Effect, linetype = Type)) +
    geom_line(color = "darkred", linewidth = 1) +
    labs(
      title = "Efeito médio nos tratados — sem suavização",
      x = "Tempo",
      y = "Efeito médio direto vs total",
      linetype = "Tipo de efeito"
    ) +
    theme_light()
}

plot_effects_treated = function(
    scenarios, 
    scenario_labels, 
    time_grid = seq(0, 10, by = 0.1)
) {
  s_val = 1
  s_index = which.min(abs(time_grid - s_val))
  
  effect_data_total = rbindlist(
    lapply(seq_along(scenarios), function(i) {
      s_dt = scenarios[[i]]$smooth
      A_mat = s_dt$A_mat
      treated_ids = which(A_mat[, s_index] == 1)
      fd_obj = s_dt$true_A_effect_fd
      effect_vals = eval.fd(time_grid, fd_obj)
      effect_vals = effect_vals[, treated_ids, drop = FALSE]
      means = rowMeans(effect_vals)
      data.table(
        t = time_grid,
        Mean_Effect = means,
        Scenario = scenario_labels[i],
        Type = "Total"
      )
    })
  )
  
  effect_data_direct = rbindlist(
    lapply(seq_along(scenarios), function(i) {
      s_dt = scenarios[[i]]$smooth
      A_mat = s_dt$A_mat
      treated_ids = which(A_mat[, s_index] == 1)
      fd_obj = s_dt$true_A_direct_effect_fd_list[[1]]
      effect_vals = eval.fd(time_grid, fd_obj)
      effect_vals = effect_vals[, treated_ids, drop = FALSE]
      means = rowMeans(effect_vals)
      data.table(
        t = time_grid,
        Mean_Effect = means,
        Scenario = scenario_labels[i],
        Type = "Direto"
      )
    })
  )
  
  combined_data = rbind(effect_data_total, effect_data_direct)
  
  ggplot(combined_data, aes(x = t, y = Mean_Effect, color = Type, linetype = Type)) +
    geom_line(size = 1) +
    facet_wrap(~ Scenario, scales = "free_y") +
    labs(
      title = "Efeito médio do tratamento nos indivíduos tratados",
      x = "Tempo",
      y = "Efeito médio suavizado",
      color = "Tipo de efeito",
      linetype = "Tipo de efeito"
    ) +
    theme_bw() +
    theme(strip.text = element_text(face = "bold"))
}

generate_covariate_comparison_df = function(
    dt, 
    time_grid = sort(unique(dt$t)), 
    s_val = 1
) {
  s_index = which.min(abs(time_grid - s_val))
  treated_ids = dt[A == 1 & t == s_val, unique(id)]
  dt_treated = dt[id %in% treated_ids]
  
  X_names = grep("^X\\d+", names(dt), value = TRUE)
  control_X_names = paste0("control_", X_names)
  K = length(X_names)
  
  results = vector("list", K)
  
  for (k in seq_len(K)) {
    X_mat        = matrix(dt_treated[[X_names[k]]], nrow = length(treated_ids), byrow = TRUE)
    control_mat  = matrix(dt_treated[[control_X_names[k]]], nrow = length(treated_ids), byrow = TRUE)
    
    df = data.table(
      t = time_grid,
      X_mean        = colMeans(X_mat),
      Control_mean  = colMeans(control_mat),
      Delta         = colMeans(X_mat - control_mat),
      Covariate     = paste0("X", k)
    )
    
    results[[k]] = df
  }
  
  final_dt = rbindlist(results)
  return(final_dt)
}


# Aplicação de função
scenarios = generate_simulation_scenarios()

stress_test = generate_simulation_scenarios(
  n = 1000,
  time_grid = seq(0, 10, by = 0.1),
  
  # 1. Efeito direto ZERO em todos os cenários
  treatment_effect = list(function(t, s) 0), 
  
  # 2. Efeitos indiretos FORTES (A → X → Y)
  effects_covariates_from_A = lapply(1:4, function(k) {
    list(function(t, s) {
      if (t >= s) {
        amplitude = 0.5 + 0.2 * k  # Efeito mais forte
        decay_rate = 2.0
        shape_param = 0.5
        time_diff = t - s
        effect = amplitude * (1 - exp(-decay_rate * time_diff)) * 
          exp(-shape_param * time_diff)
        return(effect)
      } else {
        return(0)
      }
    })
  }),
  
  # 3. Forte efeito das covariáveis em Y
  concurrent_effects_X = list(
    function(t) 1.5 * log1p(t),  # Mais forte
    function(t) 1.2 * log1p(t),  # Mais forte
    function(t) 1.8 * log1p(t),  # Mais forte
    function(t) 1.0 * log1p(t)   # Mais forte
  ),
  
  # Restante dos parâmetros mantidos
  base_y_fun = function(t) 1.5 + 0.5 * sin(0.3 * t),
  baseline_covariates = list(
    function(t) 0.5 + 0.15 * t^1.2,
    function(t) 0.3 + 0.1 * sqrt(t),
    function(t) 0.6 + 0.1 * log1p(t),
    function(t) 0.4 + 0.05 * t
  ),
  prop_func = list(function(X_list, t_index) {
    plogis(-2 + 2 * X_list[[1]][, t_index] + 1.5 * X_list[[2]][, t_index])
  }),
  noise_sd = 0.5,
  p_treat = 0.5,
  lambda = 1e-6,
  scenario_labels = c("Stress 1", "Stress 2", "Stress 3")
)

# Validação dos cenários
plot_effects_treated(scenarios, c('Cenário 1', 'Cenário 2', 'Cenário 3',
                                  'Cenário 4'))
plot_effects_treated(stress_test, c('Cenário 1', 'Cenário 2', 'Cenário 3'))

# plot_effects_raw_from_dt(scenarios[[1]]$outcome)
# plot_effects_raw_from_dt(scenarios[[2]]$outcome)
# plot_effects_raw_from_dt(scenarios[[3]]$outcome)

# generate_covariate_comparison_df(scenarios[[3]]$confounders) |> 
#   summarise(
#     'min' = min(Delta), 
#     'mean' = mean(Delta), 
#     'median' = median(Delta),
#     'max' = max(Delta),
#     .by = 'Covariate'
#   )

