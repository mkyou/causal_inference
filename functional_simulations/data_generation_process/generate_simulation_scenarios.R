#' @title Gera Múltiplos Cenários de Simulação
#' @description Orquestra a geração de cenários, chamando
#'   `generate_sequential_dgp` e `generate_functional_outcome`.
#'   Agora inclui uma etapa opcional para preparar os dados para
#'   pacotes de FDA como o `refund`.
#'
#' @param scenario_configs Lista nomeada de configurações de cenário.
#' @param global_params Lista de parâmetros globais.
#' @param prepare_for_refund Booleano. Se TRUE, a saída de cada cenário
#'   será processada por `prepare_refund_data`.
#' @param refund_nbasis,refund_lambda Parâmetros de suavização a serem
#'   passados para `prepare_refund_data`.
#'
#' @return Uma lista nomeada de cenários. O formato de cada elemento
#'   depende de `prepare_for_refund`.

generate_simulation_scenarios = function(
    scenario_configs, 
    global_params = list(),
    prepare_for_refund = FALSE,
    refund_nbasis = 20,
    refund_lambda = 1e-4
) {
  # Validação
  if (missing(scenario_configs) || !is.list(scenario_configs) || 
      length(scenario_configs) == 0) {
    stop("scenario_configs deve ser uma lista não vazia de configurações.")
  }
  
  labels = if (!is.null(names(scenario_configs))) {
    names(scenario_configs)
  } else {
    paste0("Cenario_", seq_along(scenario_configs))
  }
  
  all_scenarios_data = vector("list", length(scenario_configs))
  
  for (i in seq_along(scenario_configs)) {
    cat(paste("\n--- Gerando Cenário:", labels[i], "---\n"))
    
    # 1. Combina parâmetros
    cfg_specific = scenario_configs[[i]]
    params = utils::modifyList(global_params, cfg_specific)
    
    # 2. Etapa 1: Gerar as covariáveis e tratamentos
    cat("  - Etapa 1: Gerando covariáveis (generate_sequential_dgp)\n")
    dgp_data = generate_sequential_dgp(
      n_subjects = params$n_subjects,
      time_grid = params$time_grid,
      treatment_times = params$treatment_times,
      baseline_funcs = params$baseline_funcs,
      beta_funcs_matrix = params$beta_funcs_matrix,
      gamma_matrix = params$gamma_matrix,
      theta_matrix = params$theta_matrix,
      theta_integral_matrix = params$theta_integral_matrix,
      intercept_sd_vec = params$intercept_sd_vec,
      noise_sd_vec = params$noise_sd_vec,
      regime_a = params$regime_a,
      regime_b = params$regime_b
    )
    
    # 3. Etapa 2: Gerar o desfecho
    cat("  - Etapa 2: Gerando desfecho (generate_functional_outcome)\n")
    final_dt = generate_functional_outcome(
      dgp_data = dgp_data,
      baseline_y_func = params$baseline_y_func,
      direct_effect_funcs = params$direct_effect_funcs,
      concurrent_effect_funcs = params$concurrent_effect_funcs,
      treatment_times = params$treatment_times,
      noise_sd = params$outcome_noise_sd
    )
    
    # --- INÍCIO DA MUDANÇA ---
    # 4. Etapa 3 (Opcional): Preparar dados para refund/FDA
    if (prepare_for_refund) {
      cat("  - Etapa 3: Preparando dados para formato refund/FDA\n")
      
      # Permite que os parâmetros de suavização sejam específicos do cenário
      nbasis = params$refund_nbasis %||% refund_nbasis
      lambda = params$refund_lambda %||% refund_lambda
      
      processed_data = prepare_refund_data(
        dgp_data = final_dt,
        treatment_times = params$treatment_times,
        nbasis = nbasis,
        lambda = lambda
      )
      all_scenarios_data[[i]] = processed_data
    } else {
      all_scenarios_data[[i]] = final_dt
    }
  }
  
  names(all_scenarios_data) = labels
  cat("\n--- Geração de todos os cenários concluída. ---\n")
  return(all_scenarios_data)
}

# Operador auxiliar para definir valores padrão
`%||%` = function(a, b) {
  if (is.null(a)) b else a
}