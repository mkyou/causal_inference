library(fda)
library(data.table)

#' Otimiza Parâmetro de Suavização e Suaviza a Variável (Versão Robusta)
smooth_variable_optimally = function(time_grid, wide_matrix, nbasis) {
  
  if (length(unique(round(as.vector(wide_matrix), 6))) <= 1) {
    return(
      list(
        smoothed_matrix = as.matrix(wide_matrix),
        optimal_lambda = NA,
        min_gcv = 0,
        pseudo_r2 = 1.0
      )
    )
  }
  
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), nbasis = nbasis
  )
  
  gcv_objective_func = function(log_lambda) {
    lambda = exp(log_lambda)
    fdPar_obj = fdPar(bspline_basis, Lfdobj = 2, lambda = lambda)
    
    gcv_val = smooth.basis(
      argvals = time_grid, y = t(wide_matrix), fdPar_obj
    )$gcv
    
    # MUDANÇA 1: Usar any() para a condição sempre ter comprimento 1.
    if (any(!is.finite(gcv_val))) {
      return(.Machine$double.xmax)
    }
    
    # MUDANÇA 2: Retornar a média para garantir uma saída escalar.
    return(mean(gcv_val))
  }
  
  opt_result = optimize(gcv_objective_func, interval = c(-20, 20))
  optimal_lambda = exp(opt_result$minimum)
  min_gcv = opt_result$objective
  
  final_fdPar_obj = fdPar(
    bspline_basis, Lfdobj = 2, lambda = optimal_lambda
  )
  smoothed_fd = smooth.basis(
    argvals = time_grid,
    y = t(wide_matrix),
    final_fdPar_obj
  )$fd
  
  smoothed_matrix = t(eval.fd(time_grid, smoothed_fd))
  
  ss_error = sum((wide_matrix - smoothed_matrix)^2)
  ss_total = sum((wide_matrix - mean(as.matrix(wide_matrix)))^2)
  
  if (abs(ss_total) < .Machine$double.eps) {
    pseudo_r2 = 1.0
  } else {
    pseudo_r2 = 1 - (ss_error / ss_total)
  }
  
  return(
    list(
      smoothed_matrix = as.matrix(smoothed_matrix),
      optimal_lambda = optimal_lambda,
      min_gcv = min_gcv,
      pseudo_r2 = pseudo_r2
    )
  )
}

#' Processa e Suaviza Dados para Formato Funcional
process_and_smooth_data = function(
    dgp_data,
    treatment_times,
    nbasis = 50
) {
  cat("--- Iniciando preparação e suavização otimizada de dados ---\n")
  
  time_grid = sort(unique(dgp_data$t))
  ids = sort(unique(dgp_data$id))
  
  func_vars_names = names(dgp_data)[
    grep("^(X|a_X|b_X|y|a_y|b_y)", names(dgp_data))
  ]
  
  functional_vars = list()
  diagnostics = data.table()
  
  for (var_name in func_vars_names) {
    cat(
      paste0("  - Processando e otimizando suavização para: ", var_name, "\n")
    )
    
    wide_matrix = dcast(dgp_data, id ~ t, value.var = var_name)[, -1]
    wide_matrix = as.matrix(wide_matrix)
    
    smoothing_result = smooth_variable_optimally(
      time_grid, wide_matrix, nbasis
    )
    
    functional_vars[[var_name]] = smoothing_result$smoothed_matrix
    
    var_diagnostics = data.table(
      Variavel = var_name,
      Lambda_Otimo = smoothing_result$optimal_lambda,
      GCV_Minimo = smoothing_result$min_gcv,
      Pseudo_R2 = smoothing_result$pseudo_r2
    )
    print(var_diagnostics)
    diagnostics = rbind(diagnostics, var_diagnostics)
  }
  
  cat("  - Processando variáveis de tratamento...\n")
  treatment_data = dgp_data[t %in% treatment_times, .(id, t, A)]
  scalar_vars = dcast(treatment_data, id ~ t, value.var = "A")
  new_names = paste0("A_s", names(scalar_vars)[-1])
  setnames(scalar_vars, names(scalar_vars)[-1], new_names)
  
  cat("--- Preparação de dados concluída ---\n")
  
  return(
    list(
      functional_vars = functional_vars,
      scalar_vars = scalar_vars,
      diagnostics = diagnostics,
      id = ids,
      time_grid = time_grid
    )
  )
}