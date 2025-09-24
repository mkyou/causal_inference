library(fda)
library(data.table)

#' @title Otimiza o Parâmetro de Suavização e Suaviza a Variável
#' @description Função auxiliar que encontra o lambda ótimo via GCV e
#'   retorna os dados suavizados junto com métricas de diagnóstico.
smooth_variable_optimally = function(time_grid, wide_matrix, nbasis) {
  
  # Cria a base B-spline uma vez
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), nbasis = nbasis
  )
  
  # 1. Define a função objetivo a ser minimizada (o escore GCV)
  gcv_objective_func = function(log_lambda) {
    # Otimizamos em escala log para estabilidade numérica
    lambda = exp(log_lambda)
    fdPar_obj = fdPar(bspline_basis, Lfdobj = 2, lambda = lambda)
    
    # smooth.basis retorna o GCV na sua saída
    # t() é necessário pois a função espera (tempo x id)
    return(
      smooth.basis(argvals = time_grid, y = t(wide_matrix), fdPar_obj)$gcv
    )
  }
  
  # 2. Otimiza para encontrar o log(lambda) que minimiza o GCV
  # Procuramos em um intervalo amplo de log(lambda)
  opt_result = optimize(gcv_objective_func, interval = c(-10, 10))
  optimal_lambda = exp(opt_result$minimum)
  min_gcv = opt_result$objective
  
  # 3. Realiza a suavização final com o lambda ótimo
  final_fdPar_obj = fdPar(
    bspline_basis, Lfdobj = 2, lambda = optimal_lambda
  )
  smoothed_fd = smooth.basis(
    argvals = time_grid, 
    y = t(wide_matrix), 
    final_fdPar_obj
  )$fd
  
  # 4. Calcula as métricas de diagnóstico
  smoothed_matrix = t(eval.fd(time_grid, smoothed_fd))
  
  # Pseudo R-squared: 1 - (Soma dos Erros Quadrados / Soma Total dos Quadrados)
  ss_error = sum((wide_matrix - smoothed_matrix)^2)
  ss_total = sum((wide_matrix - mean(as.matrix(wide_matrix)))^2)
  pseudo_r2 = 1 - (ss_error / ss_total)
  
  return(
    list(
      smoothed_matrix = smoothed_matrix,
      optimal_lambda = optimal_lambda,
      min_gcv = min_gcv,
      pseudo_r2 = pseudo_r2
    )
  )
}

#' @title Processa e Suaviza Dados para Formato Funcional
#' @description Converte dados longos para o formato funcional (lista de
#'   matrizes) e otimiza a suavização para cada variável via GCV.
#'
#' @return Uma lista contendo as variáveis funcionais, escalares, metadados
#'   e uma tabela de diagnósticos da suavização.

process_and_smooth_data = function(
    dgp_data,
    treatment_times,
    nbasis = 25 # Aumentado um pouco, já que a penalidade controla a suavidade
) {
  
  # --- 1. Setup ---
  cat("--- Iniciando preparação e suavização otimizada de dados ---\n")
  
  time_grid = sort(unique(dgp_data$t))
  ids = sort(unique(dgp_data$id))
  
  func_vars_names = names(dgp_data)[
    grep("^(X|a_X|b_X|y|a_y|b_y)", names(dgp_data))
  ]
  
  # --- 2. Processamento Otimizado das Variáveis Funcionais ---
  functional_vars = list()
  diagnostics = data.table() # Tabela para guardar os diagnósticos
  
  for (var_name in func_vars_names) {
    cat(
      paste0("  - Processando e otimizando suavização para: ", var_name, "\n")
    )
    
    wide_matrix = dcast(dgp_data, id ~ t, value.var = var_name)[, -1]
    
    # Chama a nova função auxiliar para otimizar e suavizar
    smoothing_result = smooth_variable_optimally(time_grid, wide_matrix, nbasis)
    
    # Armazena os resultados
    functional_vars[[var_name]] = smoothing_result$smoothed_matrix
    
    # Armazena e imprime os diagnósticos
    var_diagnostics = data.table(
      Variavel = var_name,
      Lambda_Otimo = smoothing_result$optimal_lambda,
      GCV_Minimo = smoothing_result$min_gcv,
      Pseudo_R2 = smoothing_result$pseudo_r2
    )
    print(var_diagnostics)
    diagnostics = rbind(diagnostics, var_diagnostics)
  }
  
  # --- 3. Processamento das Variáveis de Tratamento (Escalares) ---
  cat("  - Processando variáveis de tratamento...\n")
  treatment_data = dgp_data[t %in% treatment_times, .(id, t, A)]
  scalar_vars = dcast(treatment_data, id ~ t, value.var = "A")
  new_names = paste0("A_s", names(scalar_vars)[-1])
  setnames(scalar_vars, names(scalar_vars)[-1], new_names)
  
  cat("--- Preparação de dados concluída ---\n")
  
  # --- 4. Retorno ---
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