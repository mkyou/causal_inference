# Certifique-se de que o pacote 'fda' está instalado
# install.packages("fda")
library(fda)
library(data.table)

#' @title Prepara Dados para Modelagem com Pacotes de FDA (ex: refund)
#' @description Converte o data.table longo do PGD para o formato largo
#'   (matricial) esperado por pacotes como o `refund`. A função separa
#'   variáveis funcionais (que são suavizadas) e variáveis escalares.
#'
#' @param dgp_data O data.table no formato longo gerado por
#'   `generate_functional_outcome`.
#' @param treatment_times Vetor numérico com os tempos exatos de tratamento.
#' @param nbasis Número de funções de base para a suavização.
#' @param lambda Parâmetro de penalização para a suavização.
#'
#' @return Uma lista contendo:
#'   - `functional_vars`: Uma lista nomeada de matrizes (id x tempo) para
#'     cada variável funcional suavizada.
#'   - `scalar_vars`: Um data.table (id x tratamento) com as
#'     variáveis de tratamento.
#'   - `id`: Vetor de IDs.
#'   - `time_grid`: Vetor da grade de tempo.

prepare_refund_data = function(
    dgp_data,
    treatment_times,
    nbasis = 20,
    lambda = 1e-4
) {
  
  # --- 1. Setup ---
  cat("--- Iniciando preparação de dados para formato funcional ---\n")
  
  time_grid = sort(unique(dgp_data$t))
  ids = sort(unique(dgp_data$id))
  n_subjects = length(ids)
  
  # Identifica as variáveis funcionais a serem processadas
  func_vars_names = names(dgp_data)[
    grep("^(X|a_X|b_X|y|a_y|b_y)", names(dgp_data))
  ]
  
  # --- 2. Processamento das Variáveis Funcionais ---
  
  # Cria a base de B-splines para a suavização
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), nbasis = nbasis
  )
  # Cria o objeto de penalização
  fdPar_obj = fdPar(bspline_basis, Lfdobj = 2, lambda = lambda)
  
  functional_vars = list()
  
  for (var_name in func_vars_names) {
    cat(paste0("  - Processando e suavizando: ", var_name, "\n"))
    
    # a. Transforma de longo para largo (matriz)
    wide_matrix = dcast(
      dgp_data, 
      id ~ t, 
      value.var = var_name
    )[, -1] # Remove a coluna de id
    
    # b. Suaviza os dados usando o pacote fda
    # `t()` é necessário pois smooth.basis espera (tempo x id)
    smoothed_fd = smooth.basis(
      argvals = time_grid, 
      y = t(wide_matrix), 
      fdPar = fdPar_obj
    )$fd
    
    # c. Volta para o formato matricial (id x tempo)
    # `eval.fd` avalia o objeto funcional na grade de tempo
    smoothed_matrix = t(eval.fd(time_grid, smoothed_fd))
    
    functional_vars[[var_name]] = smoothed_matrix
  }
  
  # --- 3. Processamento das Variáveis de Tratamento (Escalares) ---
  cat("  - Processando variáveis de tratamento...\n")
  
  # Filtra apenas para os tempos em que o tratamento ocorre
  treatment_data = dgp_data[t %in% treatment_times, .(id, t, A)]
  
  # Transforma de longo para largo para ter uma coluna por tratamento
  scalar_vars = dcast(
    treatment_data,
    id ~ t,
    value.var = "A"
  )
  
  # Renomeia as colunas para um formato limpo (ex: "5" -> "A_s5")
  new_names = paste0("A_s", names(scalar_vars)[-1])
  setnames(scalar_vars, names(scalar_vars)[-1], new_names)
  
  cat("--- Preparação de dados concluída ---\n")
  
  # --- 4. Retorno ---
  return(
    list(
      functional_vars = functional_vars,
      scalar_vars = scalar_vars,
      id = ids,
      time_grid = time_grid
    )
  )
}