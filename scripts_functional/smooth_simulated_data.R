#' Suavização Funcional de Dados Temporais Simulados
#'
#' Transforma trajetórias temporais discretas em objetos funcionais contínuos 
#' usando bases B-spline, preservando a estrutura de dados original e 
#' adicionando verificações de qualidade na suavização.
#'
#' @param dt_long Um `data.table` no formato longo contendo:
#' \itemize{
#'   \item \code{id}: Identificador único do indivíduo
#'   \item \code{t}: Tempo da observação
#'   \item \code{A}: Indicador de tratamento (0/1)
#'   \item \code{y}: Valores observados do desfecho
#'   \item \code{control_y}: Valores contrafactuais do desfecho
#'   \item \code{X1, X2, ...}: Covariáveis observadas
#'   \item \code{true_A_direct_effect_s1, ...}: Efeitos diretos do tratamento
#' }
#' @param nbasis Número de funções base para a expansão B-spline. 
#' Se NULL (padrão), calcula automaticamente como max(7, 
#' min(40, ceiling(n_tempos/4))).
#' @param lambda Parâmetro de suavização (default = 1e-6). Valores maiores 
#' resultam em suavizações mais intensas.
#' @param gcv_threshold Limite para Generalized Cross-Validation (GCV). Se o 
#' GCV médio exceder este valor (default = 1.5), emite um aviso.
#'
#' @return Uma lista contendo:
#' \itemize{
#'   \item \code{nbasis_used}: Número de bases efetivamente usadas
#'   \item \code{y_fd}: Objeto funcional (fd) para y observado
#'   \item \code{y_gcv}: Valores GCV para cada indivíduo (y observado)
#'   \item \code{true_A_effect_fd}: Objeto funcional para o efeito causal 
#'   (y - control_y)
#'   \item \code{true_A_effect_gcv}: GCV para o efeito causal
#'   \item \code{X_fd_list}: Lista de objetos funcionais para cada covariável
#'   \item \code{X_gcv_list}: Lista de vetores GCV para cada covariável
#'   \item \code{true_A_direct_effect_fd_list}: Lista de objetos funcionais para 
#'   efeitos diretos
#'   \item \code{true_A_direct_effect_gcv_list}: Lista de vetores GCV para 
#'   efeitos diretos
#'   \item \code{A_mat}: Matriz de tratamento (n_indivíduos × n_tempos)
#' }
#'
#' @details
#' A função realiza:
#' \enumerate{
#'   \item Verificação da estrutura de entrada
#'   \item Criação de base B-spline
#'   \item Suavização funcional para:
#'     \itemize{
#'       \item y observado e contrafactual
#'       \item Efeito causal (diferença entre y e control_y)
#'       \item Todas as covariáveis
#'       \item Todos os efeitos diretos do tratamento
#'     }
#'   \item Cálculo de métricas GCV para cada suavização
#'   \item Emissão de avisos se GCV exceder o threshold
#' }
#'
#' @examples
#' \dontrun{
#' # Gerar dados simulados
#' dt_conf = generate_confounder_pointwise(
#'   n = 100,
#'   effects_function_list = list(list(function(t,s) 0)),
#'   noise_sd = 0.1,
#'   p_treat = 0.3,
#'   time_grid = 1:10
#' )
#'
#' dt_y = generate_y_concurrent(
#'   dt = dt_conf,
#'   base_function = function(t) 2 + sin(t),
#'   effects_function_list = list(function(t,s) 0.5*exp(-0.1*(t-s))),
#'   covariates_concurrent_effect_list = list(function(t) 0.3)
#' )
#'
#' # Suavizar os dados
#' smoothed_data = smooth_simulated_data(dt_y)
#'
#' # Acessar resultados
#' plot(smoothed_data$y_fd, main = "Y Observado Suavizado")
#' summary(smoothed_data$y_gcv)
#' }
#'
#' @seealso \code{\link{generate_confounder_pointwise}}, 
#' \code{\link{generate_y_concurrent}}
#' @export

smooth_simulated_data = function(
    dt_long,
    nbasis = NULL,
    lambda = 1e-6,
    gcv_threshold = 1.5
) {
  # 1. Verificações de entrada
  required_cols = c("id", "t", "A", "y", "control_y")
  if (!all(required_cols %in% names(dt_long))) {
    missing_cols = setdiff(required_cols, names(dt_long))
    stop("O dataframe deve conter as colunas: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  if (!is.numeric(lambda) || lambda <= 0) {
    stop("lambda deve ser um número positivo")
  }
  
  # 2. Configuração inicial
  time_grid = sort(unique(dt_long$t))
  n = length(unique(dt_long$id))
  
  if (is.null(nbasis)) {
    nbasis = max(7, min(40, ceiling(length(time_grid) / 4)))
  }
  
  basis = fda::create.bspline.basis(
    rangeval = range(time_grid), 
    nbasis = nbasis
  )
  fdParobj = fda::fdPar(basis, Lfdobj = fda::int2Lfd(2), lambda = lambda)
  
  # 3. Função de suavização com verificação GCV
  smooth_with_gcv = function(data_mat, var_name = "") {
    smooth = fda::smooth.basis(
      argvals = time_grid,
      y = t(data_mat),
      fdParobj = fdParobj
    )
    
    if (mean(smooth$gcv) > gcv_threshold) {
      warning("GCV médio alto (", round(mean(smooth$gcv), 2), 
              ") para ", var_name, ". Considere ajustar nbasis ou lambda.")
    }
    
    return(list(fd = smooth$fd, gcv = smooth$gcv))
  }
  
  # 4. Processamento eficiente
  y_mat = matrix(dt_long$y, nrow = n, byrow = TRUE)
  control_y_mat = matrix(dt_long$control_y, nrow = n, byrow = TRUE)
  
  y_result = smooth_with_gcv(y_mat, "y observado")
  control_y_result = smooth_with_gcv(control_y_mat, "y contrafactual")
  effect_result = smooth_with_gcv(y_mat - control_y_mat, "efeito causal")
  
  # 5. Processar covariáveis
  covariate_names = grep("^X\\d+", names(dt_long), value = TRUE)
  X_results = lapply(covariate_names, function(v) {
    smooth_with_gcv(matrix(dt_long[[v]], nrow = n, byrow = TRUE), v)
  })
  
  # 6. Processar efeitos diretos
  treatment_effect_names = grep(
    "^true_A_direct_effect_s", 
    names(dt_long), 
    value = TRUE
  )
  
  A_effect_results = lapply(treatment_effect_names, function(v) {
    smooth_with_gcv(matrix(dt_long[[v]], nrow = n, byrow = TRUE), v)
  })
  
  # 7. Retornar resultados na estrutura original
  return(list(
    nbasis_used = nbasis,
    y_fd = y_result$fd,
    y_gcv = y_result$gcv,
    true_A_effect_fd = effect_result$fd,
    true_A_effect_gcv = effect_result$gcv,
    X_fd_list = lapply(X_results, function(x) x$fd),
    X_gcv_list = lapply(X_results, function(x) x$gcv),
    true_A_direct_effect_fd_list = lapply(A_effect_results, function(x) x$fd),
    true_A_direct_effect_gcv_list = lapply(A_effect_results, function(x) x$gcv),
    A_mat = matrix(dt_long$A, nrow = n, ncol = length(time_grid), byrow = TRUE)
  ))
}
