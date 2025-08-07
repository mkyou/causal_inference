#' Simula variável de resposta Y em função de tratamento e covariáveis temporais
#'
#' Gera dados de desfecho Y que dependem de:
#' \itemize{
#'   \item Uma trajetória basal
#'   \item Efeitos diretos do tratamento ao longo do tempo
#'   \item Efeitos concorrentes das covariáveis
#'   \item Ruído aleatório
#' }
#' Pode retornar no formato longo ou amplo.
#'
#' @param dt `data.table` gerado por `generate_confounder_pointwise` contendo:
#'            \itemize{
#'              \item `id`: Identificador do indivíduo
#'              \item `t`: Tempo da observação  
#'              \item `A`: Indicador de tratamento
#'              \item `X1...Xk`: Covariáveis observadas
#'              \item `control_X1...control_Xk`: Covariáveis contrafactuais
#'            }
#' @param base_function Função que define a trajetória basal de Y sem 
#' tratamento. Deve aceitar um vetor de tempos e retornar valores numéricos.
#' @param effects_function_list Lista de funções definindo efeitos diretos do 
#' tratamento. Cada função deve aceitar dois argumentos:
#'   \itemize{
#'     \item \code{t}: Tempo atual da observação
#'     \item \code{s}: Tempo específico de tratamento
#'   }
#' O índice da função na lista corresponde ao tempo de tratamento \code{s}.
#' 
#' @param covariates_concurrent_effect_list Lista de funções definindo efeitos 
#' concorrentes das covariáveis. Cada função deve aceitar o tempo `t`.
#' @param noise_sd Desvio padrão do ruído aditivo (numérico positivo).
#' @param return_format Formato de saída: "long", "wide" ou "both".
#'
#' @return Dependendo do `return_format`:
#' \describe{
#'   \item{long}{`data.table` com colunas adicionais:
#'     \itemize{
#'       \item `y`: Valores observados de Y
#'       \item `control_y`: Valores contrafactuais de Y
#'       \item `true_A_direct_effect_s1...`: Efeitos diretos verdadeiros por 
#'       tempo de tratamento
#'     }
#'   }
#'   \item{wide}{Lista contendo:
#'     \itemize{
#'       \item `y_mat`: Matriz n x T de Y observado
#'       \item `control_y_mat`: Matriz n x T de Y contrafactual  
#'       \item `X_mats`: Lista de matrizes n x T das covariáveis
#'       \item `t_vec`: Vetor de tempos
#'       \item `true_A_direct_effect`: Array n x T x S de efeitos diretos
#'     }
#'   }
#'   \item{both}{Lista com ambos formatos (`long` e `wide`)}
#' }
#'
#' @details
#' O modelo para Y é dado por:
#' \deqn{Y_i(t) = base(t) + \sum_{s \leq t} \beta(t,s)A_i(s) + 
#'  \sum_{k=1}^K \gamma_k(t)X_{ik}(t) + \epsilon_i(t)}
#' onde:
#' \itemize{
#'   \item \eqn{base(t)} é a função basal
#'   \item \eqn{\beta(t,s)} são os efeitos diretos do tratamento
#'   \item \eqn{\gamma_k(t)} são os efeitos concorrentes das covariáveis
#'   \item \eqn{\epsilon_i(t)} é o ruído aleatório
#' }
#'
#' @examples
#' # Exemplo mínimo
#' library(data.table)
#' 
#' # Gerar dados de covariáveis
#' dt_conf = generate_confounder_pointwise(
#'   n = 100,
#'   effects_function_list = list(list(function(t,s) 0)), 
#'   noise_sd = 0.1,
#'   p_treat = 0.3,
#'   time_grid = 1:10
#' )
#' 
#' # Gerar Y
#' dt_y = generate_y_concurrent(
#'   dt = dt_conf,
#'   base_function = function(t) 2 + sin(t),
#'   effects_function_list = list(function(t,s) 0.5 * exp(-0.1*(t-s))),
#'   covariates_concurrent_effect_list = list(function(t) 0.3),
#'   noise_sd = 0.2
#' )
#' 
#' # Exemplo com múltiplos tratamentos
#' treatment_effects = list(
#'   function(t,s) 0.8 * exp(-0.1*(t-s)),  # Efeito para tratamento em s=1
#'   function(t,s) 0.5 * exp(-0.05*(t-s))  # Efeito para tratamento em s=2
#' )
#' 
#' dt_y2 = generate_y_concurrent(
#'   dt = dt_conf,
#'   base_function = function(t) 2 + 0.5*t),
#'   effects_function_list = treatment_effects,
#'   covariates_concurrent_effect_list = list(
#'     function(t) 0.2,
#'     function(t) 0.4 * log(t+1)
#'   ),
#'   noise_sd = 0.15,
#'   return_format = "both"
#' )
#' 
#' #' # Exemplo com tratamento em múltiplos tempos
#' effects = list(
#'   # Efeito para tratamento iniciado em t=1
#'   function(t, s) 1.2 * exp(-0.2*(t-s)),
#'   # Efeito para tratamento iniciado em t=2
#'   function(t, s) 0.8 * (1 - exp(-0.3*(t-s)))
#' )
#' 
#' dt_y = generate_y_concurrent(
#'   dt = data,
#'   base_function = function(t) 2 + sin(0.5*t),
#'   effects_function_list = effects,
#'   covariates_concurrent_effect_list = list(
#'     function(t) 0.3 * t,
#'     function(t) 0.2 * log(t+1)
#'   ),
#'   noise_sd = 0.25
#' )
#'
#' @seealso \code{\link{generate_confounder_pointwise}} Para gerar os dados de 
#' entrada necessários.
#' @export
#' 

generate_y_concurrent = function(
    dt, 
    base_function,
    effects_function_list,
    covariates_concurrent_effect_list, 
    noise_sd = 1, 
    return_format = c("long", "wide", "both")
) {
  # Verificações iniciais
  requireNamespace("data.table")
  return_format = match.arg(return_format)
  
  required_cols = c("id", "t", "A")
  if (!all(required_cols %in% names(dt))) {
    stop("dt deve conter colunas: id, t, A")
  }
  
  # Identificação de covariáveis
  Xnames = grep("^X\\d+", names(dt), value = TRUE)
  control_Xnames = paste0("control_", Xnames)
  num_covariates = length(Xnames)
  
  if (num_covariates != length(covariates_concurrent_effect_list)) {
    stop("Número de covariáveis incompatível com 
         covariates_concurrent_effect_list")
  }
  
  # Preparação de dados
  dt = data.table::copy(data.table::setorder(dt, id, t))
  ids = unique(dt$id)
  times = sort(unique(dt$t))
  num_ids = length(ids)
  num_times = length(times)
  
  # Definição dos tempos de tratamento (s = 1,2,...,S)
  num_treatment_times = length(effects_function_list)
  treatment_times = 1:num_treatment_times
  
  # Conversão para formato wide
  value_vars = c("A", Xnames, control_Xnames)
  dt_wide = data.table::dcast(dt, id ~ t, value.var = value_vars)
  
  # Extração de matrizes com nomes
  A_cols = paste0("A_", times)
  A_mat = as.matrix(dt_wide[, .SD, .SDcols = A_cols])
  colnames(A_mat) = A_cols
  
  X_arr = array(NA, dim = c(num_ids, num_times, num_covariates))
  for (k in seq_len(num_covariates)) {
    X_cols = paste0(Xnames[k], "_", times)
    X_arr[, , k] = as.matrix(dt_wide[, .SD, .SDcols = X_cols])
  }
  
  control_X_arr = array(NA, dim = c(num_ids, num_times, num_covariates))
  for (k in seq_len(num_covariates)) {
    control_cols = paste0(control_Xnames[k], "_", times)
    control_X_arr[, , k] = as.matrix(dt_wide[, .SD, .SDcols = control_cols])
  }
  
  # Matriz de efeitos diretos (beta)
  beta_mat = matrix(0, nrow = num_treatment_times, ncol = num_times)
  for (s_idx in seq_len(num_treatment_times)) {
    s_val = treatment_times[s_idx]
    for (t_idx in seq_len(num_times)) {
      t_val = times[t_idx]
      if (t_val >= s_val) {
        beta_mat[s_idx, t_idx] = effects_function_list[[s_idx]](t_val, s_val)
      }
    }
  }
  
  # Identificação de colunas de tratamento
  treatment_cols = paste0("A_", treatment_times)
  treatment_col_indices = which(colnames(A_mat) %in% treatment_cols)
  
  if (length(treatment_col_indices) != num_treatment_times) {
    missing = setdiff(
      treatment_times, 
      as.numeric(sub("A_", "", treatment_cols))
    )
    stop("Tempos de tratamento ausentes: ", paste(missing, collapse = ", "))
  }
  
  # Efeitos diretos (CORRETO)
  direct_effects = A_mat[, treatment_col_indices, drop = FALSE] %*% beta_mat
  
  # Componentes de Y
  base_vals = matrix(base_function(times), num_ids, num_times, byrow = TRUE)
  covar_effects = matrix(0, num_ids, num_times)
  
  for (k in seq_len(num_covariates)) {
    gamma_k = sapply(times, covariates_concurrent_effect_list[[k]])
    gamma_mat = matrix(gamma_k, num_ids, num_times, byrow = TRUE)
    covar_effects = covar_effects + X_arr[, , k] * gamma_mat
  }
  
  # Y observado
  noise = matrix(rnorm(num_ids * num_times, 0, noise_sd), num_ids, num_times)
  y_mat = base_vals + direct_effects + covar_effects + noise
  
  # Y contrafactual
  control_covar_effects = matrix(0, num_ids, num_times)
  for (k in seq_len(num_covariates)) {
    gamma_k = sapply(times, covariates_concurrent_effect_list[[k]])
    gamma_mat = matrix(gamma_k, num_ids, num_times, byrow = TRUE)
    control_covar_effects = control_covar_effects + 
      control_X_arr[, , k] * gamma_mat
  }

  control_y_mat = base_vals + control_covar_effects + noise
  
  # Efeitos diretos verdadeiros
  true_effects = array(0, dim = c(num_ids, num_times, num_treatment_times))
  
  for (s_idx in seq_len(num_treatment_times)) {
    s_val = treatment_times[s_idx]
    s_col = paste0("A_", s_val)
    col_index = match(s_col, colnames(A_mat))
    
    for (t_idx in seq_len(num_times)) {
      t_val = times[t_idx]
      if (t_val >= s_val) {
        true_effects[, t_idx, s_idx] = A_mat[, col_index] * 
          effects_function_list[[s_idx]](t_val, s_val)
      }
    }
  }
  
  # Adicionar ao dataframe
  dt[, y := as.vector(t(y_mat))]
  dt[, control_y := as.vector(t(control_y_mat))]
  
  for (s in seq_len(num_treatment_times)) {
    effect_name = paste0("true_A_direct_effect_s", s)
    dt[[effect_name]] = as.vector(t(true_effects[,,s]))
  }
  
  # Retorno formatado
  switch(return_format,
         long = dt,
         wide = list(
           y_mat = y_mat,
           control_y_mat = control_y_mat,
           X_mats = lapply(seq_len(num_covariates), function(k) X_arr[,,k]),
           t_vec = times,
           true_A_direct_effect = true_effects
         ),
         both = list(long = dt, wide = list(
           y_mat = y_mat,
           control_y_mat = control_y_mat,
           X_mats = lapply(seq_len(num_covariates), function(k) X_arr[,,k]),
           t_vec = times,
           true_A_direct_effect = true_effects
         ))
  )
}

