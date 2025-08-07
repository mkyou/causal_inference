#' Simula covariáveis em função do tratamento ao longo do tempo
#'
#' Gera dados temporais de covariáveis afetadas por múltiplos tratamentos em 
#' diferentes tempos, considerando funções de base, efeitos causais, ruído e 
#' mecanismos de propensão. Pode retornar no formato `long` ou `wide`.
#'
#' @param n Número de indivíduos simulados (inteiro positivo).
#' @param effects_function_list Lista com listas de funções de efeito 
#'   indireto por covariável. Cada função deve aceitar (t, s) onde t é o tempo
#'   atual e s é o tempo de tratamento. Deve retornar um valor numérico.
#'   A estrutura deve ser: 
#'     effects_function_list[[k]][[s]] = função para covariável k e tratamento 
#'     no tempo s.
#' @param p_treat Probabilidade fixa de tratamento (se não usar propensão) 
#'   (numérico entre 0 e 1).
#' @param x_noise_model Modelo de ruído para covariáveis ("iid" ou "ar1")
#' @param x_noise_rho Parâmetro de autocorrelação para modelo AR(1)
#' @param noise_sd Desvio padrão do ruído (numérico positivo).
#' @param time_grid Vetor de tempos de simulação (numérico crescente). Se NULL,
#'   será criado automaticamente com base no número de funções de efeito.
#' @param baseline_function_list Lista com funções de baseline por covariável.
#'   Se NULL, assume constante zero. Cada função deve aceitar um tempo t.
#' @param propensity_function_list Lista com funções de propensão ao tratamento.
#'   Cada função deve aceitar (X_list, t_index) e retornar vetores de 
#'   probabilidade (numéricos entre 0 e 1) para cada indivíduo.
#' @param return_format Define formato de saída: "long", "wide" ou "both".
#'
#'
#' @return Retorna:
#'   - "long": \code{data.table} com colunas:
#'        \itemize{
#'          \item \code{id}: Identificador do indivíduo
#'          \item \code{t}: Tempo da observação
#'          \item \code{A}: Indicador de tratamento
#'          \item \code{X1..Xk}: Valores observados das covariáveis
#'          \item \code{control_X1..control_Xk}: Valores contrafactuais das 
#'          covariáveis (como seriam geradas na ausência de tratamento)
#'        }
#'   - "wide": Lista com elementos:
#'        \itemize{
#'          \item \code{A_mat}: Matriz n x T de tratamentos
#'          \item \code{A_control}: Matriz n x T de controles (sempre 0)
#'          \item \code{X_list}: Lista de k matrizes n x T de covariáveis 
#'          observadas
#'          \item \code{control_X_list}: Lista de k matrizes n x T de 
#'          covariáveis contrafactuais
#'          \item \code{t_vec}: Vetor de tempos
#'        }
#'   - "both": Lista contendo ambos os formatos ("long" e "wide")
#'
#' @details Cada indivíduo recebe tratamento em tempos discretos. 
#'   As covariáveis são atualizadas ao longo do tempo com base em:
#' \enumerate{
#'   \item Funções baseline (trajetória sem tratamento)
#'   \item Efeitos cumulativos de tratamentos anteriores 
#'   (via effects_function_list)
#'   \item Ruído aleatório (normal com sd = noise_sd)
#' }
#' Também são geradas versões contrafactuais (\code{control_X}) simulando A = 0 
#' para todos os indivíduos em todos os tempos.
#'
#' @examples
#' # Exemplo 1: Uso básico com efeito constante
#' eff_list = list(list(function(t, s) 0.5))
#' conf_data = generate_confounder_pointwise(
#'   n = 100, effects_function_list = eff_list, noise_sd = 0.1,
#'   p_treat = 0.2, return_format = "long"
#' )
#' head(conf_data)
#'
#' # Exemplo 2: Efeito dependente do tempo desde o tratamento
#' eff_list = list(list(function(t, s) 0.5 * (t - s)))
#' conf_data = generate_confounder_pointwise(
#'   n = 50, effects_function_list = eff_list, noise_sd = 0.2,
#'   p_treat = 0.3, time_grid = seq(0, 10, 0.5)
#' )
#'
#' # Exemplo 3: Com função de propensão e múltiplas covariáveis
#' eff_list = list(
#'   list(function(t, s) 0.3 * (t - s)),  # Efeito para X1
#'   list(function(t, s) 0.1 * (t - s)^2) # Efeito para X2
#' )
#' prop_func = function(X_list, t_index) {
#'   # Propensão depende de X1 no tempo anterior
#'   plogis(0.5 * X_list[[1]][, t_index-1])  
#' }
#' conf_data = generate_confounder_pointwise(
#'   n = 200, effects_function_list = eff_list, noise_sd = 0.15,
#'   propensity_function_list = list(prop_func), time_grid = 0:5
#' )
#' 
#' # Exemplo 4: Múltiplos tempos de tratamento
#' eff_list = list(list(
#'   function(t, s) 0.3,                 # Tratamento no tempo 1
#'   function(t, s) 0.5 * log(t - s + 1)  # Tratamento no tempo 2
#' ))
#' conf_data = generate_confounder_pointwise(
#'   n = 150, effects_function_list = eff_list, noise_sd = 0.2,
#'   p_treat = 0.4, time_grid = 0:10
#' )
#'
#' @export

generate_confounder_pointwise = function(
    n,
    effects_function_list,
    p_treat,
    x_noise_model = "iid",
    x_noise_rho = 0.6,
    noise_sd = 1,
    time_grid = NULL,
    baseline_function_list = NULL,
    propensity_function_list = NULL,
    return_format = c("long", "wide", "both")
) {
  # Validações de entrada
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) 
    stop("n deve ser um inteiro positivo")
  if (!is.numeric(noise_sd) || noise_sd <= 0) 
    stop("noise_sd deve ser um número positivo")
  if (!is.numeric(p_treat) || p_treat < 0 || p_treat > 1) 
    stop("p_treat deve estar entre 0 e 1")
  if (!is.list(effects_function_list) || length(effects_function_list) == 0) 
    stop("effects_function_list deve ser uma lista não vazia")
  if (!is.null(time_grid)) {
    if (!is.numeric(time_grid) || is.unsorted(time_grid)) 
      stop("time_grid deve ser um vetor numérico crescente")
  }
  if (!is.null(baseline_function_list) && 
      (length(baseline_function_list) != length(effects_function_list))) 
    stop("baseline_function_list e effects_function_list devem ter mesmo 
         tamanho.")
  if (!is.null(propensity_function_list)) {
    if (!is.list(propensity_function_list)) 
      stop("propensity_function_list deve ser uma lista")
    if (length(propensity_function_list) != length(effects_function_list[[1]])) 
      stop("propensity_function_list e effects_function_list[[1]] devem ter 
           mesmo tamanho.")
  }
  
  return_format = match.arg(return_format)
  
  num_covariates = length(effects_function_list)
  total_treatment_times = length(effects_function_list[[1]])
  if (is.null(time_grid)) time_grid = 0:(total_treatment_times + 1)
  total_times = length(time_grid)
  
  # Função para gerar ruído AR(1)
  gen_ar_noise = function(n, times, rho, sd) {
    n_times = length(times)
    noise = matrix(0, n, n_times)
    
    for (i in 1:n) {
      # Inicialização correta com variância teórica
      init_val = rnorm(1, 0, sd/sqrt(1 - rho^2))
      
      # Gera AR(1) usando filter()
      noise[i,] = c(
        init_val,
        stats::filter(
          x = rnorm(n_times - 1, 0, sd),
          filter = rho,
          method = "recursive",
          init = init_val
        )
      )
    }
    return(noise)
  }
  
  # --- Geração do ruído (parte modificada) ---
  Z_list = if (x_noise_model == "ar1") {
    lapply(1:num_covariates, function(k) {
      gen_ar_noise(n, time_grid, x_noise_rho, noise_sd)
    })
  } else {
    replicate(
      num_covariates,
      matrix(rnorm(n * total_times, 0, noise_sd), n, total_times),
      simplify = FALSE
    )
  }
  
  # Inicializa funções baseline
  if (is.null(baseline_function_list)) {
    baseline_function_list = replicate(
      num_covariates,
      function(t) 0,
      simplify = FALSE
    )
  }
  
  # Matriz de baseline
  baseline_matrix = vapply(
    baseline_function_list,
    function(f) f(time_grid),
    numeric(total_times)
  )
  
  # Inicialização de estruturas de dados
  A = A_control = matrix(0, n, total_times)
  X_list = control_X_list = lapply(1:num_covariates, function(k) 
    matrix(0, n, total_times))
  
  # Pré-computação de efeitos
  effects_mat_list = lapply(1:num_covariates, function(k) {
    effects_mat = matrix(0, total_treatment_times, total_times)
    for (s in 1:total_treatment_times) {
      for (t in 1:total_times) {
        if (time_grid[t] >= s) {
          effects_mat[s, t] = effects_function_list[[k]][[s]](time_grid[t], s)
        }
      }
    }
    effects_mat
  })
  
  # Valores iniciais (t=1)
  for (k in 1:num_covariates) {
    X_list[[k]][, 1] = control_X_list[[k]][, 1] = 
      baseline_matrix[1, k] + Z_list[[k]][, 1]
  }
  
  # Processo dinâmico (t >= 2)
  for (t_index in 2:total_times) {
    current_time = time_grid[t_index]
    
    # Decisão de tratamento
    if (current_time %in% 1:total_treatment_times) {
      s_index = as.integer(current_time)
      prop_func = if (!is.null(propensity_function_list)) 
        propensity_function_list[[s_index]] else NULL
      
      A[, t_index] = if (!is.null(prop_func)) {
        rbinom(n, 1, prop_func(X_list, t_index))
      } else {
        rbinom(n, 1, p_treat)
      }
    }
    
    # Atualização de covariáveis
    for (k in 1:num_covariates) {
      effects_mat = effects_mat_list[[k]]
      base_val = baseline_matrix[t_index, k]
      
      # Tempos de tratamento anteriores
      treated_cols = which(
        1:total_times < t_index &
          time_grid %in% 1:total_treatment_times
      )
      treated_cols = head(treated_cols, nrow(effects_mat))
      
      # Cálculo de efeitos cumulativos
      if (length(treated_cols) > 0) {
        effects_obs = rowSums(
          A[, treated_cols, drop = FALSE] %*%
            effects_mat[seq_along(treated_cols), t_index, drop = FALSE]
        )
        effects_ctrl = rowSums(
          A_control[, treated_cols, drop = FALSE] %*%
            effects_mat[seq_along(treated_cols), t_index, drop = FALSE]
        )
      } else {
        effects_obs = effects_ctrl = 0
      }
      
      # Atualiza com ruído pré-gerado
      X_list[[k]][, t_index] = base_val + effects_obs + Z_list[[k]][, t_index]
      control_X_list[[k]][, t_index] = base_val + effects_ctrl + 
        Z_list[[k]][, t_index]
    }
  }
  
  # Formatação de saída
  dt = data.table(
    id = rep(1:n, each = total_times),
    t = rep(time_grid, n),
    A = as.vector(t(A))
  )
  
  for (k in 1:num_covariates) {
    dt[, paste0("X", k) := as.vector(t(X_list[[k]]))]
    dt[, paste0("control_X", k) := as.vector(t(control_X_list[[k]]))]
  }
  
  # Retorno conforme formato solicitado
  switch(
    return_format,
    "long" = dt,
    "wide" = list(
      A_mat = A,
      A_control = A_control,
      X_list = X_list,
      control_X_list = control_X_list,
      t_vec = time_grid
    ),
    "both" = list(
      long = dt,
      wide = list(
        A_mat = A,
        A_control = A_control,
        X_list = X_list,
        control_X_list = control_X_list,
        t_vec = time_grid
      )
    )
  )
}
