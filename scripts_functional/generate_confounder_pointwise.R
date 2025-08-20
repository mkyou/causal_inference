generate_confounder_pointwise = function(
    n, effects_function_list, p_treat, noise_sd = 1,
    time_grid = NULL, baseline_function_list = NULL,
    propensity_function_list = NULL, a_seq = NULL,
    another_a_seq = NULL, return_format = c("long","wide","both")
) {
  if (!is.numeric(n) || n <= 0 || n != as.integer(n))
    stop("n deve ser um inteiro positivo")
  if (!is.numeric(noise_sd) || noise_sd <= 0)
    stop("noise_sd deve ser um número positivo")
  if (!is.numeric(p_treat) || p_treat < 0 || p_treat > 1)
    stop("p_treat deve estar entre 0 e 1")
  if (!is.list(effects_function_list) || length(effects_function_list)==0)
    stop("effects_function_list deve ser uma lista não vazia")
  if (!is.null(time_grid)) {
    if (!is.numeric(time_grid) || is.unsorted(time_grid))
      stop("time_grid deve ser vetor numérico crescente")
  }
  if (!is.null(baseline_function_list)) {
    if (!is.list(baseline_function_list))
      stop("baseline_function_list deve ser uma lista de listas")
    if (length(baseline_function_list) != length(effects_function_list))
      stop("baseline_function_list e effects_function_list devem ter 
           mesmo tamanho")
    if (any(sapply(baseline_function_list,length)!=n))
      stop("Cada baseline_function_list[[k]] deve ter n funções")
  }
  return_format = match.arg(return_format)
  num_covariates = length(effects_function_list)
  total_treatment_times = length(effects_function_list[[1]])
  if (is.null(time_grid)) time_grid = 0:(total_treatment_times+1)
  total_times = length(time_grid)
  if (is.null(a_seq)) a_seq = rep(1,total_treatment_times)
  if (is.null(another_a_seq)) another_a_seq = rep(0,total_treatment_times)
  if (length(a_seq) != total_treatment_times)
    stop("a_seq deve ter tamanho igual a total_treatment_times")
  if (length(another_a_seq) != total_treatment_times)
    stop("another_a_seq deve ter tamanho igual a total_treatment_times")
  
  # --- Ruído iid ---
  Z_list = replicate(
    num_covariates, matrix(rnorm(n*total_times,0,noise_sd), n, total_times),
    simplify=FALSE
  )
  
  A = matrix(0,n,total_times)
  A_a = matrix(0,n,total_times)
  A_another = matrix(0,n,total_times)
  X_list = lapply(1:num_covariates,function(k) matrix(0,n,total_times))
  a_X_list = lapply(1:num_covariates,function(k) matrix(0,n,total_times))
  another_a_X_list = lapply(
    1:num_covariates,
    function(k) matrix(0,n,total_times)
  )
  
  # Pré-computação de efeitos
  effects_mat_list = lapply(1:num_covariates,function(k){
    effects_mat = matrix(0,total_treatment_times,total_times)
    for (s in 1:total_treatment_times){
      for (t in 1:total_times){
        if (time_grid[t] >= s)
          effects_mat[s,t] = effects_function_list[[k]][[s]](time_grid[t],s)
      }
    }
    effects_mat
  })
  
  # t=1
  for (k in 1:num_covariates){
    for (i in 1:n){
      base_val = baseline_function_list[[k]][[i]](time_grid[1])
      X_list[[k]][i,1] = base_val + Z_list[[k]][i,1]
      a_X_list[[k]][i,1] = base_val + Z_list[[k]][i,1]
      another_a_X_list[[k]][i,1] = base_val + Z_list[[k]][i,1]
    }
  }
  
  # Processo dinâmico (t>=2)
  for (t_index in 2:total_times){
    current_time = time_grid[t_index]
    if (current_time %in% 1:total_treatment_times){
      s_index = as.integer(current_time)
      prop_func = if(!is.null(propensity_function_list))
        propensity_function_list[[s_index]] else NULL
      A[,t_index] = if(!is.null(prop_func))
        rbinom(n,1,prop_func(X_list,t_index)) else rbinom(n,1,p_treat)
      A_a[,t_index] = rep(a_seq[s_index],n)
      A_another[,t_index] = rep(another_a_seq[s_index],n)
    }
    for (k in 1:num_covariates){
      effects_mat = effects_mat_list[[k]]
      treated_cols_nat = which(
        1:total_times < t_index & 
          time_grid %in% 1:total_treatment_times
      )
      treated_cols_nat = head(treated_cols_nat,nrow(effects_mat))
      for (i in 1:n){
        base_val = baseline_function_list[[k]][[i]](current_time)
        
        effects_obs_nat = if(length(treated_cols_nat)>0)
          sum(
            A[i,treated_cols_nat] *
              effects_mat[seq_along(treated_cols_nat), t_index]
          ) else 0
        
        X_list[[k]][i,t_index] = base_val + effects_obs_nat + 
          Z_list[[k]][i,t_index]
        
        effects_a = if(length(treated_cols_nat)>0)
          sum(
            A_a[i,treated_cols_nat] *
              effects_mat[seq_along(treated_cols_nat),t_index]
          ) else 0
        
        a_X_list[[k]][i,t_index] = base_val + effects_a + 
          Z_list[[k]][i,t_index]
        effects_another = if(length(treated_cols_nat)>0)
          sum(
            A_another[i,treated_cols_nat] * 
              effects_mat[seq_along(treated_cols_nat), t_index]
          ) else 0
        
        another_a_X_list[[k]][i,t_index] = base_val + effects_another + 
          Z_list[[k]][i,t_index]
      }
    }
  }
  
  dt = data.table(
    id=rep(1:n,each=total_times),
    t=rep(time_grid,n),
    A=as.vector(t(A)),
    A_a=as.vector(t(A_a)),
    A_another=as.vector(t(A_another))
  )
  for (k in 1:num_covariates){
    dt[,paste0("X",k):=as.vector(t(X_list[[k]]))]
    dt[,paste0("a_X",k):=as.vector(t(a_X_list[[k]]))]
    dt[,paste0("another_a_X",k):=as.vector(t(another_a_X_list[[k]]))]
  }
  
  switch(
    return_format,
    "long" = dt,
    "wide" = list(
      A_mat = A, A_a_mat = A_a, A_another_mat = A_another,
      X_list = X_list, a_X_list = a_X_list, 
      another_a_X_list = another_a_X_list,
      t_vec = time_grid
      ),
    "both" = list(
      long = dt,
      wide = list(
        A_mat = A, A_a_mat = A_a, A_another_mat = A_another,
        X_list = X_list, a_X_list = a_X_list, 
        another_a_X_list = another_a_X_list,
        t_vec = time_grid
      )
    )
  )
}
    


