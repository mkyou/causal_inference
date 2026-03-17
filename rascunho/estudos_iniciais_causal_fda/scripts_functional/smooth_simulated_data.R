smooth_simulated_data = function(
    dt_long,
    nbasis = NULL,
    lambda = 1e-6,
    gcv_threshold = 1.5,
    lambda_max = 1e2,
    true_effect_list = NULL
) {
  requireNamespace("fda")
  
  # Time grid
  time_grid = sort(unique(dt_long$t))
  time_grid = time_grid[is.finite(time_grid)]
  if(length(time_grid) < 2) stop("time_grid vazio ou com menos de 2 pontos válidos")
  
  n = length(unique(dt_long$id))
  
  if(is.null(nbasis)) nbasis = max(7, min(40, ceiling(length(time_grid)/4)))
  
  basis = fda::create.bspline.basis(rangeval = range(time_grid), nbasis = nbasis)
  
  # Suavização adaptativa com GCV por indivíduo
  smooth_adaptive = function(data_mat, var_name = "", lambda_init = lambda) {
    l = lambda_init
    repeat {
      fdParobj = fda::fdPar(basis, Lfdobj=fda::int2Lfd(2), lambda=l)
      smooth = fda::smooth.basis(argvals=time_grid, y=t(data_mat), fdParobj=fdParobj)
      mean_gcv = mean(smooth$gcv)
      if(mean_gcv <= gcv_threshold || l >= lambda_max) break
      l = l*10
    }
    if(mean_gcv > gcv_threshold) {
      warning("GCV ainda alto (", round(mean_gcv,2), ") para ", var_name)
    }
    list(fd=smooth$fd, gcv=smooth$gcv, lambda_used=l)
  }
  
  # Matrizes y
  y_mat = matrix(dt_long$y, nrow=n, byrow=TRUE)
  a_y_mat = matrix(dt_long$a_y, nrow=n, byrow=TRUE)
  another_a_y_mat = matrix(dt_long$another_a_y, nrow=n, byrow=TRUE)
  effect_mat = matrix(dt_long$effect_a_vs_another, nrow=n, byrow=TRUE)
  
  # Suavizar
  y_result = smooth_adaptive(y_mat, "y observado")
  a_y_result = smooth_adaptive(a_y_mat, "y sob a_seq")
  another_a_y_result = smooth_adaptive(another_a_y_mat, "y sob another_a_seq")
  effect_result = smooth_adaptive(effect_mat, "efeito causal")
  
  # Covariáveis
  covariate_names = grep("^X\\d+", names(dt_long), value=TRUE)
  X_results = lapply(covariate_names, function(v) {
    smooth_adaptive(matrix(dt_long[[v]], nrow=n, byrow=TRUE), v)
  })
  
  # Efeitos diretos se fornecidos
  direct_effect_results = NULL
  if(!is.null(true_effect_list)) {
    direct_effect_results = lapply(seq_along(true_effect_list), function(s_idx) {
      mat_s = matrix(rep(true_effect_list[[s_idx]], each=n), nrow=n)
      smooth_adaptive(mat_s, paste0("efeito direto s=", s_idx))
    })
  }
  
  list(
    nbasis_used = nbasis,
    y_fd = y_result$fd, y_gcv = y_result$gcv, y_lambda = y_result$lambda_used,
    a_y_fd = a_y_result$fd, a_y_gcv = a_y_result$gcv, a_y_lambda = a_y_result$lambda_used,
    another_a_y_fd = another_a_y_result$fd, another_a_y_gcv = another_a_y_result$gcv,
    another_a_y_lambda = another_a_y_result$lambda_used,
    effect_fd = effect_result$fd, effect_gcv = effect_result$gcv, effect_lambda = effect_result$lambda_used,
    X_fd_list = lapply(X_results, function(x) x$fd),
    X_gcv_list = lapply(X_results, function(x) x$gcv),
    X_lambda_list = lapply(X_results, function(x) x$lambda_used),
    direct_effect_fd_list = if(!is.null(direct_effect_results)) lapply(direct_effect_results, function(x) x$fd) else NULL,
    direct_effect_gcv_list = if(!is.null(direct_effect_results)) lapply(direct_effect_results, function(x) x$gcv) else NULL,
    direct_effect_lambda_list = if(!is.null(direct_effect_results)) lapply(direct_effect_results, function(x) x$lambda_used) else NULL,
    A_mat = matrix(dt_long$A, nrow=n, ncol=length(time_grid), byrow=TRUE)
  )
}
