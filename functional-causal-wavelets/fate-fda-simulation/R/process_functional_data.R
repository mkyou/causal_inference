library(data.table)
library(fda)
library(wavethresh)
library(rLifting)
library(adlift)

#' Suavizador Genérico (Wrapper)
#' @param time_grid Vetor de tempo.
#' @param wide_matrix Matriz (IDs x Tempo) com dados ruidosos.
#' @param method String: "spline", "wavethresh", "rlifting_offline", "adlift".
#' @param ... Argumentos extras para os métodos.
smooth_data_matrix = function(time_grid, wide_matrix, method = "spline", ...) {

  start_time = Sys.time()

  # Bloco de proteção contra falhas
  res = tryCatch({
    switch(
      method,
      "spline" = smooth_spline_robust(time_grid, wide_matrix, ...),
      "wavethresh" = smooth_wavethresh_wrapper(wide_matrix, ...),
      "rlifting_offline" = smooth_rlifting_offline_wrapper(wide_matrix, ...),
      "adlift" = smooth_adlift_wrapper(time_grid, wide_matrix, ...),
      stop("Método desconhecido: ", method)
    )
  }, error = function(e) {
    message(sprintf("Erro no método '%s': %s", method, e$message))
    return(NULL)
  })

  end_time = Sys.time()
  if (!is.null(res)) {
    attr(res, "exec_time") = as.numeric(end_time - start_time, units = "secs")
  }
  return(res)
}

# Splines
smooth_spline_robust = function(time_grid, wide_matrix, nbasis = 50) {
  # Cria base B-Spline
  bspline_basis = create.bspline.basis(
    rangeval = range(time_grid), nbasis = nbasis
    )

  # Busca lambda ótimo numa grade logarítmica (rápido e eficaz)
  log_lambdas = seq(-10, 1, length.out = 20)
  gcvs = numeric(length(log_lambdas))

  # Testamos lambdas em uma subamostra se N for muito grande para ganhar tempo
  sample_indices = if(nrow(wide_matrix) > 50)
    sample(nrow(wide_matrix), 50) else 1:nrow(wide_matrix)
  test_data = t(wide_matrix[sample_indices, , drop = FALSE])

  for(i in seq_along(log_lambdas)) {
    fdPar_obj = fdPar(bspline_basis, Lfdobj = 2, lambda = 10^log_lambdas[i])
    # smooth.basis retorna objeto com GCV
    fd_try = smooth.basis(time_grid, test_data, fdPar_obj)
    gcvs[i] = sum(fd_try$gcv) # Soma dos GCVs das curvas
  }

  best_lambda = 10^log_lambdas[which.min(gcvs)]

  # Ajuste final com o melhor lambda
  final_fdPar = fdPar(bspline_basis, Lfdobj = 2, lambda = best_lambda)
  final_fd = smooth.basis(time_grid, t(wide_matrix), final_fdPar)$fd

  return(t(eval.fd(time_grid, final_fd)))
}

# Wavethresh (CRAN - DWT Tradicional)
smooth_wavethresh_wrapper = function(
    wide_matrix, filter.number = 4,
    family = "DaubExPhase"
    ) {

  n_rows = nrow(wide_matrix)
  n_cols = ncol(wide_matrix)
  smoothed = matrix(0, n_rows, n_cols)

  log2_n = log2(n_cols)
  target_len = 2^ceiling(log2_n)
  pad_len = target_len - n_cols

  if(family == "Haar") {
    filter.number = 1
    family = "DaubExPhase"
  }

  for(i in 1:n_rows) {
    sig = wide_matrix[i, ]
    # Padding constante
    if (pad_len > 0)
      sig_pad = c(sig, rep(tail(sig, 1), pad_len))
    else sig_pad = sig

    # Proteção contra falha do wd()
    tryCatch({
      w = wavethresh::wd(
        sig_pad, filter.number = filter.number,
        family = family, bc = "symmetric"
        )
      w_t = wavethresh::threshold(
        w, policy = "sure",
        type = "soft", dev = mad,
        verbose = FALSE
        )
      rec = wavethresh::wr(w_t)

      # Truncamento seguro
      if(length(rec) >= n_cols) {
        smoothed[i, ] = rec[1:n_cols]
      } else {
        # Fallback: se algo bizarro acontecer,
        # retorna original suavizado com média móvel simples ou original
        smoothed[i, ] = sig
      }
    }, error = function(e) {
      smoothed[i, ] = sig # Fallback silencioso em erro
    })
  }
  return(smoothed)
}

# rLifting
smooth_rlifting_offline_wrapper = function(
    wide_matrix, scheme_name = "db2",
    alpha = 0.3, beta = 1.2, threshold_method = "semisoft"
    ) {
  scheme = rLifting::lifting_scheme(scheme_name)
  n_cols = ncol(wide_matrix)
  max_levels = floor(log2(n_cols))

  smoothed = apply(wide_matrix, 1, function(x) {
    rLifting::denoise_signal_offline(
      x, scheme, levels = max_levels,
      method = threshold_method,
      extension = "symmetric",
      alpha = alpha,
      beta = beta
    )
  })
  return(t(smoothed))
}

# AdLift
smooth_adlift_wrapper = function(time_grid, wide_matrix, ...) {
  n_rows = nrow(wide_matrix)
  n_cols = ncol(wide_matrix)
  smoothed = matrix(0, n_rows, n_cols)

  for(i in 1:n_rows) {
    try({
      sig = wide_matrix[i, ]

      suppressMessages({
        decomp = adlift::fwtnp(
          time_grid, sig, LocalPred = adlift::LinearPred
          )

        denoised = adlift::denoise.fwtnp(decomp, do.plot=FALSE)
        rec = adlift::invtnp(denoised)
      })

      if(is.list(rec) && !is.null(rec$x)) {
        val_rec = as.numeric(rec$x)
      } else {
        val_rec = as.numeric(rec)
      }

      if(length(val_rec) == n_cols) {
        smoothed[i, ] = val_rec
      }
    }, silent = TRUE)
  }
  return(smoothed)
}
