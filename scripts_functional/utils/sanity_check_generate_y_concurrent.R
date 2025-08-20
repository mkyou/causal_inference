plot_effects_raw_from_dt = function(
    y_out, 
    time_grid = y_out$wide$t_vec
) {
  require(ggplot2)
  require(data.table)
  
  # Extrair dados
  dt_long = as.data.table(y_out$long)
  true_effect_list = y_out$wide$true_effect_list
  
  dt_direct_all = rbindlist(
    lapply(seq_along(true_effect_list), function(s_idx) {
      data.table(
        t = time_grid,
        Effect = true_effect_list[[s_idx]],
        s_val = time_grid[s_idx]
      )
    }
    )
  )
  
  p1 = ggplot(dt_direct_all, aes(x = t, y = Effect, color = factor(s_val))) +
    geom_line(linewidth = 1) +
    labs(
      title = "Efeitos diretos por tempo de tratamento (s)",
      x = "Tempo",
      y = "Efeito direto teórico",
      color = "s"
    ) +
    theme_light()
  
  direct_sum = Reduce("+", true_effect_list)
  mean_total_effect = dt_long[
    ,
    .(Mean_Total = mean(effect_a_vs_another)),
    by = t
  ]
  
  dt_compare = data.table(
    t = time_grid,
    Value = c(direct_sum, mean_total_effect$Mean_Total),
    Type = rep(
      c("Soma dos diretos", "Total observado"), 
      each = length(time_grid)
    )
  )
  
  p2 = ggplot(dt_compare, aes(x = t, y = Value, linetype = Type)) +
    geom_line(color = "darkred", linewidth = 1) +
    labs(
      title = "Soma dos efeitos diretos vs efeito médio total",
      x = "Tempo",
      y = "Efeito médio",
      linetype = "Tipo"
    ) +
    theme_light()
  
  return(list(plot_diretos = p1, plot_soma_vs_total = p2))
}

plot_effects_fd = function(
    smooth_data,
    time_grid
) {
  require(ggplot2)
  require(data.table)
  require(fda)
  
  # ---- 1. Efeitos diretos ----
  direct_fd_list = smooth_data$direct_effect_fd_list
  if(!is.null(direct_fd_list)) {
    dt_direct_all = rbindlist(
      lapply(seq_along(direct_fd_list), function(s_idx) {
        data.table(
          t = time_grid,
          Effect = as.numeric(fda::eval.fd(time_grid, direct_fd_list[[s_idx]])),
          s_val = time_grid[s_idx]
        )
      })
    )
    
    p1 = ggplot(dt_direct_all, aes(x = t, y = Effect, color = factor(s_val))) +
      geom_line(linewidth = 1) +
      labs(
        title = "Efeitos diretos por tempo de tratamento (s)",
        x = "Tempo",
        y = "Efeito direto suavizado",
        color = "s"
      ) +
      theme_light()
  } else {
    p1 = NULL
  }
  
  # ---- 2. Comparação soma dos diretos vs efeito médio observado ----
  if(!is.null(direct_fd_list)) {
    # Soma dos efeitos diretos
    direct_sum = Reduce("+", lapply(direct_fd_list, function(fd) {
      as.numeric(fda::eval.fd(time_grid, fd))
    }))
  } else {
    direct_sum = rep(0, length(time_grid))
  }
  
  # Efeito médio observado
  effect_mean = rowMeans(fda::eval.fd(time_grid, smooth_data$effect_fd))
  
  dt_compare = data.table(
    t = time_grid,
    Value = c(direct_sum, effect_mean),
    Type = rep(c("Soma dos diretos", "Total observado"), each = length(time_grid))
  )
  
  p2 = ggplot(dt_compare, aes(x = t, y = Value, linetype = Type)) +
    geom_line(color = "darkred", linewidth = 1) +
    labs(
      title = "Soma dos efeitos diretos vs efeito médio total",
      x = "Tempo",
      y = "Efeito médio",
      linetype = "Tipo"
    ) +
    theme_light()
  
  return(list(plot_diretos = p1, plot_soma_vs_total = p2))
}
