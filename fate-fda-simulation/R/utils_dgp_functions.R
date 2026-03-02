normalize_time = function(t, t_min, t_max) (t - t_min) / (t_max - t_min)

#' Sigmoide (Suave - Ecker)
effect_sigmoid = function(
    t, s_time, t_min = 20,
    t_max = 63, magnitude = 1, slope = 20
    ) {
  t_norm = normalize_time(t, t_min, t_max)
  s_norm = normalize_time(s_time, t_min, t_max)
  val = 1 / (1 + exp(-slope * (t_norm - s_norm)))
  val[t < (s_time - 1e-6)] = 0
  return(val * magnitude)
}

#' Doppler (Alta Frequência)
effect_doppler = function(
    t, s_time, t_min = 20,
    t_max = 63, magnitude = 1
    ) {
  val = numeric(length(t))
  idx = which(t >= s_time)
  if (length(idx) > 0) {
    tn = normalize_time(t[idx], t_min, t_max)
    tn = pmax(tn, 1e-6)
    # Frequência base em 6.0
    val[idx] = sqrt(tn * (1 - tn)) * sin((6.0 * pi) / (tn + 0.05))
  }
  return(val * magnitude)
}

#' Bumps (picos finos)
effect_bumps = function(t, s_time, t_min = 20, t_max = 63, magnitude = 1) {
  val = numeric(length(t))
  idx = which(t >= s_time)
  if (length(idx) > 0) {
    tn = normalize_time(t[idx], t_min, t_max)
    pos = c(0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81)
    h = c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 5.1, -4.2)
    w = c(
      0.005, 0.005, 0.006, 0.01, 0.01, 0.03,
      0.01, 0.01, 0.005, 0.008, 0.005
      ) * 0.5

    sum_b = numeric(length(tn))
    for (j in seq_along(pos))
      sum_b = sum_b + h[j] * (1 + abs((tn - pos[j])/w[j])^4)^(-1)
    val[idx] = sum_b
  }
  return(val * magnitude)
}

#' Heavisine (Degraus Puros)
effect_heavisine = function(t, s_time, t_min = 20, t_max = 63, magnitude = 1) {
  val = numeric(length(t))
  idx = which(t >= s_time)

  if (length(idx) > 0) {
    tn = normalize_time(t[idx], t_min, t_max)
    # Função Escada (Step Function)
    # 4 degraus em 0.2, 0.4, 0.6, 0.8
    steps = (tn > 0.2) + (tn > 0.4) + (tn > 0.6) - (tn > 0.8)*2
    val[idx] = steps
  }

  return(val * magnitude)
}
