library(data.table)

source("R/process_functional_data.R")

# Configuração do teste
N_TEST = 1000        
T_PTS = 256        
time_grid = seq(0, 1, length.out = T_PTS)

# Gera dados sintéticos rápidos (Senoide + Ruído)
set.seed(2025)
true_sig = sin(4 * pi * time_grid) * (1 + time_grid) 
noise_mat = matrix(
  rnorm(N_TEST * T_PTS, sd = 0.5), nrow = N_TEST, ncol = T_PTS
  )

wide_matrix = matrix(
  rep(true_sig, N_TEST), nrow = N_TEST, byrow = TRUE
  ) + noise_mat

cat(
  sprintf(
    "\n>>> INICIANDO BENCHMARK DE PRE-PROCESSAMENTO (N=%d, T=%d) <<<\n", 
    N_TEST, 
    T_PTS
    )
  )

cat("Comparando tempo de execução apenas da etapa de suavização.\n\n")

# 1. rLifting
cat("1. rLifting (Offline / C++ Core)...\n")
time_rlift = system.time({
  # Usando DB2 como padrão representativo
  res_rlift = smooth_data_matrix(
    time_grid, wide_matrix, 
    method = "rlifting_offline",
    scheme_name = "db2", 
    threshold_method = "semisoft"
  )
})

cat(sprintf("   Tempo: %.4f s\n", time_rlift["elapsed"]))

# 2. Wavethresh
cat("2. Wavethresh (C / DWT Tradicional)...\n")
time_wave = system.time({
  res_wave = smooth_data_matrix(
    time_grid, wide_matrix, 
    method = "wavethresh",
    family = "DaubExPhase"
  )
})

cat(sprintf("   Tempo: %.4f s\n", time_wave["elapsed"]))

# 3. AdLift
cat("3. AdLift (Adaptive Lifting - R Implementation)...\n")
cat("   (Aguarde... isso pode demorar alguns minutos)\n")
time_adlift = system.time({
  res_adlift = smooth_data_matrix(
    time_grid, wide_matrix, 
    method = "adlift"
  )
})
cat(sprintf("   Tempo: %.6f s\n", time_adlift["elapsed"]))

# Consolidação dos Resultados
results = data.table(
  Metodo = c("rLifting", "Wavethresh", "AdLift"),
  Tempo_Total_s = c(
    time_rlift["elapsed"], time_wave["elapsed"], 
    time_adlift["elapsed"]
    )
)

# Métricas derivadas
results[, Tempo_Por_Curva_ms := (Tempo_Total_s / N_TEST) * 1000]
results[, Fator_Lentidao := Tempo_Total_s / min(Tempo_Total_s)]
results[, Curvas_Por_Segundo := N_TEST / Tempo_Total_s]

cat("\n=== TABELA DE PERFORMANCE (SUAVIZAÇÃO) ===\n")
print(results)

# Salva para usar na dissertação
dir.create("output", showWarnings = FALSE)
fwrite(results, "output/adlift_viability_benchmark.csv")

# ---------------------------------------------------------------------
# Formatando artefato e salvando-o

library(data.table)
library(dplyr)
library(kableExtra)

# Carregar dados
file_path = "output/adlift_viability_benchmark.csv"
if(!file.exists(file_path)) stop("Rode o script 00_adlift_killer.R primeiro.")
res = fread(file_path)

# Formatação para publicação
table_pretty = res |>
  mutate(
    # Formata números
    Tempo_Por_Curva = sprintf("%.2f ms", Tempo_Por_Curva_ms),
    # Adiciona 'x' ao fator e trata Infinito
    Fator = ifelse(
      is.infinite(Fator_Lentidao), "> 1000x", 
      sprintf("%.1fx", Fator_Lentidao)
      ),
    # Formata throughput com separador de milhar
    Throughput = scales::comma(
      Curvas_Por_Segundo, accuracy = 1, suffix = " curves/s"
      )
  ) |>
  dplyr::select(
    `Método` = Metodo,
    `Tempo / Curva` = Tempo_Por_Curva,
    `Fator de Lentidão` = Fator,
    `Capacidade (Throughput)` = Throughput
  )

# Gerar tabela kable
kbl = kable(
  table_pretty, format = "html", booktabs = TRUE, align = "lccc",
  caption = "Comparação de performance computacional (suavização)"
  ) |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"), 
    full_width = FALSE
    ) |>
  row_spec(
    which(res$Metodo == "AdLift"), 
    color = "firebrick"
    ) |> # Destaca o problema
  row_spec(
    which(res$Metodo == "rLifting"), 
    color = "royalblue"
    )

# Salvar
save_kable(kbl, "output/table_adlift_pretty.html")
print(table_pretty)
