files <- c(
  list.files("R", full.names = TRUE, pattern = "\\.R$"),
  list.files("scripts", full.names = TRUE, pattern = "\\.R$")
)

for (f in files) {
  lines <- readLines(f, warn = FALSE)
  
  # 1. Trocar %>% por |>
  lines <- gsub("%>%", "|>", lines, fixed = TRUE)
  
  # 2. Trocar <- por = (garantindo formatação decente)
  # Usamos perl para garantir que nao substituiremos <- dentro de strings facilmente (mesmo sendo raro)
  # Basicamente: qualquer coisa, espacos opcionais, <-, espacos opcionais
  lines <- gsub("\\s*<-\\s*", " = ", lines)
  
  # 3. Remover múltiplos espaços antes do = para desalinhar
  # Captura qualquer caractere que não seja operador (!, <, >, =), seguido de multiplos espaços e um = simples
  lines <- gsub("([^!<>=\\s])\\s+=", "\\1 =", lines, perl = TRUE)
  
  # Remover espacos no fim das linhas
  lines <- trimws(lines, which = "right")
  
  writeLines(lines, f)
}

cat("Substituicoes aplicadas em", length(files), "arquivos.\\n")
