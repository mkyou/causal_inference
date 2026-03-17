library(data.table)
library(data.table)

dgp_raw = readRDS("functional_simulations/logs/CrescimentoHumano.rds")

# Encontre o índice de tempo exato para o tratamento e o passo anterior
treatment_time = 2.0
time_grid_vector = sort(unique(dgp_raw$t))
t_trat_idx = which(time_grid_vector == treatment_time)
t_antes_idx = t_trat_idx - 1

# Use os valores exatos da grade de tempo para a filtragem
t_trat_val = time_grid_vector[t_trat_idx]
t_antes_val = time_grid_vector[t_antes_idx]

# Pegue os dados relevantes: covariáveis em t_antes e tratamento em t_trat
dados_antes <- dgp_raw[t == t_antes_val, .(id, X1, X2, X3, X4)]
dados_trat <- dgp_raw[t == t_trat_val, .(id, A)]

inspecao_df <- merge(dados_antes, dados_trat, by = "id")


# AQUI ESTÁ O TESTE CRÍTICO:
# Se a atribuição em t=2 for verdadeiramente aleatória, não deve haver
# correlação entre as covariáveis em t=1.9 e o tratamento A em t=2.
# Qualquer correlação forte aqui é a prova definitiva do vazamento.
cor(inspecao_df$X1, inspecao_df$A)
cor(inspecao_df$X2, inspecao_df$A)
# ... e assim por diante para X3, X4.

# Um glm simples também revelará a associação
summary(glm(A ~ X1 + X2 + X3 + X4, data = inspecao_df, family = "binomial"))
