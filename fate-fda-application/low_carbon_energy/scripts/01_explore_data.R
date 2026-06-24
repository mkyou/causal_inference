# 01_explore_data.R
#
# Exploração e documentação estrutural do dataset Low Carbon London (LCL).
# Objetivo: entender todos os arquivos disponíveis, mapear variáveis fixas e
# variantes no tempo, verificar empiricamente a aleatorização do trial, e
# registrar as decisões de pré-processamento que orientarão os scripts
# seguintes.
#
# Fonte:  Kaggle "Smart meters in London" (jeanmidev/smart-meters-in-london)
#         https://www.kaggle.com/datasets/jeanmidev/smart-meters-in-london
# Trial:  Low Carbon London, UK Power Networks / EDF Energy, 2011–2014.
# Documentação: https://data.london.gov.uk/dataset/smartmeter-energy-use-data-in-london-households
#
# Estrutura do diretório data/raw/ após extração:
#   halfhourly_dataset/block_0.csv … block_111.csv — séries longas (LCLid, tstp, kWh/hh)
#   hhblock_dataset/hhblock_dataset/block_*.csv — formato largo (LCLid, day, hh_0…hh_47)
#   daily_dataset/daily_dataset/block_*.csv — agregados diários por domicílio
#   daily_dataset.csv — idem, arquivo único (~336 MB)
#   informations_households.csv — metadados dos domicílios (1 linha/domicílio)
#   weather_hourly_darksky.csv — clima horário, Londres, 2011–2014
#   weather_daily_darksky.csv — clima diário, Londres, 2011–2014
#   acorn_details.csv — perfil demográfico por categoria ACORN
#   uk_bank_holidays.csv — feriados bancários do Reino Unido
#   darksky_parameters_documentation.html — documentação das variáveis meteorológicas
#
# Dependências: install.packages(c("data.table", "here"))

library(data.table)
library(here)

DATA_RAW = here("low_carbon_energy", "data", "raw")
HH_DIR = file.path(DATA_RAW, "halfhourly_dataset")
HHBLK_DIR = file.path(DATA_RAW, "hhblock_dataset", "hhblock_dataset")
DAILY_DIR = file.path(DATA_RAW, "daily_dataset", "daily_dataset")
HH_INFO = file.path(DATA_RAW, "informations_households.csv")
WEATHER_H = file.path(DATA_RAW, "weather_hourly_darksky.csv")
WEATHER_D = file.path(DATA_RAW, "weather_daily_darksky.csv")
ACORN_DET = file.path(DATA_RAW, "acorn_details.csv")
HOLIDAYS = file.path(DATA_RAW, "uk_bank_holidays.csv")

if (!dir.exists(HH_DIR) || !file.exists(HH_INFO)) {
  stop(
    "Dados não encontrados em ", DATA_RAW, "\n",
    "Baixe de https://www.kaggle.com/datasets/jeanmidev/smart-meters-in-london\n",
    "e extraia em fate-fda-application/low_carbon_energy/data/raw/",
    call. = FALSE
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# ESTRUTURA CAUSAL DO PROBLEMA
# ══════════════════════════════════════════════════════════════════════════════
#
# Unidade amostral : domicílio i = 1 … n
#
# Tratamento (fixo, pré-atribuído):
#   A_i ∈ {0, 1} — 0 = tarifa padrão (Std), 1 = tarifa dinâmica (ToU)
#   Fonte: informations_households.csv, coluna stdorToU
#
# Desfecho funcional:
#   Y_i(t) — consumo de energia em kWh por meia hora no instante t
#   Fonte: halfhourly_dataset/ (longo) ou hhblock_dataset/ (largo)
#   O formato largo organiza cada curva diária como 1 linha × 48 colunas
#   (hh_0…hh_47), que é a estrutura natural para a análise FDA.
#
# Covariáveis fixas (pré-tratamento):
#   ACORN_grouped — categoria socioeconômica do CEP do domicílio, derivada do
#   Censo 2011 (CACI Ltd.). É pré-tratamento por construção: determinada pelo
#   CEP antes e independentemente do trial. Não muda ao longo do estudo.
#   A ser usada como covariável para reduzir a variância da estimativa do FATE
#   e controlar pelo desequilíbrio observado na randomização.
#   Fonte: informations_households.csv, coluna Acorn_grouped
#
# Variáveis contextuais variantes no tempo (nível Londres, não individuais):
#   Temperatura, umidade, vento, visibilidade, ponto de orvalho, pressão.
#   Fonte: weather_hourly_darksky.csv
#
#   Estas variáveis são W(t) — variam com t, não com i. O modelo de Ecker
#   aceita X_i (vetor fixo por domicílio), não covariáveis funcionais W(t).
#   Incluí-las exigiria extensão do modelo fora do escopo desta dissertação.
#   A omissão não compromete identificabilidade: como A_i é randomizado,
#   W(t) é independente de A_i por construção — sem viés, apenas variância
#   residual ligeiramente maior. Feriados (uk_bank_holidays.csv) têm o mesmo
#   estatuto: variam com t, não entram como X_i.
#
# Variáveis auxiliares:
#   uk_bank_holidays.csv — feriados bancários; não-útil para o que pretendemos.
#   acorn_details.csv — índices demográficos do Censo 2011 por categoria
#   ACORN; referência para interpretação, não entra no modelo diretamente.
#
# Estimando:
#   FATE(t) = E[Y_i^(1)(t) − Y_i^(0)(t)]
#   Efeito médio causal da tarifa dinâmica sobre a curva de consumo em t.

# ══════════════════════════════════════════════════════════════════════════════
# 1. METADADOS DOS DOMICÍLIOS (informations_households.csv)
# ══════════════════════════════════════════════════════════════════════════════

hh = fread(HH_INFO)

# Uma linha por domicílio — confirma que todas as colunas são fixas no tempo.
stopifnot(uniqueN(hh$LCLid) == nrow(hh))

print(head(hh, 5))

# Distribuição do tratamento: ~20% ToU, ~80% Std.
# O desequilíbrio (1:4) é intencional no desenho do trial — não indica problema.
print(hh[, .N, keyby = stdorToU])

# Distribuição por ACORN agrupado.
# Adversity (16% ToU) vs Affluent/Comfortable (~21–22% ToU): desequilíbrio
# documentado na seção 5 (verificação de randomização).
print(hh[, .N, keyby = Acorn_grouped])

# ══════════════════════════════════════════════════════════════════════════════
# 2. SÉRIE TEMPORAL SEMIHORORÁRIA — formato longo (halfhourly_dataset/)
# ══════════════════════════════════════════════════════════════════════════════
# Colunas: LCLid | tstp | energy(kWh/hh)
# 112 arquivos block_*.csv, ~50 domicílios cada, ~1.2M linhas por bloco.
# energy(kWh/hh) contém a string "Null" onde faltam dados — a coluna é lida
# como character; a conversão para numeric produz NA nesses slots.
# Valores negativos são possíveis em domicílios com microgeração solar.

block_files_hh = sort(
  list.files(
    HH_DIR, 
    pattern = "^block_\\d+\\.csv$", 
    full.names = TRUE
  )
)

b0 = fread(block_files_hh[1])
b0[, tstp   := as.POSIXct(tstp, tz = "Europe/London")]
b0[, energy := suppressWarnings(as.numeric(get("energy(kWh/hh)")))]

print(head(b0[!is.na(tstp), .(LCLid, tstp, energy)], 6))
print(summary(b0$energy))

# ══════════════════════════════════════════════════════════════════════════════
# 3. SÉRIE TEMPORAL — formato largo (hhblock_dataset/)
# ══════════════════════════════════════════════════════════════════════════════
# Colunas: LCLid | day | hh_0 … hh_47
# Cada linha é uma curva diária de um domicílio em formato largo:
#   hh_0 = consumo em 00:00–00:30, hh_1 = 00:30–01:00, …, hh_47 = 23:30–00:00.
# Resolução: 30 min (48 pontos por dia) — idêntica ao halfhourly_dataset,
# mas pivotada para que cada linha seja o vetor Y_i(t) de um dia inteiro.
#
# Este formato define t ∈ [0, 24h]: o FATE(t) seria o efeito médio no instante
# t do dia típico, e "em tempo real" significaria atualizar a estimativa do dia
# corrente à medida que chegam as 48 leituras daquele dia.
#
# Não é o formato que usaremos. A análise trabalhará com a série contínua
# (halfhourly_dataset), onde t percorre o calendário inteiro (jan–dez 2013,
# ~17.520 pontos por domicílio). Nesse caso o FATE(t) é o efeito no instante t
# do calendário, o lifting incremental processa cada nova leitura de 30 em 30
# min, e o CFATE acumula ao longo do ano — que é o cenário descrito na
# dissertação. O hhblock fica disponível como referência e para eventuais
# análises de padrão intradiário.

block_files_hh2 = sort(
  list.files(
    HHBLK_DIR,
    pattern = "^block_\\d+\\.csv$",
    full.names = TRUE
  )
)

b0w = fread(block_files_hh2[1])
b0w[, day := as.Date(day)]

# Inspeciona dimensões e primeiras curvas.
print(dim(b0w)) # linhas = (domicílios × dias)
print(names(b0w)[1:6]) # LCLid, day, hh_0, hh_1, …
print(b0w[1:3, 1:8, with = FALSE]) # primeiras 3 curvas, 8 colunas

# Número de dias e domicílios no primeiro bloco.
print(b0w[, .(n_dias = uniqueN(day), n_hh = uniqueN(LCLid))])

# ══════════════════════════════════════════════════════════════════════════════
# 4. AGREGADOS DIÁRIOS (daily_dataset/)
# ══════════════════════════════════════════════════════════════════════════════
# Colunas: LCLid | day | energy_median | energy_mean | energy_max |
#   energy_count | energy_std | energy_sum | energy_min
# Resumo diário por domicílio. Útil para EDA e visualização, mas não é o
# dado primário para a análise FDA (que trabalha com as curvas horárias).

block_files_daily = sort(
  list.files(
    DAILY_DIR, pattern = "^block_\\d+\\.csv$", full.names = TRUE
  )
)

bd = fread(block_files_daily[1])
bd[, day := as.Date(day)]

print(head(bd, 4))

# ══════════════════════════════════════════════════════════════════════════════
# 5. DADOS METEOROLÓGICOS
# ══════════════════════════════════════════════════════════════════════════════
# weather_hourly_darksky.csv: clima horário de Londres, nov/2011–mar/2014.
# Fonte: Dark Sky API (time machine requests). Documentação das variáveis
# em darksky_parameters_documentation.html.
#
# Variáveis numéricas disponíveis: temperature (°C), apparentTemperature,
# dewPoint, humidity (0–1), pressure (hPa), windSpeed (m/s), windBearing (°),
# visibility (km), cloudCover (0–1 — ausente neste dataset).
# Variáveis categóricas: precipType, icon, summary.
#
# Resolução horária vs. desfecho semihororário: para usar como covariável,
# a série precisa ser interpolada ou agregada para 30 min.
# weather_daily_darksky.csv oferece a versão diária, útil para modelos que
# operam na escala do dia (como o formato largo do hhblock).
#
# Não iremos usar.
# O modelo proposto não contempla covariáveis variantes no tempo.
# Ademais, essas covariáveis não são confundidoras.
# O tratamento é randomizado e as variáveis metereológicas são comuns a todos.

wh = fread(WEATHER_H)
wh[, time := as.POSIXct(time, tz = "UTC")]

print(dim(wh))
print(range(wh$time))
print(names(wh))

num_cols = c("temperature", "humidity", "windSpeed", "visibility", "pressure")
print(
  wh[
    ,
    lapply(
      .SD,
      function(x) 
        round(
          c(
            min = min(x, na.rm = TRUE), 
            median = median(x, na.rm = TRUE),
            max = max(x, na.rm = TRUE)
          ), 2
        )
    ),
    .SDcols = num_cols
  ]
)

# ══════════════════════════════════════════════════════════════════════════════
# 6. FERIADOS E DETALHES ACORN
# ══════════════════════════════════════════════════════════════════════════════

holidays = fread(HOLIDAYS)
holidays[, `Bank holidays` := as.Date(`Bank holidays`)]
# 9 feriados bancários no período 2012–2013; relevantes para filtrar dias
# atípicos no pré-processamento ou adicionar como indicador binário.
# para o que vamos fazer, não é tão útil.
print(holidays)

acorn_det = fread(ACORN_DET)
# acorn_details contém índices demográficos do Censo 2011 para cada categoria
# ACORN (A–Q): distribuição etária, composição domiciliar, ocupação, etc.
# Não entra no modelo como covariável — serve de referência para interpretação
# dos grupos ACORN na análise de heterogeneidade de efeito.
print(dim(acorn_det))
print(acorn_det[1:3, 1:6, with = FALSE])

# ══════════════════════════════════════════════════════════════════════════════
# 7. REGULARIDADE DA GRADE TEMPORAL (todos os blocos halfhourly)
# ══════════════════════════════════════════════════════════════════════════════
# Esta seção lê todos os 112 blocos (~1–2 min) para diagnosticar a regularidade
# da grade de 30 min em cada domicílio.

check_regularity = function(fpath) {
  b = fread(fpath, select = c("LCLid", "tstp"))
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b = b[!is.na(tstp)][order(LCLid, tstp)]
  b[
    ,
    diff_min := as.numeric(
      difftime(tstp, shift(tstp), units = "mins")
    ),
    by = LCLid
  ]
  b[!is.na(diff_min), .(
    n_leituras = .N + 1L,
    pct_30min = mean(diff_min == 30),
    n_non30 = sum(diff_min != 30),
    irregular = any(diff_min != 30)
  ), by = LCLid]
}

reg = rbindlist(lapply(block_files_hh, check_regularity))
reg = merge(reg, hh[, .(LCLid, stdorToU)], by = "LCLid", all.x = TRUE)

# 100% dos domicílios apresentam ao menos uma leitura fora dos 30 min,
# mas a mediana de pct_30min é 99.98%. Duas fontes de irregularidade:
# (a) Transições BST/GMT: geram intervalos ≈7 e ≈23 min (que somam 30 min)
#     no limite de cada transição — verificado empiricamente na seção abaixo.
# (b) Falhas operacionais esparsas: leituras inteiramente ausentes em instantes
#     aleatórios ao longo do ano, criando gaps de 60–150 min. A maioria dos
#     domicílios perde 1–20 leituras no ano (mediana 2). Verificado e
#     documentado na Seção 11 com diagnóstico separado.
# Único caso estruturalmente irregular: MAC005565 (pct_30min = 97.7%, n=45
# leituras totais — ~22 horas de dados). Série irrecuperavelmente curta.
# Decisão: excluir MAC005565. Com apenas 45 leituras, esse domicílio não tem
# cobertura na janela 2013 e seria removido pelo filtro de ids_2013 de qualquer
# forma. A exclusão não compromete identificabilidade: as condições causais
# (SUTVA, ignorabilidade, positividade) são propriedades da distribuição
# populacional — remover uma unidade da amostra não as viola.
#
# Sobre as transições DST (31/mar e 27/out de 2013):
# - Primavera (adiantamento): gap de 60 min entre duas leituras consecutivas.
#   Nenhuma ação necessária. O Haar em grade irregular no lifting scheme adapta
#   o coeficiente de update de 1/2 para (t_{j+1}-t_j)/(t_{j+2}-t_j), que
#   recupera 1/2 na grade regular e ajusta suavemente fora dela. rLifting
#   recebe os timestamps reais e trata o gap nativamente — em ambos os modos.
# - Outono (atraso): aparentemente geraria timestamps "duplicados" na janela
#   01:00–02:00, mas o R com tz = "Europe/London" distingue corretamente BST
#   de GMT — os dois registros das 01:xx produzem POSIXct distintos (um em
#   BST = UTC+1, outro em GMT = UTC+0). Verificado abaixo: 50 leituras no dia
#   27/10 com zero duplicatas. Nenhuma ação necessária.
#
# Verificação DST outono (27/10/2013):
b_dst = fread(block_files_hh[1], select = c("LCLid", "tstp"))
b_dst = b_dst[LCLid %in% head(ids_2013, 1)]
b_dst[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
b_dst_oct = b_dst[as.Date(tstp, tz = "Europe/London") == "2013-10-27"]
cat(
  sprintf(
  "Leituras em 27/out para %s: %d | duplicatas: %d\n",
  b_dst_oct[, unique(LCLid)], nrow(b_dst_oct), sum(duplicated(b_dst_oct$tstp))
  )
)

# Período ambíguo 00:30–02:30 (exibe os rótulos BST/GMT confirmando distinção):
print(
  b_dst_oct[
    tstp >= as.POSIXct("2013-10-27 00:30", tz = "Europe/London") &
      tstp <= as.POSIXct("2013-10-27 02:30", tz = "Europe/London"),
    tstp
  ]
)
print(quantile(reg$pct_30min, c(0, .01, .05, .25, .5, .75, .95, .99, 1)))
print(reg[pct_30min < 0.99, .(LCLid, stdorToU, n_leituras, pct_30min, n_non30)])

# Verificação crítica para identificabilidade: a irregularidade é independente
# do tratamento? Se domicílios ToU tivessem grades sistematicamente mais
# irregulares, a representação funcional seria diferencial entre grupos.
# Resultado: p = 1 — zero associação entre irregularidade e tratamento.
irr_tab = table(reg$stdorToU, reg$irregular)
print(chisq.test(irr_tab))

# ══════════════════════════════════════════════════════════════════════════════
# 8. COBERTURA TEMPORAL E JANELA DE ANÁLISE (todos os blocos)
# ══════════════════════════════════════════════════════════════════════════════

get_coverage = function(fpath) {
  b = fread(fpath, select = c("LCLid", "tstp"))
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b[!is.na(tstp), .(t_min = min(tstp), t_max = max(tstp)), by = LCLid]
}

cov = rbindlist(lapply(block_files_hh, get_coverage))
cov = merge(cov, hh[, .(LCLid, stdorToU)], by = "LCLid", all.x = TRUE)

# Distribuição dos inícios: concentrada em abr–mai/2012 (roll-out principal).
# Términos: 5.063 domicílios chegam até fev/2014 (fim do trial).
print(sort(table(format(cov$t_min, "%Y-%m"))))
print(sort(table(format(cov$t_max, "%Y-%m"))))

# Janelas candidatas: número de domicílios com cobertura completa no período.
check_window = function(start, end, dt) {
  s = as.POSIXct(start, tz = "Europe/London")
  e = as.POSIXct(end,   tz = "Europe/London")
  dt[
    t_min <= s & t_max >= e,
    .(total = .N, ToU = sum(stdorToU == "ToU"), Std = sum(stdorToU == "Std"))]
}

windows = list(
  c("2013-01-01", "2013-12-31"), # ano civil completo — escolhido
  c("2013-01-01", "2013-06-30"),
  c("2013-06-01", "2013-12-31"),
  c("2013-03-01", "2013-11-30"),
  c("2012-11-01", "2013-02-28"),
  c("2012-07-01", "2013-06-30")
)

window_results = rbindlist(
  lapply(
    windows,
    function(w) {
      r = check_window(w[1], w[2], cov)
      data.table(
        start = w[1],
        end = w[2],
        total = r$total,
        ToU = r$ToU,
        Std = r$Std
      )
    }
  )
)
print(window_results)

# Decisão: janela 2013-01-01 → 2013-12-31.
# Razões: (a) ano civil completo captura sazonalidade inteira (inverno e verão);
# (b) retém 5.096 dos 5.566 domicílios (~92%), com 1.038 ToU e 4.058 Std;
# (c) elimina o período inicial de roll-out (2012), onde havia grande variação
# no momento de entrada dos domicílios no trial.
WIN_START = as.POSIXct("2013-01-01", tz = "Europe/London")
WIN_END   = as.POSIXct("2013-12-31", tz = "Europe/London")

# Exclusões explícitas antes de aplicar o filtro de janela:
#   MAC005565 — série de apenas 45 leituras (~22h), irrecuperável e fora da janela.
#   ACORN- — 2 domicílios sem categoria válida que violam positividade (todos ToU).
EXCLUDE = c("MAC005565", hh[Acorn == "ACORN-", LCLid])

ids_2013 = cov[
  t_min <= WIN_START & t_max >= WIN_END & !LCLid %in% EXCLUDE,
  LCLid
]

# ══════════════════════════════════════════════════════════════════════════════
# 9. VERIFICAÇÃO DE RANDOMIZAÇÃO
# ══════════════════════════════════════════════════════════════════════════════
#
# O trial LCL foi desenhado como RCT com recrutamento voluntário: domicílios
# se inscreveram e foram então alocados a ToU ou Std. A proporção global de ToU
# na amostra é ~20.2%. Verificamos se essa proporção é homogênea entre os
# estratos ACORN usando ICs binomiais por categoria fina — se a randomização
# fosse perfeita dentro dos recrutados, todos os ICs deveriam conter ~ 20.2%.

overall_tou = mean(hh$stdorToU == "ToU")

strat = hh[
  !Acorn %in% c("ACORN-", "ACORN-U"),
  .(n = .N, n_tou = sum(stdorToU == "ToU")),
  by = Acorn
][order(Acorn)]

strat[
  ,
  c("lo", "hi") := {
    bts = mapply(
      function(x, n) binom.test(x, n, conf.level = 0.95)$conf.int,
      n_tou, n, SIMPLIFY = FALSE
    )
    list(
      round(100 * sapply(bts, `[`, 1), 1),
      round(100 * sapply(bts, `[`, 2), 1)
    )
  }
]

strat[, pct := round(100 * n_tou / n, 1)]
strat[, contem_global := lo <= 100 * overall_tou & 100 * overall_tou <= hi]

# 15 dos 17 estratos finos têm o IC 95% contendo a proporção global.
# Dois estratos apresentam desequilíbrio significativo:
#   ACORN-K (n=165, 13.3%, IC 8.5–19.5%): "Hard-pressed communities"
#   ACORN-Q (n=831, 14.3%, IC 12.0–16.9%): "Struggling estates" — maior grupo
#     dentro de Adversity
# Esses dois estratos agrupam domicílios de menor renda / habitação social.
# A interpretação mais provável é a seleção no recrutamento voluntário:
# domicílios de menor renda aderiram menos à tarifa dinâmica.
# Não necessariamente falha de randomização dentro dos recrutados,
# mas viés de participação.
# A randomização interna (dentro dos inscritos) não pode ser verificada sem
# acesso aos dados originais do recrutamento.
#
# Consequências para o modelo:
#   (1) O FATE estimado é condicional ao perfil de domicílios participantes —
#       a generalização externa para domicílios ACORN-K/Q é limitada.
#   (2) ACORN_grouped deve entrar como covariável de precisão para reduzir a
#       variância da estimativa e controlar pelo desequilíbrio observado.
#   (3) A identificabilidade dentro da amostra não é comprometida:
#       o desequilíbrio é de representação,
#       não de confundimento pós-recrutamento.
print(strat[, .(Acorn, n, pct, lo, hi, contem_global)])

# Balanço de ACORN_grouped na janela 2013 escolhida.
hh_2013 = hh[LCLid %in% ids_2013]
bal = dcast(
  hh_2013[, .N, .(Acorn_grouped, stdorToU)],
  Acorn_grouped ~ stdorToU, value.var = "N", fill = 0
)
bal[, pct_ToU := round(100 * ToU / (Std + ToU), 1)]
print(bal[order(Acorn_grouped)])

# ══════════════════════════════════════════════════════════════════════════════
# 10. VERIFICAÇÃO DE POSITIVIDADE
# ══════════════════════════════════════════════════════════════════════════════
#
# Positividade exige que, para todo valor das covariáveis fixas, exista
# probabilidade positiva de receber qualquer nível de tratamento:
#   0 < P(A=1 | X=x) < 1  para todo x no suporte de X.
# Como a única covariável fixa utilizada é ACORN_grouped, verificamos se
# cada estrato tem ao menos um domicílio em cada braço.

pos = dcast(
  hh[, .N, .(Acorn_grouped, stdorToU)],
  Acorn_grouped ~ stdorToU, value.var = "N", fill = 0
)
pos[, viola := Std == 0 | ToU == 0]

# Resultado: ACORN- (n = 2, ambos ToU, sem nenhum controle) viola positividade.
# ACORN-U (n=49) e todos os grupos Adversity/Affluent/Comfortable têm ambos
# os braços representados — positividade satisfeita para a análise principal.
# Decisão: excluir os 2 domicílios ACORN- do dataset de análise.
# Na janela 2013 esse estrato já some (o único ACORN- que persiste é ToU,
# e não há par de controle), então a exclusão ocorre naturalmente.
print(pos[order(Acorn_grouped)])

# Verificação fine-grained: ACORN-B tem n = 25 (4 ToU, 21 Std); ACORN-I tem
# n = 51 (8 ToU, 43 Std). Estratos pequenos mas com ambos os braços presentes.
# Não há violação de positividade nas categorias finas além de ACORN-.
pos_fine = dcast(
  hh[!Acorn %in% "ACORN-", .N, .(Acorn, stdorToU)],
  Acorn ~ stdorToU, value.var = "N", fill = 0
)
pos_fine[, viola := Std == 0 | ToU == 0]
print(pos_fine[order(Acorn)])

# ══════════════════════════════════════════════════════════════════════════════
# 11. DIAGNÓSTICO DE NAs EM ENERGY NA JANELA 2013
# ══════════════════════════════════════════════════════════════════════════════
#
# Existem dois tipos distintos de dado faltante neste dataset:
#
# TIPO 1 — NA na coluna energy: a string "Null" é gravada pelo medidor quando
#   não há leitura. Ao converter para numérico, vira NA. Detectável pela
#   contagem de is.na(energy) nas linhas existentes.
#
# TIPO 2 — Linha ausente: a leitura simplesmente não existe no arquivo (nenhuma
#   linha para aquele (LCLid, tstp)). Não detectável pelo diagnóstico de NA —
#   exige comparar a contagem de linhas com o total esperado (17.520).
#   Verificado em bloco separado abaixo.
#
# As duas perguntas relevantes para identificabilidade são as mesmas em ambos os
# casos: (a) qual a prevalência? (b) a ausência é diferencial entre ToU e Std?

cat("Lendo todos os blocos para diagnóstico de NAs (alguns minutos)...\n")

get_na_rate = function(fpath) {
  b = fread(fpath, select = c("LCLid", "tstp", "energy(kWh/hh)"))
  b[, tstp := as.POSIXct(tstp, tz = "Europe/London")]
  b = b[
    LCLid %in% ids_2013 &
    tstp >= WIN_START &
    tstp <  WIN_END + 86400
  ]
  if (nrow(b) == 0) return(NULL)
  b[, energy := suppressWarnings(as.numeric(get("energy(kWh/hh)")))]
  b[, .(
    n_obs = .N,
    n_na = sum(is.na(energy)),
    n_neg = sum(!is.na(energy) & energy < 0)
  ), by = LCLid]
}

na_diag = rbindlist(lapply(block_files_hh, get_na_rate), fill = TRUE)
na_diag = merge(
  na_diag, 
  hh[, .(LCLid, stdorToU, Acorn_grouped)],
  by = "LCLid", 
  all.x = TRUE
)
na_diag[, pct_na  := n_na  / n_obs]
na_diag[, pct_neg := n_neg / n_obs]

# TIPO 1 — NA na coluna energy:
# Resultado: pct_na = 0 em todos os quantis para ids_2013.
# Domicílios com "Null" concentrado têm cobertura incompleta e foram excluídos
# pelo filtro de ids_2013. Nenhuma ação necessária para este tipo.
print(quantile(na_diag$pct_na, c(0, .5, .75, .9, .95, .99, 1)))
print(na_diag[, .(n_total = .N, pct_na_gt5 = mean(pct_na > 0.05))])

# Teste de diferencial por grupo: todos os valores são zero, então a
# missingness é exatamente igual entre ToU e Std — sem risco para
# identificabilidade via NAs.
print(
  na_diag[
    , 
    .(mediana_na = median(pct_na), media_na = mean(pct_na)),
    by = stdorToU
  ]
)

# Valores negativos (microgeração solar): também ausentes neste subconjunto.
# Domicílios com painéis solares tendem a ter cobertura mais irregular,
# ficando fora da janela 2013. Não há ação necessária.
print(
  na_diag[
    ,
    .(
      n_com_neg = sum(n_neg > 0), 
      pct_com_neg = mean(n_neg > 0)
    ),
    by = stdorToU
  ]
)

# ══════════════════════════════════════════════════════════════════════════════
# NOTA SOBRE PERÍODO PRÉ-TRATAMENTO
# ══════════════════════════════════════════════════════════════════════════════
#
# O desenho do LCL não oferece um período pré-tratamento observável: a
# randomização (ToU vs Std) ocorre no momento da entrada de cada domicílio
# no trial, antes de qualquer leitura ser registrada para ele. Não há linha
# de base pré-randomização nos dados — quando um domicílio aparece na base,
# ele já está tratado ou controle.
#
# Um FATE estimado num suposto período "antes" seria 0 por construção: a
# randomização garante equivalência entre grupos no momento da entrada, então
# não há efeito a detectar. Esse tipo de placebo check faz sentido em estudos
# observacionais com DiD (onde o período pré serve para verificar tendências
# paralelas), mas é uninformativo num RCT. O FATE que estimamos é
# integralmente pós-randomização.
#
# A heterogeneidade no tempo de exposição pré-2013 (mediana 7.6 meses,
# intervalo 0.4–13.3, balanceado entre grupos, Wilcoxon p = 0.06) não
# compromete a identificabilidade — o domínio funcional é o calendário de
# 2013 e a comparação é feita no mesmo t para os dois grupos. Adiciona
# variância (domicílios mais antigos no trial podem ter adaptado o
# comportamento), mas não viés.
