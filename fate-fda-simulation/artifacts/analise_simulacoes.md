# Análise das Simulações da Dissertação

## 1. Como estamos gerando os dados da simulação?

O processo gerador de dados (DGP) é implementado em [generate_ecker_dgp.R](file:///c:/Users/moise/Documents/Estudos%20-%20Projetos/dissertacao/causal_inference/fate-fda-simulation/R/generate_ecker_dgp.R) e utiliza como entrada a configuração definida em [config_ecker.R](file:///c:/Users/moise/Documents/Estudos%20-%20Projetos/dissertacao/causal_inference/fate-fda-simulation/R/config_ecker.R). A geração segue o modelo:

$$Y_i(t) = \mu_Y(t) + \mathbf{X}_i^\mathsf{T} \boldsymbol{\beta}^{(A_i)}(t) + \varepsilon_i(t)$$

### Passo a passo:

**a) Covariáveis estáticas (K=9)**
São geradas de acordo com a Tabela D.1 de Ecker (perfil "Men"):

| # | Covariável | Distribuição |
|---|-----------|-------------|
| 1 | `year_first_inc` | $N(-2.3, 1.2)$ |
| 2 | `prev_income` | $N(8.3, 4.1)$ |
| 3 | `children` | $\text{Bernoulli}(0.019)$ |
| 4 | `sec_edu` | $\text{Bernoulli}(0.710)$ |
| 5 | `siblings` | $\text{Poisson}(2.3)$ |
| 6 | `parents_inc` | $N(11.8, 0.9)$ |
| 7 | `foreign_par` | $\text{Bernoulli}(0.023)$ |
| 8 | `par_sec_edu` | $\text{Bernoulli}(0.297)$ |
| 9 | `par_tert_edu` | $\text{Bernoulli}(0.112)$ |

Todas as covariáveis são **estáticas no tempo** (`static_indices = 1:9`), sem feedback do tratamento (`gamma_matrix` toda nula) e sem efeitos cruzados entre covariáveis (`beta_funcs_matrix` toda vazia).

**b) Atribuição de tratamento (A)**
O tratamento é binário e atribuído via modelo logístico (propensity score), utilizando os coeficientes da Tabela E.1 de Ecker (perfil "Men"):

$$\text{logit}\, P(A_i = 1 | \mathbf{X}_i) = \boldsymbol{\theta}^\mathsf{T} [1, \mathbf{X}_i]$$

com $\boldsymbol{\theta} = (-6.194, -0.052, -0.001, 0.020, -0.097, -0.056, 0.510, 0.568, 0.674, 0.851)$.

**c) Efeitos concorrentes $\beta^{(0)}_k(t)$ (efeito das covariáveis no controle)**
São 9 funções suaves/simples que representam como cada covariável influencia o desfecho do grupo controle ao longo do tempo. Exemplos: constante (`-0.02`), decaimento exponencial (`0.5 * exp(-2t)`), ruído cíclico (`0.01 * sin(10t)`), etc.

**d) Efeito de tratamento**
O FATE verdadeiro é distribuído **uniformemente** entre as K+1 = 10 componentes da regressão (intercepto + 9 covariáveis):

$$\delta_k(t) = \beta^{(1)}_k(t) - \beta^{(0)}_k(t) = \frac{\text{effect\_func}(t)}{10}$$

Assim o FATE populacional exato é $\theta(t) = \bar{\mathbf{X}}^\mathsf{T} \boldsymbol{\delta}(t)$, computado como `colMeans(mu_1_mat - mu_0_mat)`.

São 4 cenários de efeito, cada um com uma forma funcional distinta:

| Cenário | Função | Magnitude | Ondaleta |
|---------|--------|-----------|----------|
| Sigmoid | Suave, sigmoide | 2 | Haar |
| Bumps | Picos finos localizados | 3 | db2 |
| Doppler | Alta frequência | 5 | db2 |
| Heavisine | Degraus (escada) | 4 | Haar |

**e) Tendência global** $\mu_Y(t) = 2 \log(t + 0.5)$
Cancela na diferença $\theta(t)$ por construção.

**f) Erros $\varepsilon_i(t)$**
Seguem um **Processo Gaussiano Multivariado** (MGP) com covariância exponencial:
$$\text{Cov}(\varepsilon_i(s), \varepsilon_i(t)) = \sigma^2 \exp\!\left(-\frac{|s-t|}{\rho}\right)$$
com $\sigma = 0.5$ e $\rho = 0.8$ (no grid normalizado $[0,1]$). A geração é feita via decomposição de Cholesky, pré-computada para cada valor de $T$.

**g) Desfecho final**

$$Y_i[k] = \mu_Y(t_k) + \mathbf{X}_i^\mathsf{T} \boldsymbol{\beta}^{(A_i)}(t_k) + \varepsilon_i(t_k)$$

---

## 2. Essa geração está em conformidade com o proposto por Ecker?

**Sim, em larga medida.** O DGP implementado é fiel ao arcabouço de Ecker et al. (2024) nos seguintes aspectos:

| Aspecto | Ecker (2024) | Implementação | Conforme? |
|---------|-------------|---------------|-----------|
| Covariáveis $\mathbf{X}_i$ escalares e pré-tratamento | ✓ | 9 covariáveis estáticas (Tabela D.1 Men) | ✅ |
| Tratamento binário $A_i \in \{0,1\}$ | ✓ | Via modelo logístico (Tabela E.1 Men) | ✅ |
| Desfecho funcional $Y_i(t)$ | ✓ | Trajetória observada em grid de $T$ pontos | ✅ |
| Modelo linear funcional: $Y_i^{(a)}(t) = \mathbf{X}_i^\mathsf{T} \boldsymbol{\beta}^{(a)}(t) + \varepsilon_i^{(a)}(t)$ | ✓ | Implementado exatamente assim | ✅ |
| Erros $\varepsilon_i(t) \sim \text{MGP}(0, \sigma, \mathbb{I}_N)$ | ✓ | Cov. exponencial, $\sigma=0.5$, $\rho=0.8$ | ✅ |
| Efeito heterogêneo: $\beta^{(1)}(t) - \beta^{(0)}(t) = \boldsymbol{\delta}(t)$ | ✓ | Dividido uniformemente como $\delta_k = \text{effect}/10$ | ✅ |
| Estimador: $\hat\theta(t) = \bar{\mathbf{X}}^\mathsf{T}(\hat\beta^{(1)}(t) - \hat\beta^{(0)}(t))$ | ✓ | Implementado em `estimate_fate_fast` | ✅ |

> [!NOTE]
> **Diferenças em relação ao artigo original:**
> - Os efeitos concorrentes $\beta^{(0)}_k(t)$ não são extraídos diretamente de Ecker (que não fornece fórmulas explícitas para eles). São funções razoáveis criadas para simular efeitos realistas.
> - Ecker utiliza dados empíricos reais do mercado de trabalho dinamarquês, enquanto a simulação gera dados sintéticos com as mesmas distribuições marginais.
> - Os 4 cenários de efeito (Sigmoid, Bumps, Doppler, Heavisine) não são de Ecker — são funções-teste clássicas da literatura de ondaletas, usadas para explorar cenários que vão além do efeito suave original.

---

## 3. O que avaliamos na Simulação 1? Com quais critérios? Quais as saídas?

**Objetivo:** Comparar **três métodos de pré-processamento** na estimação do FATE, usando **todo o sinal** (abordagem retrospectiva):

| Método | Descrição | Bandas de confiança |
|--------|-----------|-------------------|
| `wavelet_offline` | DWT completa via lifting (offline) + limiarização semissuave | Bootstrap (dados ruidosos) |
| `wavelet_causal` | DWT causal via `denoise_signal_causal` com janela $W$ | Bootstrap (dados ruidosos) |
| `spline_ecker` | B-Splines com GCV (reprodução de Ecker) | Analíticas (Eq. de Ecker) |

### Grid de simulação
- **N** (amostra): 200, 500, 2000
- **T** (pontos temporais): 128, 256, 512
- **W** (janela causal): 64, 64, 128 (indexado por T)
- **Cenários**: Sigmoid, Bumps, Doppler, Heavisine
- **Repetições**: 1000 (ou 10 em modo `QUICK_TEST`)

### Critérios de avaliação (por simulação)

| Métrica | Definição | Direção desejada |
|---------|-----------|-----------------|
| **NRMSE** | $\sqrt{\text{MSE}(\hat\theta, \theta)} \,/\, \text{range}(\theta)$ | ↓ menor = melhor |
| **Bias** | $\text{mean}(\hat\theta(t) - \theta(t))$ | → próximo de 0 |
| **Cov_Simul** | 1 se $\theta(t) \in$ banda $\forall t$; 0 caso contrário | ↑ próximo de 95% |
| **Cov_Point** | Proporção dos pontos $t$ em que $\theta(t) \in$ banda | ↑ próximo de 95% |
| **Width** | Largura média da banda de confiança | ↓ menor = mais preciso |

### Saídas
- `sim1_raw.csv`: métricas por simulação, cenário, método, N, T
- `sim1_summary.csv`: médias agregadas das métricas
- `sim1_curves.csv`: estimativas pontuais, curvas verdadeiras e bandas para cada simulação

---

## 4. O que avaliamos na Simulação 2? Com quais critérios? Quais as saídas?

**Objetivo:** Avaliar o **estimador online** (*streaming*) — onde cada dado chega ponto a ponto e é processado em tempo real via `new_wavelet_stream` — medindo tanto a **qualidade da estimativa** quanto a **latência computacional** por ponto.

### Diferenças em relação à Simulação 1
- Usa `new_wavelet_stream` (processamento estritamente causal, ponto a ponto) em vez de `denoise_signal_causal` ou `denoise_signal_offline`
- **Não calcula bandas de confiança** (coerente com a limitação computacional descrita na metodologia)
- Mede **latência por indivíduo por ponto temporal** (tempo em microssegundos)

### Grid de simulação
Idêntico à Sim 1: N ∈ {200, 500, 2000}, T ∈ {128, 256, 512}, W ∈ {64, 64, 128}, 4 cenários, 1000 sims.

### Critérios de avaliação

| Métrica | Definição | Direção desejada |
|---------|-----------|-----------------|
| **NRMSE_online** | $\sqrt{\text{MSE}(\hat\theta_{\text{online}}, \theta)} \,/\, \text{range}(\theta)$ | ↓ menor = melhor |
| **Latency_mean_us** | Tempo médio (em µs) para processar um ponto, por indivíduo | ↓ menor = mais rápido |

### Saídas
- `sim2_raw.csv`: métricas por simulação, cenário, N, T
- `sim2_summary.csv`: médias agregadas de NRMSE e latência
- `sim2_curves.csv`: estimativas pontuais, curvas verdadeiras e perfil de latência ao longo de $t$

---

## 5. Divergências na implementação das bandas de Ecker

Ao verificar nossa implementação em [ecker_analytical_bands.R](file:///c:/Users/moise/Documents/Estudos%20-%20Projetos/dissertacao/causal_inference/fate-fda-simulation/R/ecker_analytical_bands.R) contra o artigo original e o pacote [`ffscb`](https://github.com/lidom/ffscb) usado por Ecker, identificamos três divergências:

### 5.1 Violação da premissa gaussiana pelo GCV

A teoria de Ecker assume que o estimador $\hat\theta(t)$ é gaussiano. A seleção do parâmetro de suavização via GCV é um procedimento não-linear aplicado aos dados, o que pode violar essa premissa. Isso não explica diretamente a sobre-cobertura de 100%, mas indica que as bandas analíticas operam sob uma hipótese que não é estritamente satisfeita.

### 5.2 Regularização Ridge (negligenciável)

Nossa implementação adiciona `ridge_lambda = 1e-7` à diagonal de $(X_a^\mathsf{T}X_a)$ antes da inversão. Ecker assume OLS puro. Com $\lambda = 10^{-7}$, o efeito é numericamente desprezível — trata-se apenas de uma salvaguarda contra matrizes quase-singulares.

### 5.3 ✅ CONFIRMADO: $u_\alpha$ constante vs. $u_\alpha(t)$ variante no tempo

> [!IMPORTANT]
> **Esta é a divergência principal e foi confirmada pela análise do código-fonte do pacote `ffscb`.**

Ecker utiliza a função `confidence_bands` do pacote `ffscb` (Liebl & Reimherr, 2023), que constrói o limiar $u(t)$ como uma **função linear por partes** (*piecewise-linear*) ao longo do domínio temporal:

- O domínio $[0,1]$ é particionado em `n_int` sub-intervalos (padrão: 4).
- Para cada sub-intervalo $j$, um coeficiente $c_j$ é otimizado via fórmula de Kac-Rice, utilizando o parâmetro de rugosidade local $\tau(t)$.
- $\tau(t)$ é estimado via `cov2tau_fun`: calcula-se a segunda derivada mista da correlação em $s=t$, resultando num vetor que varia ponto a ponto.
- A banda final é: $\hat\theta(t) \pm u(t) \cdot \sqrt{C_\theta(t,t)}$, onde **ambos** $u(t)$ e $\sqrt{C_\theta(t,t)}$ variam com $t$.

**Nossa implementação** (linhas 101–125 de `ecker_analytical_bands.R`) calcula $u_\alpha$ como um **escalar constante**: o quantil do $\sup_t |Z(t)|$ de um processo gaussiano normalizado simulado. Isso equivale ao caso especial `n_int = 1` do `ffscb` (banda de Kac-Rice constante), que é sistematicamente mais conservador do que o caso `n_int ≥ 2`.

**Consequência:** com $u_\alpha$ constante, a largura da banda é $u_\alpha \cdot \sqrt{C_\theta(t,t)}$ — modulada apenas pela variância local, mas inflada por um fator global que protege contra desvios de alta frequência em todo o domínio. Isso explica a **sobre-cobertura de 100%** no cenário Sigmoid e parte da largura constante (~0.63) observada entre cenários.
