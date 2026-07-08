# Pontos da Banca — por Capítulo

IDs: `O-x` orais (Carlos); `T-x` parecer escrito (M. Taddeo); `T-Qx` perguntas do parecer; `P-x` autor (pós-defesa).
Marcadores: `[ ]` a discutir · `[A]` acatar · `[P]` parcial · `[R]` rebater · `[X]` concluído.

---

## Geral
- Verificar uso do termo "tempo real" vs. "tempo quase real". Preferimos "tempo real". Necessário checar se está padronizado. **[]**
- **Notação \(\mathcal{T}\) vs. \(\mathcal{S}\). Verificar padronização em texto. \(\mathcal{S}\) contém \(\mathcal{T}\) unindo ao período pré-tratamento. Checar se isso segue limpo no texto. **[]**. 
- Verificar se uso de t-Student(5) está padronizado. Ocorrências de t(5) devem ser substituídas de acordo com a primeira.
- Verificar uso de ~x, em que "x" representa um número (ex: Figura ~\ref{fig}). Padronizamos o não uso do "~". 
- Verificar se termos como FATE, NRMSE, CFATE e afins estão envoltos do \mathrm ao longo de todo o texto.
- Verificar lista de siglas. Checar se cobre todas as siglas empregadas. Checar se está correto. Averiguar mudanças e decidir.
- Verificar conceitos plausíveis para "lista de conceitos" sugeridos pela banca.

## Cap. 1 — Introdução
**Status:** Concluído.
- **T-1a / T-14b'** — forward-ref a "regular" = regularidade de Hölder (Def. 2.1 no Cap. 2). **[X]**
- **T-1b** — "vazamento temporal" introduzido + conexão com causalidade operacional (Def. 2.2 no Cap. 2). **[X]**
- **O-2 / O-4** — vinheta condutora energia+UTI + CFATE online plantado. **[X]**
- **T-Q1** — CFATE introduzido: descrição textual + exemplo UTI + forward-ref ao Cap. 3. **[X]**

---

## Cap. 2 — Revisão
**Status:** Concluído (T-16/O-7 parciais por decisão).
- **T-1a / T-14b'** — "regular" = regularidade de Hölder (Def. 2.1). **[X]**
- **T-1b** — definição formal de vazamento temporal (Def. 2.2). **[X]**
- **T-1c** — "família mãe" → "família de ondaletas". **[X]**
- **T-12** (p.21) — referências ao primeiro parágrafo do capítulo. **[X]**
- **O-3** — reforçar por que splines não funcionam; já pontuado ≥3× no texto. **[R]**
- **T-14b** — argumento anti-splines: ajuste global dos coeficientes + suavidade uniformemente regular inverificável. **[X]**
- **T-2** — "a quem se atribui um tratamento". **[X]**
- **T-3** — define identificabilidade antes da lista. **[X]**
- **T-4** — forma binária da consistência; (3) e (4) mantidos separados. **[X]**
- **T-5** — exemplo de interferência (imunidade de rebanho). **[X]**
- **T-6** — ML "adaptado à inferência causal". **[X]**
- **T-7** — rodapé: ignorabilidade plena vs. média. **[X]**
- **T-7b** — derivação da eq. (2.2) em passos rotulados. **[X]**
- **T-9** — X pré-tratamento explícito + definição de confundidora. **[X]**
- **T-11** — notação \(\epsilon^{(a)}\) esclarecida (gaussianidade e exogeneidade em média aprofundadas no Cap. 3). **[X]**
- **T-14a** — \(\cos t\) removido do exemplo. **[X]**
- **O-7** — adendo técnico \(f(t)\) vs \(f[k]\). **[P]** (decisão: manter parcial)
- **T-16** — relação função de escala ↔ MRA (Mallat 1989). **[P]** (decisão: manter parcial)
- **T-17** — base do subespaço de detalhes + ref. Bachman. **[X]**
- **A11** — mapeamento notação/suavização/regularidade/tempo real com exemplos energia+UTI (§2.1.1–§2.2.5). **[X]**

---

## Cap. 3 — Metodologia
**Status:** Concluído.
- **T-11** — gaussianidade justificada pela limiarização; exogeneidade em média definida aqui; footnote robustez MAD. **[X]**
- **T-Q3** (p.41) — intervenção no meio do acompanhamento: PGD já contempla (desfecho preexiste à intervenção, como mostra o DAG no apêndice); tornar explícito no texto com ponteiro ao DAG. **[X]**
- **T-Q4** (p.43) — bootstrap vs. bandas analíticas: limiarização viola gaussianidade do estimador; justificativa escrita no Cap. 3. **[X]**
- **T-Q1** — CFATE: definição formal em seção própria (`sec-cfate`), com estimador trapezoidal, bandas bootstrap e atualização online. **[X]**
- **T-15b / P-1b / P-1c** — menção ao CRAN suficiente; instruções de instalação pertencem à documentação do pacote. **[X]**
- **rLifting** — revisar e atualizar a seção do `rLifting` para refletir o estado atual do pacote. **[X]**

---

## Cap. 4 — Estudos de Simulação
**Status:** Concluído.

- Viés(t_k) vs. Viés Global: pontual usado nas figuras do Apêndice; global é a quantidade tabulada. **[X]**
- **O-5** — explicitar se os dados são observacionais ou aleatorizados + implicação para identificabilidade. **[R]**
- **T-18** (p.55) — PGD: fórmula explícita substitui descrição informal; "Tabela 2" removida. **[X]**
- **T-19** (p.55) — desfecho preexiste à intervenção: δ(t)=0 para t < t_treat explicitado. **[X]**
- **T-20** (p.56) — escores de propensão: modelo logístico com coeficientes explícitos. **[X]**
- **T-21** (p.58) — critério de resolução máxima (ondaletas) e nº de nós (B-splines) documentados. **[X]**
- **T-10** (p.20) — confundimento não observado: três níveis discutidos e tabelados. **[X]**
- **P-2** — ruído não gaussiano: simulações t-Student(5), contraste com splines. **[X]**
- **P-2** — reforçar simulação 2 com mais métodos. **[R]**
- **T-Q4** — importância do bootstrap com ruído não gaussiano: justificativa no texto. **[X]**
- **T-Q1** — CFATE via réplicas bootstrap existentes: mencionado no encerramento do cap. **[X]**
- **P-1e** — amostragem irregular: Simulação 3 com três graus de irregularidade, lifting scheme, tabelas e análise textual. **[X]**

---

## Cap. 5 — Aplicações
**Status:** Decidido; redação pendente.
- **T-Q5** — Low Carbon London. **[A]**
- **A11** — payoff final do fio condutor da rede elétrica. **[ ]**

---

## Cap. 6 — Conclusão
**Status:** Redação pendente.
- **O-9** — substituir "correlação espúria" por confundimento não controlado / identificabilidade não garantida (UTI: confundimento por indicação). **[A]**
- **T-Q5** — reescrever de "sem aplicação" para "com aplicações e ressalvas" (ver Cap. 5). **[A]**
- **P-1d** — ressalva MAD: grade irregular comum a todos OK; por indivíduo = trabalho futuro. **[A]** não é uma questão do MAD, é do FATE. A ressalva do MAD tem a ver com a variância homoscedástica.
- **T-Q1** — CFATE online como monitor de efeito acumulado em tempo real. **[A]**
---

## Apêndice
- Revisão geral do Apêndice. **[P]** — Sim1: hierarquia seção/subseção, parágrafos analíticos (cobertura, NRMSE, viés, largura), figuras pareadas; Sim2: texto boxplots + remissão ao texto principal; Sim3: texto confirmando padrões do texto principal. Nomes de cenários padronizados (capitalizados + `\emph{}`) em todo o apêndice e nos capítulos 3 e 4. Revisão final pendente após Cap. 5 (Apêndice C de aplicações a adicionar).

---

## Transversal
- **O-1** — lista de conceitos pré-textual (aguardar estabilização da terminologia). **[A]**
- **O-6** — declaração de uso de IA (norma CCPG/UNICAMP). **[A]**
- **O-8** — tabela completa do apêndice (tendência: permanece). **[A]** ela vai crescer.
- **A11** — mover derivações pesadas ao apêndice. **[R]** não temos derivações muito pesadas.

---

## Visão geral

| Capítulo | Status |
|---|---|
| Cap. 1 Introdução | ✅ Concluído |
| Cap. 2 Revisão | ✅ Concluído (T-16/O-7 parciais por decisão) |
| Cap. 3 Metodologia | ✅ Concluído |
| Cap. 4 Simulação | ✅ Concluído |
| Cap. 5 Aplicações | 📝 Redação pendente |
| Cap. 6 Conclusão | 📝 Redação pendente |
| Transversal | ⏳ A discutir |
