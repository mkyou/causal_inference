# Pontos da Banca — por Capítulo

IDs: `O-x` orais (Carlos); `T-x` parecer escrito (M. Taddeo); `T-Qx` perguntas do parecer; `P-x` autor (pós-defesa).
Marcadores: `[ ]` a discutir · `[A]` acatar · `[P]` parcial · `[R]` rebater · `[X]` concluído.

---

## Geral
- Uso do termo "tempo real" vs. "tempo quase real". **[X]** — "tempo real" usado de forma consistente em todo o texto; nenhuma ocorrência de "tempo quase real".
- **Notação \(\mathcal{T}\) vs. \(\mathcal{S}\)** **[X]** — \(\mathcal{T}\) mantido como domínio pós-intervenção (conforme Ecker) em toda a revisão e metodologia; clarificação adicionada na revisão (linha 67) e no parágrafo introdutório do PGD do Cap. 3. \(\mathcal{S} \supset \mathcal{T}\) introduzido na seção do estimador online como relaxação explícita da hipótese de Ecker.

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
**Status:** Em andamento (DGP a reescrever; simulação em curso).

*Quatro sintomas do mesmo problema (descrição do DGP incompleta); uma reescrita detalhada com pseudocódigo resolve.*

- **O-5** — explicitar se os dados são observacionais ou aleatorizados + implicação para identificabilidade. **[R]**
- **T-18** (p.55) — esclarecer "distribuímos o efeito total uniformemente" (conflita com Tabela 2). **[A]**
- **T-19** (p.55) — geração dos dados: desfecho preexiste à intervenção. **[A]** (pronto na simulação, falta atualizar no texto)
- **T-20** (p.56) — papel dos escores de propensão na Tabela 2 (ou remover). **[A]** (acho que irei remover e escrever a fórmula de confusão explicitamente)
- **T-21** (p.58) — documentar critério de resolução máxima (ondaletas) e nº de nós (B-splines). **[A]**
- **T-10** (p.20) — sensibilidade a confundimento não observado (via `CONF_LEVELS`, 3 níveis). **[A]** (feito na simulação, falta atualizar o texto)
- **P-2** — ruído não gaussiano (caudas pesadas, heteroscedasticidade); contraste com splines. **[A]**
- **T-Q4** — reforçar a importância do bootstrap no cenário não gaussiano: com ruído de caudas pesadas, bandas analíticas não se aplicam e o bootstrap se justifica plenamente (complementa a justificativa do Cap. 3). **[A]**
- **P-1e** — amostragem irregular. **[A]** Isso aqui é a atual simulação 3, que precisa ser vista com mais carinho. Vamos ver como vamos lidar com o PGD e o que vamos estimar.
- **P-2** — reforçar simulação 2 com mais métodos e/ou cenários de dados. **[A]**
- **T-Q1** — cálculo do CFATE integrando réplicas bootstrap existentes (barato por construção). **[A]** (feito na simulação, falta texto)

*Prioridade: P-2 > P-1e > T-Q1 > T-10.*

---

## Cap. 5 — Aplicações
**Status:** Decidido; redação pendente.
- **T-Q5** — Low Carbon London (aleatorizado, principal) + UTI/MIMIC (observacional, com ressalvas). **[A]**
- **A11** — payoff final do fio condutor energia+UTI. **[ ]**

*Pendente: definição do tratamento/janela temporal no MIMIC e credenciamento PhysioNet.*

---

## Cap. 6 — Conclusão
**Status:** Redação pendente.
- **O-9** — substituir "correlação espúria" por confundimento não controlado / identificabilidade não garantida (UTI: confundimento por indicação). **[A]**
- **T-Q5** — reescrever de "sem aplicação" para "com aplicações e ressalvas" (ver Cap. 5). **[A]**
- **P-1d** — ressalva MAD: grade irregular comum a todos OK; por indivíduo = trabalho futuro. **[A]** não é uma questão do MAD, é do FATE. A ressalva do MAD tem a ver com a variância homoscedástica.
- **T-Q1** — CFATE online como monitor de efeito acumulado em tempo real. **[A]**
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
| Cap. 4 Simulação | ⏳ Em andamento |
| Cap. 5 Aplicações | 📝 Redação pendente |
| Cap. 6 Conclusão | 📝 Redação pendente |
| Transversal | ⏳ A discutir |
