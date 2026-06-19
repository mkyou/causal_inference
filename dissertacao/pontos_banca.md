# Pontos da Banca — Refinamentos da Dissertação

IDs: `O-x` orais (Carlos); `T-x` parecer escrito (M. Taddeo); `T-Qx` perguntas do parecer; `P-x` autor (pós-defesa).
Marcadores: `[ ]` a discutir · `[A]` acatar · `[P]` parcial · `[R]` rebater · `[X]` concluído.

---

## A1 — Definição de termos + lista de conceitos
**Status:** Parcial.
- **T-1a / T-14b'** — "regular" = regularidade de Hölder (Def. 2.1) + forward-ref na intro. **[X]**
- **T-1b** — "vazamento temporal" = causalidade operacional (Def. 2.2) + conexão na intro. **[X]**
- **T-1c** — "família mãe" → "família de ondaletas". **[X]**
- **O-1** — lista de conceitos pré-textual. Adiada para o final (terminologia precisa estabilizar). **[ ]**

---

## A2 — Motivação prática na introdução
**Status:** Parcial (intro concluída; capítulo de aplicações pende — ver A9).
- **O-2 / O-4** — vinheta condutora energia+UTI + CFATE online na introdução. **[X]**
- **T-Q5** — passa de "sem aplicação" para "com aplicação"; redação migra para A9. **[P]**

*Aplicações: Low Carbon London (aleatorizado, principal) + MIMIC-IV/UTI (observacional, confundimento por indicação). Variáveis verificadas. Pendente: definição do tratamento/janela no MIMIC e credenciamento PhysioNet.*

---

## A3 — Argumento contra splines (preciso, não enfático)
**Status:** Concluído.
- **O-3** (p.23) — reforçar por que splines não funcionam nesse caso. **[R]**: já é pontuado ao menos 3x ao longo do texto. E acho que com definição formal de regularidade, já está satisfeito.
- **T-14b** (p.23) — parecerista: "perturbação local afeta coeficientes de uma ampla vizinhança (de quem?)"; B-splines têm sim localização por construção. **[X]**
- *Solução (aplicada na p.23):* trocar a alegação imprecisa pelo **ajuste global dos coeficientes** + suavidade global inverificável (Def. 2.1; Donoho–Johnstone).

---

## A4 — Identificação causal (pp.17–19)
**Status:** Concluído.
- **T-2** — "a quem se atribui um tratamento". **[X]**
- **T-3** — define identificabilidade antes da lista. **[X]**
- **T-4** — forma binária da consistência; (3) e (4) **mantidos separados** (rebatido "(4) contido em (3)"). **[X]**
- **T-5** — exemplo de interferência (imunidade de rebanho). **[X]**
- **T-6** — ML "adaptado à inferência causal". **[X]**
- **T-7** — rodapé: ignorabilidade plena vs. média (rebatido "só vale em casos lineares"). **[X]**
- **T-7b** — derivação da eq. (2.2) em passos rotulados. **[X]**

---

## A5 — Modelo estrutural (2.3) (p.20)
**Status:** Concluído (pende forward-ref real à simulação P-2, que será feita).
- **T-9** — X pré-tratamento explícito + definição de confundidora. **[X]**
- **T-11** — notação \(\epsilon^{(a)}\) esclarecida; exogeneidade em média definida em Cap. 3 (metodologia l. 48); gaussianidade justificada pela limiarização; footnote robustez MAD. **[X]**

---

## A6 — Ondaletas/MRA (pp.23–30)
**Status:** Concluído.
- **T-14a** — \(\cos t\) removido do exemplo. **[X]**
- **O-7** — adendo técnico \(f(t)\) vs \(f[k]\) simplificado e reposicionado. **[P]**
- **T-16** — relação função de escala ↔ MRA (Mallat 1989) esclarecida; notação \(V_{2^j}\) mantida. **[P]**
- **T-17** — base do subespaço de detalhes e relação com função de escala (ref. Bachman et al.). **[X]**

---

## A7 — DGP da simulação (Cap. 4)
**Status:** A discutir. Quatro sintomas do mesmo problema (descrição do DGP incompleta); uma reescrita detalhada (com pseudocódigo) resolve.
- **O-5** — explicitar se os dados são observacionais ou aleatorizados (e a implicação para identificabilidade). **[ ]**
- **T-18** (p.55) — esclarecer "distribuímos o efeito total uniformemente" (conflita com a Tabela 2, onde os efeitos diferem). **[ ]**
- **T-19** (p.55) — explicar como o diagrama da Fig. 8 (p.72) é endereçado no DGP (o desfecho tem existência anterior à atribuição?). Parecerista pede detalhar a geração dos dados. **[ ]**
- **T-20** (p.56) — explicar o papel dos escores de propensão na Tabela 2 (ou removê-los, se não usados). **[ ]**

---

## A8 — Justificativas metodológicas
**Status:** A discutir. (T-Q1, T-Q2, T-10, T-21 migraram para A12/A13.)
- **T-Q3** (p.41) — onde as contribuições tratam a intervenção no "meio" do acompanhamento (Fig. 8); coerente com A7. **[ ]**
- **T-Q4** (p.43) — por que bootstrap e não bandas analíticas (Delta funcional / fórmula 2.5)? Conecta ao relaxamento Gaussiano (P-2): sem gaussianidade, o bootstrap deixa de ser ad hoc. **[ ]**

---

## A9 — Conclusão sobre dados reais
**Status:** Decidido; redação pendente.
- **O-9** — trocar "correlação espúria" (`conclusao.tex`) pela razão real: confundimento não controlado / identificabilidade não garantida sem desenho experimental (no caso da UTI, **confundimento por indicação** e não-confundimento não testável). **[A]**
- **T-Q5** — reescrever de "sem aplicação" para **aplicações com ressalvas**: Low Carbon London (aleatorizado, principal) + UTI/MIMIC (observacional, com ressalvas). **[A]**

---

## A10 — Administrativos
- **O-6** — declaração de uso de IA (correções, ortografia, código) — norma CCPG/UNICAMP. **[ ]**
- **O-8** — decidir se a tabela completa do apêndice (`tab:resultados_sim1_completo`, 108 linhas) vai para o texto (tendência: permanece). **[ ]**
- **T-12** (p.21) — adicionar referências ao primeiro parágrafo. **[ ]**
- **T-15b** (p.46) — onde encontrar/instalar/usar o `rLifting` → absorvido por A12. **[ ]**

---

## A11 — Didática transversal (P-3)
**Status:** Cap. 2 concluído. Pende aplicação real (A9).
- Intro **[X]** — dois fios energia+UTI + CFATE plantado.
- §2.1.1 notação **[X]** — mapeamento \(Y_i^{(a)}(t)\) para energia e UTI.
- §2.1.1 suavização **[X]** — DWT por curva, energia+UTI.
- §2.2 regularidade **[X]** — exemplo energia, irregularidade local.
- §2.2.5 tempo real **[X]** — FATE atualizado ponto a ponto, energia+UTI.
- Aplicação real **[ ]** — payoff final do fio condutor (ver A9).
- Mover derivações pesadas ao apêndice **[ ]**

---

## A12 — Texto ↔ novo `rLifting` (P-1)
**Status:** A discutir.
- **P-1a** — conclusão: remover "pressupõe \(\Delta_t\) fixo" (desatualizado). **[ ]**
- **P-1b / P-1c** — Cap. 3: documentar suporte a grade irregular, versão, repositório/instalação. **[ ]**
- **P-1d** — ressalva MAD: grade irregular comum a todos OK; por indivíduo = trabalho futuro. **[ ]**

*Hierarquia do caso irregular: grade comum regular/irregular → trivial ou OK com rLifting; por indivíduo offline → alinhar via interpolação; por indivíduo online → trabalho futuro.*

---

## A13 — Campanha única de simulação (P-1, P-2, T-10, T-Q1, T-21)
**Status:** A discutir. Prioridade: P-2 > P-1e > T-Q1 > T-10.
- **P-2** — ruído não gaussiano (caudas pesadas, heteroscedasticidade); contraste com splines. **[ ]**
- **P-1e** — amostragem irregular. **[ ]**
- **T-Q1** — **CFATE** \(=\int_0^t \mathrm{FATE}(s)\,ds\) (aprovado): barato (integra réplicas bootstrap existentes); a versão online vira monitor de efeito acumulado em tempo real. **[P]** — **introduzido na introdução** **[X]** (descrição textual + exemplo UTI + forward-ref ao Cap. 3); **retomar** ao longo do texto: def. formal no Cap. 3 (Metodologia, por ser contribuição nossa) **[ ]**, cálculo na simulação **[ ]**, monitor online na conclusão **[ ]**.
- **T-10** (p.20) — sensibilidade a desvios da suposição estrutural (confundimento não observado / não linearidade). **[ ]**
- **T-21** (p.58) — documentar critério de resolução máxima (ondaletas) e nº de nós (B-splines) no novo desenho. **[ ]**

---

## Visão geral

| Ação | Status |
|---|---|
| A1 Termos + glossário | Parcial (glossário adiado) |
| A2 Motivação | Parcial (capítulo de aplicações pendente) |
| A3 Anti-splines | ✅ Concluído |
| A4 Identificação causal | ✅ Concluído |
| A5 Modelo (2.3) | ✅ Concluído |
| A6 Ondaletas/MRA | ✅ Concluído (T-16/O-7 parciais por decisão) |
| A7 DGP simulação | A discutir |
| A8 Justificativas | A discutir |
| A9 Conclusão/dados reais | Decidido; redação pendente |
| A10 Administrativos | A discutir |
| A11 Didática | Parcial |
| A12 Texto ↔ rLifting | A discutir |
| A13 Simulação | A discutir |
