# Inferência causal em dados funcionais

Código, resultados e materiais da dissertação sobre estimação de efeitos causais médios funcionais com pré-processamento por ondaletas via *lifting scheme*.

## Organização

- `functional-causal-wavelets/text/dissertation/`: texto final, bibliografia, figuras e PDF da dissertação.
- `functional-causal-wavelets/text/rbras/dissertation-contest/`: fontes e PDF da apresentação para o concurso da RBras.
- `functional-causal-wavelets/text/rbras/poster/`: fontes e PDF do pôster da RBras.
- `functional-causal-wavelets/fate-fda-simulation/`: funções, scripts e resultados das três simulações.
- `functional-causal-wavelets/fate-fda-application/low_carbon_energy/`: scripts e resultados da aplicação Low Carbon London.

Os caminhos entre as pastas são usados pelos arquivos LaTeX e pelos scripts. As figuras e os resumos de simulação mantidos em `output/` fazem parte dos resultados e incluem figuras utilizadas na dissertação.

## Compilação dos documentos

É necessária uma instalação LaTeX com as classes e pacotes utilizados por cada documento. A dissertação depende da classe `abntex2`, além da classe institucional incluída no repositório. Execute os comandos a partir da respectiva pasta:

```sh
# Apresentação do concurso
cd functional-causal-wavelets/text/rbras/dissertation-contest
pdflatex -interaction=nonstopmode -halt-on-error presentation.tex
```

```sh
# Pôster, a partir da raiz do repositório
cd functional-causal-wavelets/text/rbras/poster
pdflatex -interaction=nonstopmode -halt-on-error poster.tex
bibtex poster
pdflatex -interaction=nonstopmode -halt-on-error poster.tex
pdflatex -interaction=nonstopmode -halt-on-error poster.tex
```

```sh
# Dissertação, a partir da raiz do repositório
cd functional-causal-wavelets/text/dissertation
pdflatex -interaction=nonstopmode -halt-on-error dissertation.tex
bibtex dissertation
pdflatex -interaction=nonstopmode -halt-on-error dissertation.tex
pdflatex -interaction=nonstopmode -halt-on-error dissertation.tex
```

A dissertação referencia figuras nas pastas de resultados da simulação e da aplicação; mantenha a estrutura do repositório ao compilá-la. O arquivo administrativo `chapters/auxiliary/catalog_card.pdf` não está versionado e deve ser disponibilizado localmente para a compilação completa.

## Simulações e aplicação

As simulações usam R e o pacote [rLifting](https://github.com/mkyou/rLifting), além dos pacotes carregados nos scripts. Os pontos de entrada são `functional-causal-wavelets/fate-fda-simulation/scripts/sim_1.R`, `sim_2.R` e `sim_3.R`; os scripts complementares geram figuras e resumos. As configurações completas podem demandar tempo e memória consideráveis; consulte cada script antes de executá-lo.

A aplicação usa os scripts numerados de `01_explore_data.R` a `04_plots.R`. Eles devem ser executados a partir de `functional-causal-wavelets/fate-fda-application/`. Os dados originais e arquivos processados de Low Carbon London não estão incluídos no repositório; devem ser obtidos separadamente e organizados como descrito nos próprios scripts. A fonte citada na dissertação é o [UK Data Service](https://discover.ukdataservice.ac.uk/doi?sn=7857-1).

Os resumos e figuras versionados permitem consultar os resultados, mas não constituem todos os dados intermediários: por exemplo, os scripts de geração de artefatos das simulações referenciam arquivos brutos ou curvas que não estão incluídos. Reproduzir integralmente os resultados requer gerar novamente esses arquivos, além dos dados da aplicação e das dependências R.
