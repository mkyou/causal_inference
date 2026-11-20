%% rbras-poster.cls — pôster A0 retrato, identidade visual da 70ª RBras.
%% Reproduz a geometria do Template_Poster_70RBras.pptx (slide de 841 x 1189 mm,
%% 2383,9 x 3370,4 bp). Todas as medidas estão em bp, unidade do ponto do
%% PowerPoint, de modo que 1 pt do template = 1 bp aqui.
%%
%% Compile com pdflatex, xelatex ou lualatex. Em pdflatex o texto usa o pacote
%% montserrat, da TeX Live; se ele não estiver instalado a classe recorre a
%% Helvetica e emite aviso.

\NeedsTeXFormat{LaTeX2e}[2019/10/01]
\ProvidesClass{rbras-poster}[2026/10/02 v2.0 Poster A0 retrato - 70a RBras]
\LoadClass{article}

% ---------- Página A0 retrato, sem margens ----------
\RequirePackage[paperwidth=841mm,paperheight=1189mm,margin=0pt]{geometry}
\pagestyle{empty}
\setlength\parindent{0pt}
\setlength\parskip{0pt}

% ---------- Fontes: Montserrat, com recuo seguro ----------
\RequirePackage{fix-cm}   % fontes de matemática em qualquer corpo
\RequirePackage{iftex}
\newcommand\rb@ls[1]{#1}  % espaçamento entre letras; redefinido abaixo
\ifPDFTeX
  \RequirePackage[utf8]{inputenc}
  \RequirePackage[T1]{fontenc}
  \IfFileExists{montserrat.sty}
    {\RequirePackage[defaultfam,tabular,lining]{montserrat}}
    {\RequirePackage{helvet}
     \ClassWarningNoLine{rbras-poster}{Pacote montserrat ausente; usando
       Helvetica. Instale com: tlmgr install montserrat}}
  \RequirePackage{letterspace}
  \renewcommand\rb@ls[1]{\textls[107]{#1}}
\else
  \RequirePackage{fontspec}
  \IfFontExistsTF{Montserrat}
    {\setsansfont{Montserrat}}
    {\ClassWarningNoLine{rbras-poster}{Fonte Montserrat ausente; usando a
       fonte sem serifa padrao}}
  \renewcommand\rb@ls[1]{{\addfontfeature{LetterSpace=10}#1}}
\fi
\renewcommand\familydefault{\sfdefault}
% A matemática fica na fonte matemática padrão do LaTeX, com serifa, como no
% template do PowerPoint, que usa Cambria Math no texto em Montserrat.
\RequirePackage{amsmath}
\RequirePackage{exscale}  % delimitadores grandes em qualquer corpo

% babel é opcional: entra apenas se o idioma estiver instalado, de modo que a
% classe compile também em instalações mínimas da TeX Live.
\IfFileExists{brazil.ldf}
  {\RequirePackage[brazil]{babel}
   \AtBeginDocument{\ifcsname shorthandoff\endcsname\shorthandoff*{"}\fi}}
  {\ClassWarningNoLine{rbras-poster}{babel-portuguese ausente; sem hifenizacao
     em portugues. Instale com: tlmgr install hyphen-portuguese babel-portuges}}

% ---------- Pacotes ----------
\RequirePackage[table]{xcolor}
\RequirePackage{graphicx}
\RequirePackage{tikz}
\RequirePackage{pgfplots}
\pgfplotsset{compat=1.18}
\RequirePackage{eso-pic}
\RequirePackage[absolute,overlay]{textpos}
\RequirePackage{tcolorbox}
\RequirePackage{tabularx}
\RequirePackage{booktabs}
\RequirePackage{enumitem}
\RequirePackage{ragged2e}
\RequirePackage{qrcode}

% ---------- Paleta: cores lidas do template ----------
\definecolor{rbAzul}{HTML}{073A70}        % azul institucional
\definecolor{rbAmarelo}{HTML}{F9BD38}     % amarelo das faixas
\definecolor{rbFundo}{HTML}{EAF0F7}       % fundo das caixas
\definecolor{rbTexto}{HTML}{1F2933}       % corpo de texto
\definecolor{rbCinza}{HTML}{52606D}       % rótulos e legendas
\definecolor{rbCinzaClaro}{HTML}{9AA5B1}  % série de comparação
\definecolor{rbGrade}{HTML}{E9EAEB}       % linhas de grade
\definecolor{rbEixo}{HTML}{C7CACC}        % linha do eixo
\definecolor{rbClaro}{HTML}{F4F4F4}       % texto sobre o azul do cabeçalho

% ---------- Tamanhos: pt do PowerPoint = bp ----------
\newcommand\rb@fs[2]{\fontsize{#1bp}{#2bp}\selectfont}
\newcommand\tamTitulo{\rb@fs{72}{80}}
\newcommand\tamMensagem{\rb@fs{66}{76}}
\newcommand\tamNumero{\rb@fs{80}{80}}
\newcommand\tamSecao{\rb@fs{48}{56}}
\newcommand\tamEquacao{\rb@fs{40}{50}}
\newcommand\tamFigTitulo{\rb@fs{36}{43}}
\newcommand\tamAutores{\rb@fs{36}{42}}
\newcommand\tamCorpo{\rb@fs{34}{39}}
\newcommand\tamCartao{\rb@fs{34}{41}}
\newcommand\tamQR{\rb@fs{32}{36}}
\newcommand\tamRodape{\rb@fs{31}{36}}
\newcommand\tamRotulo{\rb@fs{28}{33}}
\newcommand\tamEixo{\rb@fs{28}{32}}
\newcommand\tamTabela{\rb@fs{28}{34}}
\newcommand\tamAfil{\rb@fs{28}{33}}
\newcommand\tamKicker{\rb@fs{28}{32}}
\newcommand\tamLegenda{\rb@fs{24}{28}}
\newcommand\tamRef{\rb@fs{22}{26}}

% ---------- Metadados ----------
\newcommand\rb@titulo{Título do trabalho}
\newcommand\rb@autores{Autores}
\newcommand\rb@afiliacoes{}
\newcommand\rb@email{}
\newcommand\rb@mensagem{Escreva o achado em uma frase, com o número.}
\def\rb@qrlink{}
\newcommand\rb@qrlegenda{Artigo · código · dados}
\newcommand\rb@logoevento{logo-rbras.png}
\newcommand\rb@logoinst{logo-instituicao.png}
\newcommand\rb@rodape{70ª RBras}
\newcommand\rb@nconc{5}
\newcommand\rb@cA{}\newcommand\rb@cB{}\newcommand\rb@cC{}
\newcommand\rb@refs{}\newcommand\rb@agradec{}
\newcommand\rb@yC{2710}   % topo da seção Conclusões, em bp

\newcommand\titulo[1]{\renewcommand\rb@titulo{#1}}
\newcommand\autores[1]{\renewcommand\rb@autores{#1}}
\newcommand\afiliacoes[1]{\renewcommand\rb@afiliacoes{#1}}
\newcommand\email[1]{\renewcommand\rb@email{#1}}
\newcommand\mensagem[1]{\renewcommand\rb@mensagem{#1}}
\newcommand\qrlink[1]{\def\rb@qrlink{#1}}
\newcommand\qrlegenda[1]{\renewcommand\rb@qrlegenda{#1}}
\newcommand\logoevento[1]{\renewcommand\rb@logoevento{#1}}
\newcommand\logoinstituicao[1]{\renewcommand\rb@logoinst{#1}}
\newcommand\rodape[1]{\renewcommand\rb@rodape{#1}}
\newcommand\conclusoes[4][5]{\renewcommand\rb@nconc{#1}%
  \renewcommand\rb@cA{#2}\renewcommand\rb@cB{#3}\renewcommand\rb@cC{#4}}
\newcommand\referencias[1]{\renewcommand\rb@refs{#1}}
\newcommand\agradecimentos[1]{\renewcommand\rb@agradec{#1}}
\newcommand\posicaoconclusoes[1]{\renewcommand\rb@yC{#1}}

% ---------- Logo, ou quadro reservado se o arquivo faltar ----------
\newcommand\rb@logo[3]{%
  \begingroup
  \edef\rb@@arq{#1}%
  \IfFileExists{\rb@@arq}%
    {\includegraphics[width=#2,height=#3,keepaspectratio]{\rb@@arq}}%
    {\fcolorbox{white}{rbAzul}{\parbox[c][\dimexpr#3-2\fboxsep\relax][c]%
      {\dimexpr#2-2\fboxsep\relax}{\centering\color{white}\tamAfil
       Logo\\\expandafter\detokenize\expandafter{\rb@@arq}}}}%
  \endgroup}

% ---------- Cartão de conclusão ----------
\newcommand\rb@cartao[3]{%
  \fill[rbFundo] (#1,{\rb@yC+95}) rectangle ++(714.7,270);
  \node[anchor=west,inner sep=0,font=\tamNumero\bfseries,text=rbAzul]
    at ({#1+28},{\rb@yC+230}) {#2};
  \node[anchor=west,inner sep=0,text width=536.7bp,align=left,
        font=\tamCartao\bfseries,text=rbTexto] at ({#1+150},{\rb@yC+230}) {#3};}

% ---------- Fundo fixo: cabeçalho, mensagem, conclusões, rodapé ----------
\newcommand\rb@fundo{%
  \begin{tikzpicture}[x=1bp,y=-1bp,baseline={(0,0)}]
  \useasboundingbox (0,0) rectangle (2384,3371);
  % Cabeçalho
  \fill[rbAzul] (0,0) rectangle (2384,516);
  \fill[rbAmarelo] (0,516) rectangle (2384,535);
  \node[anchor=north west,inner sep=0] at (69,21)
    {\rb@logo{\rb@logoevento}{430bp}{400bp}};
  \node[anchor=north west,inner sep=0] at (1910,88)
    {\rb@logo{\rb@logoinst}{400bp}{330bp}};
  \node[anchor=center,inner sep=0,text width=1390bp,align=center,text=white]
    at (1190,260) {%
      {\tamTitulo\bfseries\rb@titulo\par}\vspace{24bp}%
      {\tamAutores\color{rbClaro}\rb@autores\par}\vspace{10bp}%
      {\tamAfil\color{rbClaro}\rb@afiliacoes\par}\vspace{6bp}%
      {\tamAfil\color{rbClaro}\rb@email\par}};
  % Faixa da mensagem principal
  \fill[rbFundo] (72,590) rectangle (2312,1050);
  \fill[rbAmarelo] (72,590) rectangle (96,1050);
  \node[anchor=west,inner sep=0,text width=1700bp,align=left] at (150,820) {%
    {\tamKicker\bfseries\color{rbCinza}\rb@ls{MENSAGEM PRINCIPAL}\par}%
    \vspace{16bp}%
    {\tamMensagem\bfseries\color{rbAzul}\rb@mensagem\par}};
  % QR code
  \fill[white] (1928,625) rectangle (2258,955);
  \ifx\rb@qrlink\@empty
    \draw[rbAzul,line width=4bp,dash pattern=on 14bp off 9bp]
      (1928,625) rectangle (2258,955);
    \node[font=\tamQR\bfseries,text=rbAzul] at (2093,790) {QR code};
  \else
    \node[inner sep=0] at (2093,790) {\qrcode[height=310bp]{\rb@qrlink}};
  \fi
  \node[anchor=north,align=center,font=\tamLegenda,text=rbCinza] at (2093,968) {\rb@qrlegenda};
  % Conclusões
  \fill[rbAmarelo] (72,{\rb@yC+10}) rectangle ++(16,58);
  \node[anchor=west,inner sep=0,font=\tamSecao\bfseries,text=rbAzul]
    at (108,{\rb@yC+39}) {\rb@nconc\hspace{0.45em}Conclusões};
  \rb@cartao{72}{1}{\rb@cA}
  \rb@cartao{834.7}{2}{\rb@cB}
  \rb@cartao{1597.3}{3}{\rb@cC}
  % Referências e agradecimentos
  \node[anchor=north west,inner sep=0,text width=1460bp,align=left]
    at (72,{\rb@yC+405}) {%
      {\tamLegenda\bfseries\color{rbAzul}Referências\par}\vspace{8bp}%
      {\tamRef\color{rbCinza}\rb@refs\par}};
  % Rodapé
  \fill[rbAmarelo] (0,3270) rectangle (2384,3289);
  \fill[rbAzul] (0,3289) rectangle (2384,3371);
  \node[font=\tamRodape,text=white] at (1192,3330) {\rb@rodape};
  \end{tikzpicture}}

\AtBeginDocument{%
  \AddToShipoutPictureBG*{\AtPageUpperLeft{\rb@fundo}}%
  \RaggedRight\tamCorpo\color{rbTexto}}
\AtEndDocument{\par\null}

% ---------- Colunas ----------
\newcommand\rb@inicio{\justifying\tamCorpo\color{rbTexto}}
\newenvironment{colunaesquerda}[1][900bp]
  {\begin{textblock*}{#1}(72bp,1120bp)\rb@inicio}{\end{textblock*}}
\newenvironment{colunadireita}[1][1246bp]
  {\begin{textblock*}{#1}(1066bp,1120bp)\rb@inicio}{\end{textblock*}}

% ---------- Seções numeradas ----------
\newlength\rbespacosecao \setlength\rbespacosecao{80bp}
\newenvironment{secao}[2]{%
  \par\noindent
  \raisebox{-14bp}[0pt][0pt]{\color{rbAmarelo}\rule{16bp}{58bp}}%
  \hspace{20bp}{\tamSecao\bfseries\color{rbAzul}#1\hspace{0.45em}#2}\par
  \vspace{24bp}\tamCorpo\color{rbTexto}}%
  {\par\vspace{\rbespacosecao}}

\setlist[itemize]{label={\color{rbAzul}\textbullet}, leftmargin=45bp,
  labelsep=18bp, itemsep=8bp, topsep=10bp, parsep=0pt, partopsep=0pt}

% ---------- Caixa da equação ----------
\newtcolorbox{caixaequacao}{colback=rbFundo, colframe=rbFundo, boxrule=0pt,
  arc=0pt, left=40bp, right=40bp, top=50bp, bottom=50bp, halign=center,
  before skip=20bp, after skip=0pt, fontupper=\tamEquacao,
  before upper={\setlength\abovedisplayskip{0pt}%
    \setlength\abovedisplayshortskip{0pt}\setlength\belowdisplayskip{14bp}}}
\def\tagform@#1{\maketag@@@{\textcolor{rbCinza}{\tamRotulo(\ignorespaces#1\unskip)}}}

% ---------- Figuras e tabelas ----------
%
% O bloco de figura tem geometria fixa, copiada do template: 620 bp de largura
% por 980 bp de altura, com o título em 0, o rótulo em 120, o gráfico em 205 e
% a legenda em 900. Como as duas figuras ocupam caixas idênticas, elas ficam
% alinhadas qualquer que seja o número de linhas do título.
%
% A área de plotagem também é fixa, nas mesmas proporções do manualLayout do
% gráfico do PowerPoint: margem esquerda e direita de 14% da largura, margem
% superior de 5% e inferior de 7% da altura.
\newcounter{rbfigura}\newcounter{rbtabela}
\newcommand\rbFigLarg{620}        % largura do bloco de figura, em bp
\newcommand\rbFigAlt{980}         % altura do bloco de figura, em bp
\newcommand\rbGrafAlt{680}        % altura da moldura do gráfico, em bp
\newcommand\rbGrafMargEsq{86.8}   % 0,14 da largura
\newcommand\rbGrafMargInf{47.6}   % 0,07 da altura
\newcommand\rbGrafLarg{446.4}     % 0,72 da largura: área de plotagem
\newcommand\rbGrafAltEixo{598.4}  % 0,88 da altura: área de plotagem

% \rbgrafico{ambiente axis}: moldura de 620 x 680 bp com a área de plotagem
% sempre na mesma posição. O argumento é só o ambiente axis.
\newcommand\rbgrafico[1]{%
  \begin{tikzpicture}[x=1bp,y=1bp]
  \useasboundingbox (0,0) rectangle (\rbFigLarg,\rbGrafAlt);
  \begin{scope}[shift={(\rbGrafMargEsq,\rbGrafMargInf)}]#1\end{scope}
  \end{tikzpicture}}

% \figura{título-afirmação}{rótulo: o que, unidade}{ambiente axis}{legenda}
\newcommand\figura[4]{%
  \refstepcounter{rbfigura}%
  \begin{tikzpicture}[x=1bp,y=-1bp]
  \useasboundingbox (0,0) rectangle (\rbFigLarg,\rbFigAlt);
  \node[anchor=north west,inner sep=0,text width=\rbFigLarg bp,align=left,
        font=\tamFigTitulo\bfseries,text=rbTexto] at (0,0) {#1};
  \node[anchor=north west,inner sep=0,text width=\rbFigLarg bp,align=left,
        font=\tamRotulo,text=rbCinza] at (0,120) {#2};
  \node[anchor=north west,inner sep=0] at (0,205) {\rbgrafico{#3}};
  \node[anchor=north west,inner sep=0,text width=\rbFigLarg bp,align=left,
        font=\tamLegenda,text=rbCinza] at (0,900)
    {\textcolor{rbAzul}{\bfseries Figura~\therbfigura.}~#4};
  \end{tikzpicture}}

% \tabela[largura]{título-afirmação}{tabular}{legenda}
% Bloco fixo de 355 bp: título em 0, tabela em 75, legenda em 305.
\newcommand\tabela[4][1276]{%
  \refstepcounter{rbtabela}%
  \begin{tikzpicture}[x=1bp,y=-1bp]
  \useasboundingbox (0,0) rectangle (#1,355);
  \node[anchor=north west,inner sep=0,text width=#1bp,align=left,
        font=\tamFigTitulo\bfseries,text=rbTexto] at (0,0) {#2};
  \node[anchor=north west,inner sep=0,text width=#1bp,align=left,
        font=\tamTabela] at (0,75) {#3};
  \node[anchor=north west,inner sep=0,text width=#1bp,align=left,
        font=\tamLegenda,text=rbCinza] at (0,305)
    {\textcolor{rbAzul}{\bfseries Tabela~\therbtabela.}~#4};
  \end{tikzpicture}}

% Duas figuras lado a lado, 36 bp de intervalo
\newcommand\duasfiguras[2]{\par\noindent
  \begin{minipage}[t]{\rbFigLarg bp}\rb@inicio #1\end{minipage}\hfill
  \begin{minipage}[t]{\rbFigLarg bp}\rb@inicio #2\end{minipage}\par}

% ---------- Gráficos: padrão da casa ----------
% "assume math mode" mantém os rótulos dos eixos na fonte do texto, não na
% fonte de matemática, como no template.
\pgfkeys{/pgf/number format/.cd, use comma, 1000 sep={.}, assume math mode}
\pgfplotsset{
  rbras/.style={scale only axis=true,
    width=\rbGrafLarg bp, height=\rbGrafAltEixo bp,
    axis line style={draw=none}, tick style={draw=none},
    tick label style={font=\tamEixo, text=rbCinza},
    label style={font=\tamEixo, text=rbCinza},
    grid style={draw=rbGrade, line width=1.875bp}, clip=false},
  rbras linha/.style={rbras, ymajorgrids,
    axis x line*=bottom, x axis line style={draw=rbEixo, line width=1.875bp},
    yticklabel style={/pgf/number format/.cd, fixed, fixed zerofill, precision=1}},
  rbras inclinacao/.style={rbras, xmin=-0.5, xmax=1.5, xtick={0,1}, ytick=\empty,
    axis x line*=bottom, x axis line style={draw=rbEixo, line width=1.875bp},
    xticklabel style={font=\rb@fs{31}{36}, text=rbTexto}}}
% \inclinacao[desloc]{nome}{início}{fim}: linha cinza
% \inclinacaodestaque[desloc]{nome}{início}{fim}: linha azul
% desloc afasta na vertical os rótulos que se sobrepõem (em unidades do eixo)
\newcommand\inclinacao[4][0]{%
  \addplot[rbCinzaClaro, line width=3.75bp, mark=*, mark size=5bp]
    coordinates {(0,#3) (1,#4)};
  \node[anchor=east, font=\tamEixo, text=rbCinza] at (axis cs:0,#3)
    {\pgfmathprintnumber[fixed,fixed zerofill,precision=1]{#3}\hspace*{12bp}};
  \node[anchor=west, font=\tamEixo, text=rbCinza] at (axis cs:1,{#4+#1})
    {\hspace*{12bp}\pgfmathprintnumber[fixed,fixed zerofill,precision=1]{#4}%
     \hspace*{14bp}#2};}
\newcommand\inclinacaodestaque[4][0]{%
  \addplot[rbAzul, line width=5.6bp, mark=*, mark size=6bp]
    coordinates {(0,#3) (1,#4)};
  \node[anchor=east, font=\tamEixo\bfseries, text=rbAzul] at (axis cs:0,#3)
    {\pgfmathprintnumber[fixed,fixed zerofill,precision=1]{#3}\hspace*{12bp}};
  \node[anchor=west, font=\tamEixo\bfseries, text=rbAzul] at (axis cs:1,{#4+#1})
    {\hspace*{12bp}\pgfmathprintnumber[fixed,fixed zerofill,precision=1]{#4}%
     \hspace*{14bp}#2};}

\endinput
