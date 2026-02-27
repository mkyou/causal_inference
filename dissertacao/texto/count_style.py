import re

files = {
    'Introdução': ('introducao.tex', 4.0),
    'Metodologia': ('metodologia.tex', 8.5),
    'Revisão': ('revisao.tex', 34.0)
}

patterns = {
    '\"o que [verbo]\"': r'\bo\s+que\s+\w+',
    '\"Embora\" / \"Ainda que\"': r'\b(embora|ainda que)\b',
    '\"Nesse sentido\" / \"Nesse contexto\" / \"Por outro lado\"': r'\b(nesse sentido|nesse contexto|por outro lado)\b',
    '\"Em essência\" / \"Essencialmente\" / \"Em particular\"': r'\b(em essência|essencialmente|em particular)\b',
    '\"Dessa forma\" / \"Deste modo\" / \"Desta forma\" / \"Desse modo\"': r'\b(dessa forma|desta forma|deste modo|desse modo)\b',
    '\"Ou seja\" / \"Isto é\"': r'\b(ou seja|isto é)\b',
    '\"No entanto\" / \"Contudo\" / \"Todavia\"': r'\b(no entanto|contudo|todavia)\b',
    '\"Uma vez que\" / \"Dado que\"': r'\b(uma vez que|dado que)\b',
    '\"Vale notar\" / \"Cabe ressaltar\" / \"Reforçamos\"': r'\b(vale notar|cabe ressaltar|reforçamos)\b',
    '\", [gerúndio]\"': r',\s+\w+ndo\b'
}

print('| Muleta / Estrutura Rígida | Introdução (~4 pág) | Metodologia (~8.5 pág) | Revisão (~34 pág) | Total Absoluto |')
print('| :--- | :---: | :---: | :---: | :---: |')

for name, regex in patterns.items():
    counts = []
    tot = 0
    for section, (fname, pages) in files.items():
        with open(fname, 'r', encoding='utf-8') as f:
            content = f.read()
        cnt = len(re.findall(regex, content, re.IGNORECASE))
        tot += cnt
        counts.append(f\"{cnt} ({cnt/pages:.2f}/pág)\")
    print(f'| **{name}** | {counts[0]} | {counts[1]} | {counts[2]} | **{tot}** |')
