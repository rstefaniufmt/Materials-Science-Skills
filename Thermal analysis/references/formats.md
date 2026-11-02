# Formatos de arquivo e depuração da leitura

O formato é reconhecido pelo **conteúdo** (não pela extensão). `ta.py info arquivo` mostra o
formato detectado, os canais e as unidades. Todos viram um `Run` com canais canônicos:
`time` [min], `T` [°C] (amostra), `T_ref` [°C], `mass` [mg], `mass_pct` [%], `hf`
(unidade em `r.units['hf']`), `dtg`, `segment`; o resto vai para `r.extra`.

## Sumário
1. Binário NETZSCH (.ngb-*)
2. Exportações ASCII por fabricante
3. CSV/TXT/Excel genéricos
4. Binários não suportados: como exportar
5. Depuração

## 1. Binário NETZSCH (.ngb-ss3, .ngb-bs3, .ngb-ds3)

- Lido pelo pacote **pyngb** (MIT, não oficial): `pip install pyngb`. O script avisa se faltar.
- `.ngb-ss3` = medida de amostra; `.ngb-bs3` = correção (branco); `.ngb-ds3` = DSC/DTA.
  A correção pode vir embutida ("amostra + correção"); se não vier, `info` avisa "sem
  correção" e você passa `--correction arquivo.ngb-bs3`. O pyngb recusa correções com
  programa de temperatura diferente (comportamento correto: a subtração não faz sentido).
- Massa: o arquivo guarda Δm (mg) a partir de 0; o script soma a massa da amostra dos
  metadados. DSC vem em µV; a curva de sensibilidade gravada no arquivo
  (S = (p2 + p3·z + p4·z² + p5·z³)·exp(−z²), z = (T − p0)/p1) converte para mW e, com a
  massa, para W/g. Sem sensibilidade, o sinal fica em µV (só T e caráter endo/exo).
- Metadados úteis: massa, cadinho, gases e vazões (canais de vazão de purga e proteção -
  usados para achar a troca de gás no OIT), operador, programa (segmentos).
- Direção endo: deduzida do sinal das entalpias de calibração quando presentes.

## 2. Exportações ASCII por fabricante

| Fabricante/software | Como reconhecer | Notas |
|---|---|---|
| NETZSCH Proteus (Export → ASCII) | linhas `#CHAVE:valor`, colunas `##Temp./°C;Time/min;DSC/(mW/mg);Mass/%` | `#SEPARATOR`, `#DECIMAL` respeitados; símbolo de grau corrompido (›, �) tratado |
| TA Universal Analysis (Export Data File / Signals) | `Nsig`, `Sig1 ...`, `StartOfData`; UTF-16 | linhas-sentinela (tempo negativo) ignoradas; massa em `Size` |
| TA TRIOS (.txt, Export) | blocos `[step]`, linha de nomes + linha de unidades | cada `[step]` vira segmento |
| TA TRIOS Excel (.xls/.xlsx) | planilhas por etapa | concatenadas como segmentos; `.xls` precisa de `pip install xlrd` |
| Mettler STARe (Export → ASCII/"other format") | colunas `Index t Ts Tr Value` + linha `[s] [°C] [°C] [mg]` | várias curvas no arquivo (TGA + DSC no mesmo eixo) fundidas em STA; curva de branco detectada; massa de "SampleID ..., X mg" |
| PerkinElmer Pyris (Save Data As → ASCII) | cabeçalho de 2 linhas, blocos `1) Heat from ...` | blocos viram segmentos; convenção endo ↑ |
| Setaram Calisto/SETSOFT | cabeçalho com nomes com espaço (`Sample Temperature(°C)`), UTF-16, vírgula decimal | |
| Shimadzu TA-60 (export texto) | genérico com unidades | lido pelo leitor genérico |
| Hitachi/SII, Linseis, Rigaku | exportação texto | leitor genérico; confira `info` |

Os menus mudam com a versão do software; em geral procure "Export", "Save as ASCII/Text"
ou clique-direito na curva → "Export".

## 3. CSV/TXT/Excel genéricos

- Separador (`; , tab espaço`) e decimal (ponto ou vírgula) detectados; cabeçalho e rodapé
  de texto ignorados.
- Papéis das colunas pelo nome e pela unidade: `Temperature (°C)`, `Temp/K`, `T [°C]`,
  `Time (min)`, `t/s`, `Weight (mg)`, `Mass %`, `Heat Flow (W/g)`, `DSC /(mW/mg)`, `DTA (uV)`,
  `Deriv. Weight (%/°C)`; ou linha de unidades separada (`[s] [K] [mg]`). Kelvin, segundos,
  horas, g e µg são convertidos.
- Sem cabeçalho ou nomes ambíguos: `--columns time:min,T:°C,mass:mg,hf:mW` (papéis: time, T,
  T_ref, mass, mass_pct, hf, dtg, segment, skip).
- Massa não informada: `--mass 8.5` (mg).
- Excel com várias planilhas: cada planilha é lida; `--index` escolhe a curva.

## 4. Binários não suportados: como exportar

| Arquivo | Origem | O que pedir ao usuário |
|---|---|---|
| `.tri`, `.001`/`.002` | TA TRIOS / TA Universal Analysis | exportar como texto (TRIOS: Export → Text/Excel; UA: File → Export Data File → File & Plot Signals) |
| `.stare`, `.ta` | Mettler STARe | exportar a curva como ASCII |
| `.dsc`, `.tga` (binário Pyris) | PerkinElmer | Save Data As → ASCII |
| `.tad`, `.d60` | Shimadzu | exportar como texto/CSV |
| `.ngb-*` sem pyngb | NETZSCH | `pip install pyngb` ou exportar ASCII no Proteus |

Inclua no export: tempo, temperatura da amostra, massa (mg) e/ou sinal DSC/DTA, e, se
possível, os canais de gás (OIT) e a temperatura de programa.

## 5. Depuração

- `info` mostra canais/unidades: se T saiu em K ou o tempo em s, a coluna não foi reconhecida
  - use `--columns`.
- "massa em %" sem mg: J/g e perdas em mg não podem ser calculados; peça a massa.
- Segmentos estranhos (muitos `iso` curtos): ruído na T de controle; os segmentos com menos
  de 1 min são fundidos.
- Direção endo errada: `--endo up|down` vale para todos os comandos.
- `convert -o dados.csv` grava o que foi lido (unidades canônicas) para conferir no Excel/Origin;
  `--extra` inclui os canais extras.
