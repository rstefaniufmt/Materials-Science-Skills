# Formatos de arquivo FT-IR - notas para leitura e depuração

## JCAMP-DX (.jdx, .dx, .jcamp)

Texto ASCII com rótulos `##RÓTULO=valor`. Rótulos são normalizados (maiúsculas, sem espaços,
hífens, barras ou "_": `##FIRST X` = `##FIRSTX`). Comentários após `$$`.

Rótulos que importam:
- `XUNITS`: `1/CM` (cm⁻¹), `MICROMETERS`, `NANOMETERS` → o leitor converte tudo para cm⁻¹.
- `YUNITS`: `TRANSMITTANCE` (fração 0-1, mas muitos arquivos gravam 0-100: o leitor
  detecta e trata como %T), `ABSORBANCE`, `REFLECTANCE`, `KUBELKA-MUNK`,
  `ARBITRARY UNITS`; PNNL/NIST quantitativo usa absortividade (ex. `(micromol/mol)-1m-1`).
- `XFACTOR`/`YFACTOR`: multiplicam os inteiros das linhas de dados.
- `FIRSTX`/`LASTX`/`NPOINTS`/`DELTAX`: eixo x de `XYDATA=(X++(Y..Y))`.

Formas de dados:
- `##XYDATA=(X++(Y..Y))`: cada linha = x inicial + vários y. Compressões ASDF:
  - AFFN: números separados por espaço/vírgula; PAC: separados por sinal `+`/`-`;
  - SQZ: 1º dígito vira letra com sinal (`@`=0, `A-I`=+1..9, `a-i`=-1..-9);
  - DIF: diferença em relação ao valor anterior (`%`=0, `J-R`=+1..9, `j-r`=-1..-9);
  - DUP: repetição do último valor/diferença (`S-Z`=1..8, `s`=9).
  - Em DIF, o 1º y de cada linha repete o último da linha anterior (valor de checagem) e
    é descartado; se não bater, o leitor avisa "Checagem Y falhou".
- `##XYPOINTS=(XY..XY)` e `##PEAK TABLE=(XY..XY)`: pares x,y explícitos. Peak table não é
  espectro contínuo (`meta['is_peak_table']=True`) - plote como barras, não como linha.
- Arquivos compostos (`##BLOCKS`, `##LINK`): vários `##TITLE` - `load_all` devolve cada
  bloco com dados.
- `NTUPLES` (espectros 2D, NMR complexo, MS em série) não é suportado.

Validado contra o conjunto de teste do JCAMP (Lancashire/McDonald: fix/pac/sqz/dif/dup) e
~70 espectros IR do NIST/PNNL: valores idênticos ao pacote `jcamp` onde ambos leem.

## CSV / TXT / DPT / PRN / XY / ASC

O leitor testa separadores `;`, tab, `,`, espaço e decimal `.` ou `,` e escolhe a
combinação com mais linhas numéricas consistentes. Vírgula decimal só é aceita se houver
números `123,45` de fato. Linhas não numéricas antes dos dados são cabeçalho; a última
vira nomes de coluna.

- Coluna 0 = x; cada outra coluna vira um espectro (título = nome da coluna).
- Layout "largo" (1ª linha = eixo x, cada linha seguinte = um espectro, 1ª célula = nome) é
  detectado automaticamente.
- x em µm se todos os valores < 30 (ou cabeçalho com µm/micro); nm pelo cabeçalho.
- Tipo de y: pelo cabeçalho (`Abs`, `A`, `%T`, `Transmittance`, `KM`...) ou pelo texto
  acima dos dados; senão **inferido** (`meta['ykind_guessed']=True`). A inferência erra
  em casos patológicos (espectro de absorbância com linha de base muito alta, %R): confirme.
- Encoding: UTF-8 (com/sem BOM), UTF-16 (Excel "texto Unicode"), cp1252/latin-1.
- Exportações típicas: OPUS `.dpt` (vírgula, sem cabeçalho), OMNIC `.csv` (vírgula, sem
  cabeçalho), Spectrum `.csv` (cabeçalho com `cm-1,%T` ou `cm-1,A`), Excel BR (`;` e `,`).

## PerkinElmer Spectrum .sp

Binário little-endian, assinatura `PEPE` + 40 bytes de descrição. Depois, blocos
`(uint16 id, int32 tamanho, conteúdo)`. Os membros do conjunto de dados, cada um com
código de tipo (uint16) no início do conteúdo:

| id | significado | tipo | conteúdo |
|---|---|---|---|
| 35698 | AbscissaRange | 29981 | float64 x_inicial, float64 x_final |
| 35699 | OrdinateRange | 29981 | float64 y_min, y_max |
| 35700 | Interval | 29979 | float64 passo (negativo se decrescente) |
| 35701 | NumPoints | 29995 | int32 |
| 35703 / 35704 | rótulo x / y | 29987 | int16 len + texto (`cm-1`, `A`, `%T`, `KM`...) |
| 35708 | Data | 29974 | uint32 nº de bytes + float64[] |
| 35709 / 35713 | caminho / nome | 29987 | texto |

Os mesmos ids reaparecem em blocos de histórico com tipo de string; por isso o leitor
procura `35698` **com tipo 29981** e percorre os membros a partir dali.

Armadilhas (validadas):
- O número de pontos cobre os limites **inclusivamente**: 4000→700 a 1 cm⁻¹ = 3301 pontos.
  Tratar como exclusivo desloca todas as bandas em um intervalo, sem nada parecer errado.
  Teste: `len(x) == round(x[0]-x[-1])/passo + 1`.
- Alguns arquivos têm um valor de preenchimento ~0 antes dos dados (n+1 valores): o leitor
  descarta e avisa.
- Metadados de instrumento (modelo, detector, acessório, varreduras, faixa) estão como
  strings `#u` + int16 len; ficam em `meta['sp_strings']`, e alguns campos (accessory,
  resolution, scans) são extraídos quando o par chave/valor é reconhecível.

## Não suportados (peça exportação)

- Thermo OMNIC `.spa`/`.spg`, Bruker OPUS binário (`.0`, `.1`...), Agilent `.spc`/`.asp`,
  Shimadzu `.spc`: exporte para CSV, DPT ou JCAMP-DX no software do equipamento.
  (Para OPUS existe o pacote `brukeropusreader`; para Galactic `.spc`, `spc-spectra`.)

## Depuração

```python
from ftir_io import load_all
for s in load_all("arquivo"):
    print(s.summary()); print(s.meta.get("header", {}) or s.meta.get("header_lines"))
```
Se o eixo parecer deslocado, compare `s.x[:3]` com o software do equipamento; se o tipo de
y estiver errado, ajuste `s.ykind` à mão (`s.ykind = "%T"`) antes de processar.
