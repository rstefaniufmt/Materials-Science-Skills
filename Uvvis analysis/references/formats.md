# Formatos de arquivo UV-Vis - leitura e depuração

## Texto (CSV, TXT, PRN, ASC, DAT)

O leitor testa separadores (`;`, tab, `,`, espaços) e decimal (`.`/`,`; vírgula só se houver
números `123,45` de fato), escolhe a combinação com mais linhas numéricas consistentes,
ignora cabeçalho e rodapé não numéricos e usa a última linha de cabeçalho como nomes.

Exportações reconhecidas:
- **Agilent Cary (WinUV/Cary 60/100/300/5000)**: pares de colunas por amostra
  (`Wavelength (nm),Abs,Wavelength (nm),Abs,`), nomes das amostras na 1ª linha, vírgula no
  fim da linha e metadados no rodapé → um Spectrum por par, título = nome da amostra.
- **Shimadzu UVProbe/LabSolutions**: `"Wavelength nm.","Abs."` (ou `%T`, `R%`), às vezes
  com linha de título antes; tab ou vírgula.
- **Ocean Insight/OceanView**: cabeçalho até `>>>>>Begin Spectral Data<<<<<`; integração,
  espectrômetro e médias vão para `meta`. O tipo de y vem do nome do nó (`Absorbance`,
  `Transmission`, `Reflection`); arquivos de intensidade bruta são `counts`.
- **PerkinElmer UV WinLab (ASCII)**, **Thermo VISIONlite**, **Hitachi**: x,y simples.
- **Layout largo**: 1ª linha = λ, cada linha seguinte = um espectro (1ª célula = nome).
- **Várias colunas y** com um único x: um Spectrum por coluna (título = cabeçalho).

Eixo x: nm por padrão; cm⁻¹ (> 4000), µm, eV (valores 0.3-7) detectados pelo cabeçalho ou
pela faixa. Tipo de y: pelo cabeçalho (`Abs`, `A`, `%T`, `T%`, `%R`, `R%`, `KM`, `F(R)`,
`Counts`, `log(1/R)`, `epsilon`); senão inferido (`ykind_guessed`). **%T e %R não se
distinguem pelos números** - confirme.

## Excel (.xlsx/.xls/.xlsm)

Cada planilha é convertida em linhas de texto e passa pelas mesmas regras; o nome da
planilha serve de pista para o tipo de y ("Absorbância", "Reflectância"...). Requer
`openpyxl` (já instalado na maioria dos ambientes; senão `pip install openpyxl`).

## JCAMP-DX (.jdx, .dx)

Mesmo decodificador da skill de FT-IR (AFFN, PAC, SQZ, DIF, DUP, XYPOINTS, PEAK TABLE,
blocos). `XUNITS=NANOMETERS` é o normal em UV-Vis; `1/CM` e `MICROMETERS` são convertidos.
`YUNITS` como `ABSORBANCE`, `TRANSMITTANCE`, `REFLECTANCE`, `Logarithm epsilon` (NIST:
vira `logeps`).

## Galactic GRAMS .spc

- Formato novo (byte 1 = 0x4B): cabeçalho de 512 bytes; x uniforme (ffirst→flast) ou
  explícito (flag TXVALS), multi-subarquivo (TMULTI: séries temporais, diode array, mapas),
  x por subarquivo (TXYXYS); y em float32 (expoente 0x80) ou inteiros escalados por
  2^(exp−32) (ou 16 bits com TSPREC).
- Formato antigo (0x4D): cabeçalho de 256 bytes, inteiros com as palavras de 16 bits
  trocadas - o leitor testa as duas ordens e fica com a curva suave (aviso no meta).
- Tipos de x (3 = nm, 1 = cm⁻¹, 2 = µm, 14 = eV) e de y (2 = absorbância, 3 = KM,
  128 = transmitância, 129 = reflectância, 4 = contagens) são mapeados.
Validado em arquivos de teste públicos (filtro de hólmio, UV com linha de base ruim, XY
explícito, diode array 182 subarquivos, séries 3D).

## PerkinElmer .sp (UV WinLab, Lambda; Spectrum)

Contêiner binário `PEPE`; o leitor procura o membro de faixa do eixo (35698, tipo 29981) e
percorre os membros seguintes (intervalo, nº de pontos, rótulos x/y, dados float64). Os
limites do eixo são **inclusivos** (n = |fim − início|/passo + 1); um valor de
preenchimento inicial é descartado. Rótulo x `nm` (Lambda) ou `cm-1` (IR).

## Não suportados

Cary binário (.bsw, .dsw, .vsp), Shimadzu .vspd/.uvd, Hitachi .uds, Jasco .jws: exporte
como CSV/TXT (ou JCAMP) no software do equipamento.

## Depuração

```python
from uvvis_io import load_all
for s in load_all("arquivo"):
    print(s.summary())
    print(s.meta.get("header_lines") or s.meta.get("header") or s.meta.get("comment"))
```
Se o tipo de y estiver errado: `s.ykind = "%R"` (ou `--ykind %R` na CLI) antes de processar.
Se o eixo vier invertido/em outra unidade, confira `meta['xunits_original']`.
