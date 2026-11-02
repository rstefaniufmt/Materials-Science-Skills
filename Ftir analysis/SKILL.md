---
name: ftir-analysis
description: "Lê, processa, plota, marca bandas e ajuda a interpretar espectros de infravermelho (FT-IR, ATR, KBr, IR) em JCAMP-DX (.jdx/.dx), CSV/TXT/DPT (inclusive vírgula decimal) e PerkinElmer .sp. Use sempre que aparecerem esses arquivos ou quando o usuário pedir para plotar/comparar espectros de IR, marcar ou listar picos/bandas, atribuir grupos funcionais, converter transmitância↔absorbância, corrigir linha de base, normalizar, empilhar espectros para artigo/tese, calcular razões ou índices de banda (ex.: índice de carbonila), ou interpretar um FT-IR - mesmo que não diga 'FT-IR' explicitamente (ex.: 'o espectro de IV da minha amostra', 'infravermelho do polímero')."
---

# FT-IR: leitura, figuras, bandas e interpretação

Scripts em `scripts/` (só numpy, scipy, pandas, matplotlib). Rode de qualquer pasta:

```bash
S=<caminho-desta-skill>/scripts
python $S/ftir.py info arquivo.sp                      # metadados, faixa, tipo de ordenada, avisos
python $S/ftir.py report arquivo.jdx -d saida --lang pt   # figura rotulada + bandas.csv + resumo .md
python $S/ftir.py plot a.sp b.sp c.csv -o fig.png --mode A --offset auto --label --top 12 \
       --normalize max --highlight 1760:1680:C=O --names "Puro" "5%" "10%" --lang pt
python $S/ftir.py peaks arquivo.csv --assign -o bandas.csv
python $S/ftir.py compare *.sp --band C=O:1760:1680 --band ref:1490:1440 -o indices.csv
python $S/ftir.py convert arquivo.sp -o arquivo.csv --to A --decimal ,
```

Em Python (para algo fora da CLI): `sys.path.insert(0, S)` e use
`ftir_io.load/load_all`, `ftir_process` (baseline, smooth, derivative, normalize,
atr_correction, subtract, co2_flatten), `ftir_peaks` (find_bands, band_metrics,
band_ratio, compare_series), `ftir_assign` (assign, screen_groups, candidates_for),
`ftir_plot` (plot_spectra, plot_with_derivative). As docstrings explicam parâmetros.

## Fluxo de trabalho

1. **Leia e confira** com `info` antes de qualquer coisa. Olhe três coisas:
   - `Ordenada`: A, T, %T, KM... Se aparecer `(INFERIDO - confirmar)`, o arquivo não
     declarava o tipo e ele foi deduzido pela forma (bandas para cima = A, para baixo = T).
     Diga ao usuário qual tipo foi assumido.
   - `AVISO`s: saturação (A > 2.5 ou T < 0.5%), NPOINTS inconsistente, pad no .sp.
   - Faixa e intervalo: ATR de diamante costuma cortar em ~525-400 cm⁻¹; o intervalo de
     dados não é a resolução instrumental (geralmente 4 cm⁻¹).
2. **Processe só o necessário**, sempre em absorbância (é a grandeza proporcional à
   concentração; linha de base e normalização em %T distorcem). Registre o que fez - o
   `report` já escreve o processamento no .md. Padrões razoáveis:
   - linha de base: `rubberband` para ATR/KBr de sólidos; `als` quando a linha de base
     tem curvatura suave e bandas largas; `points` (pontos-âncora sem banda) quando o
     usuário quer reproduzir uma correção manual;
   - suavização: só se o ruído atrapalha, com janela menor que as bandas mais finas;
   - normalização: `max` para comparar formas; `band` (banda interna que não muda, ex.
     CH2 1465 em poliolefinas) para comparar intensidades entre amostras.
3. **Detecte bandas** (`peaks`/`find_bands`). O limiar automático usa ruído local; a
   coluna `snr` e os `flags` (`CO2?`, `H2O(v)?`, `ruído?`, `saturada?`, `borda`) dizem o
   que merece desconfiança. Ombros vêm da 2ª derivada (`kind=shoulder`) - trate-os como
   sugestões e confirme com `plot_with_derivative` antes de afirmar que existem.
4. **Interprete** seguindo `references/interpretation.md` (leia antes de escrever a
   interpretação). A coluna `candidates` de `--assign` é uma tabela de correlação
   ranqueada, não uma resposta: ela já penaliza candidatos cuja banda parceira não está no
   espectro (`(sem parceira)`), mas a decisão é sua, por coerência entre bandas, com o que
   o usuário sabe da amostra e com a literatura do material. Para polímeros e materiais
   comuns, `references/materials.md` tem as bandas-assinatura.
5. **Figura**: `plot` / `plot_spectra` já seguem as convenções (x invertido 4000→400,
   rótulos verticais sem sobreposição com linhas-guia, ticks de y ocultos quando há
   offset, nomes ao lado dos traços). **Sempre abra a imagem gerada e olhe** antes de
   entregar: rótulos demais, legenda sobre uma banda ou um pico de ruído rotulado só se vê
   olhando. Ajuste com `--top`, `--xlim`, `--offset <número>`, `--size`.
   Atribuições na figura: `--assign-labels` escreve o 1º candidato automático e serve só
   como rascunho (ele erra - ex.: chama a combinação aromática 1774 do poliestireno de
   "anidrido"). Na figura final use as atribuições que você revisou:
   `--label-map 1715="ν(C=O)" 1240="νas(C–O–C)" ...` (só essas bandas são rotuladas, na
   banda detectada mais próxima, ±10 cm⁻¹).
6. **Entregue**: figura (PNG 300 dpi; SVG/PDF se for para artigo - basta mudar a extensão
   de `-o`), tabela de bandas (CSV) e a interpretação em texto, com as ressalvas.

## Regras que evitam erros comuns

- **Posições**: inteiros em cm⁻¹. A precisão vem do intervalo de dados e da resolução
  (tipicamente 4 cm⁻¹); "1715,37 cm⁻¹" é falsa precisão.
- **Não compare intensidades em %T** nem entre espectros com offset. Para "a banda X
  aumentou", use `band_metrics`/`compare` em absorbância, com **as mesmas janelas para
  todas as amostras**. Altura negativa em `band_metrics` = nenhuma banda resolvida na
  janela; é um resultado (ausência), reporte o sinal.
- **ATR ≠ transmissão**: no ATR as bandas de baixo número de onda ficam relativamente mais
  fortes e as posições podem deslocar alguns cm⁻¹ para menos. Ao comparar com biblioteca
  de transmissão (NIST, KBr), considere `atr_correction` e diga isso. Pergunte/registre o
  acessório e o cristal (diamante, ZnSe, Ge) - o `.sp` às vezes traz em `accessory`.
- **Artefatos antes de química**: dubleto 2360/2340 = CO₂; linhas finas em série em
  3900-3500 e 1900-1300 = vapor d'água; 3440 + 1640 largos em pastilha de KBr = água do
  KBr; 1384 fina e isolada em KBr = nitrato contaminante; ruído 2300-1900 = cristal de
  diamante. `--no-co2` só para a figura, e diga na legenda que foi feito.
- **Uma banda sozinha não prova um grupo**. Exija o conjunto (ex.: éster = C=O ~1740 +
  C–O 1300-1150 + 1150-1000; ácido = C=O ~1710 + O–H 3300-2500 muito largo + ~930) e diga
  quando a evidência é parcial. Fase cristalina, quantificação e identidade de um composto
  pedem técnica complementar (DRX, padrão, biblioteca).
- **Arquivos com vários espectros** (CSV multi-coluna, CSV "largo" com espectros em linhas,
  JCAMP composto): `load_all` devolve todos; `--index N` escolhe um na CLI.

## Formatos

| Formato | Lido | Observações |
|---|---|---|
| JCAMP-DX 4.24/5.x (.jdx .dx .jcamp) | XYDATA (AFFN, PAC, SQZ, DIF, DUP), XYPOINTS, PEAK TABLE, blocos múltiplos | XUNITS em µm/nm convertidas para cm⁻¹; TRANSMITTANCE com valores 0-100 vira %T |
| CSV/TXT/DPT/PRN/XY/ASC | separador `; , tab espaço`, vírgula decimal, cabeçalho, UTF-8/16/cp1252, várias colunas, layout largo | x em µm detectado (valores < 30); tipo de y pelo cabeçalho ou inferido |
| PerkinElmer .sp ('PEPE') | eixo, nº de pontos (limites **inclusivos**), dados float64, rótulos x/y, metadados de texto | pad inicial descartado se existir |
| Thermo .spa, Bruker OPUS (.0, .1…) | não | peça exportação para CSV/DPT/JCAMP |

Detalhes e armadilhas de cada formato: `references/formats.md`.

## Arquivos de referência

- `references/interpretation.md` - roteiro de interpretação por regiões, combinações que
  confirmam grupos, armadilhas de ATR/KBr, como redigir a interpretação. **Ler antes de
  interpretar.**
- `references/materials.md` - bandas-assinatura de polímeros, biomoléculas e inorgânicos
  comuns (PE, PP, PS, PET, PMMA, PVC, PLA, PA, PU, epóxi, silicone, celulose, quitosana,
  proteínas, carbonatos, fosfatos, sílica, óxidos...).
- `references/band_table.tsv` - tabela de correlação usada por `ftir_assign` (faixa,
  atribuição, intensidade, forma, peso, banda parceira exigida). Pode ser editada; a coluna
  `requires` aceita `hi-lo`, alternativas com `|` e qualificadores `:s :m :w :br`.
- `references/formats.md` - especificação resumida dos formatos e como depurar leitura.
