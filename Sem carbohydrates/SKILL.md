---
name: sem-carbohydrates
description: "Analisa imagens de MEV/SEM (microscopia eletrônica de varredura) de grânulos e partículas de carboidratos (amido nativo/modificado/poroso, microcápsulas, quitosana, celulose) e de filmes de carboidratos (superfície e corte transversal/fratura). Lê TIF de Thermo/FEI/Phenom, Zeiss, TESCAN, JEOL e Hitachi (escala pelos metadados ou pelo .txt que acompanha a imagem) e de qualquer origem pela barra de escala desenhada (detecção + OCR). Use sempre que o usuário trouxer imagens de MEV (.tif, .tiff, .jpg, .png) ou pedir tamanho, distribuição de tamanho (D10/D50/D90), forma, circularidade, razão de aspecto de grânulos; grânulos A e B (trigo, cevada); poros e pites de hidrólise enzimática; defeitos, poros, trincas, grânulos não gelatinizados, homogeneidade ou textura de filmes; espessura no corte transversal; painéis de micrografias com barra de escala para artigo; ou interpretação - mesmo sem dizer 'MEV' (ex.: 'micrografia do amido', 'imagem do microscópio eletrônico do filme', 'morfologia dos grânulos', 'SEM image')."
---

# MEV de carboidratos: grânulos, partículas e filmes

Scripts em `scripts/` (numpy, scipy, pandas, matplotlib, scikit-image, opencv, tifffile,
scikit-learn; `pytesseract` + programa `tesseract` só para ler a barra de escala por OCR). Rode de
qualquer pasta:

```bash
S=<caminho-desta-skill>/scripts
python $S/sem.py info *.tif                         # escala e de onde veio, kV, aumento, WD, detector, qualidade
python $S/sem.py granules amido_*.tif -o gran.png --csv gran.csv --min-d 2um     # várias imagens = 1 amostra
python $S/sem.py granules hidrolisado_*.tif --pits -o pites.png --csv hid.csv     # + pites/poros superficiais
python $S/sem.py compare nativo.csv hidrolisado.csv --names Nativo Hidrolisado --col ecd -o comp.png
python $S/sem.py film filme_superficie.tif -o sup.png --csv defeitos.csv          # poros, partículas, trincas, textura
python $S/sem.py section filme_corte.tif -o corte.png --csv espessura.csv         # espessura ao longo do filme
python $S/sem.py panel a.tif b.tif c.tif d.tif -o fig.png --titles Nativo HMT Ácido Enzima --ncols 2
python $S/sem.py scale imagem.tif --bar-value 10um -o calibrada.tif   # TIF sem databar calibrado p/ ImageJ
python $S/sem.py report amido.tif --kind granules -d saida             # análise + figura + CSV + relatorio.md
```

Opções comuns: `--pixel-size 0.05um`, `--scale 10um=235px`, `--bar-value 10um` (mede a barra desenhada
na imagem), `--databar N` (linhas a remover; `0` = não remover), `--crop x0,y0,x1,y1`, `--lang en`
(figuras em inglês), `--sep ';' --decimal ','` (CSV para Excel em português).

Em Python: `sys.path.insert(0, S)` e use `sem_io.load` (objeto com `data`, `px_um`, `px_source`,
`meta`, `scalebar`, `warnings`), `sem_process` (quality, flatten, denoise, enhance, threshold),
`sem_granules` (segment, measure, size_stats, bimodal, shape_summary, pits, compare), `sem_film`
(surface_texture, surface_defects, cross_section, measure_line) e `sem_plot` (show, scalebar,
plot_granules, plot_film, plot_section, panel, plot_compare). As docstrings explicam tudo.

## Fluxo de trabalho

1. **`info` primeiro, sempre.** Confira e diga ao usuário:
   - **Escala e de onde veio.** Metadados (FEI/Thermo, Zeiss, TESCAN, ImageJ) ou `.txt` do JEOL/Hitachi
     com o mesmo nome do TIF são confiáveis. "barra de escala + OCR - CONFIRA" precisa de
     confirmação visual do valor (abra a imagem). Sem escala, peça o valor da barra (`--bar-value`)
     - nunca invente aumento → tamanho de pixel (o aumento nominal depende da tela).
   - Aviso de que metadados e barra **discordam** = imagem redimensionada/recortada depois da
     aquisição: pergunte qual é o original.
   - **Qualidade**: saturação/riscas (carregamento), gradiente de iluminação (sombra do detector),
     SNR baixo, kV alto para carboidrato. Isso limita o que pode ser medido; diga.
2. **Escolha a análise pelo objetivo** e leia a referência antes de interpretar:
   - grânulos/partículas: `granules` + `references/starch_morphology.md` (tamanho e forma típicos
     por fonte botânica, efeito de cada modificação) e `references/methods.md` (definições);
   - filmes: `film` (superfície) e `section` (corte) + `references/films.md`;
   - preparo, aquisição, artefatos, quantos grânulos medir e como reportar →
     `references/imaging.md`; arquivos e escala → `references/formats.md`.
3. **Olhe a segmentação antes de acreditar nos números.** A figura de `granules` mostra cada objeto
   contornado pela classe de forma (cinza = tocando a borda, fora da estatística). Ajuste e refaça se:
   grânulos encostados saíram juntos (`--split-h 0.15`), um grânulo foi partido (`--split-h 0.4`
   ou `--no-split`), fundo claro/irregular (`--thresh local` ou `--flatten`), objetos escuros
   (`--dark-objects`), sujeira/ruído contado (`--min-d 2um`). Grânulos muito sobrepostos ou em
   camadas não se medem por imagem 2D - diga isso e peça imagens com o pó mais disperso.
4. **Estatística honesta.** Uma amostra = várias imagens/campos (passe todas juntas no mesmo comando
   `granules`); mire ≥ 300-500 grânulos válidos. O script avisa abaixo disso. Reporte média ± DP,
   mediana, D10/D50/D90 (número) e Dv50/D[4,3] (volume, para comparar com difração laser) e diga
   quantos grânulos e imagens. Comparação entre tratamentos: `compare` (Kruskal-Wallis/Mann-Whitney
   com Holm); com réplicas de preparo, a unidade estatística é a réplica, não o grânulo.
5. **Figura**: as figuras já saem com barra de escala limpa (a databar original é removida), letras
   (a), (b)... e cores acessíveis. **Abra a imagem e olhe** antes de entregar.
6. **Entregue** número + método + n: "D50 = 13,2 µm (ECD; n = 412 grânulos em 8 campos; Dv50 =
   15,8 µm)", "grânulos poliédricos (61%) e esféricos (29%), AR 1,18 ± 0,10", "espessura 86,5 ± 3,0
   µm (1266 posições em 1 fratura)", "35 poros/0,049 mm² (712 mm⁻²), ECD 2,6 ± 1,1 µm".

## Regras que evitam erros comuns

- **MEV é projeção 2D da face visível.** Lenticulares (trigo A) de perfil parecem alongados;
  ECD/volume assumem esfera. Diga isso ao comparar com difração laser (que dá Dv maior, em volume).
- **Objetos cortados pela borda** não entram na estatística (tamanho truncado) - por isso aumento
  alto demais reduz o n. Grânulo pequeno demais (< ~10 px de diâmetro) não é medido com confiança.
- **Brilho de MEV não é altura nem composição** (em SE): "rugosidade" e textura de superfície são
  índices comparativos entre imagens adquiridas nas mesmas condições. Rugosidade real → AFM.
- **Pites/poros superficiais** (`--pits`): contagem da face visível; use sempre o amido nativo,
  nas mesmas condições, como controle - o centro escuro de grânulos pequenos e o hilo podem gerar
  falsos positivos.
- **Defeitos de filme** são candidatos: o overlay precisa ser conferido (bordas de partículas,
  sujeira, riscos de manuseio e trincas do próprio feixe parecem defeitos).
- **Espessura por MEV** é de um campo, filme seco, sob vácuo, com fratura possivelmente inclinada:
  confira com micrômetro e reporte os dois. Resolução ≈ 2 px por borda.
- **Não diagnostique fonte botânica só pela imagem** quando as faixas se sobrepõem: use
  `starch_morphology.md` como plausibilidade, não como identificação.
- **Números**: diâmetros com 3 algarismos significativos (ou 1 casa em µm), forma com 2 casas,
  porcentagens inteiras, sempre com n.

## Formatos e escala

| Origem | Escala lida de | Observações |
|---|---|---|
| Thermo Fisher/FEI/Phenom (Quanta, Inspect, Helios, Prisma, Apreo, Phenom) | tag 34682 (`[Scan] PixelWidth`) | databar removida pela altura gravada |
| Zeiss (SmartSEM: EVO, Sigma, Gemini, Merlin, Auriga) | tag 34118 (`Image Pixel Size`) | data zone detectada e removida |
| TESCAN (Vega, Mira, Maia, Clara) | tag 50431 (`[MAIN] PixelSizeX`) | |
| JEOL (JSM-6xxx, IT/IT-x00, InTouchScope) | `.txt` com o mesmo nome (`$$SM_MICRON_BAR` + `MARKER`) | sem `.txt` → barra + OCR |
| Hitachi (S-, SU-, TM-, Regulus, FlexSEM) | `.txt` com o mesmo nome (`[SemImageFile] PixelSize`) | sem `.txt` → barra + OCR |
| TIF/PNG/JPG de qualquer origem | barra desenhada (+ OCR) ou `--bar-value`/`--scale`/`--pixel-size` | TIF calibrado no ImageJ também é lido |

Detalhes e depuração: `references/formats.md`.

## Arquivos de referência

- `references/starch_morphology.md` - tamanho/forma típicos por fonte (milho, mandioca, batata,
  trigo, arroz, leguminosas, banana, inhame, pinhão...), grânulos A/B, poros naturais, efeito de
  gelatinização, hidrólise ácida/enzimática, HMT, annealing, modificações químicas, moagem,
  ultrassom; partículas de spray-drying/liofilização, quitosana, celulose.
- `references/films.md` - leitura de superfície e fratura de filmes de amido, quitosana, alginato,
  pectina, CMC, blendas, compósitos e emulsionados; o que cada defeito indica.
- `references/imaging.md` - preparo (dispersão do pó, recobrimento, criofratura), kV/detector/baixo
  vácuo, artefatos (carregamento, dano do feixe), aumento certo, quantas imagens e grânulos, como
  reportar no artigo.
- `references/methods.md` - definições (ECD, Feret, AR, circularidade, solidez), algoritmos,
  limites, distribuições em número e volume, bimodalidade, textura (GLCM, PSD), espessura.
- `references/formats.md` - metadados por fabricante, barra de escala/OCR, databar, depuração.
