# Arquivos de MEV, escala e depuração

## Sumário
1. Onde está a escala, por fabricante
2. Databar e barra de escala desenhada
3. OCR
4. Imagens exportadas (JPG/PNG/BMP) e editadas
5. Depuração

## 1. Onde está a escala, por fabricante

| Fabricante | Onde | Chave | Unidade |
|---|---|---|---|
| Thermo Fisher / FEI / Phenom (xT, Quanta, Inspect, Nova, Helios, Scios, Apreo, Prisma, Phenom) | TIF, tag 34682 (texto INI) | `[Scan] PixelWidth` (+ `[EBeam] HV, WD, HFW`, `[Detectors] Name`, `[Vacuum] UserMode`) | m |
| Zeiss (SmartSEM: EVO, Sigma, Gemini, Merlin, Auriga, Crossbeam) | TIF, tag 34118 (CZ_SEM) | `ap_image_pixel_size` ou `ap_pixel_size` (+ `ap_actualkv`, `ap_wd`, `ap_mag`, detector) | nm/µm conforme o campo |
| TESCAN (Vega, Mira, Lyra, Maia, Clara, Amber) | TIF, tag 50431 (texto `[MAIN]`) | `PixelSizeX` (+ `HV`, `WD`, `Detector`) | m |
| JEOL (JSM-6xxx, 7xxx, IT-100/200/300/500/700, NeoScope) | **arquivo `.txt` com o mesmo nome** | `$$SM_MICRON_BAR` (px da barra) + `$$SM_MICRON_MARKER` (ex.: `10um`); `$CM_FULL_SIZE` (resolução de referência), `$CM_MAG`, `$CM_ACCEL_VOLT` | - |
| Hitachi (S-3400/4800, SU-xxxx, TM-3000/4000, FlexSEM, Regulus) | **arquivo `.txt` com o mesmo nome** (`[SemImageFile]`) | `PixelSize` (+ `DataSize`, `Magnification`, `AcceleratingVoltage`, `WorkingDistance`) | nm |
| ImageJ/Fiji (TIF salvo calibrado) | tags ImageJ | `unit` + XResolution | unidade indicada |
| OME-TIFF, Olympus SIS | metadados próprios | PhysicalSizeX / pixelsizex | |

Copie o `.txt` junto do `.tif` (mesmo nome) para JEOL e Hitachi - sem ele a escala sai só da barra.
Se a imagem foi salva em resolução diferente da de aquisição, o script corrige a escala do JEOL por
`$CM_FULL_SIZE` e a do Hitachi por `DataSize`.

## 2. Databar e barra de escala desenhada

- **Databar** = faixa com texto (kV, aumento, barra) na base da imagem. FEI/Thermo: altura vem dos
  metadados. Demais: detectada como linhas quase uniformes na base (ou no topo) com texto. Sempre é
  removida antes da análise (não é amostra!). Se a detecção falhar: `--databar N` (linhas) ou
  `--crop`. Barra/legenda **sobreposta à imagem** (Phenom, alguns Hitachi TM, imagens editadas):
  recorte com `--crop` antes de segmentar.
- **Barra de escala**: procurada dentro da databar como linha horizontal fina (inteira, dividida
  pelo rótulo como no FEI, com ticks nas pontas como no Zeiss, ou tracejada como em alguns
  Hitachi/JEOL). O comprimento em px vem da detecção; o valor vem do OCR ou de `--bar-value`.
- Com metadados **e** barra, o script compara os dois e avisa se diferirem > 3% (> 10% para barra
  tracejada) - sinal de imagem redimensionada ou de OCR errado.

## 3. OCR

- Usa `pytesseract` + o programa `tesseract` (Linux: `apt install tesseract-ocr`; Windows: instalador
  UB-Mannheim; `pip install pytesseract`). Sem eles, a skill pede `--bar-value`.
- Isola o rótulo mais próximo da barra, amplia em 3 escalas e vota; corrige leituras típicas (`pm`,
  `prn`, `urn`, `y` → µm; `S` → 5; `O` → 0). Fontes muito pequenas (< 8 px) são o caso difícil.
- O resultado diz "CONFIRA o valor lido": abra a imagem e confirme o número; se estiver errado,
  `--bar-value 20um`.

## 4. Imagens exportadas (JPG/PNG/BMP) e editadas

- Exportações com "data bar" queimada funcionam (barra + OCR), mas perdem os metadados (kV, WD...).
- **Nunca** use imagens redimensionadas, recortadas em programa de edição ou capturas de tela sem
  recalibrar pela barra visível nelas.
- JPG comprime e cria artefatos em bordas finas: prefira TIF para medir.
- TIF com compressão JPEG precisa de `pip install imagecodecs` (o script tenta o Pillow antes).
- RGB/colorido (falsa cor, overlay de EDS) é convertido para cinza por luminância.

## 5. Depuração

- `info` mostra: dimensões antes/depois da databar, a escala e a fonte, a barra detectada e o texto
  lido, kV/aumento/WD/detector, qualidade (SNR, saturação, gradiente) e avisos.
- `scale arquivo -o calibrado.tif` grava o TIF sem databar já calibrado para conferir no ImageJ
  (Analyze → Set Scale mostra a escala; meça a barra).
- `measure arquivo --points x1,y1,x2,y2` mede uma distância (px da imagem sem databar) - útil para
  conferir um grânulo ou a espessura em um ponto.
