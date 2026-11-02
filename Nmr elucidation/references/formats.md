# Arquivos de RMN

## JCAMP-DX (preferido)

- Versões 4.24, 5.x e 6.0; `##DATA TYPE= NMR SPECTRUM` (processado) ou `NMR FID` (bruto).
- `##DATA CLASS= XYDATA` (uma coluna real) ou `NTUPLES` (real + imaginário em páginas; Bruker/MestReNova).
- Eixo: `##XUNITS= HZ` ou `PPM`. Em Hz, o deslocamento sai de `##$OFFSET`/`##$SF` (Bruker),
  `##.SHIFT REFERENCE= INTERNAL, solvente, ponto, ppm` ou Hz/frequência observada (JEOL Delta,
  MestReNova). Blocos `LINK` (MestReNova): usa o espectro; lê `PEAK TABLE` e ignora atribuições.
- Compressões ASDF (SQZ/DIF/DUP) via nmrglue; tabela de picos densa (XY..XY) também é lida.

## Como exportar

| Software | Caminho |
|---|---|
| Bruker TopSpin | comando `tojdx` (escolha "processed" para o espectro ou "FID"); ou File → Save → JCAMP-DX. Também pode enviar a pasta do experimento (fid + acqus + pdata) |
| MestReNova | File → Save As → JCAMP-DX (*.jdx) (espectro processado; inclui real e imaginário) |
| JEOL Delta | File → Export → JCAMP-DX |
| Agilent/Varian VnmrJ | enviar a pasta .fid (fid + procpar) ou exportar JCAMP |
| Magritek Spinsolve / Nanalysis / Oxford X-Pulse | o software grava JCAMP-DX (FID ou espectro) |
| NMRium | Export → JCAMP-DX |
| ACD/Labs | Export → JCAMP-DX |

Para ¹³C, envie também os DEPT-135/90 (ou APT) do mesmo tubo/sessão.

## CSV/TXT

Duas colunas (ppm, intensidade), separador vírgula/ponto e vírgula/tab; lista curta (< 400 linhas) é
tratada como lista de picos. Informe `--sf` (MHz) e `--solvent` (o arquivo não diz).

## Depuração

- `info` mostra eixo usado, referência, fase, ruído e resolução digital.
- Espectro "espelhado" (aromáticos em 1-2 ppm): eixo em Hz interpretado errado → informe o campo/`--ref`.
- Sinais todos negativos/dispersivos: fase; a skill tenta refasar, senão processe no software.
- `process fid.jdx -o espectro.jdx` grava o espectro como JCAMP-DX 5.01 (AFFN, eixo em ppm) para abrir em
  qualquer programa.
