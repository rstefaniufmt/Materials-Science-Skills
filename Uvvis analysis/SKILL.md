---
name: uvvis-analysis
description: "Lê, processa, plota e interpreta espectros UV-Vis(-NIR) e de reflectância difusa (DRS): JCAMP-DX, CSV/TXT (Cary, Shimadzu, OceanView, vírgula decimal), Excel, Galactic .spc e PerkinElmer .sp. Use sempre que o usuário trouxer esses arquivos ou pedir λmax, curva de calibração/Lei de Beer (LOD, LOQ, ε, concentração de amostra), cinética de degradação/fotocatálise (C/C0, k, t½), band gap óptico (Tauc, Kubelka-Munk, Urbach), pontos isosbésticos, ressonância plasmônica de nanopartículas, verificação de comprimento de onda com filtro de hólmio, ou figuras de espectros UV-Vis para artigo/tese - mesmo sem dizer 'UV-Vis' (ex.: 'espectro de absorção', 'curva analítica do azul de metileno', 'gap do TiO2', 'espectrofotômetro')."
---

# UV-Vis: leitura, quantificação, cinética, band gap e figuras

Scripts em `scripts/` (numpy, scipy, pandas, matplotlib; openpyxl só para Excel). Rode de
qualquer pasta:

```bash
S=<caminho-desta-skill>/scripts
python $S/uvvis.py info arquivo.csv                  # formato, faixa, ordenada, avisos, degraus de detector
python $S/uvvis.py report arquivo.spc -d saida --lang pt        # λmax + figura + tabela + resumo
python $S/uvvis.py plot *.txt -o fig.png --series 0 15 30 60 --series-label "t (min)" --label --top 2
python $S/uvvis.py calib padroes.csv --conc 0.5 1 2 3 4 5 --wl 664 --base 800 --unit mg/L \
       --unknown amostra.csv --dilution 10 -o curva.png --lang pt
python $S/uvvis.py kinetics t*.txt --times -30 0 15 30 60 90 --base 800 --model all -o cin.png
python $S/uvvis.py tauc DRS.csv --transition indirect -o tauc.png      # %R -> Kubelka-Munk automático
python $S/uvvis.py isosbestic serie.xlsx --region 350 750
python $S/uvvis.py wlcheck holmio.spc                    # desvio de calibração do comprimento de onda
```

Em Python: `sys.path.insert(0, S)` e use `uvvis_io.load/load_all`, `uvvis_process`
(offset_correct, scatter_correct, baseline, detect_steps, splice_correct, smooth, derivative,
normalize, subtract, average), `uvvis_peaks` (find_bands, band_metrics, isosbestic,
wavelength_check), `uvvis_quant` (calibrate, predict, kinetics, compare_kinetic_models),
`uvvis_bandgap` (tauc, urbach, derivative_edge, onset_wavelength) e `uvvis_plot`
(plot_spectra, plot_tauc, plot_calibration, plot_kinetics). As docstrings explicam tudo.

## Fluxo de trabalho

1. **`info` primeiro.** Confira:
   - `Ordenada`: A, %T, %R, KM, counts... `(INFERIDO - confirmar)` = o arquivo não dizia.
     %T e %R têm a mesma faixa numérica: se a amostra é pó/filme opaco medido com esfera
     integradora, é **%R** - force com `--ykind %R` (Python: `s.ykind = "%R"`). Errar isso
     invalida Kubelka-Munk e Tauc.
   - Avisos de A > ~2.5-3 (luz espúria/saturação: dilua ou use caminho menor) e linhas
     `DEGRAU` (troca de detector ~800-870 nm ou lâmpada ~320-350 nm) - corrija com
     `splice_correct(s, at=λ)` antes de qualquer análise que cruze esse ponto.
   - Faixa útil: abaixo do corte do solvente/cubeta o sinal é lixo (tabela em
     `references/reference_data.md`).
2. **Corrija a linha de base do jeito certo para UV-Vis** (não é como no IR):
   - solução límpida: `--offset-correct 750 800` (média de uma região sem absorção vira 0);
   - turbidez/nanopartículas/suspensão: `--scatter 600 800` ajusta e remove A ∝ λ⁻ⁿ;
   - para medir uma banda numa série: `--base λ` ou `--base λ1 λ2` em `calib`/`kinetics`
     (linha de base local, mesma para todas as amostras).
3. **Escolha a ferramenta pelo objetivo** e leia a referência correspondente antes de
   concluir:
   - concentração/Beer-Lambert/LOD → `calib` + `references/quantitation.md`;
   - degradação/cinética → `kinetics` + `references/quantitation.md` (seção cinética);
   - band gap/Urbach → `tauc` + `references/bandgap.md` (**obrigatório** - as escolhas de
     transição, de F e de região mudam o resultado e precisam ser justificadas);
   - identificar/interpretar bandas → `peaks`/`report` + `references/interpretation.md`
     e `references/reference_data.md` (λmax de corantes, fármacos, ensaios, complexos,
     nanopartículas, biomoléculas, band gaps típicos).
4. **Figura**: `plot` já segue as convenções (x em nm crescente; `--ev-axis` adiciona eV no
   topo; séries com gradiente + barra de cores via `--series`; `--isosbestic` marca os
   pontos; λmax rotulados sem sobreposição - em séries, só o 1º e o último espectro, use
   `--label-all` para todos). **Abra a imagem e olhe** antes de entregar.
5. **Entregue** número + incerteza + método: "Eg = 3.19 ± 0.01 eV (Tauc indireto, F(R),
   ajuste 3.63-3.96 eV)", "c = 24.9 ± 0.1 mg/L (IC95%, n = 1, curva de 6 pontos, 664 nm)",
   "k = 0.0179 min⁻¹ (pseudo-1ª ordem, R² = 0.999)".

## Regras que evitam erros comuns

- **Quantifique em λ fixo** (`--wl`), não no "máximo da janela": o máximo desloca com
  ruído e interferentes. Use o mesmo λ e a mesma linha de base para padrões e amostras.
- **Amostra fora da faixa da curva** é extrapolação: o script marca `FORA da faixa`;
  dilua e meça de novo em vez de reportar. Abaixo do LOQ, reporte "< LOQ".
- **R² alto não prova linearidade**: veja resíduos (figura de `calib`), o teste de Mandel e
  o teste do intercepto que o script imprime.
- **Cinética de fotocatálise**: C₀ é o valor em t = 0 *após* o equilíbrio de adsorção no
  escuro (tempos negativos = fase escura; o script separa e reporta a % adsorvida). Diga
  se houve controle de fotólise (sem catalisador) e de adsorção (sem luz). A queda de A
  pode ser descoloração sem mineralização - não chame de "degradação completa".
- **Band gap**: Tauc só vale para a borda de um semicondutor; o tipo de transição deve vir
  da literatura do material (não escolha pelo R² maior). Pó → reflectância difusa →
  Kubelka-Munk (nunca "absorbância" de pó em cubeta). Se o script avisar "outra borda",
  há absorção sub-gap (dopante, defeito, sensibilizador): reporte as duas e use também
  `Eg_baseline`. Reporte Eg com 2 casas e a faixa de ajuste.
- **Posições**: λmax inteiro (nm) ou uma casa se o intervalo de dados < 1 nm. Bandas na
  borda do espectro (`flag borda`) podem ser só a subida de absorção, não um máximo.
- **Nada acima de A ≈ 3** deve ser interpretado quantitativamente; bandas "achatadas" no
  topo são saturação.

## Formatos

| Formato | Lido | Observações |
|---|---|---|
| CSV/TXT/PRN/ASC/DAT | sim | `; , tab espaço`, vírgula decimal, cabeçalho/rodapé, Cary pareado (x,y,x,y,…) com nomes na 1ª linha, Shimadzu UVProbe, OceanView (`>>>>>Begin Spectral Data<<<<<`), layout largo |
| Excel .xlsx/.xls | sim | cada planilha como tabela; nome da planilha serve de pista para o tipo de y |
| JCAMP-DX (.jdx/.dx) | sim | todas as compressões ASDF; XUNITS em nm, cm⁻¹, µm |
| Galactic .spc | sim | novo (0x4B) e antigo (0x4D); multi-subarquivo (séries, diode array) |
| PerkinElmer .sp | sim | UV WinLab/Lambda e Spectrum (limites do eixo inclusivos) |
| Cary .bsw/.dsw/.vsp, Shimadzu .vspd/.uvd | não | exporte como CSV/TXT no software |

Detalhes: `references/formats.md`.

## Arquivos de referência

- `references/bandgap.md` - teoria mínima, escolha da transição e de F, Kubelka-Munk,
  espessura, Urbach, armadilhas e como reportar. **Ler antes de qualquer band gap.**
- `references/quantitation.md` - Beer-Lambert e seus desvios, desenho da curva, LOD/LOQ,
  ε, incerteza, adição de padrão, misturas de dois componentes, cinética e controles.
- `references/interpretation.md` - transições eletrônicas, cromóforos, efeitos de solvente
  e pH, complexos (d-d, transferência de carga), agregação, SPR, controle de qualidade.
- `references/reference_data.md` - tabelas: λmax de corantes/poluentes/fármacos, ensaios
  colorimétricos, biomoléculas, complexos, nanopartículas, band gaps, cortes de solventes
  e cubetas.
- `references/formats.md` - formatos de arquivo e depuração da leitura.
