---
name: thermal-analysis
description: "Lê, processa, plota e interpreta análise térmica: TG/TGA, DTG, DSC, DTA e STA (TG-DSC/TG-DTA simultâneos). Lê binário NETZSCH .ngb-ss3/.ngb-bs3 e exportações ASCII de NETZSCH Proteus, TA Instruments (Universal Analysis, TRIOS, Excel), Mettler Toledo STARe, PerkinElmer Pyris, Setaram, além de CSV/TXT/Excel genéricos (vírgula decimal). Use sempre que o usuário trouxer esses arquivos ou pedir etapas de perda de massa, onset/endset, T5%/Td, resíduo/teor de cinzas/carga, estequiometria de decomposição ou nº de águas de hidratação, Tg, Tm, Tc, cristalização a frio, ΔH, grau de cristalinidade (Xc), OIT, Cp com safira, cinética de decomposição (Friedman, Ozawa-Flynn-Wall, KAS, Starink, Vyazovkin, Kissinger, Ea), ou figuras TG/DTG/DSC para artigo/tese - mesmo sem dizer 'análise térmica' (ex.: 'curva termogravimétrica', 'estabilidade térmica do polímero', 'fusão da amostra no DSC', 'termograma', 'calorimetria')."
---

# Análise térmica: TG/DTG, DSC, DTA e STA

Scripts em `scripts/` (numpy, scipy, pandas, matplotlib; `pyngb` só para binário NETZSCH,
`xlrd`/`openpyxl` só para Excel). Rode de qualquer pasta:

```bash
S=<caminho-desta-skill>/scripts
python $S/ta.py info arquivo.txt                       # formato, tipo, massa, segmentos, direção endo
python $S/ta.py report arquivo.ngb-ss3 -d saida        # análise completa conforme o tipo (TGA/DSC/STA)
python $S/ta.py tga oxalato.txt -o tg.png --csv etapas.csv --theory "CaC2O4·H2O:H2O,CO,CO2"
python $S/ta.py tga amostra.txt --basis 150            # base seca: 100% = massa a 150 °C
python $S/ta.py dsc PET.txt --segment heat1 heat2 -o dsc.png --polymer PET   # automático
python $S/ta.py dsc PET.txt --segment heat2 --tg 60:100 --peak 110:150 --peak 220:270 --polymer PET
python $S/ta.py sta arquivo.ngb-ss3 -o sta.png         # TG/DTG + DSC: caráter endo/exo de cada etapa
python $S/ta.py kinetics b05.txt b10.txt b20.txt b40.txt --window 250:450 -o ea.png --csv ea.csv
python $S/ta.py oit oit.txt -o oit.png                 # troca N2→O2 detectada no canal de gás
python $S/ta.py cp --sample s.txt --blank b.txt --sapphire saf.txt --range 50:300 -o cp.csv
python $S/ta.py stoich "CuSO4·5H2O" --release 2H2O 2H2O H2O       # perdas teóricas por etapa
python $S/ta.py hydrate --loss 14.4 --anhydrous CuSO4             # nº de águas pela perda
python $S/ta.py plot a.txt b.txt c.txt -o comp.png --names A B C  # sobreposição TG ou DSC
python $S/ta.py convert arquivo.ngb-ss3 -o dados.csv [--decimal ,] # exporta em unidades canônicas
```

Opções comuns: `--endo up|down` (direção endotérmica do sinal no arquivo), `--display
exo_up|endo_up` (como desenhar), `--mass mg` (se o arquivo não traz), `--columns
time:min,T:°C,mass:mg` (arquivo sem cabeçalho), `--correction corr.ngb-bs3` (NETZSCH),
`--index n` (arquivo com várias curvas), `--lang en` (figuras/saída em inglês).

Em Python: `sys.path.insert(0, S)` e use `ta_io` (load, load_all, save_csv - objeto `Run`
com `r.time` [min], `r.T` [°C], `r.mass` [mg], `r.mass_pct`, `r.hf`, `r.units`, `r.extra`,
`r.meta`), `ta_process` (segments, select, heating_rate, hf_specific, endo_direction,
endo_signal, mass_pct, dtg, subtract_blank, smooth), `ta_tga` (find_steps, char_temps,
residues, theoretical_losses, hydration_number, compare_with_theory, molar_mass),
`ta_dsc` (find_events, glass_transition, integrate_peak, crystallinity, DH0, oit,
cp_sapphire), `ta_kinetics` (conversion, isoconversional, kissinger) e `ta_plot`
(plot_tga, plot_dsc, plot_sta, plot_kinetics, plot_oit). As docstrings explicam tudo.

## Fluxo de trabalho

1. **`info` primeiro.** Confira e diga ao usuário o que não bate:
   - **Tipo** (TGA, DSC, DTA, STA), **massa da amostra** e **unidade do sinal** (W/g, mW,
     µV, µV/mg, K). DSC em mW sem massa não dá J/g → peça a massa (`--mass`).
   - **Segmentos do programa** (`heat1`, `cool1`, `iso1`, `heat2`...) com taxa medida.
     Polímeros: o 1º aquecimento mostra a história térmica/processamento, o 2º a
     propriedade do material - reporte o que o objetivo pede e diga qual foi.
   - **Direção endo** e de onde veio. A ordem de confiança é: `--endo` do usuário → rótulo
     explícito no arquivo → dedução pelos dados (aquecimento vs resfriamento) → sinal da
     calibração NETZSCH → convenção do fabricante (TA, Mettler, NETZSCH: endo ↓;
     PerkinElmer: endo ↑). Se terminar em "(confirmar)" ou "desconhecida", **confira com
     um evento conhecido** (fusão é sempre endo; cristalização no resfriamento é exo) e
     force com `--endo`. Errar isso troca fusão por cristalização.
   - NETZSCH `.ngb-ss3` "sem correção": o sinal DSC/massa não teve a linha de base
     subtraída - peça o `.ngb-bs3` correspondente e use `--correction` (o pyngb recusa
     correções com programa de temperatura diferente, e está certo).
2. **Escolha a análise pelo objetivo** e leia a referência antes de concluir:
   - perda de massa, estabilidade, composição, hidratos → `tga` + `references/tga.md`;
   - Tg, fusão, cristalização, ΔH, Xc, OIT, Cp, pureza, cura → `dsc`/`oit`/`cp` +
     `references/dsc.md`;
   - TG+DSC/DTA simultâneos → `sta` (cada etapa de massa com seu evento térmico);
   - Ea, cinética, previsão → `kinetics` + `references/kinetics.md` (**obrigatório**: as
     recomendações ICTAC mudam o que pode ser afirmado);
   - valores de referência (ΔH⁰ 100% cristalino, Tg/Tm de polímeros, decomposições
     clássicas, padrões de calibração) → `references/reference_data.md`.
3. **No DSC, os eventos automáticos são candidatos.** O modo automático acha picos e
   degraus, sugere o que podem ser e integra com janelas limitadas pelos vizinhos. **Olhe
   a figura**; se a janela cortar o pico, pegar um transiente ou juntar eventos, refaça
   com `--tg T1:T2` / `--peak T1:T2` (janelas com linha de base dos dois lados). Linhas de
   base muito curvas (amostras que perdem massa, transientes de início) geram falsos Tg e
   picos: eventos marcados `borda (transiente?)` não são quantificados.
4. **Figura**: as funções já seguem convenções (TG em % com DTG no eixo direito; DSC com
   seta `exo ↑` ou `endo ↑`, segmentos empilhados, áreas integradas sombreadas, Tg
   marcado; transientes do início/fim de segmento escondidos). **Abra a imagem e olhe**
   antes de entregar.
5. **Entregue** número + método + condições: "Etapa 1: 12,2% (teórico 12,3% para 1 H₂O),
   T_onset 176 °C, T_pico DTG 206 °C, N₂, 10 K/min", "Tg = 78 °C (ponto médio, 2º
   aquecimento, 10 K/min, ΔCp 0,32 J g⁻¹ K⁻¹)", "Xc = 27,8% (ΔHm − ΔHcc, ΔH⁰ = 140 J/g)",
   "Ea = 180 ± 3 kJ/mol (Friedman, 0,1 ≤ α ≤ 0,9, 4 taxas)".

## Regras que evitam erros comuns

- **Atmosfera, taxa, massa e cadinho** fazem parte do resultado: T_onset e T_pico sobem
  com β; reporte sempre com eles. Comparação entre amostras só nas mesmas condições.
- **TG**: Onset é por tangentes (ISO 11358) e depende da escolha; T5% é mais reprodutível
  para "estabilidade". Perda até ~150 °C geralmente é umidade/solvente: use `--basis 150`
  para expressar em base seca. Ganhos de massa pequenos no início costumam ser empuxo
  - subtraia um branco (`--blank`, cadinho vazio, mesmo programa).
- **Etapas sobrepostas** (`sobreposta = True`): a divisão é pelo mínimo da DTG e a
  quantificação de cada uma é aproximada. Diga isso; não force estequiometria.
- **Resíduo ≠ carga inorgânica** sem cuidado: em N₂ polímeros deixam carvão; em ar a
  carga pode oxidar/decompor (CaCO₃ → CaO perde 44% acima de ~600 °C).
- **DSC**: ΔH só com linha de base dos dois lados do pico; Tg pelo ponto médio (meia
  altura do ΔCp) com a taxa; pico de relaxação entálpica sobre o Tg (1º aquecimento) não
  é fusão. Cristalinidade: Xc = (ΔHm − ΔHcc)/(ΔH⁰·w) - desconte a cristalização a frio e
  divida pela fração do polímero w em compósitos/blendas; cite o ΔH⁰ usado (varia na
  literatura).
- **Massa perdida durante o evento invalida ΔH em J/g** (desidratação, decomposição):
  diga isso, ou use STA. Picos de decomposição no DSC de cadinho fechado podem ser
  artefatos de pressão.
- **Cinética**: ≥ 3 taxas (ideal 4-5), mesma massa/cadinho/gás; nunca Ea de uma única
  taxa (Coats-Redfern etc.). Se Ea varia > 20% com α, o processo é multietapa: reporte
  Ea(α), não um valor só. Kissinger é estimativa.
- **Números**: temperaturas com 1 casa no máximo (inteiro costuma bastar), perdas com
  1 casa (%), ΔH com 3 algarismos significativos, Ea inteiro (kJ/mol) com o IC.

## Formatos

| Formato | Lido | Observações |
|---|---|---|
| NETZSCH `.ngb-ss3/.ngb-bs3/.ngb-ds3` | sim (`pip install pyngb`) | massa, DSC (µV→mW pela sensibilidade do arquivo), gases, segmentos; `--correction` para a correção |
| NETZSCH Proteus ASCII (`#` cabeçalho, `##` colunas) | sim | SEPARATOR/DECIMAL lidos do cabeçalho |
| TA Universal Analysis (`Nsig`/`StartOfData`) | sim | UTF-16; linhas-sentinela ignoradas |
| TA TRIOS `.txt` e Excel `.xls/.xlsx` | sim | `[step]` e planilhas viram segmentos (`.xls` precisa de `xlrd`) |
| Mettler STARe ASCII | sim | várias curvas (TGA + DSC) fundidas em STA; branco detectado |
| PerkinElmer Pyris ASCII | sim | blocos "N) etapa" viram segmentos; endo ↑ |
| Setaram (Calisto/SETSOFT) | sim | UTF-16, vírgula decimal, nomes com espaço |
| CSV/TXT/Excel genérico | sim | unidades no nome `(mg)`, `[°C]`, `/min` ou em linha própria; sem cabeçalho → `--columns` |
| TA `.tri`, `.001` binários, Mettler `.stare`, Pyris `.dsc/.tga` binários, Shimadzu `.tad` | não | exporte como texto no software (ver `references/formats.md`) |

## Arquivos de referência

- `references/tga.md` - leitura da curva TG/DTG, onset/endset, T5%, resíduo, base seca,
  branco/empuxo, estequiometria, hidratos, polímeros, biomassa, compósitos, atmosfera.
- `references/dsc.md` - Tg (métodos, ΔCp, relaxação), fusão/cristalização, Xc, cristalização
  a frio, polimorfismo, pureza, cura, OIT, Cp com safira, DTA, calibração, armadilhas.
- `references/kinetics.md` - recomendações ICTAC, métodos isoconversionais, Kissinger,
  dependência Ea(α), o que não fazer, como reportar.
- `references/reference_data.md` - ΔH⁰ de polímeros, Tg/Tm/Td típicos, decomposições
  clássicas com perdas teóricas, padrões de calibração (T, ΔH, Curie).
- `references/formats.md` - formatos, como exportar de cada software, depuração da leitura.
