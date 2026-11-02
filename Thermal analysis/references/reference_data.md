# Dados de referência

Valores típicos de literatura, para orientar a interpretação. Dependem de grau, massa molar,
aditivos, taxa e atmosfera: cite uma fonte primária quando o número entrar no trabalho.

## Sumário
1. ΔH⁰ de fusão (100% cristalino) para Xc
2. Tg, Tm, Tc e Td de polímeros
3. Decomposições clássicas (perdas teóricas)
4. Padrões de calibração de T e ΔH (DSC/DTA)
5. Pontos de Curie (calibração de T em TG)
6. Faixas típicas em outros materiais

## 1. ΔH⁰ de fusão (100% cristalino) para Xc

Valores embutidos em `ta_dsc.DH0` (use `--dh0` para outro valor; a literatura varia).

| Polímero | ΔH⁰ (J/g) | Faixa na literatura | | Polímero | ΔH⁰ (J/g) | Faixa |
|---|---|---|---|---|---|---|
| PE (HDPE, LDPE, UHMWPE) | 293 | 277-293 | | PEEK | 130 | 130-161 |
| PP (iPP) | 207 | 165-209 | | PEKK | 130 | |
| PET | 140 | 117-140 | | PPS | 112 | 77-112 |
| PBT | 140 | 140-145 | | PVDF | 104,7 | 93-105 |
| PEN | 103 | 103-190 | | PTFE | 82 | 69-82 |
| PLA (PLLA) | 93 | 93-106 | | POM | 326 | 190-326 |
| PCL | 139,5 | 136-156 | | PEO/PEG | 196,8 | 196-213 |
| PHB | 146 | 146 | | iPS / sPS | 86 / 82 | |
| PA6 | 230 | 190-240 | | PBS | 110,3 | 110-200 |
| PA66 | 196 | 196-255 | | PVA | 138,6 | 138-163 |
| PA11 / PA12 | 206 / 209 | 189-245 | | PBAT | 114 | |

## 2. Tg, Tm, Tc e Td de polímeros

DSC a 10 K/min; Td5% em N₂ a 10 K/min (ordens de grandeza).

| Polímero | Tg (°C) | Tm (°C) | Observações | Td5% N₂ (°C) |
|---|---|---|---|---|
| HDPE / LDPE | ~ −120 (difícil de ver) | 130-137 / 105-115 | Tc ~115 / ~95 | 420-450 |
| PP | −10 a 0 | 160-170 | forma β 145-155; Tc 110-125 | 380-420 |
| PS (atático) | 95-105 | - | amorfo | 370-400 |
| PMMA | 100-120 | - | despolimeriza | 280-350 |
| PVC (rígido) | 80-85 | - | HCl a 250-350 °C | 250-280 |
| PC | 145-150 | - | | 450-480 |
| ABS | ~105 (SAN) e ~ −80 (PB) | - | | 360-390 |
| PET | 70-80 | 245-260 | cristalização a frio 120-140 (amorfo) | 390-420 |
| PBT | 40-50 | 220-225 | | 360-380 |
| PLA | 55-62 | 150-180 | cristalização a frio 90-125; pré-fusão α' | 310-340 |
| PCL | ~ −60 | 58-63 | | 350-380 |
| PHB | 0-5 | 170-180 | decompõe logo após a fusão | 250-280 |
| PA6 / PA66 | 45-60 / 50-70 | 215-225 / 255-265 | absorve água (Tg cai) | 380-420 |
| PEEK | ~143 | ~343 | | 520-560 |
| PTFE | ~115 (contestado) | ~327 | transições 19 e 30 °C | 500-550 |
| PVDF | ~ −35 | 165-175 | polimorfos α/β/γ | 420-460 |
| POM | ~ −60 | 165-178 | | 280-320 |
| PEO / PEG (alto M) | ~ −60 | 60-68 | | 340-380 |
| PVA | 75-85 | 220-230 | perde água antes | 250-300 |
| PU termoplástico | segmento mole −50 a −20 | domínios rígidos 150-220 (largo) | | 280-320 |
| Epóxi curado | 60-200 (depende da cura) | - | exo de cura 100-250 | 330-380 |
| Borracha natural (NR) | ~ −65 | - | | 330-360 |
| SBR | −50 a −55 | - | | 380-420 |

## 3. Decomposições clássicas (perdas teóricas)

Perdas relativas à massa inicial (calculadas com `theoretical_losses`). Faixas de T
aproximadas em N₂/ar, 10 K/min.

| Composto | Etapa (espécie) | Perda teórica (%) | T típica (°C) |
|---|---|---|---|
| CaC₂O₄·H₂O (padrão clássico) | H₂O / CO / CO₂ | 12,33 / 19,17 / 30,12 | 100-250 / 400-550 / 600-850 |
| CuSO₄·5H₂O | 2H₂O / 2H₂O / H₂O / SO₃ | 14,43 / 14,43 / 7,22 / 32,06 | ~90 / ~120 / ~250 / >750 |
| CaSO₄·2H₂O (gesso) | 1,5H₂O / 0,5H₂O | 15,70 / 5,23 | 100-160 / 160-220 (frequentemente sobrepostas) |
| CaCO₃ | CO₂ | 43,97 | 600-850 (sobe sob CO₂) |
| MgCO₃ | CO₂ | 52,20 | 400-600 |
| SrCO₃ | CO₂ | 29,81 | 900-1150 |
| Ca(OH)₂ | H₂O | 24,31 | 380-480 |
| Mg(OH)₂ | H₂O | 30,89 | 300-420 |
| Al(OH)₃ (gibbsita) | 1,5H₂O (total) | 34,64 | 200-350 (+ boehmita ~500) |
| Caulinita Al₂Si₂O₅(OH)₄ | 2H₂O (desidroxilação) | 13,96 | 450-650; exo ~980 (mulita/espinélio) sem massa |
| NaHCO₃ | ½H₂O + ½CO₂ (simultâneas) | 36,91 total | 100-200 |
| BaCl₂·2H₂O | 2H₂O | 14,75 | 100-170 (duas etapas de 1 H₂O) |
| FeC₂O₄·2H₂O / CoC₂O₄·2H₂O | 2H₂O | 20,03 / 19,69 | 150-250 |
| ZnC₂O₄·2H₂O | 2H₂O / CO + CO₂ | 19,02 / 38,02 | 100-180 / 350-420 → ZnO |
| NiC₂O₄·2H₂O | 2H₂O / CO + CO₂ | 19,72 / 39,41 | 200-280 / ~350 (em N₂ dá Ni, não NiO: confira o resíduo) |
| KHC₈H₄O₄ (biftalato) | - | - | funde/decompõe ~295 |

Polímero → resíduo: PE/PP/PS/PMMA deixam ~0% em N₂; PVC ~10%; PC ~20-25%; PA ~0-2%;
fenólicas e PI 40-60%; celulose 10-20%; lignina 30-45%.

## 4. Padrões de calibração de T e ΔH (DSC/DTA)

Use o valor do certificado do lote; os abaixo são de referência (ITS-90; GEFTA/ICTAC/NIST).

| Padrão | T_fus/trans (°C) | ΔH (J/g) | Uso |
|---|---|---|---|
| Ciclo-hexano | −87,06 (sólido-sólido); 6,54 (fusão) | 79,6; 31,3 | baixa T |
| Mercúrio | −38,83 | 11,4 | |
| Água | 0,00 | 333,6 | |
| Gálio | 29,76 | 80,1 | |
| Bifenila | 68,93 | 120,4 | orgânico |
| Ácido benzoico | 122,35 | 147,3 | orgânico |
| Índio | 156,60 | 28,62 | **verificação diária** |
| Estanho | 231,93 | 60,4 | |
| Bismuto | 271,40 | 53,1 | |
| Chumbo | 327,46 | 23,0 | |
| Zinco | 419,53 | 107,5 | |
| KClO₄ | ~299,5 (sólido-sólido) | ~104 | DTA |
| Ag₂SO₄ | ~426 | ~51 | DTA |
| Quartzo α→β | 573 | pequeno | DTA |
| K₂SO₄ | ~583 (sólido-sólido) | ~33 | DTA |
| K₂CrO₄ | ~668 | ~36 | DTA |
| Alumínio | 660,32 | 397 | alta T (cadinho de alumina) |
| BaCO₃ | ~808 | ~95 | DTA |
| SrCO₃ | ~925 | ~130 | DTA |
| Prata | 961,78 | 104,6 | alta T |
| Ouro | 1064,18 | 63,7 | alta T |

## 5. Pontos de Curie (calibração de T em TG)

Ímã sob o forno; a massa aparente muda bruscamente na T de Curie. Valores aproximados (ver
certificado; materiais ICTAC/NIST GM-761 e similares): Alumel ~154-163 °C; Níquel ~354-358 °C;
Mumetal ~382 °C; Perkalloy ~596 °C; Ferro ~770 °C; Hisat-50 ~750 °C; Trafoperm ~754 °C.

## 6. Faixas típicas em outros materiais

- **Água**: adsorvida/superficial < 120-150 °C; água de cristalização 50-300 °C (em etapas);
  hidroxilas estruturais (desidroxilação) 300-700 °C.
- **Orgânicos em sólidos inorgânicos** (MOFs, híbridos, catalisadores com coque): solvente
  nos poros < 200 °C; ligante/molécula orgânica 300-550 °C; coque em ar 400-700 °C.
- **Fármacos**: fusão aguda (paracetamol forma I ~169 °C; ibuprofeno ~76 °C; ácido
  acetilsalicílico ~135-140 °C, decompõe na fusão); hidratos com endo + perda de massa
  < 150 °C.
- **Cimento/argamassa**: água livre e C-S-H 50-200 °C; etringita ~100-130 °C; gesso
  ~140-160 °C; portlandita (Ca(OH)₂) 400-500 °C; carbonatos 600-800 °C.
- **Biomassa**: ver tga.md §6.
- **TiO₂**: anatase → rutilo exo largo 600-900 °C (sem massa).
- **Metais em ar**: ganho de massa por oxidação; Fe→Fe₂O₃ +43,0%; Cu→CuO +25,2%;
  Al oxida pouco (camada passiva); Mg→MgO +65,8%.
