# Interpretação de espectros FT-IR

Roteiro para transformar uma lista de bandas numa interpretação defensável. A ideia
central: cada conclusão deve se apoiar num **conjunto coerente** de bandas, não numa
banda isolada, e o texto deve deixar claro o grau de certeza.

## Sumário
1. Antes de interpretar
2. Roteiro por perguntas (ordem recomendada)
3. Combinações que confirmam grupos
4. Região de impressão digital (< 1500 cm⁻¹)
5. Amostras reais: misturas, materiais, comparação de séries
6. Armadilhas (artefatos, ATR, KBr)
7. Como redigir

## 1. Antes de interpretar

- **Contexto**: o que o usuário espera da amostra (composto sintetizado, polímero,
  compósito, biomaterial, mineral)? Interpretação "às cegas" é muito mais fraca que
  verificar uma hipótese. Se ele disse o que é, a tarefa vira "o espectro é consistente
  com X? o que sobra/falta?".
- **Técnica de amostragem**: ATR (cristal?), pastilha de KBr, filme, solução, DRIFTS.
  Muda intensidades relativas, posições e artefatos esperados.
- **Qualidade**: saturação (bandas achatadas no fundo, T ≈ 0, A > 2.5), ruído, linha de
  base inclinada/ondulada (franjas de interferência em filmes finos), CO₂/vapor.
  Bandas saturadas têm posição e forma pouco confiáveis - diga isso.
- **Rode `screen_groups`** (vem no `report` e em `peaks --assign`) para ter o checklist
  das perguntas abaixo já preenchido; depois confira cada item no espectro.

## 2. Roteiro por perguntas (ordem recomendada)

1. **Há C=O (1870-1650, banda forte)?** É a pergunta mais informativa. Se sim, a posição
   restringe muito: ~1800+ (anidrido, haleto de acila, carbonato cíclico) ›
   ~1740 (éster) › ~1725-1705 (aldeído, cetona, ácido) › ~1690-1630 (amida, ureia,
   cetona conjugada). Conjugação e ligação de H baixam ~20-40 cm⁻¹; anel tenso sobe.
2. **Se há C=O, qual o tipo?** Procure as parceiras (seção 3): O–H muito largo (ácido),
   C–O forte 1300-1000 (éster/anidrido), C–H 2820/2720 (aldeído), N–H + amida II (amida).
3. **Há O–H ou N–H (3650-3100)?** Larga e arredondada → O–H em ligação de H (álcool,
   fenol, água, polissacarídeo). Mais estreita, 1 ou 2 bandas → N–H (2 = NH₂). Muito
   larga 3300-2500 com "dentes" do C–H por cima → ácido carboxílico.
4. **C–H: acima ou abaixo de 3000?** Acima (3100-3000) = sp² (aromático/alceno); abaixo
   (2990-2840) = sp³. Razão das bandas CH₃ (2960/2870) e CH₂ (2925/2855) dá ideia de
   cadeias longas vs ramificadas. ~3300 fina = ≡C–H.
5. **Triplas/cumuladas (2300-2000)?** Região quase vazia; qualquer banda real aqui é muito
   informativa (C≡N ~2250, C≡C, N=C=O ~2270, N₃ ~2100, SCN⁻, Si–H, M–CO). Antes, descarte
   CO₂ (2360/2340) e o ruído do diamante (ATR).
6. **Aromático?** Precisa de ~1600 e ~1500 (às vezes 1580 e 1450) + C–H > 3000 + oop
   forte 900-680 (padrão de substituição) + sobretons fracos 2000-1650.
7. **Bandas fortes 1300-1000?** C–O (éster, éter, álcool), Si–O, P–O, S=O, C–F. Quase
   sempre as mais fortes do espectro quando presentes - a posição exata e a forma ajudam.
8. **Abaixo de 1000**: oop de C–H (aromático/alceno), C–Cl/Br, e para inorgânicos as
   bandas de M–O, carbonato (875/712), fosfato (600/560), sulfato (670/600), silicato (460).

## 3. Combinações que confirmam grupos

| Grupo | Bandas que precisam aparecer juntas |
|---|---|
| Ácido carboxílico | C=O ~1710 (dímero) + O–H 3300-2500 muito largo + C–O ~1300-1210 + O–H oop ~930 largo |
| Carboxilato (sal) | sem C=O ~1710; COO⁻ ass 1610-1550 + sim 1440-1335 (Δ indica coordenação) |
| Éster | C=O ~1740 + C–O–C 1300-1150 (forte) + O–C–C 1150-1000 ("regra dos três") |
| Cetona | C=O ~1715 sem O–H largo, sem C–O forte, sem 2720 |
| Aldeído | C=O 1740-1685 + C–H 2820 e **2720** (dubleto de Fermi) |
| Anidrido | duas C=O (~1820 e ~1760, Δ≈60) + C–O 1300-1000 |
| Amida 2ª / proteína | N–H ~3300 + amida I ~1650 + amida II ~1540 (+ amida III 1300-1230) |
| Amida 1ª | N–H duplo 3350/3180 + amida I ~1670 + δNH₂ ~1630 |
| Amina 1ª | N–H duplo 3500-3300 + δNH₂ 1650-1580 + C–N |
| Uretano (PU) | N–H ~3330 + C=O 1730/1700 + amida II ~1530 + C–O–C ~1220/1070 |
| Álcool | O–H largo 3550-3200 + C–O 1260-1000 (1° ~1050, 2° ~1100, 3° ~1150) |
| Éter | C–O–C forte 1150-1085 (alifático) ou 1275-1200 + 1075-1020 (aril), sem O–H e sem C=O |
| Nitro | 1570-1500 + 1390-1330, ambas fortes |
| Nitrila | ~2250 fina; aromática ~2230 |
| Epóxido (terminal) | ~915 (e ~830, ~1250); some na cura |
| Aromático monossubst. | 3100-3000 + 1600/1500 + **~750 e ~700 fortes** |
| Aromático para | 3100-3000 + 1600/1510 + **~830** forte |
| Água | O–H 3600-3000 largo + δ H–O–H ~1640 |
| Carbonato | 1500-1400 muito forte e larga + ~875 fina (+ ~712) |
| Fosfato (apatita) | 1100-1000 muito forte + dubleto 603/565 (+ ~960) |
| Sílica/silicato | 1100-1000 muito forte + ~800 + ~460 |

Se só parte do conjunto aparece, diga "compatível com", não "confirma".

## 4. Região de impressão digital (< 1500 cm⁻¹)

Bandas aqui são acopladas e dependem da molécula inteira. Use-a para:
- **comparar** com um padrão/biblioteca medido na mesma técnica (é aqui que dois
  compostos parecidos se distinguem);
- confirmar grupos já sugeridos acima (C–O, oop de aromáticos/alcenos);
- **não** para atribuir cada banda individualmente a um grupo sem suporte.

## 5. Amostras reais

- **Misturas/compósitos**: o espectro é ≈ soma (em absorbância). Compare com os
  componentes puros na mesma técnica; `subtract(a, b, region=...)` remove a matriz
  ajustando o fator numa banda exclusiva dela. Bandas novas ou deslocadas indicam
  interação (ex.: C=O deslocando para menos = ligação de H; COO⁻ aparecendo = sal/complexo).
- **Deslocamentos pequenos (2-10 cm⁻¹)** só são significativos se a resolução, a técnica e
  o processamento forem idênticos entre amostras, e se maiores que o intervalo de dados.
- **Séries (cura, degradação, dopagem, concentração)**: use `compare` com janelas fixas e
  normalização por banda interna estável. Índices clássicos:
  - carbonila (oxidação de poliolefinas): A(1760-1680)/A(1490-1440) (ou ref. 1375/2020);
  - conversão de epóxi: A(915)/A(1510 aromático);
  - grau de desacetilação da quitosana: A(1655)/A(3450) ou A(1320)/A(1420) (método-dependente);
  - cristalinidade do PE: dubleto 730/719; do PET: 1340/1410, 973;
  - razão amida I/amida II para conformação/hidratação de proteínas.
  Diga sempre qual definição foi usada.
- **Proteínas (estrutura secundária)**: amida I deconvoluída (2ª derivada + ajuste de
  bandas). As sub-bandas da tabela (`components=True`) são faixas típicas, com
  sobreposição - trate como estimativa.

## 6. Armadilhas

- **CO₂** 2360/2340 (+ 667): diferença entre fundo e amostra. Positivo ou negativo.
- **Vapor d'água**: dezenas de linhas finas 3900-3500 e 1900-1300. Distorcem amida I e
  qualquer banda fraca nessas regiões. Não atribua linhas finas aí.
- **KBr**: absorve umidade (3440 e 1640 largos); KNO₃ dá 1384 fina; pastilha grossa ou mal
  moída dá efeito Christiansen (bandas assimétricas, linha de base inclinada); troca iônica
  com sais da amostra.
- **ATR**: profundidade ∝ λ (bandas < 1500 relativamente mais fortes que em transmissão);
  dispersão anômala desloca bandas fortes alguns cm⁻¹ para menos; contato ruim = espectro
  fraco e ruidoso; diamante tem fônons 2300-1900 (região ruidosa) e corta abaixo de ~525
  em alguns modelos. Contaminação de limpeza (isopropanol: 950, 817; acetona: 1715) se o
  cristal não secou.
- **Contaminantes comuns**: graxa de silicone (1260, 1100-1020, 800), ftalatos de
  plastificante (1725, 1280, 1120, 1070, 745), óleo/digital (2920, 2850, 1460), talco
  (3675, 1015) e CaCO₃ (1420, 875) como cargas em plásticos.
- **Franjas de interferência** (filmes finos/lisos): ondulação senoidal da linha de base;
  não são bandas.
- **Saturação**: bandas com fundo plano em %T ≈ 0; posição incerta, intensidades relativas
  inválidas. Sugira amostra mais fina/diluída ou ATR.

## 7. Como redigir

- Comece pelo que é inequívoco (grupos com conjunto completo), depois o provável, depois
  o que é ambíguo e o que **não** está presente (ausências também informam: "ausência de
  banda em 1700-1750 descarta carbonila em quantidade detectável").
- Tabela: posição (cm⁻¹, inteiro) | intensidade (vs/s/m/w, br/sh) | atribuição (modo
  vibracional: ν estiramento, δ deformação, ρ rocking, ω wagging, oop fora do plano;
  "as"/"s" assimétrico/simétrico) | comentário.
- Cite a técnica, o acessório, a resolução e qualquer processamento (linha de base,
  normalização, CO₂ removido, correção ATR).
- Separe o que vem do espectro do que vem do conhecimento prévio sobre a amostra.
- Para conclusões fortes (identidade, pureza, fase), recomende confirmação: padrão de
  referência medido no mesmo equipamento, biblioteca, DRX, RMN, TGA etc.
