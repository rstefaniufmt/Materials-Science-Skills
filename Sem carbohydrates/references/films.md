# Filmes de carboidratos no MEV: superfície e fratura

Filmes de amido, quitosana, alginato, pectina, CMC/HPMC, carragena, ágar, gelatina+polissacarídeo,
blendas (amido/PVA, amido/quitosana, amido/PLA, amido/PBAT), compósitos (nanocelulose, argila,
ZnO, TiO₂, Ag) e filmes ativos/emulsionados (óleos essenciais, lipídios, ceras).

## Sumário
1. Superfície: o que cada característica indica
2. Fratura (corte transversal)
3. Medidas que a skill faz
4. Blendas, compósitos e emulsões
5. Lado do ar × lado da placa
6. Como escrever

## 1. Superfície: o que cada característica indica

| Observação | Interpretação provável | Confirme com |
|---|---|---|
| Superfície lisa, contínua, homogênea | boa miscibilidade/dispersão, gelatinização completa | DSC (sem endoterma de gelatinização residual), DRX |
| **Grânulos ou "fantasmas"** (protuberâncias arredondadas do tamanho dos grânulos) | gelatinização incompleta (temperatura/tempo/água insuficientes) | microscopia de luz polarizada (cruz de malta), DSC |
| Pontos/manchas claras pequenas, cristais facetados | **exsudação de plastificante** (glicerol, sorbitol), sais (CaCl₂, NaCl), cristais de aditivo | armazenamento em UR alta acelera; EDS para sais |
| Poros/furos circulares | bolhas de ar (sem desgaseificação), evaporação rápida do solvente, gotas de óleo removidas | desgaseificar a solução; secar mais devagar |
| **Trincas** | filme frágil (pouco plastificante, retrogradação), tensão de secagem, **dano do feixe/vácuo** (trincas que crescem durante a observação) | propriedades mecânicas (alongamento); imagem a kV menor |
| Rugosidade/textura granular fina | separação de fases, cristalização (retrogradação), aditivos/cargas | AFM (rugosidade real), DRX |
| Domínios/gotículas (ilhas) | **separação de fases** em blendas, gotas de emulsão | medir tamanho de domínio (partículas/poros), FTIR, DSC (duas Tg) |
| Aglomerados claros (µm) | carga mal dispersa (argila, nanocelulose, óxidos) | EDS (para inorgânicos), DRX |
| Riscos lineares paralelos | marcas da placa/suporte ou de manuseio - não é propriedade do material | lado do ar × lado da placa |
| Regiões brilhantes que mudam durante a varredura | carregamento (filme isolante) | recobrimento melhor, kV menor, baixo vácuo |

## 2. Fratura (corte transversal)

Fratura criogênica: imersão em N₂ líquido (~1-3 min) e quebra; corte com lâmina deforma e esconde
a estrutura interna.

| Observação | Interpretação |
|---|---|
| Fratura lisa, plana, "vítrea" | matriz homogênea, frágil (abaixo da Tg); boa compatibilidade |
| Fratura rugosa, com degraus/escamas | heterogeneidade, deformação plástica (filme mais dúctil/plastificado), cargas |
| Estrutura estratificada/lamelar | secagem em camadas, multicamadas, orientação de cargas lamelares (argila) |
| Vazios esféricos | bolhas ou gotas de óleo (filmes emulsionados) - tamanho de gota indica estabilidade da emulsão |
| Partículas arrancadas (*pull-out*) e cavidades ao redor de cargas | adesão interfacial ruim carga/matriz |
| Fibras quebradas rente à matriz | boa adesão (fibras/celulose) |
| Interface visível entre camadas | multicamadas, revestimento; adesão ruim se houver vão |
| Faixa com espessura muito irregular | casting irregular, filme encolhido; medir em várias posições |

## 3. Medidas que a skill faz

- **`film`** (superfície): conta e mede poros escuros, partículas/protuberâncias claras e trincas;
  densidade por mm², fração de área, ECD médio; textura (CV de brilho, GLCM: contraste,
  homogeneidade, energia, correlação, entropia; comprimento de correlação; dimensão fractal pela
  PSD). Use para **comparar formulações nas mesmas condições de aquisição** (mesmo kV, detector,
  aumento, contraste e brilho, mesmo recobrimento). Não é rugosidade física.
- **`section`** (corte): espessura média ± DP ao longo da faixa (corrigida pela inclinação),
  mín/máx/p5-p95, porosidade aparente (fração de área de poros na seção) e tamanho de poros.
- Para domínios de separação de fases ou gotas de emulsão claros/escuros sobre a matriz, `granules`
  (com `--dark-objects` se os domínios forem escuros) dá distribuição de tamanho e forma dos domínios.

Parâmetros: `--pore-max` (maior poro esperado), `--particle-max` (maior partícula), `--crack-min`
(menor trinca), `-k` (sensibilidade: maior = só os mais contrastados). Comece com o padrão, olhe o
overlay e ajuste uma coisa por vez.

## 4. Blendas, compósitos e emulsões

- **Blendas**: fase contínua + fase dispersa em gotas/ilhas = imiscível; tamanho de domínio
  diminui com compatibilizante. Superfície homogênea não prova miscibilidade em escala molecular
  (confirme com DSC/DMA: uma Tg).
- **Nanocompósitos**: a dispersão "boa" no MEV de baixo aumento só mostra ausência de aglomerados
  grandes (> ~1 µm). Dispersão nanométrica → MET/AFM/DRX (argila esfoliada).
- **Filmes emulsionados**: gotas de óleo (vazios na fratura) e sua distribuição; coalescência e
  cremeação aparecem como gotas grandes e concentração perto de uma face.
- **Filmes ativos com óleos essenciais**: canais/poros na fratura por volatilização do óleo; efeito
  sobre a permeabilidade ao vapor d'água.

## 5. Lado do ar × lado da placa

Em filmes por casting, a face em contato com a placa (lado de baixo) é mais lisa e replica a
superfície da placa; a face do ar mostra exsudação, cristais e partículas que migram/sedimentam.
Sempre indique qual face foi observada e compare faces iguais.

## 6. Como escrever

"A superfície do filme controle (amido de mandioca, 30% glicerol) apresentou-se contínua e
homogênea, sem grânulos residuais, indicando gelatinização completa. A adição de 5% de óleo
essencial gerou poros superficiais (712 poros/mm², ECD = 2,6 ± 1,1 µm) e, na fratura criogênica,
vazios esféricos associados às gotas de óleo (porosidade aparente de 0,7%). A espessura no corte
foi de 86,5 ± 3,0 µm (MEV), concordante com o micrômetro (88 ± 4 µm, n = 10)."
