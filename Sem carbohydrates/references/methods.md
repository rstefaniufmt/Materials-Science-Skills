# Métodos: definições, algoritmos e limites

## Sumário
1. Segmentação de grânulos
2. Medidas por objeto
3. Classes de forma
4. Estatística de tamanho (número × volume)
5. Bimodalidade
6. Pites/poros superficiais
7. Defeitos e textura de filmes
8. Espessura em corte
9. Comparação entre amostras
10. Validação desta skill

## 1. Segmentação de grânulos (`segment`)

1. Suavização gaussiana (σ = 1 px; `--smooth`), correção opcional de iluminação (`--flatten`:
   divide por uma versão muito suavizada da imagem).
2. Limiar global (Otsu padrão; Li, triângulo, Yen, valor fixo) ou local (`--thresh local`, bloco
   gaussiano) - objetos claros sobre fundo escuro (SE típico sobre fita de carbono); `--dark-objects`
   inverte.
3. Abertura morfológica (1 px) e preenchimento de buracos (o centro escuro dos grânulos em SE e os
   pites ficam dentro do objeto). Buracos cercados por **dois ou mais** grânulos (vazio entre grânulos
   encostados) voltam a ser fundo depois da separação.
4. **Separação de grânulos encostados**: watershed na transformada de distância, com marcadores =
   máximos regionais cuja proeminência é ≥ `split_h` × raio do objeto (h-máxima adaptativo). Menor
   `--split-h` separa mais (risco de partir grânulos alongados); maior separa menos.
5. Remove objetos com ECD < `--min-d` (padrão 6 px) ou > `--max-d`.
6. Objetos que tocam a borda da imagem - e fragmentos do watershed com contato longo com eles - ficam
   marcados `borda` e saem da estatística.

Limites: grânulos sobrepostos em várias camadas, fundo com o mesmo brilho dos grânulos, grânulos
transparentes ao feixe (bordas claras e interior igual ao fundo em baixo kV) e imagens com
carregamento forte. Nesses casos ajuste o preparo/aquisição; o software não recupera informação.

## 2. Medidas por objeto (`measure`) - ISO 9276-6 / ISO 13322-1

| Coluna | Definição |
|---|---|
| `area` | área projetada (µm²) |
| `ecd` | diâmetro do círculo de mesma área: √(4A/π) |
| `feret_max`, `feret_min` | maior e menor distância entre retas paralelas tangentes (calibre; Feret mín por pinças rotativas na envoltória convexa) |
| `eixo_maior`, `eixo_menor` | eixos da elipse de mesmos momentos de 2ª ordem |
| `perimetro` | perímetro de Crofton (menos sensível à pixelização que contar pixels) |
| `AR` | razão de aspecto = Feret máx / Feret mín (1 = isométrico) |
| `circularidade` | 4πA/P² (1 = círculo; hexágono 0,91, quadrado 0,79; cai com contorno rugoso e com a resolução) |
| `arredondamento` | 4A/(π Feret_máx²) (insensível à rugosidade do contorno) |
| `solidez` | A / A da envoltória convexa (reentrâncias, fragmentos, agregados reduzem) |
| `convexidade` | perímetro da envoltória convexa / perímetro |
| `vertices` | vértices do polígono que aproxima a envoltória (tolerância 4% do perímetro) - poliédrico × arredondado |
| `brilho` | intensidade média (comparativa) |

Resolução: abaixo de ~10 px de diâmetro, perímetro e circularidade perdem significado; abaixo de
~5 px, até o ECD tem erro > 10%.

## 3. Classes de forma (`classify`, descritivas)

Em ordem: solidez < 0,80 → **agregado/sobreposto**; < 0,90 → **irregular/fragmentado**; AR ≥ 2 →
**alongado/lenticular (de perfil)**; 3-7 vértices, circularidade < 0,92 e AR < 1,5 →
**poliédrico/angular**; AR < 1,2 → **esférico/arredondado**; demais → **oval/elíptico**. São
categorias para descrever e comparar amostras, não taxonomia; reporte também AR e circularidade.

## 4. Estatística de tamanho (`size_stats`)

- **Em número**: média, DP, CV, mediana, D10/D50/D90 (percentis do ECD), span = (D90 − D10)/D50.
- **Em volume** (cada grânulo pesa ECD³ - esfera equivalente): Dv10/Dv50/Dv90, D[4,3] = Σd⁴/Σd³
  (média em volume, a "média" da difração laser), D[3,2] = Σd³/Σd² (Sauter, relacionado à área
  específica).
- A distribuição em volume é dominada pelos poucos grânulos grandes: com n pequeno ela é instável.
- Comparação com difração laser: MEV dá distribuição **em número de área projetada**; laser dá
  volume de esfera equivalente por espalhamento - Dv50 do MEV e D50 do laser costumam concordar só
  aproximadamente (forma, orientação, aglomerados no laser).

## 5. Bimodalidade (`bimodal`)

Mistura gaussiana de 1 e 2 componentes em ln(ECD); bimodal se ΔBIC > 10, separação de Ashman D > 2
e cada população com ≥ 3%. Corte = ponto onde as duas densidades ponderadas se igualam. Exige ≥ 40
objetos (prefira centenas). Mistura de amidos de fontes diferentes também aparece como bimodal.

## 6. Pites/poros superficiais (`pits`)

Top-hat preto (fechamento em tons de cinza − imagem) com disco de diâmetro `--pit-max`, dentro dos
grânulos erodidos (longe da borda); limiar = mediana + 5 × desvio robusto; mantém só manchas
compactas (AR ≤ 2,5, solidez ≥ 0,75), menores que 35% do ECD do grânulo hospedeiro e com **degrau
nítido** (anel ao redor ≥ 0,35 × contraste grânulo-fundo mais claro que o interior). Saída:
pites/grânulo, % de grânulos com pites, diâmetro médio, % de área. É contagem da face visível
(metade do grânulo): use para comparar com o nativo, não como porosidade absoluta.

## 7. Defeitos e textura de filmes (`surface_defects`, `surface_texture`)

- Partículas claras: top-hat branco; anéis claros com interior escuro (borda de poro) são descartados.
- Trincas: filtro de cristas escuras de Sato (σ = 1-3 px), limiar robusto, esqueleto; mantém
  estruturas com comprimento ≥ `--crack-min` e esbeltez ≥ 8; ignora o halo das partículas.
- Poros: top-hat preto (diâmetro ≤ `--pore-max`), fora das trincas e dos halos; compactos.
- Densidades por mm² usam a área total analisada.
- Textura: CV de brilho; GLCM (64 níveis, 4 direções, distâncias em µm); comprimento de correlação
  (autocorrelação cai a 1/e); inclinação β da densidade espectral radial e D = (8 − β)/2 quando
  2 < β < 4. São índices comparativos (mesma aquisição).

## 8. Espessura em corte (`cross_section`)

Faixa do filme = maior região contínua mais clara/texturizada (brilho e desvio local combinados,
Otsu) que atravessa a imagem; bordas por coluna refinadas no máximo do gradiente, suavizadas por
mediana; espessura perpendicular = vertical × cos(inclinação local). Colunas em que a faixa toca a
borda da imagem são descartadas. Exatidão ~±2 px por borda (o brilho de borda em SE desloca a
borda para fora em ~1-2 px): prefira aumento em que a espessura tenha ≥ 200 px. Porosidade: manchas
escuras compactas dentro da faixa (longe das bordas), limiar robusto (6 × MAD).

## 9. Comparação entre amostras (`compare`)

Kruskal-Wallis (3+ amostras) e Mann-Whitney par a par com correção de Holm, efeito por correlação
rank-biserial (|r| ≈ 0,1 pequeno, 0,3 médio, 0,5 grande). Com milhares de grânulos, diferenças
mínimas ficam "significativas": olhe o efeito e a diferença de medianas. Com réplicas de preparo,
compare as médias por réplica (ANOVA/Kruskal com n = réplicas).

## 10. Validação desta skill

Imagens sintéticas tipo SE com verdade conhecida: grânulos poliédricos (milho), mistura A+B
(trigo), amido poroso, filme com poros/partículas/trincas e corte transversal. ECD dos grânulos
detectados: erro médio 0,2-1% (IoU > 0,8 em 95-100% dos grânulos inteiros); populações A/B
recuperadas (médias 18,0 e 4,4 µm × verdade 18,6 e 4,4 µm); pites 6,7/grânulo × 9 verdadeiros
(pites encostados se fundem; 0 e 1 falso positivo em 21 e 88 grânulos sem pites); espessura 86,5 µm ×
85,4 µm verdadeiros. Metadados conferidos em TIFs reais de FEI/Thermo e Zeiss; barra + OCR conferidos
em databars FEI e Zeiss reais (erro de escala ≤ 1,2%).
