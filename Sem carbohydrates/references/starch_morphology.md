# Morfologia de grânulos de amido e outras partículas de carboidrato

Faixas típicas de literatura (MEV e microscopia óptica; ECD ou eixo maior). Variam com cultivar,
maturidade, extração e medida: use como **plausibilidade**, não como identificação. Fontes de base:
Jane et al., *Starch/Stärke* 46 (1994) 121 (atlas de MEV de 54 amidos); Lindeboom, Chang & Tyler,
*Starch/Stärke* 56 (2004) 89 (tamanho de grânulos); Fannon, Hauber & BeMiller, *Cereal Chem.* 69
(1992) 284 (poros superficiais); Tester, Karkalas & Qi, *J. Cereal Sci.* 39 (2004) 151.

## Sumário
1. Amidos nativos por fonte
2. Grânulos A e B (trigo, cevada, centeio, triticale)
3. Superfície: poros, canais, hilo, fissuras
4. Efeito de tratamentos e modificações
5. Outras partículas de carboidrato
6. Como descrever a morfologia no texto

## 1. Amidos nativos por fonte

| Fonte | Tamanho típico (µm) | Forma | Distribuição / observações | Cristalinidade |
|---|---|---|---|---|
| Milho normal / ceroso | 5-25 (média ~10-15) | poliédrico/angular (endosperma vítreo) e arredondado (farináceo) | unimodal; **poros superficiais naturais** e canais | A |
| Milho alta amilose | 3-20 | irregular, alongado, "salsicha"/filamentoso, alguns agregados | heterogêneo | B |
| Sorgo, milheto | 5-25 | poliédrico, como milho | poros superficiais | A |
| Arroz | 2-8 | poliédrico, angular | muitos em **grânulos compostos** (agregados); difícil individualizar | A |
| Aveia | 2-15 (individuais); compostos 20-150 | poliédrico irregular | grânulos compostos que se desfazem na extração | A |
| Quinoa, amaranto | 0,5-3 | poliédrico muito pequeno | agregados esféricos de muitos grânulos; exige FEG/aumento alto | A |
| Trigo, cevada, centeio, triticale | A: 10-35; B: 2-10 | A lenticular/disco (com sulco equatorial); B esférico/poliédrico | **bimodal** (ver §2) | A |
| Batata | 5-100 (média ~30-40) | oval/elipsoide (grandes), esférico (pequenos); superfície lisa | assimétrica, cauda para grandes | B |
| Mandioca (fécula/polvilho) | 4-35 (média ~12-15) | arredondado/oval com **face truncada** (forma de sino/"chapéu") | unimodal | A ou Cₐ |
| Batata-doce | 3-30 | arredondado, poliédrico, truncado | | C/A |
| Inhame (*Dioscorea*) | 10-70 | oval, triangular arredondado, achatado | | B/C |
| Taro/cará (*Colocasia*) | 1-5 | poliédrico muito pequeno | agregados | A |
| Araruta (*Maranta*) | 5-50 | oval a elíptico, alguns truncados | | A/C |
| Banana verde | 10-60 | alongado, elipsoide irregular, às vezes com estrias/cristas | grande, irregular | B/C |
| Ervilha, feijão, lentilha, grão-de-bico (leguminosas) | 5-45 | oval, reniforme ("rim"); às vezes **fissura central** (hilo) | unimodal ou levemente bimodal | C |
| Pinhão (*Araucaria angustifolia*) | ≈ 5-25 | arredondado, oval, alguns truncados | | A/C |
| Mandioquinha/batata-baroa, gengibre, outras raízes | variável | variável | consulte literatura do material | |

Leituras práticas:
- Superfície lisa e bordas nítidas = grânulo íntegro. Fissuras e reentrâncias podem ser naturais
  (leguminosas, hilo) ou de dano mecânico/secagem - compare com o nativo.
- Faces planas truncadas (mandioca, batata-doce) e poliedros (milho, arroz) vêm do empacotamento
  dos grânulos no amiloplasto; não são dano.
- Grânulos compostos (arroz, aveia, taro, quinoa) aparecem como agregados de baixa solidez (`classe =
  agregado/sobreposto`): medir os grânulos individuais exige aumento maior e é limitado em MEV.

## 2. Grânulos A e B (trigo, cevada, centeio, triticale)

- **A**: lenticulares (disco), 10-35 µm, ~70-80% do volume e < 10% do número;
  **B**: esféricos/poliédricos 2-10 µm, ~90% do número. Alguns autores separam **C** < 2-3 µm
  (cevada).
- `bimodal()` testa 1 vs 2 populações em log(ECD) (mistura gaussiana, ΔBIC > 10) e dá o corte, a
  fração em número e em volume de cada população. O corte típico fica entre 8 e 12 µm.
- Lenticulares deitados mostram a face (ECD grande, AR ~1,1-1,3); de perfil parecem alongados (AR
  > 2, classe "alongado/lenticular") - a mistura de orientações alarga a distribuição de A.
- A fração de B influencia a viscosidade, a gelatinização e a qualidade de panificação/malte; é a
  informação mais útil desta análise para cereais de inverno.

## 3. Superfície: poros, canais, hilo, fissuras

- **Poros naturais** (0,1-0,3 µm, abertura de canais que vão ao hilo) em milho, sorgo, milheto
  (e em trigo/cevada sobretudo na região equatorial); ausentes ou raros em batata, mandioca e
  leguminosas. São visíveis só em FEG/aumento alto (≥ 5-10 k×) - não espere detectá-los a 1000×.
- **Hilo**: ponto central (às vezes fendido em leguminosas e em grânulos secos rapidamente).
- **Fissuras** radiais ou "estrela" no hilo: secagem rápida, congelamento, ultrassom, dano do feixe.
- Detecção automática (`--pits`): pites escuros compactos e nítidos. O diâmetro mínimo detectável
  é ~5 px - a 1000× (≈ 0,1 µm/px) isso é ~0,5 µm. Para poros naturais finos, use aumento alto.

## 4. Efeito de tratamentos e modificações (o que procurar)

| Tratamento | Aparência típica no MEV |
|---|---|
| Gelatinização completa | perda da estrutura granular; matriz contínua, esponjosa ou em placas |
| Gelatinização parcial / inchamento | grânulos inchados, deformados, "fantasmas" (envoltórios colapsados), fusão entre grânulos |
| Retrogradação / amido resistente tipo 3 | massas compactas, superfície rugosa, sem grânulos |
| Hidrólise ácida (lintnerização, HCl/H₂SO₄) | erosão e rugosidade superficial, fissuras; hilo exposto; fragmentação nos casos severos; forma geral preservada em hidrólises brandas |
| Hidrólise enzimática (α-amilase, amiloglucosidase, pululanase) | **pites e poros** (pinholes), canais até o hilo (endocorrosão), grânulos "esponja"/ocos; fragmentos. Milho, trigo e mandioca são atacados por poros; batata (B) resiste e mostra erosão superficial |
| Amido poroso (produção) | muitos poros 0,5-2 µm distribuídos na superfície; reporte pites/grânulo e % de grânulos com pites |
| HMT (calor-umidade) | em geral pouca mudança; às vezes depressões/cavidades perto do hilo, superfície levemente rugosa, agregação |
| Annealing | praticamente sem mudança visível |
| Ligação cruzada, fosfatação, hidroxipropilação, acetilação (GS baixo) | pouca ou nenhuma mudança visível - **ausência de mudança é o resultado esperado** |
| Acetilação/esterificação com GS alto, succinilação (OSA) intensa | superfície rugosa, fusão e agregação de grânulos, perda de definição |
| Oxidação (hipoclorito, H₂O₂, ozônio) | superfície levemente rugosa, fissuras, poros discretos |
| Ultrassom | fissuras, sulcos, pequenas depressões e fragmentos |
| Moagem/ball milling | fragmentação, grânulos quebrados com faces de fratura, "amido danificado", agregados |
| Extrusão, pré-gelatinização em tambor | estrutura granular destruída, placas ou partículas porosas irregulares |
| Plasma frio, irradiação (γ, UV) | mudanças sutis de superfície; muitas vezes invisíveis |
| Nanopartículas/nanocristais de amido | agregados de partículas 20-200 nm; exige FEG, recobrimento fino; não usar o detector de grânulos sem rever escala |
| Complexos amilose-lipídio, inclusão | partículas irregulares, agregados |

Comparação quantitativa útil: ECD/D50 (fragmentação reduz; inchamento aumenta), solidez (cai com
erosão, fusão e agregação), fração de grânulos com pites e pites/grânulo (hidrólise enzimática),
% de "irregular/fragmentado" e "agregado" (moagem, modificação química intensa).

## 5. Outras partículas de carboidrato

| Material | Aparência típica |
|---|---|
| Microcápsulas por **spray-drying** (maltodextrina, goma arábica, amido OSA, inulina, proteína+polissacarídeo) | esferas 1-50 µm, com **depressões/rugas** (encolhimento na secagem) ou lisas; ocas; trincas e poros = proteção ruim do núcleo; aglomerados = higroscopicidade/umidade |
| **Liofilização** | partículas irregulares, flocos/placas porosas, aspecto de vidro quebrado |
| Quitosana em pó / micropartículas (gelificação iônica com TPP) | flocos irregulares (pó comercial); esferas/aglomerados (micro/nanopartículas) |
| Esferas de alginato (gotejamento em CaCl₂) | esferas de 0,5-3 mm quando úmidas; secas: superfície rugosa, enrugada, colapsada |
| Celulose microcristalina (MCC) | partículas irregulares/alongadas 20-200 µm com fibrilas |
| Nanocelulose (CNF/CNC), quitina nanofibrilar | rede de fibrilas (CNF), bastonetes (CNC) - exige FEG; secagem cria filmes/agregados |
| Fibras vegetais (bagaço, casca) | feixes de fibras, células, poros, sílica; medir comprimento/diâmetro em vez de ECD |
| Açúcares cristalinos, sais | cristais facetados com arestas retas |

Para fibras e bastonetes, use Feret máx/mín (comprimento/largura) e AR em vez do ECD.

## 6. Como descrever a morfologia no texto

Uma frase de forma + uma de tamanho + uma de superfície, sempre com números:
"Os grânulos de amido de mandioca nativo apresentaram forma arredondada a oval com faces
truncadas (esférico/arredondado 58%, oval 31%), D50 = 13,4 µm (D10-D90 = 7,9-19,6 µm; n = 486) e
superfície lisa, sem poros visíveis a 1500×. Após hidrólise com α-amilase, 74% dos grânulos
exibiram pites (3,1 pites/grânulo; d = 0,9 ± 0,3 µm), sem alteração significativa do D50
(Mann-Whitney, p = 0,21)."
