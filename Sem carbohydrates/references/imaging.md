# Preparo, aquisição, artefatos e como reportar

## Sumário
1. Preparo de pós e grânulos
2. Preparo de filmes
3. Condições de aquisição
4. Aumento certo e quantas imagens
5. Artefatos e como reconhecê-los
6. Como reportar (métodos e legendas)

## 1. Preparo de pós e grânulos

- Fita de carbono dupla face no stub; espalhe pouco pó e **remova o excesso** (bater o stub de lado
  ou jato de ar seco). Camada única, grânulos separados: é o que torna a segmentação possível.
- Alternativa para pós que aglomeram: suspensão diluída em etanol absoluto ou acetona (amido não
  incha), gota no stub/lamínula e secagem; ultrassom curto para dispersar (cuidado: ultrassom
  longo fissura grânulos).
- Pós úmidos (amido modificado, hidrolisado) devem ser secos antes (estufa branda ou liofilização)
  - água sob vácuo e feixe gera bolhas e trincas.
- Recobrimento: ouro ou ouro-paládio 10-20 nm (sputter, ~60-120 s); carbono se for fazer EDS.
  Recobrimento grosso esconde poros finos; fino demais carrega.
- **Baixo vácuo/ESEM** (≈ 30-130 Pa, detector de baixo vácuo/BSE) permite amostras sem recobrimento;
  resolução e contraste de superfície menores.

## 2. Preparo de filmes

- Condicione (UR e T conhecidas) e seque sob vácuo/dessecador antes; anote a face (ar/placa).
- **Fratura criogênica**: tira do filme mergulhada em N₂ líquido 1-3 min e quebrada com pinças;
  monte de pé (suporte de seção transversal ou fita na lateral do stub) com a fratura para cima.
- Superfície: pedaço plano colado na fita, sem dobrar; evite tocar a face a observar.
- Recobrir dos dois tipos de amostra; a fratura precisa de recobrimento também nas laterais
  (incline o stub no sputter ou faça duas deposições).

## 3. Condições de aquisição

| Parâmetro | Recomendação para carboidratos | Por quê |
|---|---|---|
| Tensão (kV) | **1-5 kV** (FEG) ou 5-10 kV (tungstênio) | menos dano e carregamento; mais detalhe de superfície. 15-20 kV "atravessa" a superfície e queima a amostra |
| Detector | SE (ETD/Everhart-Thornley, SE2, InLens/TLD) para topografia | BSE para contraste de composição (cargas inorgânicas) ou em baixo vácuo |
| Corrente/spot | baixa | dano do feixe |
| WD | 5-10 mm (resolução) ou maior (profundidade de foco para pós) | |
| Varredura | lenta ou com média de quadros (integração) | reduz ruído; mas tempo longo no mesmo campo danifica - foque em área vizinha e mova |
| Inclinação | 0° para medidas de tamanho; inclinado só para ilustrar 3D | inclinação encurta dimensões em uma direção |
| Contraste/brilho | sem saturar (histograma sem picos em 0 e 255) | saturação destrói bordas e textura |
| Resolução de gravação | ≥ 1024 × 768; prefira 2048 × 1536 para estatística | mais grânulos por imagem com resolução suficiente |
| Formato | **TIF original do equipamento** com metadados (e o .txt do JEOL/Hitachi) | JPG/PNG exportado perde metadados e às vezes é redimensionado |

## 4. Aumento certo e quantas imagens

- O menor grânulo de interesse deve ter **≥ 10 px de diâmetro** (≥ 20 px para forma confiável); o
  maior deve caber bem dentro do campo (≤ ~1/5 da largura), senão a borda corta os grandes e a
  distribuição fica enviesada para pequenos (o script exclui os cortados e avisa quando são muitos).
- Distribuição larga (batata, trigo A+B): pode exigir dois aumentos - um para os pequenos, outro
  para os grandes - analisados juntos (em µm) cuidando para não contar o mesmo campo duas vezes.
- **Número**: ≥ 300-500 grânulos válidos por amostra, de **≥ 5-10 campos aleatórios** (não escolha
  "os mais bonitos"), idealmente de 2-3 preparos independentes. Para D10/D90 e populações, mais.
  (ISO 13322-1 trata o número de partículas necessário para a precisão desejada.)
- Filmes: 3-5 campos por face e por amostra para defeitos; espessura em ≥ 3 fraturas.

## 5. Artefatos e como reconhecê-los

| Artefato | Aparência | O que fazer |
|---|---|---|
| **Carregamento** | regiões brancas saturadas, riscas horizontais, imagem que "desliza", bordas superbrilhantes | recobrimento melhor, aterramento com fita de carbono/tinta de prata, kV menor, baixo vácuo, varredura rápida + integração |
| **Dano do feixe** | bolhas, buracos, trincas que aparecem/crescem durante a observação, encolhimento, área retangular escurecida da varredura anterior | kV e corrente menores, menos tempo por campo, focar ao lado |
| Trincas de secagem/vácuo | trincas retas em filmes muito secos | comparar com amostra condicionada; baixo vácuo |
| Recobrimento | grãos de ouro (aumentos > 50 k×), "craquelado" do metal | recobrimento mais fino, Au-Pd, Pt |
| Fita de carbono | textura do adesivo como fundo, grânulos afundados | menos pó, pressionar menos, outro suporte |
| Sobreposição/aglomeração | grânulos empilhados, contornos ocultos | preparo mais disperso; não medir agregados |
| Contaminação | partículas estranhas (poeira, fibras, cristais de sal) | limpeza, EDS |
| Sombra do detector | gradiente de iluminação de um lado para outro | `--flatten` ou `--thresh local`; posição/tipo do detector |
| Imagem redimensionada | metadados ≠ barra desenhada | usar o TIF original |

## 6. Como reportar (métodos e legendas)

**Métodos**: equipamento e fabricante; preparo (dispersão, fita, recobrimento: metal, espessura ou
tempo); tensão, detector, WD, vácuo; aumentos; número de campos/imagens e de grânulos medidos;
software (esta skill: Python, scikit-image) e critérios: limiar (Otsu/local), separação de grânulos
encostados (watershed na transformada de distância), exclusão dos objetos que tocam a borda,
diâmetro mínimo; definições (ECD, Feret, AR, circularidade). Ex.:

"As amostras foram dispersas sobre fita de carbono, recobertas com ouro (15 nm) e observadas em MEV
(JEOL JSM-6610LV) a 5 kV com detector de elétrons secundários. O tamanho e a forma dos grânulos
foram determinados por análise de imagens (Python/scikit-image) em 8 campos por amostra (≥ 400
grânulos): segmentação por limiar de Otsu, separação de grânulos em contato por watershed e
exclusão dos grânulos que tocavam as bordas. Reportam-se o diâmetro de área projetada equivalente
(ECD; D10, D50, D90 em número e Dv50 em volume), a razão de aspecto (Feret máx/mín) e a
circularidade (4πA/P²)."

**Legenda**: material, tratamento, aumento **e** barra de escala (a barra vale mais que o aumento
nominal), detector/kV se relevante. "Figura 3. Micrografias de MEV (5 kV, SE) dos grânulos de
amido de mandioca: (a) nativo; (b) hidrolisado com α-amilase por 4 h. Barras: 10 µm."
