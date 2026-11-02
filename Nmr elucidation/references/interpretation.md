# Da lista de sinais à estrutura: roteiro com RMN 1D (e quando pedir mais)

## Sumário
1. Antes de interpretar
2. Inventário (o que os dados dizem sozinhos)
3. Fragmentos e sistemas de spin
4. Juntando os fragmentos
5. Quando 1D não basta: quais 2D pedir
6. Desreplicação e candidatas
7. Como apresentar

## 1. Antes de interpretar

- **Fórmula molecular** (EM de alta resolução, ± análise elementar). Sem ela, o espaço de estruturas é
  aberto demais. Com ela: IDH = C − H/2 − X/2 + N/2 + 1 (anéis + ligações π). Benzeno consome 4,
  C=O consome 1, cada anel 1.
- **Pureza**: sinais pequenos sem integral inteira, solvente residual, água, graxa → `impurities`.
  Mistura de isômeros/rotâmeros: conjuntos paralelos de sinais com integrais fracionárias constantes.
- **Solvente certo**: H trocáveis só são vistos em DMSO-d6 (bem) ou CDCl3 (às vezes); somem em CD3OD/D₂O.

## 2. Inventário (o que os dados dizem sozinhos)

`report` gera um resumo. Leia assim:
- **¹³C**: nº de sinais vs nº de C (simetria); C=O > 160 (ácido/éster/amida 165-180, cetona > 195,
  aldeído 190-205); sp² 100-160 (aromáticos com H 115-135; C–O aromático 145-165); C–O sp³ 60-90;
  anoméricos 90-110; sp³ < 60.
- **DEPT**: nº de CH₃, CH₂, CH e C. Σ H em carbono vs H da fórmula → nº de H em O/N.
- **¹H**: H aromáticos (6,5-8,5: padrão de substituição pela multiplicidade - ver abaixo), vinílicos
  (4,5-6,5; J trans 12-18, cis 6-12), O–CH (3,3-4,5), CH₃ singletos (OCH₃ 3,7-3,9; NCH₃ 2,2-3,0; ArCH₃ 2,3;
  CH₃CO 2,0-2,2; CH₃ de acetato 2,0-2,1), CHO 9,5-10,1, COOH 10-13.

Padrões aromáticos (¹H, 1ª ordem aproximada):
- monossubstituído: 5H (2H d/m, 2H t/m, 1H t/m), ou multipleto 5H;
- 1,4 (para): dois "dubletos" de 2H (AA'BB', J ≈ 8-9 Hz);
- 1,2 (orto): quatro sinais 1H: dd, td, td, dd (J orto 7-8, meta 1-2);
- 1,3 (meta): s largo (t, J ≈ 2), d, t, d;
- 1,2,4: d (J ≈ 2), dd (8, 2), d (8);
- 1,2,3: d, t, d (ou dd, t, dd);
- 1,3,5: singleto(s) 1:2 ou três s largos.

## 3. Fragmentos e sistemas de spin

- `parceiros_J` lista multipletos com o mesmo J: candidatos a acoplados entre si. Construa cadeias
  (A–B–C) pela compatibilidade de J e nº de vizinhos (t → 2 H vizinhos equivalentes; q → 3; dd → 2 H
  diferentes...). Confirme com COSY.
- Fragmentos clássicos: etila (t 3H + q 2H, J 7; OCH₂ 4,1-4,4 = éster; NCH₂ 2,5-3,5), isopropila
  (d 6H + sept 1H), terc-butila (s 9H, 0,9-1,5), n-propila (t, sext, t), metoxila (s 3H 3,7-3,9 / δC 55-56
  aromática, 51-52 éster), acetila (s 3H 2,0-2,6 / δC 20-30 + C=O), vinila (dd 1H 5,5-6,8 + 2 dd 1H 5,0-5,9).
- Açúcares (glicosídeos): H anomérico 4,3-5,6 (d; J 7-8 Hz = β-gluco/galacto diaxial; 1-4 Hz = α),
  δC anomérico 95-110, demais CH–O 60-85, CH₂OH 60-62; ramnose: CH₃ d (J ≈ 6) em 0,9-1,3 / δC 17-18.
- Flavonoides: OH-5 quelado ~12-13 ppm (s, DMSO); C-4 (C=O) 175-182; anel A 5,7-dihidroxi: H-6/H-8
  d (J ≈ 2) em 6,1-6,5; anel B 3',4'-di-OH: 1,2,4 (d, dd, d).
- Terpenos/esteroides: muitos CH₃ singletos/dubletos 0,7-1,3; C quaternários 30-50; olefinas trissubstituídas.

## 4. Juntando os fragmentos

1. Some átomos e IDH dos fragmentos e compare com a fórmula: o que sobrou (heteroátomos, C quaternários,
   anéis) precisa ser encaixado.
2. Use deslocamentos para posicionar heteroátomos (C–O vs C–N vs C–C) e conjugação (enonas: H β 6,8-7,7,
   C β 140-160; C=O conjugada 165-200).
3. Gere TODAS as montagens compatíveis (não só a primeira que parece certa) e passe-as em `candidates`.
4. Estereoquímica relativa: J (axial-axial 8-13 Hz em anéis de 6; cis/trans em alcenos) e NOE; absoluta:
   não sai de RMN sem derivatização (Mosher) ou cálculo + DC/ROA/raio X.

## 5. Quando 1D não basta: quais 2D pedir

| Pergunta | Experimento |
|---|---|
| Qual H está em qual C? | HSQC (editado: CH/CH₃ × CH₂ com sinais opostos) |
| Quem acopla com quem (H–H, 2-3 ligações)? | COSY (e TOCSY para sistemas inteiros, ex.: açúcares) |
| Ligações entre fragmentos, posição de C quaternário/heteroátomo | HMBC (²J/³J C–H) |
| Proximidade espacial (estereoquímica, conformação) | NOESY/ROESY (ROESY para MM 700-1500) |
| N | ¹H-¹⁵N HMBC |
| Sobreposição em ¹H | HSQC ou ¹H em outro solvente (C6D6 desloca), campo maior, pure-shift |
| Pureza/mistura | DOSY |

Diga ao usuário qual pergunta concreta cada 2D responderia no caso dele.

## 6. Desreplicação e candidatas

- **Composto conhecido?** Fórmula + ¹³C são as melhores chaves. Busque na literatura do organismo/síntese,
  PubChem (pela fórmula), COCONUT, LOTUS, NPAtlas (produtos naturais), SDBS e nmrshiftdb2 (espectros
  experimentais). Compare a lista de δC (diferenças ≤ 0,5-1 ppm em todos os C, mesmo solvente, = mesmo
  composto, salvo estereoisômeros).
- **Candidatas para `candidates`**: o produto esperado + regioisômeros + produtos de rearranjo/
  epimerização + compostos da literatura com a mesma fórmula. Para glicosídeos: posição de ligação do
  açúcar e tipo de ligação interglicosídica; para aromáticos: posições de substituição.
- **Inédito de verdade**: a skill ajuda a montar e testar hipóteses (fragmentos, contagens, J,
  predição), mas a prova vem de 2D + EM; não declare estrutura nova só com 1D.

## 7. Como apresentar

Tabela de atribuição (posição | δC, tipo | δH (mult., J em Hz, nH)), com numeração da estrutura; nas
notas: solvente, frequência, referência (solvente residual: CDCl3 7,26/77,16), e o que foi confirmado por
2D. A lista "δ 7.27 (dd, J = 9.0, 6.8 Hz, 1H), ..." vai para a parte experimental.
