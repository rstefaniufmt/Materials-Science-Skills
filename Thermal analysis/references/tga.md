# TG/TGA e DTG: interpretação

## Sumário
1. Leitura da curva
2. Temperaturas características
3. Etapas, sobreposição e DTG
4. Branco, empuxo e artefatos
5. Estequiometria e hidratos
6. Polímeros, compósitos e biomassa
7. Atmosfera e condições
8. Como reportar

## 1. Leitura da curva

- Eixo y: massa (%) relativa à massa inicial (ou à massa numa T de referência com `--basis`).
  DTG = dm/dT (%/°C) ou dm/dt (%/min); `ta_tga.find_steps` dá os dois. A DTG é negativa
  para perda; o pico (mínimo) marca a taxa máxima de perda.
- Patamar = ausência de reação. Cada degrau entre patamares é uma etapa; a DTG separa
  degraus que a TG mostra como uma rampa contínua.
- Ganho de massa real: oxidação de metais/ligas, adsorção, carbonatação, formação de
  óxidos de valência maior. Ganho de poucos décimos de % ao longo de toda a corrida é quase
  sempre empuxo (ver §4).

## 2. Temperaturas características

| Grandeza | Definição | Uso |
|---|---|---|
| T_onset (extrapolado) | interseção da tangente ao patamar anterior com a tangente no ponto de inflexão (DTG máx) - ISO 11358 | início "termodinâmico" de uma etapa isolada |
| T_endset | idem, com o patamar posterior | fim da etapa |
| T_pico (DTG) | máximo da taxa de perda | identificar/comparar etapas; Kissinger |
| T1%, T2%, T5%, T10%, T50% | T em que a perda acumulada atinge x% | estabilidade térmica comparativa (T5% é o mais usado; Td5%) |
| T_inicio / T_fim (script) | bordas da etapa no patamar (mínimos da DTG) | integração da perda |

O onset depende da escolha das tangentes; com etapas sobrepostas ou patamar inclinado ele
perde o sentido - o script não reporta onset fora do intervalo da etapa. Para "estabilidade
térmica" de polímeros, use T5% (ou T2% se há umidade) **em base seca** se houver perda de
água/solvente antes.

## 3. Etapas, sobreposição e DTG

- `find_steps` acha os picos da −DTG, separa as etapas pelos mínimos da DTG entre eles e
  integra a perda entre essas bordas. `min_loss` (padrão 0,5%) descarta etapas pequenas;
  `rel_prominence` e `min_sep_K` controlam a separação.
- `sobreposta = True`: a DTG não volta perto de zero entre duas etapas. A divisão pelo mínimo
  é convencional - a perda de cada uma é aproximada (erro de alguns %). Para separar
  melhor: taxa menor (2-5 K/min), massa menor, ou deconvolução da DTG (não implementada:
  Fraser-Suzuki/assimétricas; ajuste com cuidado e diga que é modelo).
- Ombro na DTG = processo sobreposto não resolvido; mencione.
- Etapas em isotermas (`tipo isoterma`): perda durante um patamar de T (secagem, programas
  de cinzas). Em programas com troca de gás (N₂ → ar para queimar carbono), a etapa após a
  troca é a combustão do carbono/negro de fumo.
- Programas com resfriamento: só aquecimentos/isotermas entram em `find_steps`; o
  segmento usado aparece na coluna `segmento`.

## 4. Branco, empuxo e artefatos

- **Empuxo/flutuação**: a densidade do gás cai com T, a massa aparente sobe (dezenas de µg).
  Correção: corrida de branco (cadinho vazio, mesmo programa, gás e vazão) subtraída
  (`--blank branco.txt`; `--blank-time` se o programa tem isotermas). Arquivos NETZSCH
  `.ngb-ss3` com correção já aplicada não precisam disso.
- **Picos espúrios na DTG**: troca de gás, vibração, bolhas/projeção de amostra (perda
  súbita, DTG muito estreita) - o script filtra picos de DTG muito estreitos.
- **Amostra que funde, espuma ou projeta** distorce a curva; cadinho com tampa furada ajuda.
- **Primeiro minuto**: transiente de estabilização (vazão/temperatura); o script ignora.
- **Termopar**: a temperatura medida é perto da amostra, não nela; com β alto o atraso cresce
  (T_pico sobe). Calibração de temperatura em TG: pontos de Curie (ver reference_data.md).

## 5. Estequiometria e hidratos

- `--theory "Reagente:esp1,esp2,..."` calcula as perdas teóricas de cada espécie liberada e
  compara com as etapas medidas (diferença e `ok` com tolerância 1,5 pontos percentuais).
  Coeficientes e hidratos: `"CuSO4·5H2O:2H2O,2H2O,H2O"`, `"Mg(OH)2:H2O"`.
- `stoich` faz só o cálculo teórico; `hydrate --loss x --anhydrous Sal` dá o nº de águas:
  n = (perda/100 · M_anidro) / ((1 − perda/100) · M_H₂O).
- A perda teórica é relativa à massa **inicial** do reagente indicado; se a amostra já perdeu
  umidade antes, use `--basis` para recalcular.
- Discrepâncias típicas: água adsorvida (perda maior na 1ª etapa), hidratação incompleta,
  carbonatação da amostra (CO₂ extra), impureza, etapas sobrepostas.
- Interprete a identidade do gás com cautela: TG só dá massa; TG-MS/TG-FTIR ou STA confirmam.

## 6. Polímeros, compósitos e biomassa

- **Umidade/voláteis** até ~150 °C; **aditivos/plastificantes** 150-300 °C; **cadeia
  polimérica** em uma ou mais etapas; **carvão** (N₂) e **cinzas/carga** no resíduo.
- **Teor de carga** (vidro, talco, sílica, TiO₂) = resíduo em ar a 600-800 °C, corrigido pelo
  resíduo do polímero puro na mesma condição. Cargas que decompõem: CaCO₃ (−44% CO₂ acima
  de ~600 °C), Mg(OH)₂ e Al(OH)₃ (retardantes, perdem água a 300-400 / 200-350 °C),
  negro de fumo (só queima em ar/O₂).
- **Blendas**: etapas separadas permitem estimar a composição se as Td forem distintas
  (cada componente puro como referência).
- **Borracha** (método tipo ASTM E1131 / ISO 9924): voláteis/óleo em N₂, polímero em N₂,
  troca para ar → negro de fumo, resíduo = cinzas.
- **Biomassa/lignocelulósicos** (N₂, 10 K/min): umidade < 150 °C; hemicelulose ~220-315 °C
  (ombro), celulose ~315-400 °C (pico DTG ~340-360), lignina larga 200-700 °C. A análise
  imediata (proximate analysis) segue umidade → voláteis (N₂ até ~900 °C) → carbono fixo
  (troca para ar) → cinzas.
- **PVC**: 1ª etapa ~250-350 °C = HCl (58,3% teórico do polímero puro) + plastificante;
  2ª ~400-500 °C = cadeia poliênica.
- **Polímeros oxidam antes em ar** (Td menor, etapas extras, resíduo menor). Diga a
  atmosfera sempre.

## 7. Atmosfera e condições

- N₂/Ar: pirólise (carvão residual). Ar/O₂: oxidação/combustão (resíduo = cinzas, óxidos).
  Gás reativo (CO₂, H₂, vapor) muda as reações (ex.: CaCO₃ decompõe mais tarde sob CO₂).
- Taxa β: 10 K/min é o padrão; 20 K/min desloca T_pico ~10-20 °C para cima. Massa 5-15 mg;
  amostras maiores pioram a resolução e o gradiente térmico.
- Cadinho: alumina (até ~1600 °C, inerte para a maioria), platina (não usar com metais/
  redutores - forma liga), alumínio (≤ 600 °C).

## 8. Como reportar

Tabela por etapa: faixa de T (início-fim), T_onset, T_pico (DTG), perda (%), atribuição
(e teórico, se houver), mais T5% e resíduo a T definida. Na legenda/método: massa, cadinho,
atmosfera e vazão, β, branco subtraído ou não, base (inicial ou seca).
