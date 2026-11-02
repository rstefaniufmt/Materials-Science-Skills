# DSC e DTA: interpretação

## Sumário
1. Sinal, unidades e direção endo/exo
2. Eventos: como reconhecer
3. Transição vítrea (Tg)
4. Fusão, cristalização e cristalinidade
5. Integração de picos (ΔH)
6. OIT (tempo de indução oxidativa)
7. Cp com safira
8. Cura, pureza, polimorfismo, proteínas
9. DTA
10. Calibração e verificação
11. Armadilhas
12. Como reportar

## 1. Sinal, unidades e direção endo/exo

- DSC de fluxo de calor mede a diferença de fluxo amostra-referência: mW, ou W/g
  normalizado pela massa. DTA mede ΔT (K, µV). NETZSCH bruto vem em µV; com a curva de
  sensibilidade do arquivo o script converte para mW (e W/g com a massa).
- **Direção**: cada fabricante desenha de um jeito (TA/Mettler/NETZSCH: endo ↓ por padrão;
  PerkinElmer: endo ↑). Internamente o script orienta tudo com endo positivo e desenha com
  `--display exo_up` (padrão, convenção ICTAC/ASTM moderna) ou `endo_up`. A figura sempre
  traz a seta. Confirme com um evento de natureza conhecida.
- Linha de base: no aquecimento o sinal é proporcional a Cp·β; degraus de linha de base
  = mudança de Cp (Tg, fusão, reação); inclinação = Cp variando com T e assimetria do
  sensor.

## 2. Eventos: como reconhecer

| Aparência (exo ↑) | Evento provável | Confirmação |
|---|---|---|
| degrau para baixo (endo) | Tg | some/permanece no 2º aquecimento; sem massa perdida |
| degrau + pico endo pequeno | Tg com relaxação entálpica (envelhecimento) | pico some no 2º aquecimento |
| pico exo acima de Tg, no aquecimento | cristalização a frio | amostra amorfa/temperada; some após recozimento |
| pico endo agudo | fusão (cristal), transição sólido-sólido | reaparece no 2º aquecimento; exo correspondente no resfriamento |
| pico endo largo < 150 °C | evaporação de água/solvente, desidratação | perda de massa (TG/STA); muda com cadinho aberto/fechado |
| pico exo no resfriamento | cristalização (Tc) | |
| pico exo largo, irreversível | cura, reação, oxidação, decomposição exo | não volta no 2º aquecimento |
| pico endo + perda de massa | decomposição/desidratação | STA |
| degrau pequeno sem pico, só no 1º aquecimento | relaxação de tensões, transiente | |
| subida/queda forte no início/fim do segmento | transiente instrumental | ignore (o script não quantifica `borda`) |

O detector automático (`find_events`) dá **candidatos** com uma sugestão; a palavra final
vem do contexto (material, programa, TG).

## 3. Transição vítrea (Tg)

- `glass_transition(r, (T1, T2))` / `--tg T1:T2`: ajusta retas na linha de base antes e
  depois e a tangente no ponto de inflexão. Retorna Tg_onset (interseção), Tg_midpoint
  (meia altura entre as linhas de base - **padrão para reportar**, ISO 11357-2/ASTM E1356),
  Tg_inflexão, Tg_endset e ΔCp (J g⁻¹ K⁻¹, se o sinal está em W/g e β é conhecido).
- A janela deve incluir uns 10-20 K de linha de base reta de cada lado.
- Tg depende de β (+~3 K por década), da história térmica e de umidade/plastificante. Use o
  2º aquecimento (após resfriamento controlado) para propriedade do material e o 1º para
  o estado da amostra como recebida.
- Relaxação entálpica (pico endo no fim do degrau) desloca o onset; o midpoint é mais robusto.
- Tg de blendas: uma Tg intermediária = miscível; duas Tg próximas às dos componentes =
  imiscível. Fox: 1/Tg = w₁/Tg₁ + w₂/Tg₂ (K).
- Tg fraca (compósitos muito carregados, termofixos curados): aumente a massa e β, ou use
  DSC modulado (reversível).

## 4. Fusão, cristalização e cristalinidade

- T_pico de fusão (Tm) é o reportado normalmente; para substâncias puras e metais, o
  **onset** é a temperatura termodinâmica (calibração, pureza).
- **Cristalinidade**: Xc = (ΔHm − ΔHcc) / (ΔH⁰m · w) × 100
  - ΔHm: entalpia de fusão (positiva); ΔHcc: cristalização a frio no mesmo aquecimento
    (positiva, se houver); ΔH⁰m: fusão do polímero 100% cristalino (`ta_dsc.DH0`, tabela
    em reference_data.md); w: fração mássica do polímero (compósito/blenda).
  - `--polymer PET` ou `--dh0 140`; `--w 0.7`.
  - No modo automático, o script soma as cristalizações a frio do segmento e calcula Xc
    "inicial" (estado como recebido, com ΔHcc) e a fusão total.
- Fusão múltipla (dois picos): populações de cristais, recristalização durante o
  aquecimento, polimorfos (PLA α/α', PP α/β, PA6 α/γ). Integre o conjunto para Xc e
  descreva os picos.
- Tc e o super-resfriamento (Tm − Tc) indicam nucleação (agentes nucleantes, cargas sobem Tc).
- Pré-fusão de PET/PLA logo após a cristalização a frio: janelas mal escolhidas somam
  eventos de sinais opostos - mantenha a cristalização e a fusão em janelas separadas.

## 5. Integração de picos (ΔH)

- `integrate_peak(r, (T1, T2), baseline='linear'|'sigmoidal')` / `--peak T1:T2`: integra no
  tempo (ΔH = ∫(q − base) dt / m) e dá T_onset/T_endset por tangente, T_pico, altura,
  largura a meia altura. ΔH > 0 = endo.
- Linha de base **linear** entre as bordas para picos com pouca mudança de Cp; **sigmoidal**
  quando a linha de base é diferente antes e depois (fusão com ΔCp, cura).
- As bordas da janela precisam estar na linha de base: janela curta subestima ΔH; janela que
  invade outro evento mistura áreas.
- ΔH em J/g exige massa constante durante o evento - com perda de massa (desidratação,
  decomposição) o número é aparente.
- `cal_factor`: multiplica pelo fator de calibração de entalpia se o arquivo não estiver
  calibrado (ex.: 28,62/ΔH_medido do índio).

## 6. OIT (tempo de indução oxidativa)

- ISO 11357-6 / ASTM D3895: aquecer em N₂ até a T de ensaio (200 °C para PE/PP típico),
  estabilizar, trocar para O₂ (ou ar) e medir o tempo até o início da oxidação (exo).
- `oit(r)` / `ta.py oit`: o instante da troca vem do canal de vazão de gás (NETZSCH, TA) ou
  `--t-switch min`. OIT = t_onset (interseção da linha de base isotérmica com a tangente no
  ponto de maior inclinação antes do máximo exo) − t_switch.
- OIT dinâmico (OOT, onset em rampa sob O₂): use `--peak` no exo e o T_onset.
- Reporte T, gás, vazão, cadinho (Al aberto; Cu catalisa), massa.

## 7. Cp com safira

- Três corridas idênticas (mesmo cadinho, programa isot.-rampa-isot.): branco (cadinhos
  vazios), safira (padrão), amostra. Cp = Cp_safira · (DSC_amostra − DSC_branco) /
  (DSC_safira − DSC_branco) · m_safira / m_amostra (ASTM E1269, ISO 11357-4).
- `ta.py cp --sample --blank --sapphire --m-sample --m-sapphire --range 50:300`. A tabela de
  Cp da safira embutida é aproximada (tipo JANAF/NIST α-Al₂O₃); para trabalho de
  referência, use a do certificado.
- Erros típicos: cadinhos/posições diferentes entre corridas, linha de base sem
  estabilização, massa de safira muito diferente da amostra (em capacidade térmica).

## 8. Cura, pureza, polimorfismo, proteínas

- **Cura de termofixos (epóxi, poliéster, fenólicas)**: integre o exo de cura (ΔH_total na
  resina não curada; ΔH_residual na peça). Grau de cura = 1 − ΔH_residual/ΔH_total. A Tg
  sobe com o grau de cura (DiBenedetto). Cinética de cura: `kinetics --signal hf`.
- **Pureza (van 't Hoff, ASTM E928)**: 1/F vs Ts para o pico de fusão de substância ≥ ~98,5%
  pura; não implementado - pode ser feito com `integrate_peak` + área parcial. Diga que é
  estimativa e que não vale para impurezas que formam solução sólida.
- **Polimorfismo (fármacos)**: fusão de forma metaestável seguida de recristalização exo e
  fusão da forma estável; varie β para confirmar (a recristalização some em β alto).
  Pseudopolimorfos (hidratos/solvatos): endo com perda de massa.
- **Proteínas (microDSC)**: Tm de desnaturação, ΔH calorimétrico; requer subtração de
  tampão e linha de base progressiva - não é o foco dos scripts.

## 9. DTA

- Mesmo raciocínio qualitativo do DSC (picos endo/exo, Tg como degrau), mas o sinal ΔT não
  é diretamente proporcional ao fluxo de calor: ΔH só com calibração de área contra
  padrões na mesma faixa de T. Use DTA para T de eventos e caráter endo/exo.
- Em STA (TG-DTA/TG-DSC) cruze cada etapa de massa com seu evento térmico: perda + endo =
  desidratação, descarboxilação, decomposição; perda + exo = combustão/oxidação; sem perda
  + endo = fusão/transição sólido-sólido; sem perda + exo = cristalização, transição de
  fase de óxido (ex.: anatase→rutilo), cura; ganho + exo = oxidação de metal.

## 10. Calibração e verificação

- Temperatura e entalpia com padrões de fusão (In, Sn, Zn, Bi, Pb, Al; tabela em
  reference_data.md), na mesma β e cadinho das amostras. Índio é o teste rápido: onset
  156,6 °C e 28,6 J/g. Desvios > 0,5 K ou > 2% pedem recalibração.
- Use o **onset** do padrão (não o pico) para calibrar temperatura.

## 11. Armadilhas

- Transientes de início de rampa, troca de segmento e fim de corrida não são eventos.
- Cadinho deformado/tampa mal fechada → linha de base ruim; evaporação em cadinho aberto
  gera endo enorme.
- Amostra grande + β alto → picos largos e deslocados; amostra pequena → Tg fraco.
- Comparar 1º aquecimento de amostras com histórias diferentes leva a conclusões erradas.
- "Pico" no DSC acima da Td vista na TG é decomposição - não é fusão.
- Sinal em mW sem massa → não há J/g; diga e peça a massa.

## 12. Como reportar

Massa, cadinho (Al selado/furado, aberto), gás e vazão, programa (taxa, faixa, segmentos,
isotermas), segmento analisado. Por evento: tipo, T_onset/T_pico (Tg midpoint), ΔH (J/g) ou
ΔCp, Xc com ΔH⁰ e w. Figura com seta exo/endo e segmentos identificados.
