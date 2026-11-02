# Processamento, referência, integração e armadilhas

## Referência

- Recomendação IUPAC: TMS = 0 (¹H e ¹³C) ou, na prática, sinal residual do solvente (Fulmer 2010):
  CDCl3 7,26 / 77,16; DMSO-d6 2,50 / 39,52; CD3OD 3,31 / 49,00; acetona-d6 2,05 / 29,84 e 206,26;
  C6D6 7,16 / 128,06; CD3CN 1,94 / 1,32 e 118,26; CD2Cl2 5,32 / 53,84; D₂O 4,79 (HDO, varia com T;
  prefira DSS/TSP em água).
- A skill reconhece o padrão do solvente (acoplamento com ²H: CDCl3 ¹³C 1:1:1, J 32 Hz; CD2Cl2 1:2:3:2:1,
  27 Hz; DMSO-d6 septeto, 21 Hz; CHD₂ ¹H quintetos, ~2 Hz) e só corrige desvios pequenos automaticamente.
  Se o solvente residual estiver sob um sinal da amostra, a referência do arquivo é mantida.
- `--ref tms` (pico em ~0 → 0,00), `--ref 7.31=7.26` (pico observado → valor), `--ref file`.
- Espectros de ¹³C e DEPT do mesmo experimento recebem o mesmo deslocamento (o DEPT não mostra solvente).

## Fase e linha de base

- Espectro processado com fase ruim (área negativa > 12% dos pontos significativos) é refasado: parte
  imaginária do arquivo ou reconstruída por transformada de Hilbert + ACME (entropia + penalidade de
  negativos). Confira a figura; DEPT/APT não são refasados (têm sinais negativos legítimos).
- Linha de base: suavização de Whittaker ponderada com exclusão iterativa dos sinais (funciona com sinais
  positivos e negativos). Corcundas largas (sinais muito largos, supressão de água, filtro) podem ser
  parcialmente removidas - cuidado com integrais de sinais largos.

## FID

- Bruker: filtro digital removido (GRPDLY ou tabela DECIM/DSPFVS) - sem isso, fase de 1ª ordem enorme.
- LB 0,3 Hz (¹H) / 1 Hz (¹³C), zero-filling ×2, FFT, fase automática; `--lb`, `--zf` em `process`.
- Eixo: portadora ($SFO1 − $SF), $OFFSET, ou .SHIFT REFERENCE; depois o solvente confere.

## Integração

- Regiões contíguas acima de k·σ, separadas em vales profundos entre sinais distantes (> J máx), unidas
  quando juntas formam um multipleto de 1ª ordem; fragmentos pequenos de multipletos complexos são unidos.
- Normalização: pela fórmula (total de H; testa também com 1-3 H trocáveis ausentes), por um sinal
  (`--norm 2.29=3`) ou automática (menor desvio de inteiros). Sinais < 0,03 H relativos são descartados;
  satélites de ¹³C (0,55% de cada lado, a ±60-85 Hz) e impurezas comuns são rotulados.
- Erros de integral: relaxação (d1 curto), saturação, linha de base, sobreposição, supressão de água,
  satélites de ¹³C dentro da janela. Integrais de ¹³C comum não são quantitativas (NOE e T1).

## Multipletos

- Linhas do multipleto (≥ 4% da maior) → unidades inteiras de intensidade (Pascal) → desconvolução em
  dubletos sucessivos (árvore de J). Simetrização antes (remove o efeito de telhado).
- `m` = sem solução de 1ª ordem (2ª ordem, sobreposição, acoplamentos não resolvidos). `AB` = quatro
  linhas com internas mais altas (δA, δB e J calculados). `sl` = sinal largo (> 3,5 Hz, trocáveis).
- Resolução digital: J < 2× (Hz/ponto) não é resolvido; use zero-filling e LB menor.

## DEPT e APT

- DEPT-135: CH e CH₃ positivos, CH₂ negativo, C ausente. DEPT-90: só CH. → CH₃ = positivo no 135 e
  ausente no 90. Fase dos DEPT às vezes vem invertida: confira com um CH₂ óbvio.
- APT/JMOD: C e CH₂ de um sinal, CH e CH₃ do outro (a convenção varia entre equipamentos: a skill assume
  C/CH₂ positivos; inverta mentalmente se for o caso).
- Σ H (DEPT) vs fórmula → H em heteroátomos.

## Armadilhas de aquisição (para comentar com o usuário)

- d1 curto e pulso de 90° em ¹H quantitativo; 30° + d1 1 s é ok para rotina, não para qRMN.
- Amostra concentrada demais: linhas largas, desvio de δ (OH/NH), ¹³C demorado mas ok.
- Shimming ruim: linhas assimétricas → multipletos "m" falsos.
- Rotâmeros (amidas terciárias, carbamatos): sinais duplicados/largos à temperatura ambiente; aqueça.
- Troca química (OH/NH, equilíbrios conformacionais): sinais largos ou ausentes.
