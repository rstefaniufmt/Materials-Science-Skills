# Ranqueamento de estruturas candidatas

## Sumário
1. O que é calculado
2. Predição (CSP5)
3. Atribuição ótima
4. Escores e probabilidade relativa
5. Como ler o resultado
6. Limitações
7. Como reportar

## 1. O que é calculado

Para cada candidata (SMILES): fórmula, IDH, classes de simetria (RDKit, ranking canônico sem desempate)
→ nº de sinais de ¹³C esperados, multiplicidades (CH₃/CH₂/CH/C), grupos de H equivalentes (CH₂
diastereotópicos são separados quando o C está em anel ou a molécula é quiral), H trocáveis. Isso já
permite descartar candidatas por contagem (DEPT, nº de sinais, integrais) mesmo sem predição.

## 2. Predição (CSP5)

`pip install csp5` (MIT; modelos inclusos, roda offline em CPU; gera conformação 3D com RDKit e prediz
δ ¹³C e ¹H por átomo, condicionado ao solvente: cdcl3, dmso-d6, cd3od, acetone-d6, cd2cl2, c6d6, cd3cn,
d2o, thf-d8, toluene-d8, dmf-d7). Erro típico de modelos desse tipo em moléculas orgânicas comuns:
~1-2 ppm (¹³C) e ~0,1-0,2 ppm (¹H); maior para metais, íons, tautômeros, moléculas muito flexíveis e
classes raras. `--conformers 5` usa média de conformações (mais lento).

## 3. Atribuição ótima

- **¹³C**: classes previstas → picos experimentais por atribuição linear (Hungarian). Um pico pode receber
  até 2 classes (sobreposição) com penalidade; classe sem pico e pico sem dono são penalizados. Com DEPT, o
  tipo do C precisa bater (penalidade alta se não bater).
- **¹H**: grupos de H equivalentes → multipletos, por programação linear inteira: cada grupo vai inteiro
  para UM multipleto (H equivalentes têm um só δ), a integral do multipleto é a capacidade (penalidade por H
  faltando/sobrando); dentro de um "m" (faixa), erro zero se o previsto cai na faixa. Trocáveis (OH/NH)
  vão para sinais largos ou ficam sem atribuição sem custo.

## 4. Escores e probabilidade relativa

- MAE e erro máximo de ¹³C e ¹H; C/H sem sinal; sinais experimentais sobrando; desajuste de integrais.
- log-verossimilhança com t de Student (ν = 4; σ = 2,0 ppm ¹³C e 0,22 ppm ¹H) somada sobre as classes,
  com penalidades para faltas → **prob_rel_pct** = exp(logL) normalizada entre as candidatas (ideia do
  DP4 de Smith & Goodman, J. Am. Chem. Soc. 132 (2010) 12946; aqui com predição por ML, não DFT).
- Candidatas com fórmula diferente da informada (`--formula`) recebem probabilidade 0.

## 5. Como ler o resultado

- **Compatível**: MAE ¹³C ≲ 1,5-2 ppm, nenhum erro > 5 ppm, MAE ¹H ≲ 0,15-0,2 ppm, contagens (DEPT,
  sinais, integrais) coerentes.
- **Prob. alta + MAE ruim** → nenhuma candidata serve; amplie a lista.
- **Duas candidatas próximas** (prob. 40-60% e MAE parecidos) → 1D não decide; diga qual 2D decidiria
  (HMBC para posição de ligação; NOESY para estereo).
- Diferenças concentradas em 1-2 átomos (ex.: C da ligação glicosídica) são as mais informativas: olhe a
  tabela de atribuição (`*_atribuicao.csv`) e a figura (laranja = erro grande).
- Estereoisômeros (epímeros, E/Z, diastereoisômeros): predições 2D/1 conformação distinguem mal; use
  `--conformers`, J e NOE, ou DP4+ (DFT) fora desta skill.

## 6. Limitações

- Só escolhe entre as candidatas dadas.
- Não modela efeitos de concentração, pH, agregação, ligação de hidrogênio intramolecular forte além do
  que o modelo aprendeu; sais e zwitteríons devem ser desenhados no estado em que estão no solvente.
- Sinais sobrepostos/omitidos no experimental reduzem a separação entre candidatas.
- Probabilidade relativa ≠ certeza: valide com 2D, EM, comparação com padrão autêntico/literatura.

## 7. Como reportar

"Entre as N estruturas consideradas (lista no material suplementar), X apresentou o melhor ajuste aos
deslocamentos experimentais (MAE ¹³C = 0,7 ppm, ¹H = 0,07 ppm; predição CSP5; probabilidade relativa 99%),
com atribuição consistente com o DEPT. A estrutura foi confirmada por HSQC/HMBC (...)". Nunca apresente a
probabilidade relativa como prova isolada.
