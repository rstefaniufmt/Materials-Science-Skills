# Band gap óptico a partir de UV-Vis

## Sumário
1. O que o método mede (e o que não mede)
2. Escolher F: absorbância, α ou Kubelka-Munk
3. Escolher a transição (expoente)
4. Escolher a região de ajuste
5. Absorção sub-gap e o método da linha de base
6. Energia de Urbach
7. Métodos alternativos (derivada, onset)
8. Armadilhas
9. Como reportar

## 1. O que o método mede

Perto da borda de absorção de um semicondutor, α(hν)·hν = B (hν − Eg)ⁿ. Linearizando,
(α·hν)^(1/n) é reta em hν e corta o eixo em Eg (gráfico de Tauc). O resultado é o **gap
óptico**, que pode diferir do gap eletrônico (excitons, efeito Burstein–Moss em materiais
muito dopados, confinamento quântico em partículas pequenas, desordem). Tauc foi deduzido
para semicondutores amorfos; para cristalinos é uma aproximação aceita, não uma lei.

## 2. Escolher F (a grandeza ∝ α)

| Amostra | Medida | F a usar | `uvvis.py tauc` |
|---|---|---|---|
| Filme fino transparente em substrato | transmitância/absorbância | α = 2.303·A/d (ou A se d desconhecida) | `--thickness d_cm` opcional |
| Solução/dispersão coloidal estável (pontos quânticos, nanopartículas) | absorbância | A (corrigir espalhamento antes: `--scatter`) | padrão |
| Pó, pastilha, filme opaco/espesso | **reflectância difusa** com esfera integradora vs. BaSO₄/Spectralon | F(R) = (1−R)²/2R (Kubelka-Munk) | `--ykind %R` se o arquivo não disser |

Eg não depende de fatores de escala (espessura, concentração, constante de espalhamento
de KM): só a forma da borda importa. Por isso A serve no lugar de α para Eg - mas não para
comparar coeficientes de absorção entre amostras.

**Nunca** use "absorbância" de pó suspenso em cubeta para band gap: o espalhamento domina.
Kubelka-Munk supõe amostra "infinitamente" espessa (> 2-3 mm de pó), partículas pequenas e
absorção fraca; para materiais muito absorventes, diluir em BaSO₄ (5-10%) melhora.

## 3. Escolher a transição

| Transição | n (Tauc) | expoente em (F·hν) | `--transition` |
|---|---|---|---|
| direta permitida | 1/2 | 2 | `direct` |
| indireta permitida | 2 | 1/2 | `indirect` |
| direta proibida | 3/2 | 2/3 | `direct_forbidden` |
| indireta proibida | 3 | 1/3 | `indirect_forbidden` |

A transição é propriedade do material e deve vir da literatura/cálculo de estrutura de
bandas, **não** do expoente que dá o melhor R² (os dois costumam dar retas razoáveis e Eg
diferentes em 0.1-0.3 eV). Exemplos usuais: TiO₂ anatase indireta (~3.2 eV), rutilo
direta/indireta (~3.0), ZnO direta (~3.3), CdS direta (2.4), WO₃ indireta (2.6-2.8),
α-Fe₂O₃ indireta (~2.1), BiVO₄ monoclínico direta (~2.4), g-C₃N₄ (~2.7, frequentemente
tratada como indireta), Cu₂O direta (~2.1). Quando a literatura diverge, mostre os dois
(`--transition both`) e diga qual adota e por quê.

## 4. Região de ajuste

O script escolhe automaticamente a região: identifica as bordas (picos de d(F·hν)ᵖ/dE),
pega a de maior subida (borda principal), e dentro dela, entre o início da subida e 85% do
platô local, o trecho linear (R² ≥ 0.995) mais íngreme. Isso imita a "tangente na parte
mais íngreme e linear da borda" feita à mão.

Sempre:
- **olhe o gráfico** (`-o tauc.png`): a reta tracejada deve cobrir a parte reta da borda,
  não a curvatura da cauda nem a região saturada/ruidosa acima;
- confira `variação com a janela`: se > ±0.03 eV, a escolha da região importa - reporte
  essa incerteza ou fixe a região com `--fit E1 E2` justificando;
- use `--window E1 E2` para excluir picos excitônicos, outra fase ou a região onde o
  substrato/vidro/solvente absorve (vidro comum corta < ~3.6 eV / 340 nm!).

## 5. Absorção sub-gap e o método da linha de base

Dopantes, defeitos (vacâncias de O), sensibilizadores, carbono residual ou espalhamento
produzem absorção abaixo do gap. Então a reta de Tauc cruzaria y = 0 em energia "errada".
Makuła, Pacia & Macyk (J. Phys. Chem. Lett. 2018, 9, 6814) recomendam intersectar a reta
da borda com a reta da região sub-gap (linha de base), não com o eixo. O script reporta os
dois:

- `Eg (interseção y=0)`: padrão quando a região abaixo da borda é ~0 (amostra pura);
- `Eg (interseção c/ sub-gap)`: preferível quando há absorção sub-gap significativa
  (platô visível acima de zero antes da borda). Em amostras limpas ele fica um pouco acima
  do anterior (a cauda de Urbach puxa a reta sub-gap) - não use nesse caso.

Se o script avisar **"outra borda em X eV"**, há uma segunda absorção (ex.: TiO₂ dopado com
N: ombro no visível ~2.4-2.8 eV). Reporte: "Eg = 3.16 eV (borda do TiO₂) e absorção
adicional com início em ~2.4 eV atribuída a ..." - não chame o ombro de "gap reduzido"
sem evidência de que a borda principal se deslocou.

## 6. Energia de Urbach

Abaixo do gap, α ∝ exp(hν/E_U). ln F vs hν é reta e E_U = 1/inclinação (meV). Mede a
desordem/caudas de estados. Valores típicos: dezenas de meV (cristais bons ~10-30, filmes
nanocristalinos/amorfos 50-300+).

Só tem significado se a cauda é exponencial e medida acima do ruído e do espalhamento. O
script marca como **duvidoso** quando há outra borda (absorção sub-gap não exponencial) ou
E_U > 300 meV. Em pós por DRS e em gaps indiretos é semi-quantitativo.

## 7. Métodos alternativos (conferência)

- **Inflexão (derivada)**: E de máximo dF/dE - em borda direta fica perto de Eg (um pouco
  acima); em indireta não tem significado claro.
- **Onset/corte**: tangente em A vs λ intercepta a linha de base; Eg = 1240/λ_onset. Comum
  em química de materiais, dá valores geralmente um pouco menores que Tauc.
Não misture métodos numa mesma série; diga qual usou.

## 8. Armadilhas

- Saturação (A > 2.5 ou R < ~2%) achata a borda e desloca Eg para baixo.
- %T lido como %R (ou vice-versa) - confira `info`.
- Degrau de troca de lâmpada/detector dentro da borda (ex.: ~350 nm) cria "bordas" falsas.
- Partículas pequenas (< ~10 nm) têm Eg maior por confinamento; partículas grandes/
  agregadas espalham e "esticam" a borda.
- Fases misturadas (anatase/rutilo, P25) dão bordas sobrepostas: o Eg é efetivo.
- Comparar Eg de amostras com F diferentes (KM vs A) não é válido.
- Diferenças < ~0.03-0.05 eV entre amostras geralmente estão dentro da incerteza do
  método - não interprete como efeito do dopante sem replicatas.

## 9. Como reportar

"O gap óptico foi estimado por reflectância difusa (esfera integradora, referência BaSO₄),
convertida pela função de Kubelka-Munk, e gráfico de Tauc para transição indireta
permitida [(F(R)·hν)^1/2 vs hν]. A extrapolação da região linear (3.63-3.96 eV, R² =
0.9997) forneceu Eg = 3.19 ± 0.01 eV."
Inclua a figura de Tauc com a reta e a região de ajuste (o `plot_tauc` já faz).
