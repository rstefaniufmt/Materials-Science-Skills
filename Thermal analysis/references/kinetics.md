# Cinética não isotérmica (recomendações ICTAC)

Base: Vyazovkin et al., *ICTAC Kinetics Committee recommendations for performing kinetic
computations on thermal analysis data*, Thermochim. Acta 520 (2011) 1-19; e
*... for collecting experimental thermal analysis data for kinetic computations*,
Thermochim. Acta 590 (2014) 1-23.

## Sumário
1. Dados necessários
2. Conversão α
3. Métodos implementados
4. Lendo o resultado
5. O que não fazer
6. Como reportar

## 1. Dados necessários

- **≥ 3 taxas de aquecimento, idealmente 4-5**, cobrindo um fator ≥ 5-10 (ex.: 2, 5, 10, 20 K/min).
  Com 3 taxas, o IC95% de cada Ea é largo; o script avisa.
- Mesma massa (pequena, 2-10 mg), mesmo cadinho, mesmo gás e vazão. Verifique que a massa
  não muda o resultado (gradientes térmicos e de difusão em β alto).
- A **taxa medida** é usada (dT/dt na janela), não a nominal.
- Um **processo** por análise: escolha a janela (`--window T1:T2`) que contém uma etapa (ou um
  pico DSC) com patamares dos dois lados em todas as taxas. Etapas sobrepostas dão Ea(α)
  variável, que é real mas deve ser interpretada como processo complexo.

## 2. Conversão α

- TG: α = (m₀ − m)/(m₀ − m∞), com m₀ e m∞ medianas nas bordas da janela.
- DSC: α = área parcial / área total do pico (linha de base linear na janela; `--signal hf`).
- dα/dt por Savitzky-Golay no tempo (Friedman). Ruído em dα/dt afeta Friedman nas pontas
  (α < 0,1 e > 0,9) - por isso o resumo usa 0,1 ≤ α ≤ 0,9.

## 3. Métodos implementados

| Método | Relação linearizada (em cada α) | Observação |
|---|---|---|
| Friedman (diferencial) | ln(β dα/dT) = ln(dα/dt) = ln[A f(α)] − Ea/RT | sem aproximações da integral; sensível a ruído e linha de base |
| Ozawa-Flynn-Wall (OFW) | ln β = const − 1,052 Ea/RT | aproximação de Doyle; erro de alguns % quando Ea/RT < 20 |
| Kissinger-Akahira-Sunose (KAS) | ln(β/T²) = const − Ea/RT | mais exato que OFW |
| Starink | ln(β/T¹·⁹²) = const − 1,0008 Ea/RT | ainda mais exato |
| Vyazovkin (não linear) | minimiza Σᵢ Σⱼ≠ᵢ I(Ea,Tα,ᵢ)βⱼ / I(Ea,Tα,ⱼ)βᵢ | integral p(x) de Senum-Yang; o mais preciso dos integrais |
| Kissinger | ln(β/Tp²) = ln(AR/Ea) − Ea/(R Tp) | um Ea só; rigoroso só para 1ª ordem; dá ln A (1ª ordem) |

Os integrais (OFW, KAS, Starink, Vyazovkin "integral") assumem Ea constante de 0 até α;
quando Ea varia muito, Friedman (ou Vyazovkin modificado/incremental) é mais fiel. A
diferença Friedman vs integrais é diagnóstico disso.

## 4. Lendo o resultado

- `Ea(α)` quase constante (variação < 10-20% em 0,1-0,9): processo de etapa única efetiva -
  um Ea médio pode ser reportado, e faz sentido buscar o modelo de reação (master plots,
  não implementados) e A.
- `Ea(α)` crescente: reações paralelas/competitivas ou mudança de mecanismo; decrescente:
  reação reversível, mudança de controle cinético para difusão, autocatalítico. Reporte a
  curva, não a média (o script avisa quando a variação passa de 20%).
- R² < 0,98 em boa parte de α: pontos de alguma taxa fora (massa diferente, linha de base),
  janela errada ou sobreposição.
- O IC95% vem da regressão em cada α com n = nº de taxas (t de Student com n − 2 g.l.); com
  3 taxas ele é grande - isso é a incerteza real.
- Kissinger: compare com a média isoconversional; se α no pico varia entre as taxas, o
  script avisa.

## 5. O que não fazer

- **Ea de uma única taxa** (Coats-Redfern, Horowitz-Metzger, Broido, ajuste de modelo a uma
  curva): a ICTAC desaconselha - qualquer modelo ajusta bem com Ea/A errados (efeito de
  compensação).
- Usar o ln A do Kissinger como "o" fator pré-exponencial de um processo complexo.
- Escolher o modelo f(α) pelo maior R² do ajuste.
- Extrapolar vida útil (lifetime) muito além da faixa medida sem validar com isotermas.
- Comparar Ea de amostras medidas em condições diferentes (massa, gás).
- Misturar OFW e KAS como "confirmação independente": usam os mesmos T_α - concordarem é
  esperado.

## 6. Como reportar

Taxas medidas, massa, gás, janela, definição de α; gráfico Ea(α) com barras de erro (o
`plot_kinetics` faz) e tabela por α; Ea médio ± desvio (ou faixa) em 0,1-0,9 por método;
Kissinger como comparação. Ex.: "Ea = 180 ± 3 kJ/mol (Friedman, 0,1 ≤ α ≤ 0,9; β = 2,
5, 10 e 20 K/min, N₂); Ea praticamente constante com α, sugerindo processo de etapa única".
