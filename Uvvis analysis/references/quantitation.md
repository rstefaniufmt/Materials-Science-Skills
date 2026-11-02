# Quantificação e cinética por UV-Vis

## Sumário
1. Beer-Lambert e quando ele falha
2. Desenho da curva de calibração
3. O que o `calib` calcula e como ler
4. Amostras: predição, diluição, incerteza
5. Situações especiais (matriz, adição de padrão, dois componentes)
6. Cinética (degradação, fotocatálise, reações)
7. Como reportar

## 1. Beer-Lambert

A = ε·b·c (A adimensional, b em cm, c em mol/L → ε em L mol⁻¹ cm⁻¹). Vale para luz
monocromática, soluções diluídas e sem interações. Desvios:
- **Instrumentais**: luz espúria (curvatura para baixo em A alto - mantenha A ≲ 1.5,
  nunca > ~2.5-3), largura de banda espectral grande frente à banda (medir no máximo,
  onde a banda é plana, minimiza), absorbância do branco mal subtraída (intercepto ≠ 0).
- **Químicos**: equilíbrios dependentes da concentração - dimerização/agregação (azul de
  metileno forma dímero com banda ~610 nm em concentrações altas), pH (indicadores,
  fenóis), complexação; fluorescência da amostra; reação/fotodegradação durante a medida.
- **Físicos**: turbidez/espalhamento (filtre ou corrija com `--scatter`), bolhas, cubeta
  suja ou mal posicionada.

## 2. Desenho da curva

- ≥ 5-6 níveis (ideal 6-8), espaçados uniformemente cobrindo a faixa esperada das
  amostras (amostras perto do centro da curva têm menor incerteza).
- Mesmo solvente/matriz e mesmo branco das amostras; mesma cubeta (ou cubetas casadas).
- Faixa de A ~0.1-1.0 dá a melhor precisão fotométrica.
- Réplicas independentes (preparos separados) estimam melhor a incerteza que leituras
  repetidas da mesma solução.
- λ de medida: máximo da banda (sensibilidade máxima, menos sensível a erro de λ), longe
  de interferentes; mantenha o mesmo λ para tudo.

## 3. Leitura do `calib`

| Saída | Significado | O que fazer |
|---|---|---|
| A = b·c + a | inclinação = sensibilidade | b/caminho = ε (se c em mol/L) |
| R², R²aj | ajuste | R² ≥ 0.995 é usual, mas não prova linearidade |
| s(y/x) | desvio residual | base de LOD/LOQ e da incerteza |
| LOD = 3.3·s(y/x)/b, LOQ = 10·s(y/x)/b | ICH Q2(R2), pela curva | Alternativa: 3.3σ/10σ de ≥ 10 brancos |
| teste de Mandel (p < 0.05) | modelo quadrático é significativamente melhor | reduza a faixa |
| intercepto p < 0.05 | a ≠ 0 | branco/interferente/desvio; investigue antes de forçar pela origem |
| resíduos com tendência | curvatura | reduza a faixa ou retire o ponto saturado |

Forçar pela origem (`--origin`) só se o intercepto não for significativo e houver razão
física; não é um jeito de "melhorar" a curva.

## 4. Amostras

`--unknown arquivo(s)` ou `--unknown-A valores`; `--dilution` multiplica o resultado.
O script dá c ± s e IC95% (fórmula da incerteza de predição, considera n padrões e m
réplicas da amostra) e marca fora da faixa / abaixo do LOQ. Meça amostras em triplicata
quando possível (m = 3 reduz a incerteza).

## 5. Situações especiais

- **Efeito de matriz**: compare as inclinações de curva em solvente e em matriz
  fortificada; se diferirem, use curva na matriz ou **adição de padrão** (curva A vs
  c adicionada na própria amostra; c_amostra = |intercepto em x| = a/b).
- **Mistura de dois componentes** com espectros sobrepostos: meça em dois λ e resolva
  A₁ = ε₁ₓ c_x + ε₁ᵧ c_y, A₂ = ε₂ₓ c_x + ε₂ᵧ c_y (ε de cada componente puro em cada λ).
  Escolha λ onde as diferenças de ε são grandes.
- **Turbidez**: correção de três pontos (A_λ − média dos ombros) ou `--scatter`.
- **Proteínas/ácidos nucleicos**: A280/A260 com ε tabelados; razão A260/A280 indica pureza
  (ver reference_data.md).

## 6. Cinética

Com Beer-Lambert válido, C/C₀ = (A − A∞)/(A₀ − A∞) (A∞ = 0 se o produto não absorve
nesse λ; senão informe A∞).

| Modelo | Forma linear | k | t½ |
|---|---|---|---|
| ordem zero | C₀ − C = k₀ t | conc./tempo | C₀/2k₀ |
| pseudo-1ª ordem (L-H a baixa conc.) | ln(C₀/C) = k t | tempo⁻¹ | ln2/k |
| pseudo-2ª ordem | 1/C − 1/C₀ = k t | conc.⁻¹ tempo⁻¹ | 1/(k C₀) |

`--model all` compara os três por R²; prefira o modelo coerente com o mecanismo e com
resíduos sem tendência, não só o maior R². Langmuir-Hinshelwood: a baixa concentração vira
pseudo-1ª ordem; k aparente depende de C₀, massa de catalisador, intensidade de luz e
geometria - só compare k obtidos nas mesmas condições.

**Fotocatálise - checklist**:
- fase escura até o equilíbrio de adsorção (tempos negativos); C₀ no fim do escuro;
- controles: fotólise (luz sem catalisador) e adsorção (catalisador sem luz);
- centrifugar/filtrar antes de medir (partículas espalham: aumento aparente de A);
- verificar se o espectro inteiro cai junto (degradação) ou se λmax desloca (ex.: rodamina
  B desetila e desloca 554 → ~498 nm; azul de metileno pode ser reduzido à forma leuco,
  incolor, reversível) - então "descoloração" ≠ degradação;
- mineralização exige TOC/COD; UV-Vis não mede.
- Para comparar catalisadores, reporte também k normalizado por massa/área se fizer sentido.

**Redução de 4-nitrofenol (teste catalítico de nanopartículas)**: com NaBH₄, monitore o
4-nitrofenolato em 400 nm; aparece 4-aminofenol ~300 nm; ponto isosbéstico ~315 nm indica
conversão limpa. Há tempo de indução (descarte da regressão).

## 7. Como reportar

"Curva analítica em 664 nm (seis níveis, 0.5-5.0 mg/L; A = 0.2384c − 0.0015; R² =
0.99999; s(y/x) = 0.0015; LOD = 0.021 mg/L e LOQ = 0.063 mg/L, ICH Q2). A amostra (diluição
1:10) apresentou 24.9 ± 0.1 mg/L (IC95%)."
"A descoloração seguiu cinética de pseudo-primeira ordem (k = 0.0179 min⁻¹, R² = 0.999,
t½ = 39 min), após 30 min de adsorção no escuro (15% de remoção)."
