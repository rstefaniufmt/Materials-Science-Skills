# Interpretação de espectros UV-Vis

## Sumário
1. Controle de qualidade antes de interpretar
2. Transições eletrônicas e intensidades
3. Cromóforos orgânicos, conjugação e efeitos de solvente/pH
4. Complexos de metais de transição
5. Agregação, interações e pontos isosbésticos
6. Nanopartículas metálicas (SPR)
7. Semicondutores e sólidos
8. Como redigir

## 1. Controle de qualidade

- **Faixa de A**: confiável ~0.05-2 (até ~3 em equipamentos de duplo feixe bons). Bandas
  com topo achatado/ruidoso = saturação.
- **Corte do solvente e da cubeta** (tabela em reference_data.md): abaixo deles o
  "espectro" é ruído; não atribua bandas ali.
- **Degraus** em ~320-350 nm (troca de lâmpada D₂/W) e ~800-870 nm (troca PMT/InGaAs) são
  instrumentais (`info` detecta; `splice_correct` corrige).
- **Linha de base**: A ≠ 0 onde nada absorve indica branco errado, cubeta diferente, ou
  espalhamento (subida suave para o UV ∝ λ⁻ⁿ). Espalhamento não é banda.
- **Picos muito finos** (< 2-3 nm) em espectros de solução raramente são da amostra: lâmpada
  (linha do D₂ em 486/656 nm), ruído, ou padrões de lantanídeos/hólmio.
- **Fluorescência** da amostra pode reduzir a absorbância aparente.
- Diferenças de λmax < ~1-2 nm geralmente não são significativas (verifique a calibração
  com hólmio: `wlcheck`).

## 2. Transições e intensidades (ε, L mol⁻¹ cm⁻¹)

| Transição | Região típica | ε | Observação |
|---|---|---|---|
| σ→σ* | < 200 nm | - | fora do UV comum |
| n→σ* | 150-250 | 100-3 000 | álcoois, aminas, haletos |
| π→π* (isolado) | < 200 | ~10⁴ | alcenos |
| π→π* (conjugado/aromático) | 200-500+ | 10³-10⁵ | desloca para o vermelho com conjugação |
| n→π* | 270-350 (carbonila), maior para C=N, N=N | 10-100 | fraca, desloca para o azul em solvente polar/prótico |
| d-d (octaédrico) | visível/NIR | 1-100 | proibida por Laporte; tetraédrico 100-1 000 |
| transferência de carga (LMCT/MLCT/IVCT) | UV-visível | 10³-5×10⁴ | intensa; cores fortes (MnO₄⁻, [Fe(phen)₃]²⁺) |
| SPR (Au, Ag, Cu) | 380-1100 | muito alta (10⁸-10¹⁰ por partícula) | depende de tamanho/forma/meio |

Intensidade (ε) distingue o tipo de transição melhor que a posição.

## 3. Cromóforos orgânicos

- **Conjugação**: cada dupla conjugada adicional desloca λmax ~30-40 nm para o vermelho
  (batocrômico) e aumenta ε (hipercrômico). Regras de Woodward-Fieser estimam λmax de
  dienos e enonas; para sistemas extensos use comparação com análogos/literatura.
- **Benzeno**: 184, 204 e banda fina estruturada ~255 nm (B); substituintes doadores/
  aceptores deslocam e intensificam (anilina 230/280; fenol 210/270; fenolato 235/287;
  nitrobenzeno ~268).
- **Solvatocromismo**: π→π* desloca para o vermelho com polaridade; n→π* para o azul
  (solvente prótico estabiliza o par n). Estrutura vibracional some em solventes polares.
- **pH**: fenóis, aminas, indicadores e muitos corantes mudam de espectro com a protonação
  (ex.: 4-nitrofenol 317 nm → 4-nitrofenolato 400 nm; alaranjado de metila 505 (ácido) →
  464 nm (básico)). Controle o pH em quantificação.
- **Corantes azo/xanteno/fenotiazina**: banda visível intensa (ε 10⁴-10⁵) + bandas UV dos
  anéis; perda da banda visível = quebra do cromóforo (ou redução à forma leuco).

## 4. Complexos de metais de transição

- Bandas d-d fracas e largas; número e posição dependem da configuração dⁿ, geometria e
  campo ligante (série espectroquímica: I⁻ < Br⁻ < Cl⁻ < F⁻ < OH⁻ < H₂O < NH₃ < en < CN⁻ < CO).
  Ligante de campo mais forte → bandas para energia maior (azul).
- 10Dq (Δo) para d¹ (Ti³⁺) e d⁹ (Cu²⁺) ≈ energia da banda: Δ (cm⁻¹) = 10⁷/λ(nm).
- Jahn-Teller (d⁹ Cu²⁺, d⁴ alto spin, d¹ excitado) alarga/desdobra bandas.
- Tetraédricos têm bandas ~10-100× mais intensas que octaédricos (ex.: [CoCl₄]²⁻ azul
  intenso vs [Co(H₂O)₆]²⁺ rosa pálido).
- Bandas intensas (ε > 1 000) no UV-visível em complexos = transferência de carga, não d-d.
- Use Tanabe-Sugano para atribuir e estimar Δ e B (parâmetro de Racah) quando há ≥ 2 bandas
  d-d.

## 5. Agregação e isosbésticos

- **Agregados H** (empilhamento face a face): nova banda para o azul (MB: 664 → ~610 nm do
  dímero); **agregados J** (cabeça-cauda): banda estreita para o vermelho.
- Desvio da Lei de Beer com a concentração e mudança da razão entre bandas indicam
  agregação/equilíbrio.
- **Ponto isosbéstico**: λ onde todos os espectros de uma série se cruzam → interconversão
  entre DUAS espécies absorventes com estequiometria fixa (ácido/base, livre/complexado,
  reagente/produto). Perda do ponto isosbéstico ao longo da série → terceira espécie,
  intermediário, precipitação/espalhamento, ou deriva instrumental.
- Titulação espectrofotométrica: A(λ) vs pH ajustado por Henderson-Hasselbalch dá o pKa;
  razão molar/Job dão estequiometria de complexos.

## 6. Nanopartículas metálicas (SPR)

- Au esférico 10-20 nm em água: ~520 nm; desloca para o vermelho e alarga com o tamanho
  (~530 nm em ~50 nm, ~550+ em ~80-100 nm) e com o índice de refração do meio/
  recobrimento. Ag esférica ~400 nm; Cu ~570-590 nm (sensível a oxidação).
- Bastões (nanorods): banda transversal ~520 nm (Au) + longitudinal 600-1100 nm,
  proporcional à razão de aspecto.
- **Agregação**: banda alarga, nova absorção 600-800 nm (Au: vermelho → roxo/azul); base de
  sensores colorimétricos.
- Largura (FWHM) estreita e simétrica sugere distribuição de tamanho estreita; cauda no
  vermelho sugere agregados/polidispersão. A razão A(SPR)/A(450) é usada para estimar
  tamanho e concentração de AuNP (Haiss et al., Anal. Chem. 2007) - valores dependem do
  modelo; cite a fonte.
- Interbanda do Au abaixo de ~500 nm (A400 é usada para estimar concentração de Au(0)).
- Espalhamento contribui muito para partículas > ~40 nm: considere `--scatter` antes de
  medir λmax/FWHM.

## 7. Semicondutores e sólidos

- Borda de absorção (não um pico): use `tauc` e `references/bandgap.md`.
- Picos excitônicos logo acima da borda (pontos quânticos: posição dá o tamanho via curvas
  de calibração da literatura, ex. CdSe de Yu et al., Chem. Mater. 2003).
- Absorção no visível em óxidos brancos dopados (ombro/cauda) = estados no gap.
- Em DRS, compare sempre em F(R) (não em %R) e mesma referência branca.

## 8. Como redigir

- Comece pelos λmax principais com ε (se conhecido) e atribua o tipo de transição pela
  intensidade e pelo contexto químico.
- Para séries, descreva tendências com números (deslocamento de X nm, razão entre bandas,
  ponto isosbéstico em Y nm) e o que elas implicam (equilíbrio de duas espécies, agregação,
  crescimento de partícula).
- Deixe claro o que é inferência: UV-Vis raramente identifica um composto sozinho; combine
  com padrão de referência, HPLC, IV, RMN, DRX conforme o caso.
