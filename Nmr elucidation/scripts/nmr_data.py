"""
nmr_data.py - tabelas de referência: solventes deuterados, impurezas comuns, faixas de deslocamento
químico e constantes de acoplamento típicas, fórmula molecular / IDH (DBE)

Fontes: H. E. Gottlieb, V. Kotlyar, A. Nudelman, J. Org. Chem. 62 (1997) 7512; G. R. Fulmer et al.,
Organometallics 29 (2010) 2176 (valores a ~298 K; variam ±0.02-0.1 ppm com T, concentração, pH).
Faixas de deslocamento e J: Pretsch, Bühlmann, Badertscher, "Structure Determination of Organic
Compounds"; Silverstein, Webster & Kiemle.
"""
from __future__ import annotations

import re

# sinal residual do solvente: (¹H [ppm, multiplicidade], ¹³C [ppm...], água ¹H)
SOLVENTS = {
    "cdcl3": dict(h=[(7.26, "s")], c=[77.16], water=1.56, name="CDCl3"),
    "dmso-d6": dict(h=[(2.50, "quint")], c=[39.52], water=3.33, name="DMSO-d6"),
    "cd3od": dict(h=[(3.31, "quint")], c=[49.00], water=4.87, name="CD3OD"),
    "d2o": dict(h=[(4.79, "s")], c=[], water=4.79, name="D2O"),
    "acetone-d6": dict(h=[(2.05, "quint")], c=[29.84, 206.26], water=2.84, name="acetona-d6"),
    "c6d6": dict(h=[(7.16, "s")], c=[128.06], water=0.40, name="C6D6"),
    "cd3cn": dict(h=[(1.94, "quint")], c=[1.32, 118.26], water=2.13, name="CD3CN"),
    "cd2cl2": dict(h=[(5.32, "t")], c=[53.84], water=1.52, name="CD2Cl2"),
    "thf-d8": dict(h=[(1.72, "m"), (3.58, "m")], c=[25.31, 67.21], water=2.46, name="THF-d8"),
    "toluene-d8": dict(h=[(2.08, "quint"), (6.97, "m"), (7.01, "m"), (7.09, "m")],
                       c=[20.43, 125.13, 127.96, 128.87, 137.48], water=0.43, name="tolueno-d8"),
    "pyridine-d5": dict(h=[(7.22, "m"), (7.58, "m"), (8.74, "m")], c=[123.87, 135.91, 150.35], water=4.96,
                        name="piridina-d5"),
    "dmf-d7": dict(h=[(2.75, "quint"), (2.92, "quint"), (8.03, "s")], c=[29.76, 34.89, 163.15], water=3.51,
                   name="DMF-d7"),
    "tfa-d": dict(h=[(11.50, "s")], c=[116.6, 164.2], water=11.5, name="TFA-d"),
}

# padrão do sinal residual (acoplamento com ²H): (δ, J_Hz, intensidades relativas)
SOLVENT_PATTERNS = {
    "1H": {
        "cdcl3": [(7.26, 0.0, [1])], "c6d6": [(7.16, 0.0, [1])], "d2o": [(4.79, 0.0, [1])],
        "dmso-d6": [(2.50, 1.9, [1, 2, 3, 2, 1])], "cd3od": [(3.31, 1.7, [1, 2, 3, 2, 1])],
        "acetone-d6": [(2.05, 2.2, [1, 2, 3, 2, 1])], "cd3cn": [(1.94, 2.5, [1, 2, 3, 2, 1])],
        "cd2cl2": [(5.32, 1.0, [1, 1, 1])],
    },
    "13C": {
        "cdcl3": [(77.16, 32.0, [1, 1, 1])], "cd2cl2": [(53.84, 27.2, [1, 2, 3, 2, 1])],
        "dmso-d6": [(39.52, 21.0, [1, 3, 6, 7, 6, 3, 1])], "cd3od": [(49.00, 21.4, [1, 3, 6, 7, 6, 3, 1])],
        "acetone-d6": [(29.84, 19.4, [1, 3, 6, 7, 6, 3, 1]), (206.26, 0.0, [1])],
        "c6d6": [(128.06, 24.3, [1, 1, 1])], "cd3cn": [(1.32, 21.0, [1, 3, 6, 7, 6, 3, 1]), (118.26, 0.0, [1])],
        "thf-d8": [(67.21, 22.0, [1, 2, 3, 2, 1]), (25.31, 20.5, [1, 2, 3, 2, 1])],
    },
}

# impurezas comuns: {solvente: {composto: ([δH...], [δC...])}}
IMPURITIES = {
    "cdcl3": {
        "ácido acético": ([2.10, 11.4], [20.81, 175.99]), "acetona": ([2.17], [30.92, 207.07]),
        "acetonitrila": ([2.10], [1.89, 116.43]), "benzeno": ([7.36], [128.37]),
        "terc-butanol": ([1.28], [31.25, 69.15]), "MTBE": ([1.19, 3.22], [27.0, 49.45, 72.87]),
        "BHT": ([6.98, 5.01, 2.27, 1.43], [21.2, 30.3, 34.25, 125.55, 128.27, 135.87, 151.55]),
        "clorofórmio (CHCl3)": ([7.26], [77.36]), "ciclo-hexano": ([1.43], [26.94]),
        "1,2-dicloroetano": ([3.73], [43.5]), "diclorometano": ([5.30], [53.52]),
        "éter etílico": ([1.21, 3.48], [15.20, 65.91]), "DMF": ([8.02, 2.96, 2.88], [162.62, 36.50, 31.39]),
        "DMSO": ([2.62], [40.76]), "dioxano": ([3.71], [67.14]),
        "etanol": ([1.25, 3.72], [18.41, 58.28]), "acetato de etila": ([2.05, 4.12, 1.26], [21.04, 60.49, 171.36, 14.19]),
        "etilenoglicol": ([3.76], [63.79]), "graxa (hidrocarboneto)": ([0.86, 1.26], [29.76]),
        "n-hexano": ([0.88, 1.26], [14.14, 22.70, 31.64]), "metanol": ([3.49], [50.41]),
        "nitrometano": ([4.33], [62.50]), "n-pentano": ([0.88, 1.27], [14.08, 22.38, 34.16]),
        "2-propanol": ([1.22, 4.04], [25.14, 64.50]), "piridina": ([8.62, 7.29, 7.68], [149.90, 123.75, 135.96]),
        "graxa de silicone": ([0.07], [1.19]), "THF": ([1.85, 3.76], [25.62, 67.97]),
        "tolueno": ([2.36, 7.17, 7.25], [21.46, 125.33, 128.26, 129.07, 137.89]),
        "trietilamina": ([1.03, 2.53], [11.61, 46.25]), "água": ([1.56], []), "TMS": ([0.00], [0.00]),
    },
    "dmso-d6": {
        "ácido acético": ([1.91], [20.95, 171.93]), "acetona": ([2.09], [30.56, 206.31]),
        "acetonitrila": ([2.07], [1.03, 117.91]), "benzeno": ([7.37], [128.30]),
        "terc-butanol": ([1.11, 4.19], [30.38, 66.88]), "BHT": ([6.87, 6.65, 2.18, 1.36], []),
        "clorofórmio": ([8.32], [79.16]), "ciclo-hexano": ([1.40], [26.33]), "diclorometano": ([5.76], [54.84]),
        "éter etílico": ([1.09, 3.38], [15.12, 62.05]), "DMF": ([7.95, 2.89, 2.73], [162.29, 35.73, 30.73]),
        "dioxano": ([3.57], [66.36]), "etanol": ([1.06, 3.44, 4.35], [18.51, 56.07]),
        "acetato de etila": ([1.99, 4.03, 1.17], [20.68, 59.74, 170.31, 14.40]),
        "n-hexano": ([0.86, 1.25], [13.88, 22.05, 30.95]), "metanol": ([3.16, 4.01], [48.59]),
        "2-propanol": ([1.04, 3.78, 4.34], [25.43, 64.92]), "piridina": ([8.58, 7.39, 7.79], [149.58, 123.84, 136.05]),
        "THF": ([1.76, 3.60], [25.14, 67.03]), "tolueno": ([2.30, 7.18, 7.25], [20.99, 125.29, 128.18, 128.88, 137.35]),
        "trietilamina": ([0.93, 2.43], [11.74, 45.74]), "água": ([3.33], []), "graxa": ([0.85, 1.24], [29.0]),
    },
    "cd3od": {
        "ácido acético": ([1.99], [20.56, 175.11]), "acetona": ([2.15], [30.67, 209.67]),
        "acetonitrila": ([2.03], [0.85, 118.06]), "benzeno": ([7.33], [129.34]), "clorofórmio": ([7.90], [79.44]),
        "diclorometano": ([5.49], [54.78]), "éter etílico": ([1.18, 3.49], [15.46, 66.88]),
        "DMF": ([7.97, 2.99, 2.86], [164.73, 36.89, 31.61]), "DMSO": ([2.65], [40.45]), "dioxano": ([3.66], [68.11]),
        "etanol": ([1.19, 3.60], [18.40, 58.26]), "acetato de etila": ([2.01, 4.09, 1.24], [20.88, 61.50, 172.89, 14.49]),
        "n-hexano": ([0.90, 1.29], [14.45, 23.68, 32.73]), "metanol (CH3OH)": ([3.34], [49.86]),
        "2-propanol": ([1.15, 3.92], [25.27, 64.71]), "piridina": ([8.53, 7.44, 7.85], [150.07, 125.53, 138.35]),
        "THF": ([1.87, 3.71], [26.48, 68.83]), "tolueno": ([2.32, 7.16], [21.50]), "trietilamina": ([1.05, 2.58], [9.73, 47.02]),
        "água (HDO)": ([4.87], []),
    },
    "d2o": {
        "ácido acético": ([2.08], [21.03, 177.21]), "acetona": ([2.22], [30.89, 215.94]), "acetonitrila": ([2.06], [1.47, 119.68]),
        "terc-butanol": ([1.24], [30.29, 70.36]), "DMF": ([7.92, 3.01, 2.85], [165.53, 37.54, 32.03]), "DMSO": ([2.71], [39.39]),
        "dioxano": ([3.75], [67.19]), "etanol": ([1.17, 3.65], [17.47, 58.05]), "éter etílico": ([1.17, 3.56], [14.77, 66.42]),
        "metanol": ([3.34], [49.50]), "2-propanol": ([1.17, 4.02], [24.38, 64.88]), "THF": ([1.88, 3.74], [25.67, 68.68]),
        "trietilamina": ([0.99, 2.57], [9.07, 47.19]), "HDO": ([4.79], []),
    },
    "acetone-d6": {
        "ácido acético": ([1.96], [20.51, 172.31]), "acetonitrila": ([2.05], [1.12, 117.60]), "benzeno": ([7.36], [129.15]),
        "clorofórmio": ([8.02], [79.19]), "diclorometano": ([5.63], [54.95]), "éter etílico": ([1.11, 3.41], [15.78, 66.12]),
        "DMF": ([7.96, 2.94, 2.78], [162.79, 36.15, 31.03]), "DMSO": ([2.52], [41.23]),
        "etanol": ([1.12, 3.57], [18.89, 57.72]), "acetato de etila": ([1.97, 4.05, 1.20], [20.83, 60.56, 170.96, 14.50]),
        "n-hexano": ([0.88, 1.28], [14.29, 23.28, 32.30]), "metanol": ([3.31], [49.77]), "THF": ([1.79, 3.63], [26.15, 68.07]),
        "tolueno": ([2.32, 7.1, 7.2], [21.46]), "água": ([2.84], []),
    },
    "c6d6": {
        "acetona": ([1.55], [30.14, 204.43]), "acetonitrila": ([0.58], [0.63, 116.02]), "clorofórmio": ([6.15], [77.79]),
        "diclorometano": ([4.27], [53.46]), "éter etílico": ([1.11, 3.26], [15.46, 65.94]), "DMSO": ([1.68], [40.03]),
        "acetato de etila": ([1.65, 3.89, 0.92], [20.56, 60.21, 170.44, 14.19]), "n-hexano": ([0.89, 1.24], [14.32, 23.04, 31.96]),
        "metanol": ([3.07], [49.97]), "THF": ([1.40, 3.57], [25.72, 67.80]), "tolueno": ([2.11, 7.02, 7.13], [21.10]),
        "água": ([0.40], []), "graxa de silicone": ([0.29], [1.38]),
    },
    "cd3cn": {
        "acetona": ([2.08], [30.91, 207.43]), "clorofórmio": ([7.58], [79.21]), "diclorometano": ([5.44], [55.32]),
        "éter etílico": ([1.12, 3.42], [15.63, 66.32]), "DMF": ([7.92, 2.89, 2.77], [163.31, 36.57, 31.32]),
        "DMSO": ([2.50], [41.31]), "acetato de etila": ([1.97, 4.06, 1.20], [21.16, 60.98, 171.68, 14.54]),
        "metanol": ([3.28], [50.41]), "THF": ([1.80, 3.64], [26.27, 68.33]), "tolueno": ([2.33], [21.50]), "água": ([2.13], []),
    },
    "cd2cl2": {
        "acetona": ([2.12], [31.00, 206.78]), "acetonitrila": ([1.97], [2.03, 116.92]), "clorofórmio": ([7.32], [77.99]),
        "éter etílico": ([1.15, 3.43], [15.44, 66.11]), "DMSO": ([2.59], [40.57]),
        "acetato de etila": ([2.00, 4.08, 1.24], [21.15, 60.63, 171.03, 14.38]), "n-hexano": ([0.89, 1.27], [14.28, 23.07, 31.96]),
        "metanol": ([3.42], [50.45]), "THF": ([1.82, 3.69], [25.98, 68.16]), "tolueno": ([2.34], [21.53]), "água": ([1.52], []),
        "graxa": ([0.86, 1.27], [30.10]),
    },
}

# faixas típicas (lo, hi, descrição) - várias podem se aplicar a um mesmo sinal
RANGES_1H = [
    (-1.0, 0.4, "Si–CH₃, ciclopropano, M–H; TMS (0,00)"),
    (0.7, 1.1, "CH₃ alifático terminal (t ou d)"),
    (0.9, 1.5, "CH₃ geminal/terc-butila (s), CH₂ de cadeia"),
    (1.2, 2.1, "CH/CH₂ alifático, CH₂ de anel; CH₃ alílico (1,6-1,8)"),
    (1.9, 2.7, "H α a C=O (CH₃CO 2,0-2,2; CH₂CO 2,2-2,6); Ar–CH₃ (2,2-2,4); H alílico/benzílico"),
    (2.2, 3.1, "N–CH (aminas), S–CH, N–CH₃ de amida (2,8-3,1), H acetilênico (1,8-3,1)"),
    (3.2, 4.5, "O–CH (álcool, éter, éster O–CH 3,6-4,3), OCH₃ (3,3-3,9), CH–N⁺, CH–halogênio, CH₂ de açúcares"),
    (4.3, 5.6, "H anomérico de açúcares; O–CH–O (acetal); CH–O de éster secundário; =CH₂ terminal (4,6-5,2)"),
    (4.6, 6.6, "H vinílico (alceno); enol éter (=CH–O 6,2-6,6)"),
    (6.0, 7.0, "H aromático orto/para a doador (OH, OR, NH₂); H β de enona/acrilato (5,8-6,5)"),
    (6.8, 8.3, "H aromático (benzeno 7,27; orto a C=O/NO₂ 7,8-8,3); heteroaromático (furano, tiofeno, pirrol 6,2-7,5)"),
    (7.0, 9.3, "piridina/heteroaromático deficiente (H-2 8,5-9,0); H α de enona (6,8-7,7); formato HCO₂R (8,0-8,1)"),
    (5.0, 8.5, "N–H de amida (largo, troca com D₂O)"),
    (9.3, 10.2, "aldeído CHO (alifático 9,5-9,8, d/t J 1-3 Hz; aromático 9,8-10,1, s)"),
    (10.0, 13.5, "COOH (largo); OH fenólico quelado (11-14); N–H de indol/pirrol/imida"),
    (0.5, 5.5, "OH/NH trocáveis (largos, posição varia com concentração e solvente)"),
]
RANGES_13C = [
    (-5, 8, "Si–CH₃, ciclopropano, CH₃–M"),
    (8, 25, "CH₃ alifático (CH₃–C)"),
    (15, 45, "CH₂/CH alifático, C quaternário alifático (30-45)"),
    (20, 32, "CH₃ de acetila/CH₃CO (20-30), Ar–CH₃ (19-22)"),
    (28, 50, "N–CH₃ (30-45), CH₂–N, CH₂–S, CH₂ α a C=O (30-45)"),
    (50, 58, "OCH₃ de éster (51-53), OCH₃ aromático (55-56), N–CH (45-60)"),
    (58, 90, "C–O sp³ (álcool, éter, éster: CH₂O 60-70; CHO 65-80; C quaternário–O 70-90); C de açúcar (60-85)"),
    (65, 92, "C≡C alquino (70-90)"),
    (90, 112, "C anomérico (O–C–O, 90-110); C=C β de enol éter/enamina; =CH₂ terminal (105-120)"),
    (110, 125, "C≡N nitrila (115-125); C aromático orto a OR/NR₂ (110-118)"),
    (110, 150, "C=C alceno e C aromático/heteroaromático (CH 115-135; C quaternário 125-150)"),
    (140, 168, "C aromático–O/N (fenol, éter arílico 145-165); C=N imina/heteroaromático (140-165)"),
    (150, 160, "carbonato, ureia, carbamato (C=O 150-158)"),
    (160, 185, "C=O de ácido (170-185), éster/lactona (165-175), amida/lactama (165-175), anidrido, conjugado (160-170)"),
    (180, 192, "quinona (180-190), C=O de enona/conjugada (185-200)"),
    (190, 206, "aldeído CHO (190-205)"),
    (195, 222, "cetona (195-220; cíclica de 5 anéis 210-220)"),
]
J_TYPICAL = [
    ("²J geminal sp³ (H–C–H)", "10-15 Hz (sinal negativo; só visível se H diastereotópicos)"),
    ("³J vicinal sp³ com rotação livre", "6-8 Hz (≈ 7 Hz em CH₃CH₂)"),
    ("³J axial-axial (ciclo-hexano, piranose)", "8-13 Hz; ax-eq e eq-eq 2-5 Hz"),
    ("³J alceno trans", "12-18 Hz"),
    ("³J alceno cis", "6-12 Hz"),
    ("²J =CH₂ terminal", "0-3 Hz"),
    ("³J aromático orto", "6-9 Hz (benzeno 7-8,5)"),
    ("⁴J aromático meta", "1-3 Hz"),
    ("⁵J aromático para", "0-1 Hz"),
    ("³J H-2/H-3 piridina", "4-6 Hz; H-3/H-4 7-9 Hz"),
    ("³J furano H-2/H-3", "≈ 1,8 Hz; H-3/H-4 ≈ 3,4 Hz"),
    ("³J CHO–CH", "1-3 Hz (aldeído alifático)"),
    ("⁴J alílico / W", "0-3 Hz"),
    ("³J H-1/H-2 anomérico β (gluco)", "7-8 Hz; α 3-4 Hz"),
    ("³J CH–NH amida", "6-9 Hz"),
    ("¹J C–H (¹³C acoplado)", "125 (sp³), 160 (sp²), 250 (sp)"),
    ("²J/³J H–F", "45-50 Hz (geminal), 5-25 Hz (vicinal); ¹J C–F ≈ 245 Hz"),
    ("²J/³J H–P", "7-25 Hz"),
]

MULT_NAMES = {1: "d", 2: "t", 3: "q", 4: "quint", 5: "sext", 6: "sept", 7: "oct", 8: "non"}

ATOMIC = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "S": 32.065, "P": 30.973762, "F": 18.9984032,
          "Cl": 35.453, "Br": 79.904, "I": 126.90447, "Si": 28.0855, "B": 10.811, "Na": 22.98977, "K": 39.0983,
          "Li": 6.941, "Se": 78.96, "D": 2.014101}


def parse_formula(f: str) -> dict:
    """'C9H8O4' -> {'C': 9, 'H': 8, 'O': 4}."""
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", f.replace(" ", "")):
        if el not in ATOMIC:
            raise ValueError(f"elemento desconhecido: {el}")
        out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def dbe(formula) -> float:
    """Índice de deficiência de hidrogênio (insaturações + anéis): C − H/2 − X/2 + N/2 + 1."""
    c = parse_formula(formula) if isinstance(formula, str) else formula
    C = c.get("C", 0) + c.get("Si", 0)
    H = c.get("H", 0) + c.get("D", 0)
    X = sum(c.get(x, 0) for x in ("F", "Cl", "Br", "I")) + c.get("Na", 0) + c.get("K", 0) + c.get("Li", 0)
    N = c.get("N", 0) + c.get("P", 0) + c.get("B", 0)
    return C - H / 2 - X / 2 + N / 2 + 1


def mass(formula) -> float:
    c = parse_formula(formula) if isinstance(formula, str) else formula
    return sum(ATOMIC[e] * n for e, n in c.items())


def ranges_for(nucleus: str, d: float):
    tab = RANGES_1H if nucleus == "1H" else RANGES_13C
    return [desc for lo, hi, desc in tab if lo <= d <= hi]


def impurity_lines(solvent: str | None, nucleus: str):
    """Lista (δ, nome) de impurezas comuns no solvente."""
    tab = IMPURITIES.get(solvent or "", {})
    k = 0 if nucleus == "1H" else 1
    out = []
    for name, vals in tab.items():
        for d in vals[k]:
            out.append((d, name))
    return sorted(out)
