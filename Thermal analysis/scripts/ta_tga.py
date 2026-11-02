"""
ta_tga.py - análise de TG/DTG: etapas de perda de massa, temperaturas características, resíduo
e estequiometria (perdas teóricas, número de hidratação)

    from ta_tga import find_steps, char_temps, residues, theoretical_losses, hydration_number
    st = find_steps(r)                 # DataFrame: uma linha por etapa (T_inicio, T_onset, T_pico DTG, Δm %...)
    tc = char_temps(r, (1, 5, 10, 50)) # T1%, T5%, T10%, T50%
    theoretical_losses("CaC2O4·H2O", ["H2O", "CO", "CO2"])

Convenções: Δm em % da massa inicial (ou da base escolhida); T_onset/T_endset extrapolados pelo
método das tangentes (ISO 11358-1): reta do patamar anterior × tangente no ponto de inflexão
(pico da DTG). T_pico = máximo de −dm/dT.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, savgol_filter

from ta_io import Run
from ta_process import dtg, mass_pct, segments, heating_rate


def _dyn(r: Run):
    """Pontos dos segmentos de aquecimento (onde etapas vs T fazem sentido)."""
    segs = segments(r)
    m = np.zeros(r.n, bool)
    for s in segs:
        if s["type"] == "heat":
            m[s["i0"]:s["i1"]] = True
    if not m.any():
        m[:] = True
    return m, segs


def find_steps(r: Run, min_loss: float = 0.5, rel_prominence: float = 0.03, min_sep_K: float = 20.0,
               window_K: float = 6.0, basis=None) -> pd.DataFrame:
    """Etapas de perda (ou ganho) de massa nos trechos de aquecimento, delimitadas pelos mínimos da
    −DTG entre picos. Também lista perdas em isotermas (> min_loss) como etapas 'isoterma'.
    Colunas: etapa, tipo, T_inicio, T_fim, T_onset, T_endset, T_pico, DTG_max (%/°C e %/min),
    massa_inicio, massa_fim, perda_pct, perda_mg, sobreposta (bool), segmento."""
    m_all, unit = mass_pct(r, basis)
    if m_all is None:
        raise ValueError("sem canal de massa")
    d_all = dtg(r, "T", window_K=window_K, y=m_all)
    d_time = dtg(r, "time", window_K=window_K, y=m_all) if r.has("time") else None
    dyn, segs = _dyn(r)
    rows = []
    sm = r.sample_mass
    beta = heating_rate(r)
    for s in segs:
        a, b = s["i0"], s["i1"]
        T, m = r.T[a:b], m_all[a:b]
        ok = np.isfinite(T) & np.isfinite(m)
        if ok.sum() < 20:
            continue
        small_ramp = s["type"] == "heat" and abs(float(np.nanmax(T) - np.nanmin(T))) < 20
        if s["type"] == "iso" or small_ramp or (s["type"] != "heat" and not any(x["type"] == "heat" for x in segs)):
            dm = m[ok][0] - m[ok][-1]
            if abs(dm) >= min_loss:
                lab = "isoterma" if s["type"] == "iso" or small_ramp else ("resfriamento" if s["type"] == "cool" else s["type"])
                rows.append(dict(tipo=f"{lab} ({'perda' if dm > 0 else 'ganho'})", segmento=s["n"],
                                 T_inicio=float(T[ok][0]), T_fim=float(T[ok][-1]), T_onset=np.nan, T_endset=np.nan,
                                 T_pico=np.nan, DTG_max_pct_C=np.nan, DTG_max_pct_min=np.nan,
                                 massa_inicio=float(m[ok][0]), massa_fim=float(m[ok][-1]), perda_pct=float(dm),
                                 duracao_min=float(r.time[a:b][ok][-1] - r.time[a:b][ok][0]) if r.has("time") else np.nan,
                                 sobreposta=False))
            continue
        if s["type"] != "heat":
            continue
        d = -d_all[a:b]                       # -dm/dT: picos positivos para perda
        Tm = T.copy()
        good = ok & np.isfinite(d)
        if good.sum() < 20:
            continue
        idx = np.where(good)[0]
        Tg, dg, mg = Tm[idx], d[idx], m[idx]
        step = float(np.median(np.abs(np.diff(Tg)))) or 0.5
        dist = max(1, int(round(min_sep_K / step)))
        span = float(np.nanpercentile(np.abs(dg), 99.5))
        noise = float(1.4826 * np.median(np.abs(np.diff(dg) - np.median(np.diff(dg))))) + 1e-12
        prom = max(rel_prominence * span, 6 * noise)
        pk, props = find_peaks(dg, prominence=prom, distance=dist)
        # perdas de massa (picos positivos); ganhos (oxidação de metal) como picos negativos
        pk_gain, _ = find_peaks(-dg, prominence=prom * 2, distance=dist)
        events = sorted([(p, "perda") for p in pk] + [(p, "ganho") for p in pk_gain])
        if not events:
            continue
        # limites entre eventos = mínimo de |dg| entre picos; nas pontas, onde |dg| cai a 3% do pico
        bounds = []
        for k, (p, kind) in enumerate(events):
            h = abs(dg[p])
            if k == 0:
                L = p
                while L > 0 and abs(dg[L]) > 0.03 * h:
                    L -= 1
            else:
                L = bounds[-1][1]
            if k == len(events) - 1:
                R = p
                while R < len(dg) - 1 and abs(dg[R]) > 0.03 * h:
                    R += 1
            else:
                q = events[k + 1][0]
                R = p + int(np.argmin(np.abs(dg[p:q + 1])))
            bounds.append((L, R))
        edge = max(3, int(round(window_K / step)))
        for (p, kind), (L, R) in zip(events, bounds):
            dm = float(mg[L] - mg[R])
            if abs(dm) < min_loss:
                continue
            try:
                fw = float(peak_widths(np.abs(dg), [p], rel_height=0.5)[0][0] * step)
            except Exception:
                fw = np.nan
            if np.isfinite(fw) and fw < max(2.5, 0.8 * window_K) and abs(dm) < 3:
                continue                       # pico estreito demais: salto/ruído da balança, não etapa
            if (p < edge or p > len(dg) - 1 - edge) and abs(dm) < 5:
                continue                       # artefato de borda da derivada
            kind = "perda" if dm > 0 else "ganho"
            # tangente no ponto de inflexão (pico DTG) cruzando os patamares
            Tp, mp, slope = float(Tg[p]), float(mg[p]), -float(dg[p])     # dm/dT no pico
            w = max(3, int(round(5.0 / step)))
            pre = slice(max(0, L - w), L + 1)
            post = slice(R, min(len(Tg), R + w + 1))
            def level(sl):
                if sl.stop - sl.start >= 3:
                    c = np.polyfit(Tg[sl], mg[sl], 1)
                    return c
                return np.array([0.0, mg[sl.start]])
            c1, c2 = level(pre), level(post)
            def cross(c):
                den = slope - c[0]
                return float((c[1] - (mp - slope * Tp)) / den) if abs(den) > 1e-12 else np.nan
            T_on, T_end = cross(c1), cross(c2)
            if not (Tg[L] - 30 <= T_on <= Tp):
                T_on = np.nan
            if not (Tp <= T_end <= Tg[R] + 30):
                T_end = np.nan
            sobre = bool(L > 0 and abs(dg[L]) > 0.15 * abs(dg[p])) or bool(R < len(dg) - 1 and abs(dg[R]) > 0.15 * abs(dg[p]))
            ii = idx[p]
            rows.append(dict(tipo="perda" if kind == "perda" else "ganho", segmento=s["n"],
                             T_inicio=float(Tg[L]), T_fim=float(Tg[R]), T_onset=T_on, T_endset=T_end, T_pico=Tp,
                             DTG_max_pct_C=float(-dg[p]),
                             DTG_max_pct_min=float(d_time[a + ii]) if d_time is not None else np.nan,
                             massa_inicio=float(mg[L]), massa_fim=float(mg[R]), perda_pct=dm, duracao_min=np.nan,
                             largura_DTG_K=round(fw, 1) if np.isfinite(fw) else np.nan, sobreposta=sobre))
    df = pd.DataFrame(rows)
    if len(df):
        df.insert(0, "etapa", range(1, len(df) + 1))
        if sm:
            df["perda_mg"] = df["perda_pct"] / 100 * sm
        for c in ("T_inicio", "T_fim", "T_onset", "T_endset", "T_pico"):
            df[c] = df[c].round(1)
        for c in ("perda_pct", "massa_inicio", "massa_fim"):
            df[c] = df[c].round(2)
        df["DTG_max_pct_C"] = df["DTG_max_pct_C"].round(4)
        df["DTG_max_pct_min"] = df["DTG_max_pct_min"].round(3)
    df.attrs.update(unit=unit, beta=beta, title=r.title, basis=basis)
    return df


def char_temps(r: Run, losses=(1, 2, 5, 10, 20, 50), basis=None, segment_type: str = "heat") -> dict:
    """Temperatura em que a perda acumulada atinge x% (T1%, T5% = 'Td5%' ...), só nos trechos de
    aquecimento. basis=150 -> perdas relativas à massa a 150 °C (exclui umidade)."""
    m, _ = mass_pct(r, basis)
    dyn, segs = _dyn(r)
    T = np.where(dyn, r.T, np.nan)
    out = {}
    start = 0
    if basis is not None:
        start = int(np.where(r.T >= float(basis))[0][0])
    loss = 100.0 - m
    for x in losses:
        idx = np.where((np.arange(r.n) >= start) & dyn & (loss >= x))[0]
        out[f"T{x:g}%"] = float(np.round(T[idx[0]], 1)) if idx.size else np.nan
    return out


def residues(r: Run, at=(200, 400, 600, 800, 900, 1000), basis=None) -> dict:
    """Massa residual (%) em temperaturas fixas (1º cruzamento em aquecimento) e no fim do ensaio."""
    m, _ = mass_pct(r, basis)
    dyn, _ = _dyn(r)
    out = {}
    for T0 in at:
        idx = np.where(dyn & (r.T >= T0))[0]
        if idx.size:
            out[f"m@{T0:g}°C"] = float(np.round(m[idx[0]], 2))
    fin = m[np.isfinite(m)]
    out["final"] = float(np.round(fin[-1], 2)) if fin.size else np.nan
    out["T_final"] = float(np.round(r.T[np.isfinite(m)][-1], 1)) if fin.size else np.nan
    return out


# ============================================================================ estequiometria
ATOMIC = {"H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
          "F": 18.998, "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
          "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867,
          "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
          "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904, "Kr": 83.798,
          "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.95, "Ru": 101.07,
          "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71, "Sb": 121.76,
          "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
          "Pr": 140.91, "Nd": 144.24, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93, "Dy": 162.50,
          "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95,
          "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
          "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Th": 232.04, "U": 238.03}


def _parse_group(s: str) -> dict:
    tokens = re.findall(r"([A-Z][a-z]?|\(|\)|\[|\]|\d+(?:\.\d+)?)", s)
    stack = [{}]
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t in "([":
            stack.append({})
        elif t in ")]":
            grp = stack.pop()
            mult = 1.0
            if i + 1 < len(tokens) and re.match(r"^\d", tokens[i + 1]):
                mult = float(tokens[i + 1])
                i += 1
            for k, v in grp.items():
                stack[-1][k] = stack[-1].get(k, 0) + v * mult
        elif re.match(r"^[A-Z]", t):
            if t not in ATOMIC:
                raise ValueError(f"elemento desconhecido: {t}")
            n = 1.0
            if i + 1 < len(tokens) and re.match(r"^\d", tokens[i + 1]):
                n = float(tokens[i + 1])
                i += 1
            stack[-1][t] = stack[-1].get(t, 0) + n
        i += 1
    return stack[0]


def composition(formula: str) -> dict:
    """'CuSO4·5H2O', 'CaC2O4.H2O', 'Al2(SO4)3*18H2O', 'Mg(OH)2' -> {elemento: n}."""
    parts = re.split(r"[·•∙*]|\.(?=\d*[A-Z(])", formula.replace(" ", ""))
    total = {}
    for p in parts:
        m = re.match(r"^(\d+(?:\.\d+)?)(.*)$", p)
        k, body = (float(m.group(1)), m.group(2)) if m else (1.0, p)
        for el, n in _parse_group(body).items():
            total[el] = total.get(el, 0) + k * n
    return total


def molar_mass(formula: str) -> float:
    return float(sum(ATOMIC[e] * n for e, n in composition(formula).items()))


def theoretical_losses(reactant: str, released) -> pd.DataFrame:
    """Perdas teóricas sucessivas. released: lista de espécies liberadas em cada etapa, ex.
    ['H2O', 'CO', 'CO2'] ou ['2H2O', '2H2O', 'H2O'] (coeficiente opcional) ou [['H2O','CO2'], ...].
    Retorna etapa, espécie, perda_pct (da massa inicial), acumulada, resíduo e fórmula do resíduo."""
    M0 = molar_mass(reactant)
    rows, acc = [], 0.0
    for k, sp in enumerate(released, 1):
        items = sp if isinstance(sp, (list, tuple)) else [sp]
        dm = 0.0
        for it in items:
            mm = re.match(r"^(\d+(?:\.\d+)?)\s*(.*)$", it)
            coef, f = (float(mm.group(1)), mm.group(2)) if mm else (1.0, it)
            dm += coef * molar_mass(f)
        acc += dm
        rows.append(dict(etapa=k, especie=" + ".join(items), perda_pct=round(100 * dm / M0, 2),
                         acumulada_pct=round(100 * acc / M0, 2), residuo_pct=round(100 * (M0 - acc) / M0, 2)))
    df = pd.DataFrame(rows)
    df.attrs["M"] = M0
    return df


def hydration_number(loss_pct: float, anhydrous: str, water: str = "H2O") -> float:
    """n em X·nH2O a partir da perda de água (% da massa do hidrato):
    loss = n·M(H2O) / (M(X) + n·M(H2O))  ->  n = loss·M(X) / (M(H2O)·(100 − loss))."""
    return float(loss_pct * molar_mass(anhydrous) / (molar_mass(water) * (100.0 - loss_pct)))


def compare_with_theory(steps: pd.DataFrame, theory: pd.DataFrame, tol: float = 1.5) -> pd.DataFrame:
    """Casa etapas experimentais (em ordem) com as teóricas e calcula o desvio."""
    exp = steps[steps["tipo"].isin(["perda", "isoterma"])].reset_index(drop=True)
    rows = []
    for i, th in theory.iterrows():
        e = exp.iloc[i] if i < len(exp) else None
        rows.append(dict(etapa=th["etapa"], especie=th["especie"], teorica_pct=th["perda_pct"],
                         experimental_pct=round(float(e["perda_pct"]), 2) if e is not None else np.nan,
                         T_pico=e["T_pico"] if e is not None else np.nan,
                         desvio=round(float(e["perda_pct"]) - th["perda_pct"], 2) if e is not None else np.nan))
    df = pd.DataFrame(rows)
    df["ok"] = df["desvio"].abs() <= tol
    return df
