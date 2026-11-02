"""
ftir_assign.py - candidatos de atribuição para bandas e triagem de grupos funcionais

    from ftir_assign import assign, screen_groups, load_table
    tab = find_bands(s)                  # de ftir_peaks
    tab = assign(tab)                    # adiciona colunas 'assignment' e 'candidates'
    print(screen_groups(tab))            # checklist de presença/ausência por região

A tabela (references/band_table.tsv) é uma tabela de correlação de grupos: diz o que
PODE estar em cada faixa. Ela não decide sozinha - a atribuição final exige coerência
entre bandas (ex.: C=O 1710 só é ácido se houver O–H 3300-2500 muito largo).
Use os candidatos como ponto de partida para o raciocínio, não como resposta.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE = os.path.join(HERE, "..", "references", "band_table.tsv")

_INT_RANK = {"vw": 0, "w": 1, "m": 2, "s": 3, "vs": 4}


def load_table(path: str | None = None) -> pd.DataFrame:
    t = pd.read_csv(path or DEFAULT_TABLE, sep="\t", dtype={"note": str}).fillna("")
    t["hi"], t["lo"] = t[["hi", "lo"]].max(axis=1), t[["hi", "lo"]].min(axis=1)
    return t


def _row_int(code: str):
    """'m-s' -> (2,3); 'v' (variável) -> (0,4)."""
    code = str(code).strip()
    if code in ("v", ""):
        return 0, 4
    parts = [p for p in code.split("-") if p in _INT_RANK]
    if not parts:
        return 0, 4
    r = [_INT_RANK[p] for p in parts]
    return min(r), max(r)


_QUAL = {"s": 0.40, "m": 0.15, "w": 0.05}


def _requirement_met(req: str, bands: pd.DataFrame | None, self_pos: float, min_rel: float = 0.03):
    """req = 'hi-lo[:q]|hi-lo[:q]' (qualquer uma basta). q = s|m|w (intensidade relativa
    mínima) ou br (banda larga). None = não avaliado (sem lista de bandas)."""
    if not req or bands is None or not len(bands):
        return None
    for part in str(req).split("|"):
        rng, _, q = part.partition(":")
        try:
            hi, lo = (float(v) for v in rng.split("-"))
        except ValueError:
            continue
        hi, lo = max(hi, lo), min(hi, lo)
        thr = _QUAL.get(q, min_rel)
        m = (bands.position <= hi) & (bands.position >= lo) & (bands.rel_intensity >= thr)
        if q == "br" and "shape" in bands:
            m &= bands["shape"] == "br"
        m &= (bands.position - self_pos).abs() > 5  # não vale a própria banda
        if m.any():
            return True
    return False


def candidates_for(pos: float, intensity: str = "", shape: str = "", table: pd.DataFrame | None = None,
                   tol: float = 10.0, n: int = 3, bands: pd.DataFrame | None = None,
                   components: bool = False, fwhm: float | None = None) -> list[dict]:
    """Lista ranqueada de linhas da tabela compatíveis com uma banda em `pos` cm-1.
    `bands` (a tabela completa de find_bands) ativa a checagem de bandas parceiras
    (coluna 'requires'): candidatos sem parceiro no espectro caem no ranking.
    `components=True` inclui sub-bandas de amida I de proteína (só faz sentido com
    2ª derivada/deconvolução)."""
    t = load_table() if table is None else table
    t = t[t["diagnostic"].astype(int) > 0]
    if not components:
        t = t[t["shape"] != "component"]
    m = (t["lo"] - tol <= pos) & (t["hi"] + tol >= pos)
    out = []
    for _, r in t[m].iterrows():
        width = r.hi - r.lo
        center = (r.hi + r.lo) / 2
        inside = r.lo <= pos <= r.hi
        # especificidade: faixas estreitas valem um pouco mais
        score = (1.0 + float(r.diagnostic)) / (1.0 + width / 400.0)
        ok = _requirement_met(r.get("requires", ""), bands, pos)
        if ok is True:
            score *= 1.25
        elif ok is False:
            score *= 0.35
        # proximidade do centro (e penalidade se só entrou pela tolerância)
        score *= 1.0 - 0.3 * min(1.0, abs(pos - center) / (width / 2 + tol))
        if not inside:
            score *= 0.5
        # intensidade compatível?
        if intensity in _INT_RANK:
            lo_i, hi_i = _row_int(r.intensity)
            k = _INT_RANK[intensity]
            gap = (lo_i - k) if k < lo_i else (k - hi_i if k > hi_i else 0)
            # banda observada bem mais fraca que o esperado pesa muito (ex.: 'C=O' fraca
            # quase nunca é carbonila principal); mais forte que o esperado pesa menos
            if k < lo_i:
                score *= {1: 0.8, 2: 0.5}.get(gap, 0.2)
            elif gap:
                score *= 0.8 if gap == 1 else 0.6
        # forma compatível? (banda larga observada vs. esperada)
        if shape == "br":
            score *= 1.3 if "broad" in r["shape"] else 0.7
        elif "broad" in r["shape"] and "medium" not in r["shape"]:
            score *= 0.35 if (fwhm is not None and np.isfinite(fwhm) and fwhm < 25) else 0.9
        out.append(dict(assignment=r.assignment, range=f"{int(r.hi)}-{int(r.lo)}", mode=r["mode"],
                        expected=r.intensity, category=r.category, score=round(score, 3), note=r.note,
                        partner={True: "parceira ok", False: "SEM parceira", None: ""}[ok]))
    out.sort(key=lambda d: -d["score"])
    return out[:n]


def assign(df: pd.DataFrame, table: pd.DataFrame | None = None, tol: float = 10.0, n: int = 3) -> pd.DataFrame:
    """Acrescenta ao DataFrame de find_bands:
       'assignment'  melhor candidato (texto curto)
       'candidates'  até n candidatos 'atribuição [faixa]' separados por ' | '"""
    t = load_table() if table is None else table
    best, cands = [], []
    for _, r in df.iterrows():
        c = candidates_for(float(r.position), r.get("intensity", ""), r.get("shape", ""), t, tol, n, bands=df,
                           fwhm=r.get("fwhm", None))
        best.append(c[0]["assignment"] if c else "")
        cands.append(" | ".join(f"{d['assignment']} [{d['range']}]" + (" (sem parceira)" if d["partner"] == "SEM parceira" else "")
                                for d in c))
    out = df.copy()
    out["assignment"] = best
    out["candidates"] = cands
    return out


# ------------------------------------------------------------------ triagem
def _has(df, hi, lo, min_rel=0.0, shape=None, kinds=("peak", "shoulder")):
    m = (df.position <= hi) & (df.position >= lo) & (df.rel_intensity >= min_rel) & (df.kind.isin(kinds))
    if shape:
        m &= df["shape"] == shape
    return df[m]


def screen_groups(df: pd.DataFrame) -> str:
    """Checklist das perguntas clássicas de interpretação (ordem de Pavia/Silverstein).
    Baseado só na lista de bandas: confirme olhando o espectro."""
    L = []

    def line(q, hit, detail):
        L.append(f"[{'SIM' if hit else 'não'}] {q}: {detail}")

    def fmt(sub):
        return ", ".join(f"{p}({i}{'/' + s if s else ''})" for p, i, s in
                         zip(sub.position, sub.intensity, sub["shape"])) or "—"

    co = _has(df, 1870, 1650, 0.25, kinds=("peak",))
    line("C=O (banda forte 1870-1650)", len(co) > 0, fmt(co))
    oh = _has(df, 3650, 3150, 0.10)
    oh_br = oh[oh["shape"] == "br"]
    line("O–H/N–H (3650-3150)", len(oh) > 0, fmt(oh) + (" -> larga: O–H/água em ligação de H" if len(oh_br) else ""))
    acid = _has(df, 3300, 2500, 0.05, shape="br")
    line("O–H de ácido muito largo (3300-2500)", len(acid) > 0 and len(co) > 0, fmt(acid))
    nh2 = _has(df, 3500, 3300, 0.05)
    line("Dupla N–H (NH2) 3500-3300", len(nh2) >= 2, fmt(nh2))
    sp2 = _has(df, 3100, 3000, 0.02)
    sp3 = _has(df, 2990, 2840, 0.05)
    line("C–H sp2 (>3000) / sp3 (<3000)", len(sp2) > 0 or len(sp3) > 0, f"sp2: {fmt(sp2)}; sp3: {fmt(sp3)}")
    ald = _has(df, 2760, 2690, 0.02)
    line("C–H de aldeído ~2720", len(ald) > 0 and len(co) > 0, fmt(ald))
    trip = _has(df, 2280, 2000, 0.05)
    trip = trip[~((trip.position <= 2400) & (trip.position >= 2280))]
    line("Tripla/cumulada (2280-2000)", len(trip) > 0, fmt(trip) + "  (cuidado: ruído do diamante no ATR)")
    aro = _has(df, 1625, 1575, 0.05)
    aro2 = _has(df, 1520, 1470, 0.05)
    oop = _has(df, 900, 680, 0.15)
    line("Anel aromático (~1600 + ~1500 + oop 900-680)", len(aro) > 0 and len(aro2) > 0, f"{fmt(aro)}; {fmt(aro2)}; oop: {fmt(oop)}")
    ce = _has(df, 1680, 1620, 0.05)
    line("C=C / C=N (1680-1620, fraca-média)", len(ce) > 0, fmt(ce) + "  (água δ ~1640 e amida I também caem aqui)")
    amide2 = _has(df, 1570, 1510, 0.15)
    line("Amida II (1570-1510) junto com C=O 1690-1630", len(amide2) > 0 and len(_has(df, 1690, 1630, 0.25)) > 0, fmt(amide2))
    cof = _has(df, 1300, 1000, 0.30)
    line("C–O / Si–O / P–O / S=O forte (1300-1000)", len(cof) > 0, fmt(cof))
    no2 = _has(df, 1570, 1500, 0.3)
    no2b = _has(df, 1390, 1330, 0.3)
    line("NO2 (1570-1500 + 1390-1330, ambas fortes)", len(no2) > 0 and len(no2b) > 0, f"{fmt(no2)}; {fmt(no2b)}")
    ch2 = _has(df, 730, 715, 0.03)
    line("(CH2)n n≥4 (~720)", len(ch2) > 0, fmt(ch2))
    low = _has(df, 700, 400, 0.30)
    line("Bandas fortes < 700 (M–O, C–X, inorgânicos)", len(low) > 0, fmt(low))
    art = df[df["flags"].astype(str).str.contains("CO2|H2O", regex=True)]
    if len(art):
        L.append(f"[!] Bandas em regiões atmosféricas (CO2/vapor): {fmt(art)}")
    return "\n".join(L)
