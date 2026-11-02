"""
nmr_peaks.py - picos, integração, análise de multipletos (J) e listas no formato de artigo

    from nmr_peaks import pick_peaks, analyze_1h, analyze_13c, acs_1h, acs_13c
    res = analyze_1h(s, formula="C9H8O4")      # multipletos com δ, multiplicidade, J, integral (nH)
    print(acs_1h(res))                         # "¹H NMR (300 MHz, CDCl3) δ 8.12 (dd, J = 7.9, 1.7 Hz, 1H), ..."
    r13 = analyze_13c(s13, dept135=d135, dept90=d90)

Multipletos de 1ª ordem: a lista de linhas (com intensidades em unidades inteiras) é "desconvoluída"
em dubletos sucessivos (árvore de J; Hoye et al., J. Org. Chem. 59 (1994) 4096 e 67 (2002) 4014). Se
não houver solução de 1ª ordem, o sinal é reportado como m (faixa) e, se for o caso, como AB.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths

from nmr_data import MULT_NAMES, SOLVENTS, SOLVENT_PATTERNS, impurity_lines, parse_formula, ranges_for
from nmr_process import noise


# ----------------------------------------------------------------------------- picos

def pick_peaks(s, k: float = 5.0, rel: float = 0.0, lo=None, hi=None, negative=False, sigma=None) -> pd.DataFrame:
    """Picos acima de k·σ (e de rel × máximo). Posição por interpolação parabólica; largura a meia
    altura em Hz. negative=True também procura picos negativos (DEPT/APT)."""
    sigma = sigma or noise(s.y)
    y = np.asarray(s.y, float)
    m = np.ones_like(y, bool)
    if lo is not None:
        m &= s.ppm >= lo
    if hi is not None:
        m &= s.ppm <= hi
    rows = []
    signs = (1, -1) if negative else (1,)
    for sg in signs:
        yy = sg * y
        thr = max(k * sigma, rel * np.max(np.abs(y[m])) if m.any() else 0)
        idx, prop = find_peaks(yy, height=thr, prominence=0.5 * k * sigma)
        idx = idx[m[idx]]
        if not idx.size:
            continue
        w = peak_widths(yy, idx, rel_height=0.5)[0]
        for i, wi in zip(idx, w):
            if 0 < i < y.size - 1:
                a, b, c = yy[i - 1], yy[i], yy[i + 1]
                den = a - 2 * b + c
                off = 0.5 * (a - c) / den if den != 0 else 0.0
                off = max(-0.5, min(0.5, off))
            else:
                off = 0.0
            step = s.ppm[1] - s.ppm[0] if y.size > 1 else 0
            p = s.ppm[i] + off * step
            rows.append(dict(ppm=float(p), hz=float(p * s.sf), altura=float(sg * yy[i]), snr=float(yy[i] / sigma),
                             fwhm_hz=float(wi * abs(step) * s.sf), sinal=sg, idx=int(i)))
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("ppm", ascending=False).reset_index(drop=True)
    return df


# ----------------------------------------------------------------------------- solvente e impurezas

def solvent_mask(s, extra_hz=3.0):
    """Faixas (ppm) ocupadas pelo sinal residual do solvente (e água em ¹H)."""
    out = []
    pats = SOLVENT_PATTERNS.get(s.nucleus, {}).get(s.solvent or "", [])
    for c, J, rel in pats:
        half = (len(rel) - 1) / 2 * J + extra_hz + 2.0
        out.append((c - half / s.sf, c + half / s.sf, "solvente"))
    sol = SOLVENTS.get(s.solvent or "")
    if sol and s.nucleus == "1H" and sol.get("water") is not None:
        out.append((sol["water"] - 0.04, sol["water"] + 0.04, "água"))
    return out


def label_impurity(s, ppm, tol=None):
    """Nome do solvente/impureza comum mais próxima (ou None)."""
    tol = tol or (0.02 if s.nucleus == "1H" else 0.3)
    for lo, hi, name in solvent_mask(s):
        if lo <= ppm <= hi:
            return name
    best = None
    for d, name in impurity_lines(s.solvent, s.nucleus):
        if abs(d - ppm) <= tol and (best is None or abs(d - ppm) < best[0]):
            best = (abs(d - ppm), name)
    return f"possível {best[1]}" if best else None


# ----------------------------------------------------------------------------- multipletos

def _units(heights, max_units=128):
    """Converte alturas em 'unidades' inteiras (contagem de linhas coincidentes). Testa várias escalas e
    prefere total potência de 2 com menor erro de arredondamento."""
    h = np.asarray(heights, float)
    base_opts = [min(h[0], h[-1]), 0.5 * (h[0] + h[-1]), h.min()]
    best = None
    for b in base_opts:
        for f in (1.0, 0.9, 1.1, 0.8, 1.25, 0.75, 1.33, 0.67, 0.5, 2.0):
            u = b * f
            if u <= 0:
                continue
            r = h / u
            n = np.maximum(1, np.round(r)).astype(int)
            tot = n.sum()
            if tot > max_units:
                continue
            err = np.mean(np.abs(r - n) / n)
            pow2 = (tot & (tot - 1)) == 0
            score = err + (0 if pow2 else 0.35) + 0.01 * tot / 8
            if best is None or score < best[0]:
                best = (score, n, err, pow2)
    return best


def _deconvolve(P, tol, depth=0):
    """P: posições (Hz) ordenadas, uma por unidade de linha. Retorna lista de J ou None."""
    if len(P) == 1:
        return []
    if len(P) % 2:
        return None
    seen = set()
    for k in range(1, len(P)):
        J = P[k] - P[0]
        if J < tol * 0.8:
            continue
        key = round(J / tol)
        if key in seen:
            continue
        seen.add(key)
        rem = list(P)
        Q = []
        ok = True
        while rem:
            p = rem.pop(0)
            target = p + J
            j = None
            for ii, q in enumerate(rem):
                if abs(q - target) <= tol:
                    j = ii
                    break
                if q > target + tol:
                    break
            if j is None:
                ok = False
                break
            q = rem.pop(j)
            Q.append(0.5 * (p + q))
        if not ok:
            continue
        sub = _deconvolve(sorted(Q), tol * 1.15, depth + 1)
        if sub is not None:
            return [J] + sub
    return None


def _label(Js, tolJ):
    """Agrupa J iguais: [7.0, 7.1, 1.5] -> ('td', [7.05, 1.5])."""
    Js = sorted(Js, reverse=True)
    groups = []
    for J in Js:
        for g in groups:
            if abs(np.mean(g) - J) <= tolJ:
                g.append(J)
                break
        else:
            groups.append([J])
    letters = "".join(MULT_NAMES.get(len(g), f"({len(g) + 1})") if len(g) <= 8 else "m" for g in groups)
    if len(letters) > 12:
        return "m", []
    return letters, [float(np.mean(g)) for g in groups]


def _simulate(center, Js, n_equiv):
    pos = np.array([center])
    for J, k in zip(Js, n_equiv):
        for _ in range(k):
            pos = np.concatenate([pos - J / 2, pos + J / 2])
    return np.sort(pos)


def analyze_lines(hz, heights, fwhm, tol, sf, max_j=20.0):
    """Analisa uma lista de linhas (Hz, alturas) de um multipleto. Retorna dict com mult, J, center."""
    o = np.argsort(hz)
    hz, heights = np.asarray(hz, float)[o], np.asarray(heights, float)[o]
    n = hz.size
    cen_w = float(np.sum(hz * heights) / np.sum(heights))
    if n == 1:
        br = fwhm is not None and fwhm > 3.5
        return dict(mult="sl" if br else "s", J=[], center_hz=float(hz[0]), ok=True, note="largo" if br else "")
    # simetria
    mid = 0.5 * (hz[0] + hz[-1])
    sym = np.allclose(hz - mid, -(hz[::-1] - mid), atol=2 * tol)
    hr = heights / heights[::-1]
    sym_h = np.all((hr > 0.55) & (hr < 1.8))
    if sym:
        # simetriza (remove o efeito de telhado nas posições e intensidades) antes de desconvoluir
        d = hz - mid
        dsym = 0.5 * (d - d[::-1])
        hsym = 0.5 * (heights + heights[::-1])
        u = _units(hsym)
        if u is not None:
            _, units, err, pow2 = u
            P = sorted(np.repeat(mid + dsym, units))
            Js = _deconvolve(P, tol) if pow2 else None
            if Js is not None and len(Js) >= 1 and max(Js) <= max_j:
                lab, Jg = _label(Js, max(0.35, tol))
                # pede intensidades coerentes (o padrão simulado tem que reproduzir as linhas)
                note = "" if sym_h else "intensidades assimétricas (efeito de telhado/2ª ordem?)"
                if lab != "m":
                    return dict(mult=lab, J=Jg, center_hz=float(mid), ok=True, note=note, units=units.tolist())
    # AB (4 linhas, internas mais altas)
    if n == 4:
        l = hz
        J1, J2 = l[1] - l[0], l[3] - l[2]
        if abs(J1 - J2) <= 2 * tol and heights[1] > heights[0] and heights[2] > heights[3]:
            J = 0.5 * (J1 + J2)
            dnu = np.sqrt(max((l[3] - l[0]) * (l[2] - l[1]), 0))
            c = 0.5 * (l[0] + l[3])
            return dict(mult="AB", J=[float(J)], center_hz=float(c), ok=True, AB=(float(c + dnu / 2), float(c - dnu / 2)),
                        note=f"sistema AB: Δν = {dnu:.1f} Hz, Δν/J = {dnu / J:.1f}")
    # assimétrico com 2 linhas iguais-ish = d com telhado
    if n == 2 and (hz[1] - hz[0]) <= max_j:
        return dict(mult="d", J=[float(hz[1] - hz[0])], center_hz=float(np.mean(hz)), ok=True,
                    note="dubleto com intensidades diferentes (telhado/sobreposição)" if not sym_h else "")
    return dict(mult="m", J=[], center_hz=cen_w, ok=False, note="")


def _regions(s, sigma, k, gap_hz, exclude=None):
    """Regiões contíguas acima de k·σ (unidas se a lacuna < gap_hz), expandidas até a linha de base
    (y < 2σ ou subida do vale). exclude: máscara booleana de pontos que não entram (solvente)."""
    y = s.y.copy()
    above = y > k * sigma
    if exclude is not None:
        above &= ~exclude
    if not above.any():
        return []
    idx = np.where(above)[0]
    step_hz = abs(s.ppm[1] - s.ppm[0]) * s.sf
    gap = max(1, int(gap_hz / step_hz))
    regs = []
    a = prev = idx[0]
    for i in idx[1:]:
        if i - prev > gap:
            regs.append([a, prev])
            a = i
        prev = i
    regs.append([a, prev])
    lim = 2.0 * sigma
    maxpad = max(2, int(6.0 / step_hz))
    for r in regs:
        a, b = r
        for _ in range(maxpad):
            if a > 0 and y[a - 1] > lim and y[a - 1] <= y[a] * 1.02 and not (exclude is not None and exclude[a - 1]):
                a -= 1
            else:
                break
        for _ in range(maxpad):
            if b < y.size - 1 and y[b + 1] > lim and y[b + 1] <= y[b] * 1.02 and not (exclude is not None and exclude[b + 1]):
                b += 1
            else:
                break
        r[0], r[1] = a, b
    out = []
    for a, b in regs:
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(r) for r in out]


def _split_valleys(y, a, b, thr, frac=0.2, min_sep_pts=0):
    """Divide [a, b] nos vales cujo mínimo < frac × o menor dos dois máximos vizinhos, quando os máximos
    estão a mais de min_sep_pts pontos (linhas de um mesmo multipleto ficam juntas)."""
    seg = y[a:b + 1]
    if seg.size < 5:
        return [(a, b)]
    pk, _ = find_peaks(seg, height=thr)
    if pk.size < 2:
        return [(a, b)]
    # só máximos relevantes (picos minúsculos entre dois sinais não impedem o corte)
    pk = pk[seg[pk] >= 0.08 * seg[pk].max()]
    if pk.size < 2:
        return [(a, b)]
    cuts = []
    for p, q in zip(pk[:-1], pk[1:]):
        if q - p <= min_sep_pts:
            continue
        j = p + int(np.argmin(seg[p:q + 1]))
        if seg[j] < frac * min(seg[p], seg[q]):
            cuts.append(j)
    if not cuts:
        return [(a, b)]
    out = []
    start = a
    for c in cuts:
        out.append((start, a + c))
        start = a + c + 1
    out.append((start, b))
    return out


def _broad_regions(s, sigma, taken, width_hz=12.0, k=6.0):
    """Sinais largos (OH, NH, COOH) que escapam do limiar de picos: busca no espectro suavizado."""
    from scipy.ndimage import gaussian_filter1d
    step_hz = abs(s.ppm[1] - s.ppm[0]) * s.sf
    ys = gaussian_filter1d(s.y, max(1.0, width_hz / step_hz))
    sig_s = noise(ys)
    m = ys > k * sig_s
    for a, b in taken:
        m[max(0, a - int(40 / step_hz)):b + int(40 / step_hz) + 1] = False
    out = []
    if not m.any():
        return out
    idx = np.where(m)[0]
    a = prev = idx[0]
    for i in list(idx[1:]) + [None]:
        if i is None or i - prev > 1:
            if (prev - a) * step_hz >= 2 * width_hz:
                out.append((a, prev))
            if i is not None:
                a = i
        if i is not None:
            prev = i
    return out


def _integrate(s, a, b):
    return float(np.sum(s.y[a:b + 1]) * abs(s.ppm[1] - s.ppm[0]))


def _scale_integrals(I, nH=None, fixed=None):
    """Escala as integrais em número de H. fixed=(índice, nH) fixa um sinal; nH = total (fórmula)."""
    I = np.asarray(I, float)
    if fixed is not None:
        j, h = fixed
        return h / I[j], f"sinal {j + 1} = {h} H"
    cands = []
    big = I >= 0.05 * I.max()
    if nH:
        for e in range(0, 4):
            if nH - e <= 0:
                continue
            k = (nH - e) / I.sum()
            r = I * k
            dev = np.sum(np.abs(r - np.round(r))[big]) / max(1, big.sum())
            cands.append((dev + 0.08 * e, k, f"total = {nH - e} H" + (f" (fórmula {nH} H; {e} H trocável(is) ausente(s)?)" if e else " (fórmula)")))
    else:
        for j in np.where(big)[0]:
            for m in (1, 2, 3, 4, 6, 9):
                k = m / I[j]
                r = I * k
                if np.any(r[big] < 0.6):
                    continue
                dev = np.sum(np.abs(r - np.round(r))[big]) / max(1, big.sum())
                cands.append((dev + 0.004 * np.round(r).sum(), k, "escala automática (menor desvio de inteiros)"))
    if not cands:
        k = 1.0 / I[big].min()
        return k, "menor sinal = 1 H"
    bestdev = min(c[0] for c in cands)
    ok = [c for c in cands if c[0] <= bestdev + 0.04]
    best = min(ok, key=lambda c: c[1]) if not nH else min(cands, key=lambda c: c[0])
    return best[1], best[2]


def analyze_1h(s, k: float = 6.0, gap_hz: float = 2.5, formula=None, nH=None, norm=None, sigma=None,
               min_rel=0.004, max_j=18.0) -> dict:
    """Análise completa do ¹H: regiões (multipletos), integrais em nH, multiplicidade e J, rótulos de
    solvente/impurezas, acoplamentos compartilhados. norm='2.36=3' fixa a integral do sinal em 2,36 ppm."""
    sigma = sigma or noise(s.y)
    if nH is None and formula:
        nH = parse_formula(formula).get("H")
    # solvente fora das regiões da amostra (vira região própria)
    excl = np.zeros(s.y.size, bool)
    sol_regs = []
    for lo, hi, name in solvent_mask(s):
        mm = (s.ppm >= lo) & (s.ppm <= hi)
        if mm.any():
            excl |= mm
            ii = np.where(mm)[0]
            sol_regs.append((ii[0], ii[-1], name))
    # referência de altura: mediana dos 10 maiores picos (um sinal gigante isolado não esconde os demais)
    pk_all = pick_peaks(s, k=k, sigma=sigma)
    if len(pk_all):
        keep = np.ones(len(pk_all), bool)
        for lo, hi, _ in solvent_mask(s):
            keep &= ~((pk_all.ppm >= lo) & (pk_all.ppm <= hi)).values
        hs = np.sort(pk_all.altura.values[keep])[::-1][:10]
        ymax = float(np.median(hs)) if hs.size else float(np.max(s.y))
    else:
        ymax = float(np.max(s.y))
    k_eff = max(k, min_rel * ymax / sigma)
    regs = _regions(s, sigma, k_eff, gap_hz, exclude=excl)
    step_hz = abs(s.ppm[1] - s.ppm[0]) * s.sf
    tol = max(0.3, 1.6 * step_hz)

    def one(a, b, forced=None):
        ppm_hi, ppm_lo = float(s.ppm[a]), float(s.ppm[b])
        I = _integrate(s, a, b)
        seg = s.copy(ppm=s.ppm[a:b + 1], y=s.y[a:b + 1])
        pk0 = pick_peaks(seg, k=max(4.0, k_eff * 0.6), sigma=sigma)
        pk = pk0[pk0.altura >= 0.04 * pk0.altura.max()] if len(pk0) else pk0
        if not len(pk):
            if forced:
                cen = float(np.sum(s.ppm[a:b + 1] * s.y[a:b + 1]) / max(np.sum(s.y[a:b + 1]), 1e-30))
                return dict(ppm=round(cen, 4), ppm_max=ppm_hi, ppm_min=ppm_lo, integral=I, mult=forced, J=[], n_linhas=0,
                            fwhm_hz=float((b - a) * step_hz / 2), nota="sinal largo (trocável?)", rotulo="", AB=None,
                            linhas_hz=[], alturas=[], i0=a, i1=b)
            return None
        res = analyze_lines(pk.hz.values, pk.altura.values, float(pk.fwhm_hz.median()), tol, s.sf)
        if res["mult"] == "m":
            # linhas fracas extras (satélites, impurezas, 2ª ordem) podem impedir a análise: tenta sem elas
            for frac in (0.08, 0.12):
                pk2 = pk0[pk0.altura >= frac * pk0.altura.max()]
                lost = 1 - pk2.altura.sum() / pk.altura.sum()
                if 1 < len(pk2) < len(pk) and lost < 0.08:
                    r2 = analyze_lines(pk2.hz.values, pk2.altura.values, float(pk2.fwhm_hz.median()), tol, s.sf)
                    if r2["mult"] != "m":
                        res, pk = r2, pk2
                        res["note"] = (res.get("note", "") + " linhas fracas ignoradas").strip()
                        break
        cen_ppm = res["center_hz"] / s.sf
        return dict(ppm=round(cen_ppm, 4), ppm_max=ppm_hi, ppm_min=ppm_lo, integral=I, mult=forced or res["mult"],
                    J=[round(j, 1) for j in res["J"]], n_linhas=len(pk), fwhm_hz=float(pk.fwhm_hz.median()),
                    nota=res.get("note", ""), rotulo=label_impurity(s, cen_ppm) or "", AB=res.get("AB"),
                    linhas_hz=pk.hz.values.tolist(), alturas=pk.altura.values.tolist(), i0=a, i1=b)
    # divide regiões em vales profundos (caudas de Lorentz unem sinais distintos)
    regs2 = []
    for a, b in regs:
        regs2.extend(_split_valleys(s.y, a, b, k_eff * sigma, min_sep_pts=int(max_j / step_hz)))
    items = [one(a, b) for a, b in regs2]
    items = [it for it in items if it]
    # une regiões vizinhas (< max_J) quando juntas formam um multipleto de 1ª ordem e separadas não
    merged = True
    while merged and len(items) > 1:
        merged = False
        for i in range(len(items) - 1):
            A, B = items[i], items[i + 1]
            gap = (B["ppm_max"] - A["ppm_min"]) * -s.sf if B["ppm_max"] < A["ppm_min"] else 0
            sep = abs(A["ppm"] - B["ppm"]) * s.sf
            if sep > max_j or A["rotulo"] or B["rotulo"]:
                continue
            ratio = A["integral"] / B["integral"] if B["integral"] > 0 else 9
            if not (0.2 < ratio < 5.0) or B["i0"] - A["i1"] > 3:
                continue
            C = one(A["i0"], B["i1"])
            if C and C["mult"] not in ("m", "s", "sl") and C["n_linhas"] == A["n_linhas"] + B["n_linhas"] \
                    and (A["mult"] in ("s", "d", "sl") or B["mult"] in ("s", "d", "sl")):
                items[i:i + 2] = [C]
                merged = True
                break
    # fragmentos de multipletos de 2ª ordem/sobrepostos: junta pedaços pequenos encostados
    if len(items) > 2:
        typical = np.median([it["integral"] for it in items if not it["rotulo"]] or [1.0])
        merged = True
        while merged and len(items) > 1:
            merged = False
            for i in range(len(items) - 1):
                A, B = items[i], items[i + 1]
                if A["rotulo"] or B["rotulo"] or B["i0"] - A["i1"] > 3:
                    continue
                if abs(A["ppm"] - B["ppm"]) * s.sf > max_j:
                    continue
                if min(A["integral"], B["integral"]) < 0.45 * typical or "m" in (A["mult"], B["mult"]) and \
                        min(A["integral"], B["integral"]) < 0.8 * typical:
                    C = one(A["i0"], B["i1"])
                    if C:
                        items[i:i + 2] = [C]
                        merged = True
                        break
    # sinais largos
    for a, b in _broad_regions(s, sigma, [(it["i0"], it["i1"]) for it in items] + [(x[0], x[1]) for x in sol_regs]):
        it = one(a, b, forced="sl")
        if it:
            it["nota"] = "sinal largo (OH/NH/COOH trocável? confirme com D₂O)"
            items.append(it)
    for a, b, name in sol_regs:
        it = one(a, b)
        if it:
            it["rotulo"] = name
            items.append(it)
    for it in items:
        if it["rotulo"]:
            continue
        for a, b, name in sol_regs:
            lo_, hi_ = s.ppm[b], s.ppm[a]
            if it["ppm_min"] - 5 / s.sf <= hi_ and it["ppm_max"] + 5 / s.sf >= lo_:
                it["nota"] = (it["nota"] + f" sobreposto ao sinal de {name}: multiplicidade/integral incertas").strip()
    rows = sorted(items, key=lambda r: -r["ppm"])
    df = pd.DataFrame(rows)
    if not len(df):
        return dict(table=df, sigma=sigma, norm_desc="", nH=nH)
    sample = df.rotulo.eq("") | df.rotulo.str.startswith("possível")
    # sinais minúsculos (satélites de ¹³C, ruído) não entram na normalização
    real = sample & (df.integral > 0.02 * df.integral[sample].max())
    fixed = None
    if norm:
        a_, _, b_ = str(norm).partition("=")
        j = int(np.argmin(np.abs(df.ppm.values - float(a_))))
        fixed = (int(np.where(df.index == j)[0][0]), float(b_))
    Ireal = df.integral[real].values
    if fixed is not None:
        kscale = fixed[1] / df.integral.iloc[fixed[0]]
        ndesc = f"δ {df.ppm.iloc[fixed[0]]:.2f} = {fixed[1]:g} H"
    else:
        kscale, ndesc = _scale_integrals(Ireal, nH)
    df["nH"] = df.integral * kscale
    df["nH_int"] = np.where(df.nH >= 0.5, np.round(df.nH), np.round(df.nH, 1))
    df["amostra"] = real
    # rótulo de impureza só faz sentido para sinais pequenos
    big = df.nH >= 0.5
    df.loc[big & df.rotulo.str.startswith("possível"), "rotulo"] = ""
    # satélites de ¹³C (≈0,55% de cada lado, a ±¹J_CH/2 ≈ 60-85 Hz de um sinal grande)
    for i, r in df.iterrows():
        if r.nH < 0.1 and not r.rotulo:
            for j, q in df.iterrows():
                d = abs(r.ppm - q.ppm) * s.sf
                if q.nH >= 1 and 55 <= d <= 90 and r.integral < 0.012 * q.integral:
                    df.at[i, "rotulo"] = f"satélite ¹³C de {q.ppm:.2f}"
                    break
    # descarta lixo: largos minúsculos e picos insignificantes sem rótulo
    junk = (df.nH < 0.15) & (df.rotulo == "") & df.nota.str.contains("largo")
    junk |= df.nota.str.contains("largo") & (df.fwhm_hz > 150)        # deriva de linha de base
    junk |= (df.nH < 0.03) & (df.rotulo == "")
    df = df[~junk].reset_index(drop=True)
    real = df.amostra.values
    # ambientes prováveis
    small = (df.nH < 0.5) & (df.rotulo == "")
    df.loc[small, "nota"] = (df.loc[small, "nota"].astype(str) + " sinal pequeno (impureza?)").str.strip()
    df["amostra"] = df.amostra & (df.rotulo == "")
    df["ambientes"] = [", ".join(ranges_for("1H", p)[:2]) for p in df.ppm]
    # acoplamentos compartilhados
    partners = []
    for i, r in df.iterrows():
        ps = []
        for j, q in df.iterrows():
            if i == j or not r.J or not q.J or not r.amostra or not q.amostra:
                continue
            for J in r.J:
                for Jq in q.J:
                    if abs(J - Jq) <= max(0.4, tol):
                        ps.append(f"{q.ppm:.2f} ({J:.1f})")
        partners.append("; ".join(sorted(set(ps))))
    df["parceiros_J"] = partners
    tot = float(df.nH_int[df.amostra & (df.nH >= 0.5)].sum())
    warn = []
    if nH and abs(tot - nH) >= 1:
        warn.append(f"soma das integrais ({tot:.0f} H) ≠ fórmula ({nH} H): H trocáveis (OH/NH/COOH) largos/ausentes, "
                    "sobreposição com solvente/água ou integrais imprecisas")
    frac = df.nH[df.amostra] - np.round(df.nH[df.amostra])
    if (np.abs(frac[df.nH[df.amostra] >= 0.8]) > 0.25).any():
        warn.append("integrais longe de inteiros: mistura, impureza, relaxação incompleta (d1 curto) ou linha de base")
    return dict(table=df, sigma=sigma, norm_desc=ndesc, nH=nH, total_H=tot, warnings=warn, tol_hz=tol)


def fmt_j(J):
    return ", ".join(f"{j:.1f}" for j in J)


def acs_1h(res, s, decimal=".") -> str:
    df = res["table"]
    sol = SOLVENTS.get(s.solvent or "", {}).get("name", s.solvent or "?")
    parts = []
    for r in df[df.amostra & (df.nH >= 0.5)].itertuples():
        if r.mult in ("m",) or (r.mult not in ("s", "sl") and not r.J):
            d = f"{r.ppm_max:.2f}–{r.ppm_min:.2f}"
            body = f"m, {r.nH_int:.0f}H"
        elif r.mult == "AB" and r.AB:
            d = f"{r.AB[0] / s.sf:.2f}, {r.AB[1] / s.sf:.2f}"
            body = f"AB, J = {fmt_j(r.J)} Hz, {r.nH_int:.0f}H"
        else:
            d = f"{r.ppm:.2f}"
            body = f"{r.mult}" + (f", J = {fmt_j(r.J)} Hz" if r.J else "") + f", {r.nH_int:.0f}H"
        parts.append(f"{d} ({body})")
    txt = f"¹H NMR ({s.sf:.0f} MHz, {sol}) δ " + ", ".join(parts) + "."
    if decimal == ",":
        txt = txt.replace(".", "#").replace(",", ";").replace("#", ",").replace(";", ",")
    return txt


# ----------------------------------------------------------------------------- ¹³C e DEPT

def analyze_13c(s, k: float = 5.0, dept135=None, dept90=None, apt=None, sigma=None, tol=0.12) -> dict:
    """Picos de ¹³C (sem solvente), rótulos de impurezas e multiplicidade por DEPT/APT:
    DEPT-90 positivo → CH; DEPT-135 positivo e ausente no DEPT-90 → CH₃; DEPT-135 negativo → CH₂;
    ausente nos DEPT → C (quaternário). APT: C/CH₂ de um sinal, CH/CH₃ do outro."""
    sigma = sigma or noise(s.y)
    pk = pick_peaks(s, k=k, sigma=sigma)
    if len(pk) > 1:
        # funde máximos duplicados da mesma linha (< 2,5 Hz)
        pk = pk.sort_values("ppm", ascending=False).reset_index(drop=True)
        keep = np.ones(len(pk), bool)
        for i in range(1, len(pk)):
            j = i - 1
            while j >= 0 and not keep[j]:
                j -= 1
            if j >= 0 and abs(pk.ppm[i] - pk.ppm[j]) * s.sf < 2.5:
                if pk.altura[i] > pk.altura[j]:
                    keep[j] = False
                else:
                    keep[i] = False
        pk = pk[keep].reset_index(drop=True)
    if len(pk):
        # remove "ombros"/oscilações: pico < 4% de um vizinho a menos de 0,6 ppm
        keep = np.ones(len(pk), bool)
        h = pk.altura.values
        for i in range(len(pk)):
            near = np.abs(pk.ppm.values - pk.ppm.values[i]) < 0.6
            if h[i] < 0.04 * h[near].max():
                keep[i] = False
        pk = pk[keep].reset_index(drop=True)
    if len(pk):
        # remove picos de solvente
        mask = solvent_mask(s)
        sol = np.zeros(len(pk), bool)
        for lo, hi, _ in mask:
            sol |= (pk.ppm >= lo) & (pk.ppm <= hi)
        pk["rotulo"] = [label_impurity(s, p) or "" for p in pk.ppm]
        pk.loc[sol, "rotulo"] = "solvente"
    pk["mult"] = ""
    info = {}
    if (dept135 is not None or dept90 is not None) and len(pk):
        def dpk(sp):
            if sp is None:
                return None
            d = pick_peaks(sp, k=max(k, 8.0), negative=True)
            if len(d):
                big = np.sort(np.abs(d.altura.values))[::-1][:8]
                d = d[np.abs(d.altura) >= 0.08 * np.median(big)].reset_index(drop=True)
            return d
        d135, d90 = dpk(dept135), dpk(dept90)

        def near(df, p, sign=None):
            if df is None or not len(df):
                return None
            d = np.abs(df.ppm.values - p)
            j = int(np.argmin(d))
            if d[j] > tol:
                return None
            if sign is not None and np.sign(df.altura.values[j]) != sign:
                return None
            return df.iloc[j]
        thr90 = None
        if d90 is not None and len(d90):
            thr90 = 0.15 * d90.altura.abs().max()
        mults = []
        for p in pk.ppm:
            m135 = near(d135, p)
            m90 = near(d90, p)
            is90 = m90 is not None and m90.altura > 0 and (thr90 is None or m90.altura > thr90)
            if m135 is None and not is90:
                mults.append("C")
            elif is90:
                mults.append("CH")
            elif m135 is not None and m135.altura < 0:
                mults.append("CH2")
            elif m135 is not None and m135.altura > 0:
                mults.append("CH3" if d90 is not None else "CH/CH3")
            else:
                mults.append("?")
        pk["mult"] = mults
        info["dept"] = True
        # pico dentro da janela do solvente mas com sinal no DEPT é carbono da amostra
        sol_dept = (pk.rotulo == "solvente") & pk.mult.isin(["CH", "CH2", "CH3", "CH/CH3"])
        pk.loc[sol_dept, "rotulo"] = ""
        pk.loc[(pk.rotulo == "solvente"), "mult"] = ""
    elif apt is not None and len(pk):
        da = pick_peaks(apt, k=k, negative=True)
        mults = []
        for p in pk.ppm:
            d = np.abs(da.ppm.values - p) if len(da) else np.array([9e9])
            j = int(np.argmin(d))
            if d[j] > tol:
                mults.append("?")
            else:
                mults.append("C/CH2" if da.altura.values[j] > 0 else "CH/CH3")
        pk["mult"] = mults
        info["apt"] = "convenção: C/CH₂ positivo (inverta se o seu APT usa a outra fase)"
    if len(pk):
        pk["ambientes"] = [", ".join(ranges_for("13C", p)[:2]) for p in pk.ppm]
        samp0 = pk.rotulo.eq("") | pk.rotulo.str.startswith("possível")
        if samp0.sum() > 3:
            ref_h = np.median(np.sort(np.abs(pk.altura[samp0].values))[::-1][:max(3, samp0.sum() // 2)])
            tiny = samp0 & (np.abs(pk.altura) < 0.03 * ref_h)
            pk.loc[tiny & pk.rotulo.eq(""), "rotulo"] = "pequeno (impureza/ruído?)"
    samp = pk[pk.rotulo.eq("")] if len(pk) else pk
    counts = samp["mult"].value_counts().to_dict() if len(samp) else {}
    nH_dept = None
    if info.get("dept"):
        nH_dept = int(sum({"CH3": 3, "CH2": 2, "CH": 1}.get(m, 0) for m in samp["mult"]))
    return dict(table=pk, sigma=sigma, counts=counts, nH_CH=nH_dept, **info)


def acs_13c(res, s) -> str:
    df = res["table"]
    sol = SOLVENTS.get(s.solvent or "", {}).get("name", s.solvent or "?")
    df = df[df.rotulo.eq("")]
    if df.mult.astype(bool).any():
        items = [f"{p:.1f} ({m})" for p, m in zip(df.ppm, df.mult)]
    else:
        items = [f"{p:.1f}" for p in df.ppm]
    return f"¹³C NMR ({s.sf:.0f} MHz, {sol}) δ " + ", ".join(items) + "."


def inventory(r1=None, r13=None, formula=None) -> list:
    """Resumo estrutural a partir dos dados (sem propor estrutura)."""
    from nmr_data import dbe
    L = []
    if formula:
        L.append(f"Fórmula {formula}: IDH (anéis + ligações π) = {dbe(formula):g}")
    if r1 is not None and len(r1["table"]):
        t = r1["table"]
        t = t[t.amostra & (t.nH >= 0.5)]
        def h(lo, hi):
            return float(t.nH_int[(t.ppm >= lo) & (t.ppm < hi)].sum())
        L.append(f"¹H: {h(6.0, 9.5):.0f} H aromáticos/heteroaromáticos/vinílicos conjugados (6,0-9,5); "
                 f"{h(4.5, 6.0):.0f} H em 4,5-6,0 (vinílicos, anoméricos, O–CH–O); {h(3.0, 4.5):.0f} H em 3,0-4,5 "
                 f"(X–CH: O, N, halogênio); {h(0.0, 3.0):.0f} H alifáticos (0-3,0); {h(9.5, 16):.0f} H > 9,5 "
                 "(CHO, COOH, OH quelado, NH)")
        sing = t[(t.mult == "s") & (t.nH_int == 3)]
        if len(sing):
            L.append("CH₃ isolados (s, 3H): " + ", ".join(f"{p:.2f}" for p in sing.ppm) +
                     " → OCH₃ (3,3-4,0), NCH₃ (2,2-3,1), CH₃CO/ArCH₃ (2,0-2,6), C–CH₃ quaternário (0,8-1,6)")
        trip = t[(t.mult == "t") & (t.nH_int == 3)]
        quad = t[(t.mult == "q") & (t.nH_int == 2)]
        for a in trip.itertuples():
            for b in quad.itertuples():
                if a.J and b.J and abs(a.J[0] - b.J[0]) < 0.5:
                    L.append(f"etila: t {a.ppm:.2f} (3H) + q {b.ppm:.2f} (2H), J = {a.J[0]:.1f} Hz"
                             + (" → OCH₂CH₃ (éster/éter)" if b.ppm > 3.3 else " → CH₂CH₃ ligado a C/N"))
    if r13 is not None and len(r13["table"]):
        t = r13["table"]
        t = t[t.rotulo.eq("") | t.rotulo.str.startswith("possível")]
        n = len(t)
        def c(lo, hi):
            return int(((t.ppm >= lo) & (t.ppm < hi)).sum())
        L.append(f"¹³C: {n} sinais; {c(160, 230)} C=O/C=N (>160; cetona/aldeído > 190: {c(190, 230)}); "
                 f"{c(100, 160)} sp² (100-160); {c(60, 100)} C–O sp³/anoméricos (60-100); {c(0, 60)} sp³ (< 60)")
        if r13.get("counts"):
            L.append("Multiplicidades (DEPT): " + ", ".join(f"{k} {v}" for k, v in r13["counts"].items()))
        if r13.get("nH_CH") is not None:
            L.append(f"H ligados a C (DEPT): {r13['nH_CH']}")
    return L
