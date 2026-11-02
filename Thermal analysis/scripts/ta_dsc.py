"""
ta_dsc.py - DSC e DTA: eventos, transição vítrea, integração de picos, cristalinidade, OIT, Cp

    from ta_dsc import find_events, glass_transition, integrate_peak, crystallinity, oit
    seg = select(r, "heat2")                         # (ta_process) segmento a analisar
    ev = find_events(seg)                            # DataFrame de picos endo/exo e degraus (candidatos)
    tg = glass_transition(seg, (60, 100))            # Tg onset/midpoint/inflexão/endset, ΔCp
    pk = integrate_peak(seg, (220, 265))             # ΔH (J/g), T_onset, T_pico, T_endset
    crystallinity(pk["dH"], "PET", dHcc=...)         # Xc (%)

Convenções internas: o sinal é reorientado para endo POSITIVO (ta_process.endo_signal); os
resultados dizem 'endo'/'exo' explicitamente e ΔH é reportado em módulo com o tipo.
Integração em TEMPO (∫ q dt): q em W/g -> ΔH em J/g. Sem canal de tempo, usa dt = dT/β.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter

from ta_io import Run
from ta_process import endo_signal, segments, heating_rate, _noise

# ΔH de fusão do polímero 100% cristalino (J/g) - valores mais citados; a literatura varia, cite a fonte
DH0 = {"PE": 293.0, "HDPE": 293.0, "LDPE": 293.0, "UHMWPE": 293.0, "PP": 207.0, "iPP": 207.0,
       "PET": 140.0, "PBT": 140.0, "PEN": 103.0, "PLA": 93.0, "PLLA": 93.0, "PCL": 139.5, "PHB": 146.0,
       "PA6": 230.0, "PA66": 196.0, "PA11": 206.0, "PA12": 209.0, "PEEK": 130.0, "PPS": 112.0,
       "PVDF": 104.7, "PTFE": 82.0, "POM": 326.0, "PEO": 196.8, "PEG": 196.8, "iPS": 86.0, "sPS": 82.0,
       "PBS": 110.3, "PVA": 138.6, "PEKK": 130.0, "PBAT": 114.0, "PPC": None}


def _seg_info(r: Run):
    s = r.meta.get("segment_info")
    if s:
        return s
    segs = segments(r)
    return segs[0] if len(segs) == 1 else {"type": "heat" if np.nanmedian(np.diff(r.T)) >= 0 else "cool"}


def _xy(r: Run, override=None):
    T, t, y, u, d, src = endo_signal(r, override)
    ok = np.isfinite(T) & np.isfinite(y)
    if t is not None:
        ok &= np.isfinite(t)
    T, y = T[ok], y[ok]
    t = t[ok] if t is not None else None
    return T, t, y, u, d, src


def _dt_seconds(r: Run, T, t):
    if t is not None:
        return t * 60.0
    beta = heating_rate(r) or 10.0
    return (T - T[0]) / beta * 60.0


def _uniform(x, y, step=None):
    """Reamostra em grade uniforme crescente (remove x repetidos - comum em T de DSC)."""
    o = np.argsort(x, kind="stable")
    x, y = np.asarray(x, float)[o], np.asarray(y, float)[o]
    ux, inv = np.unique(x, return_inverse=True)
    uy = np.bincount(inv, weights=y) / np.bincount(inv)
    if ux.size < 3:
        return ux, uy
    st = step or float(np.median(np.diff(ux)))
    g = np.arange(ux[0], ux[-1] + st / 2, st)
    return g, np.interp(g, ux, uy)


def _lin(x, y):
    if len(x) < 2:
        return np.array([0.0, float(np.mean(y)) if len(y) else 0.0])
    return np.polyfit(x, y, 1)


# ============================================================================ eventos (candidatos)
def _robust_poly(x, y, deg=3, it=6, k=2.5):
    """Linha de base polinomial robusta (exclui iterativamente pontos fora de k·σ)."""
    m = np.ones_like(x, bool)
    xc = (x - x.mean()) / (x.std() or 1)
    c = np.polyfit(xc, y, deg)
    for _ in range(it):
        r = y - np.polyval(c, xc)
        s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))) + 1e-15
        m_new = np.abs(r - np.median(r[m])) < k * s
        if m_new.sum() < deg + 5 or (m_new == m).all():
            break
        m = m_new
        c = np.polyfit(xc[m], y[m], deg)
    return np.polyval(c, xc)


def find_events(r: Run, endo: str | None = None, rel_prominence: float = 0.05, min_width_K: float = 1.0,
                smooth_K: float = 2.0, edge_K: float | None = None) -> pd.DataFrame:
    """Varre um segmento e lista candidatos: picos endo/exo (resíduo sobre linha de base
    polinomial robusta) e degraus (Tg) pela derivada. São SUGESTÕES de janela; quantifique com
    integrate_peak / glass_transition passando a janela conferida no gráfico."""
    T, t, y, u, d, src = _xy(r, endo)
    if T.size < 30:
        return pd.DataFrame()
    info = _seg_info(r)
    typ = info.get("type", "heat")
    oo = np.argsort(T)
    Ts, ys = T[oo], y[oo]
    step = float(np.median(np.diff(Ts))) or 0.1
    span_T = float(Ts[-1] - Ts[0])
    edge = edge_K if edge_K is not None else max(5.0, 0.03 * span_T)
    inner = (Ts > Ts[0] + edge) & (Ts < Ts[-1] - edge)
    w = max(5, int(round(smooth_K / step)) | 1)
    w = min(w, (Ts.size - 2) | 1)
    ysm = savgol_filter(ys, w, 2)
    base = np.full_like(ysm, np.nan)
    base[inner] = _robust_poly(Ts[inner], ysm[inner], deg=3)
    res = ysm - base
    nz = _noise(ys)
    rmax = float(np.nanmax(np.abs(res[inner]))) if inner.any() else 0.0
    prom = max(rel_prominence * rmax, 8 * nz)
    rows = []
    rr = np.nan_to_num(res)
    for sign, kind in ((1, "endo"), (-1, "exo")):
        pk, pr = find_peaks(sign * rr, prominence=prom, width=max(1, int(min_width_K / step)))
        for i, p in enumerate(pk):
            if not inner[p]:
                continue
            h = sign * rr[p]
            a = p
            while a > 0 and sign * rr[a] > 0.05 * h and inner[a]:
                a -= 1
            b = p
            while b < len(rr) - 1 and sign * rr[b] > 0.05 * h and inner[b]:
                b += 1
            rows.append(dict(tipo=kind, evento="pico", T_pico=float(Ts[p]), T_ini=float(Ts[a]), T_fim=float(Ts[b]),
                             altura=float(h), unidade=u))
    # degraus (Tg): pico da derivada SEM lobo de sinal oposto ao lado (um pico gera um par +/−
    # na derivada; um degrau, um lobo só), com níveis por retas locais dos dois lados
    wd = min(max(w, int(round(4 / step)) | 1), (Ts.size - 2) | 1)
    dy = savgol_filter(ys, wd, 2, deriv=1, delta=step)
    trend = np.gradient(np.nan_to_num(base, nan=np.nanmedian(ysm)), Ts)
    dyb = np.where(inner, dy - trend, 0.0)
    k = max(4, int(round(8 / step)))
    dn = _noise(dyb[inner]) if inner.sum() > 10 else 1e-15
    cand, cpr = find_peaks(dyb, distance=k, width=int(1.5 / step), prominence=6 * dn * np.sqrt(wd))
    for ci, p in enumerate(cand):
        if not inner[p] or p - 3 * k < 0 or p + 3 * k >= len(Ts):
            continue
        if cpr["widths"][ci] * step > 25:
            continue                      # derivada larga demais: curvatura de linha de base, não Tg
        Tp = float(Ts[p])
        win = dyb[max(0, p - 3 * k):p + 3 * k]
        if np.min(win) < -0.35 * dyb[p]:
            continue                      # lobo oposto: flanco de pico, não degrau
        L = slice(p - 3 * k, p - k)
        R = slice(p + k, p + 3 * k)
        cl, cr = np.polyfit(Ts[L], ysm[L], 1), np.polyfit(Ts[R], ysm[R], 1)
        jump = float(np.polyval(cr, Tp) - np.polyval(cl, Tp))
        if jump <= 0 or jump < max(10 * nz, 0.01 * rmax):
            continue
        if any(e["evento"] == "pico" and e["T_ini"] <= Tp <= e["T_fim"] and e["altura"] > jump for e in rows):
            continue                      # dentro de um pico maior
        rows = [e for e in rows if not (e["evento"] == "pico" and abs(e["T_pico"] - Tp) < 2 * k * step
                                        and e["altura"] < 1.5 * jump)]
        rows.append(dict(tipo="endo", evento="degrau (Tg?)", T_pico=Tp,
                         T_ini=float(Ts[p - 2 * k]), T_fim=float(Ts[p + 2 * k]), altura=jump, unidade=u))
    # picos pequenos encostados na borda útil = curvatura de linha de base; grandes = transiente
    lo_in, hi_in = Ts[0] + edge, Ts[-1] - edge
    keep = []
    for e in rows:
        if e["evento"] == "pico" and (e["T_fim"] - e["T_ini"]) > 0.35 * span_T:
            continue                      # "pico" mais largo que 1/3 do segmento = linha de base curva
        touches = e["T_ini"] <= lo_in + step or e["T_fim"] >= hi_in - step
        if e["evento"] == "pico" and touches:
            if e["altura"] < 0.3 * rmax:
                continue
            e["evento"] = "borda (transiente?)"
        keep.append(e)
    rows = keep
    # pico endo pequeno com mudança de nível através dele = Tg com relaxação entálpica
    for e in rows:
        if e["evento"] == "pico" and e["tipo"] == "endo" and e["altura"] < 0.35 * rmax:
            i0 = int(np.searchsorted(Ts, e["T_ini"]))
            i1 = int(np.searchsorted(Ts, e["T_fim"]))
            kk = max(4, int(round(6 / step)))
            if i0 - kk >= 0 and i1 + kk < len(Ts):
                lvl = float(np.median(ysm[i1:i1 + kk] - base[i1:i1 + kk]) - np.median(ysm[i0 - kk:i0] - base[i0 - kk:i0])) \
                    if np.isfinite(base[i0 - kk:i1 + kk]).all() else 0.0
                raw = float(np.median(ysm[i1:i1 + kk]) - np.median(ysm[i0 - kk:i0]))
                if raw > 0.25 * e["altura"]:
                    e["evento"] = "pico+degrau"
                    e["_tg"] = True
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values("T_pico").reset_index(drop=True)
        def guess(row):
            if row["evento"].startswith("borda"):
                return "início/fim do segmento: provável transiente instrumental - não quantificado"
            if row["evento"] == "pico+degrau":
                return "transição vítrea com relaxação entálpica (ou fusão sobre degrau)"
            if row["evento"] != "pico":
                return "transição vítrea" if row["tipo"] == "endo" else "degrau de linha de base (exo?)"
            if typ == "cool":
                return "cristalização" if row["tipo"] == "exo" else "endo no resfriamento (raro: confira)"
            return "fusão / evaporação / desidratação / desnaturação" if row["tipo"] == "endo" else \
                "cristalização a frio / cura / oxidação / decomposição"
        df["sugestao"] = df.apply(guess, axis=1)
        for c in ("T_pico", "T_ini", "T_fim"):
            df[c] = df[c].round(1)
        df["altura"] = df["altura"].round(4)
        df = df.drop(columns=[c for c in df.columns if c.startswith("_")])
    df.attrs.update(endo=d, endo_source=src, unit=u, noise=nz)
    return df


# ============================================================================ transição vítrea
def glass_transition(r: Run, window, endo: str | None = None, fit_K: float | None = None,
                     smooth_K: float = 2.0) -> dict:
    """Tg pelo método das tangentes (ISO 11357-2 / ASTM E1356) numa janela (T1, T2) que contenha
    linha de base dos dois lados do degrau. Retorna T_onset, T_mid (meia altura entre as linhas
    de base extrapoladas), T_inflexao, T_endset, ΔCp (J/(g·K), se q em W/g e β conhecido) e flag de
    relaxação entálpica (pico de overshoot)."""
    T, t, y, u, d, src = _xy(r, endo)
    lo, hi = min(window), max(window)
    m = (T >= lo) & (T <= hi)
    if m.sum() < 20:
        raise ValueError(f"poucos pontos na janela {window}")
    o = np.argsort(T[m])
    x, yy = _uniform(T[m][o], y[m][o])
    step = float(np.median(np.diff(x))) or 0.1
    w = max(5, int(round(smooth_K / step)) | 1)
    ys = savgol_filter(yy, min(w, (x.size - 2) | 1), 2)
    fk = fit_K or max(3.0, 0.15 * (hi - lo))
    pre, post = x <= lo + fk, x >= hi - fk
    c1, c2 = _lin(x[pre], ys[pre]), _lin(x[post], ys[post])
    info = _seg_info(r)
    sgn = 1.0 if info.get("type") != "cool" else -1.0      # aquecimento: degrau endo (+)
    dy = np.gradient(ys, x)
    core = (x > lo + fk) & (x < hi - fk)
    if core.sum() < 3:
        core = np.ones_like(x, bool)
    ii = np.where(core)[0]
    j = ii[int(np.argmax(sgn * dy[ii]))]
    Ti, yi, si = float(x[j]), float(ys[j]), float(dy[j])
    def inter(c):
        den = si - c[0]
        return float((c[1] - (yi - si * Ti)) / den) if abs(den) > 1e-15 else np.nan
    T_on, T_end = inter(c1), inter(c2)
    mid = (np.polyval(c1, x) + np.polyval(c2, x)) / 2
    diff = ys - mid
    T_mid = np.nan
    for k in range(1, len(x)):
        if (x[k] >= min(T_on, T_end) - 5) and np.sign(diff[k]) != np.sign(diff[k - 1]) and diff[k - 1] != 0:
            T_mid = float(x[k - 1] + (x[k] - x[k - 1]) * (-diff[k - 1]) / (diff[k] - diff[k - 1]))
            break
    step_h = float(np.polyval(c2, T_mid if np.isfinite(T_mid) else Ti) - np.polyval(c1, T_mid if np.isfinite(T_mid) else Ti))
    beta = abs(info.get("rate") or heating_rate(r) or np.nan)
    dcp = step_h * sgn / (beta / 60.0) if u == "W/g" and np.isfinite(beta) else np.nan
    # overshoot de relaxação entálpica: sinal acima da linha de base pós-Tg logo após o degrau
    after = (x > Ti) & (x < hi - fk)
    over = float(np.nanmax(sgn * (ys[after] - np.polyval(c2, x[after])))) if after.any() else 0.0
    relax = bool(over > 0.3 * abs(step_h)) if step_h else False
    warn = []
    if not (lo < T_on < T_end < hi):
        warn.append("Onset/endset fora da janela: amplie a janela para incluir linha de base dos dois lados.")
    if step_h * sgn <= 0:
        warn.append("Degrau no sentido exo: confira a direção endo ou se isto é mesmo uma Tg.")
    if relax:
        warn.append("Pico de relaxação entálpica sobre a Tg (envelhecimento físico): prefira o 2º aquecimento "
                    "ou reporte Tg pelo método de Moynihan/área.")
    return dict(T_onset=round(T_on, 1), T_mid=round(T_mid, 1) if np.isfinite(T_mid) else np.nan,
                T_inflexao=round(Ti, 1), T_endset=round(T_end, 1), degrau=step_h, unidade=u,
                dCp=round(dcp, 3) if np.isfinite(dcp) else np.nan, beta=beta, relaxacao=relax,
                janela=(lo, hi), endo=d, endo_source=src, avisos=warn,
                _plot=dict(c1=c1, c2=c2, Ti=Ti, yi=yi, si=si, x=x, sgn=sgn))


# ============================================================================ integração de picos
def integrate_peak(r: Run, window, endo: str | None = None, baseline: str = "linear", fit_K: float = 3.0,
                   smooth_K: float = 1.0, cal_factor: float | None = None) -> dict:
    """Integra um pico na janela (T1, T2). baseline: 'linear' (reta entre as bordas) ou
    'sigmoidal' (ISO 11357-1: linha de base que muda de nível em proporção à área parcial -
    use quando Cp muda no evento, p.ex. fusão com degrau). ΔH = ∫(q − b) dt (J/g se q em W/g).
    Para DTA (q em µV ou K) informe cal_factor (J/g por unidade·s) para obter J/g.
    Retorna tipo (endo/exo), dH (módulo), T_onset, T_pico, T_endset, altura, FWHM."""
    T, t, y, u, d, src = _xy(r, endo)
    lo, hi = min(window), max(window)
    m = (T >= lo) & (T <= hi)
    if m.sum() < 10:
        raise ValueError(f"poucos pontos na janela {window}")
    idx = np.where(m)[0]
    x, yy = T[idx], y[idx]
    ts = _dt_seconds(r, T, t)[idx]
    order = np.argsort(ts, kind="stable")
    x, yy, ts = x[order], yy[order], ts[order]
    keep = np.concatenate([[True], np.diff(ts) > 0])
    x, yy, ts = x[keep], yy[keep], ts[keep]
    step = float(np.median(np.abs(np.diff(x)))) or 0.1
    w = max(5, int(round(smooth_K / step)) | 1)
    ys = savgol_filter(yy, min(w, (x.size - 2) | 1), 2) if x.size > 7 else yy
    k = max(2, int(round(fit_K / step)))
    # direção dominante do pico em relação à corda
    chord = np.interp(ts, [ts[0], ts[-1]], [np.mean(ys[:k]), np.mean(ys[-k:])])
    sgn = 1.0 if np.nanmax(ys - chord) >= -np.nanmin(ys - chord) else -1.0
    kind = "endo" if sgn > 0 else "exo"
    if baseline == "linear":
        b = chord
    elif baseline == "sigmoidal":
        c1 = _lin(ts[:k], ys[:k])
        c2 = _lin(ts[-k:], ys[-k:])
        b = chord.copy()
        for _ in range(8):
            diff = sgn * (ys - b)
            diff = np.clip(diff, 0, None)
            cum = np.concatenate([[0], np.cumsum(0.5 * (diff[1:] + diff[:-1]) * np.diff(ts))])
            a = cum / cum[-1] if cum[-1] > 0 else np.linspace(0, 1, len(ts))
            b = (1 - a) * np.polyval(c1, ts) + a * np.polyval(c2, ts)
    else:
        raise ValueError(baseline)
    net = sgn * (ys - b)
    area = float(np.trapezoid(net, ts))                     # (unidade)·s
    p = int(np.argmax(net))
    Tp = float(x[p])
    # tangentes nos pontos de inflexão das bordas (no tempo, para valer em aquecimento e resfriamento)
    dn = np.gradient(net, ts)
    lead = slice(0, p + 1)
    trail = slice(p, len(ts))
    def tangent_cross(sl, use_max):
        seg = dn[sl]
        if seg.size < 2:
            return np.nan
        j = sl.start + (int(np.argmax(seg)) if use_max else int(np.argmin(seg)))
        s_, y0, t0 = dn[j], net[j], ts[j]
        if s_ == 0:
            return np.nan
        tc = t0 - y0 / s_                                   # cruza net=0 (a linha de base)
        return float(np.interp(tc, ts, x, left=np.nan, right=np.nan)) if ts[0] - 60 <= tc <= ts[-1] + 60 else np.nan
    T_on = tangent_cross(lead, True)
    T_end = tangent_cross(trail, False)
    half = net >= 0.5 * net[p]
    hw = np.where(half)[0]
    fwhm = float(abs(x[hw[-1]] - x[hw[0]])) if hw.size else np.nan
    dH = area
    unit_H = {"W/g": "J/g", "mW": "mJ", "µV/mg": "µV·s/mg", "µV": "µV·s", "K": "K·s"}.get(u, u + "·s")
    if cal_factor:
        dH, unit_H = area * cal_factor, "J/g"
    warn = []
    if net[0] > 0.15 * net[p] or net[-1] > 0.15 * net[p]:
        warn.append("O pico não volta à linha de base nas bordas da janela: amplie a janela.")
    if (sgn * (ys - b) < -0.1 * net[p]).sum() > 0.1 * len(net):
        warn.append("Sinal cruza a linha de base dentro da janela: eventos sobrepostos ou linha de base inadequada.")
    return dict(tipo=kind, dH=round(dH, 3), unidade_dH=unit_H, T_onset=round(T_on, 1), T_pico=round(Tp, 1),
                T_endset=round(T_end, 1), altura=float(net[p]), unidade=u, FWHM=round(fwhm, 1),
                janela=(lo, hi), linha_base=baseline, endo=d, endo_source=src, avisos=warn,
                _plot=dict(x=x, y=ys * (1 if d == "up" else -1), b=b * (1 if d == "up" else -1)))


def crystallinity(dHm: float, polymer: str | float, dHcc: float = 0.0, w: float = 1.0) -> dict:
    """Xc (%) = (ΔHm − ΔHcc) / (ΔH0 · w) × 100. polymer: nome ('PET', 'PLA', 'PP'...) ou ΔH0 (J/g).
    w = fração mássica do polímero (compósitos/blendas). ΔHcc = cristalização a frio no mesmo aquecimento."""
    dH0 = DH0.get(str(polymer).upper().replace(" ", ""), None) if not isinstance(polymer, (int, float)) else float(polymer)
    if not dH0:
        raise ValueError(f"ΔH0 desconhecido para {polymer!r}; informe o valor (J/g). Tabela: "
                         + ", ".join(f"{k} {v}" for k, v in DH0.items() if v))
    xc = (abs(dHm) - abs(dHcc)) / (dH0 * w) * 100
    return dict(Xc=round(xc, 1), dH0=dH0, dHm=abs(dHm), dHcc=abs(dHcc), w=w)


# ============================================================================ OIT
def gas_switch_time(r: Run):
    """Tempo (min) da troca de gás: maior variação em canais de vazão/gás (NETZSCH purge, TA gas)."""
    best = None
    for k, v in r.extra.items():
        if not any(s in k for s in ("flow", "purge", "gas", "mfc")) or not np.isfinite(v).any():
            continue
        dv = np.abs(np.diff(np.nan_to_num(v)))
        if dv.size and dv.max() > 2:
            i = int(np.argmax(dv)) + 1
            if best is None or dv.max() > best[1]:
                best = (float(r.time[i]), float(dv.max()), k)
    return best


def oit(r: Run, t_switch: float | None = None, endo: str | None = None, baseline_min: float = 2.0,
        threshold: float | None = None) -> dict:
    """Tempo de indução oxidativa (ISO 11357-6 / ASTM D3895): t_onset(exotérmico) − t_troca_de_gás,
    onset pela tangente no ponto de maior inclinação cruzando a linha de base isotérmica.
    t_switch: tempo (min) da troca N2->O2; se None, tenta achar nos canais de vazão; senão usa o
    início do último segmento isotérmico."""
    T, t, y, u, d, src = _xy(r, endo)
    if t is None:
        raise ValueError("OIT exige canal de tempo")
    src_sw = "informado"
    if t_switch is None:
        g = gas_switch_time(r)
        if g:
            t_switch, src_sw = g[0], f"canal {g[2]}"
        else:
            isos = [s for s in segments(r) if s["type"] == "iso"]
            if not isos:
                raise ValueError("sem segmento isotérmico; informe t_switch")
            t_switch, src_sw = isos[-1]["t0"], "início da isoterma (sem canal de gás)"
    m = t >= t_switch
    tt, yy = t[m], -y[m]                        # exo positivo
    w = max(5, int(round(0.2 / max(np.median(np.diff(tt)), 1e-6))) | 1)
    ys = savgol_filter(yy, min(w, (tt.size - 2) | 1), 2)
    b_m = (tt >= t_switch + 0.3) & (tt <= t_switch + 0.3 + baseline_min)
    base = float(np.median(ys[b_m])) if b_m.any() else float(ys[0])
    nz = _noise(yy)
    rise = ys - base
    big = np.where(rise > max(20 * nz, 0.1 * np.nanmax(rise)))[0]
    if big.size == 0:
        return dict(OIT=np.nan, t_switch=t_switch, fonte_troca=src_sw, aviso="nenhuma exoterma detectada")
    j0 = big[0]
    # até o primeiro máximo da exoterma (ou o máximo global)
    pk, _ = find_peaks(ys[j0:], prominence=max(10 * nz, 0.2 * np.nanmax(rise)))
    jend = j0 + (int(pk[0]) if pk.size else int(np.argmax(ys[j0:])))
    while j0 > 0 and rise[j0] > 3 * nz:
        j0 -= 1
    seg = slice(j0, max(jend, j0 + 3))
    ds = np.gradient(ys, tt)
    j = seg.start + int(np.argmax(ds[seg]))
    t_on = float(tt[j] - (ys[j] - base) / ds[j]) if ds[j] > 0 else float(tt[j0])
    out = dict(OIT=round(t_on - t_switch, 2), t_onset=round(t_on, 2), t_switch=round(t_switch, 2),
               fonte_troca=src_sw, T_iso=round(float(np.nanmedian(T[m])), 1), endo=d, endo_source=src)
    if threshold is not None:
        k = np.where(rise > threshold)[0]
        out["OIT_limiar"] = round(float(tt[k[0]] - t_switch), 2) if k.size else np.nan
    return out


# ============================================================================ Cp (safira, ASTM E1269)
# Cp da safira (α-Al2O3), J/(g·K) - aproximado (tabelas JANAF); use o certificado do seu padrão
SAPPHIRE_CP = np.array([[0, 0.718], [25, 0.775], [50, 0.825], [100, 0.880], [150, 0.928], [200, 0.966],
                        [250, 0.997], [300, 1.023], [350, 1.044], [400, 1.061], [450, 1.076], [500, 1.088],
                        [600, 1.108], [700, 1.123], [800, 1.136], [900, 1.147], [1000, 1.157]])


def cp_sapphire(sample: Run, blank: Run, sapphire: Run, m_sample: float | None = None,
                m_sapphire: float | None = None, T_range=None, step: float = 1.0) -> pd.DataFrame:
    """Cp(T) pelo método de três corridas (branco, safira, amostra; mesmo programa e cadinhos):
    Cp = Cp_safira · (q_amostra − q_branco)/(q_safira − q_branco) · m_safira/m_amostra.
    Use o sinal em mW (não normalizado) dos três arquivos."""
    ms = m_sample or sample.sample_mass
    mz = m_sapphire or sapphire.sample_mass
    if not (ms and mz):
        raise ValueError("informe as massas da amostra e da safira (mg)")
    def sig(r):
        y = r.hf.astype(float)
        if r.units.get("hf") == "W/g" and r.sample_mass:
            y = y * r.sample_mass
        o = np.argsort(r.T)
        return r.T[o], y[o]
    Ts, ys = sig(sample)
    Tb, yb = sig(blank)
    Tz, yz = sig(sapphire)
    lo = max(Ts.min(), Tb.min(), Tz.min())
    hi = min(Ts.max(), Tb.max(), Tz.max())
    if T_range:
        lo, hi = max(lo, min(T_range)), min(hi, max(T_range))
    g = np.arange(np.ceil(lo), hi, step)
    qs, qb, qz = np.interp(g, Ts, ys), np.interp(g, Tb, yb), np.interp(g, Tz, yz)
    cpz = np.interp(g, SAPPHIRE_CP[:, 0], SAPPHIRE_CP[:, 1])
    cp = cpz * (qs - qb) / (qz - qb) * mz / ms
    return pd.DataFrame({"T": g, "Cp_J_gK": np.round(cp, 4), "Cp_safira": np.round(cpz, 4)})
