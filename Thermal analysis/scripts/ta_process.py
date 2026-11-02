"""
ta_process.py - segmentos do programa térmico, direção endo/exo, DTG, normalizações, branco

    from ta_process import segments, select, dtg, endo_direction, hf_specific, subtract_blank
    segs = segments(r)                    # [{'n':1,'type':'heat','T0':..,'T1':..,'rate':10.0,...}, ...]
    h2 = select(r, "heat2")               # Run só com o 2º aquecimento ('heat1','cool1','iso1','last_heat', 3, ...)
    d = dtg(r, per="T")                   # dm/dT em %/°C (negativo = perda), alinhado com r
    endo, fonte = endo_direction(r)       # 'up' | 'down'
"""
from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter

from ta_io import Run

TYPE_PT = {"heat": "aquecimento", "cool": "resfriamento", "iso": "isoterma"}


# ============================================================================ segmentos
def _rate(r: Run, smooth_min: float = 0.5):
    """dT/dt (K/min) suavizado, alinhado com r. NaN se faltar tempo."""
    if not (r.has("time") and r.has("T")):
        return np.full(r.n, np.nan)
    t, T = r.time.astype(float), r.T.astype(float)
    ok = np.isfinite(t) & np.isfinite(T)
    if ok.sum() < 7:
        return np.full(r.n, np.nan)
    tt, TT = t[ok], T[ok]
    dt = np.median(np.diff(tt)) if tt.size > 1 else 1.0
    if not np.all(np.diff(tt) > 0) or dt <= 0:
        o = np.argsort(tt, kind="stable")
        tt, TT = tt[o], TT[o]
        u = np.concatenate([[True], np.diff(tt) > 0])
        tt, TT = tt[u], TT[u]
        dt = np.median(np.diff(tt))
    grid = np.arange(tt[0], tt[-1] + dt / 2, dt)
    Tg = np.interp(grid, tt, TT)
    w = max(7, int(round(smooth_min / dt)) | 1)
    w = min(w, (grid.size - 2) | 1) if grid.size > 9 else 5
    if grid.size <= w:
        d = np.gradient(Tg, grid)
    else:
        d = savgol_filter(Tg, w, 2, deriv=1, delta=dt)
    out = np.full(r.n, np.nan)
    out[ok] = np.interp(t[ok], grid, d)
    return out


def segments(r: Run, iso_tol: float = 0.4, min_minutes: float = 1.0):
    """Segmentos do programa: usa o canal 'segment' se existir; senão detecta pelo sinal de dT/dt
    (|dT/dt| < iso_tol K/min = isoterma). Devolve lista de dicts com índices [i0, i1)."""
    n = r.n
    rate = _rate(r)
    if r.has("segment"):
        lab = r.segment.astype(float)
        cuts = [0] + [i for i in range(1, n) if lab[i] != lab[i - 1]] + [n]
    else:
        cls = np.where(np.isnan(rate), 0, np.where(rate > iso_tol, 1, np.where(rate < -iso_tol, -1, 0)))
        cuts = [0] + [i for i in range(1, n) if cls[i] != cls[i - 1]] + [n]
        # fundir pedaços curtos (transições, ruído) com o vizinho
        t = r.time if r.has("time") else np.arange(n, dtype=float)
        changed = True
        while changed and len(cuts) > 2:
            changed = False
            for k in range(len(cuts) - 1):
                a, b = cuts[k], cuts[k + 1]
                dur = t[min(b, n) - 1] - t[a]
                if dur < min_minutes and len(cuts) > 2:
                    del cuts[k + 1 if k + 1 < len(cuts) - 1 else k]
                    changed = True
                    break
        # regras: após fundir, reclassificar por mediana
    out, counts = [], {"heat": 0, "cool": 0, "iso": 0}
    for k in range(len(cuts) - 1):
        a, b = cuts[k], cuts[k + 1]
        if b - a < 2:
            continue
        rr = np.nanmedian(rate[a:b]) if np.isfinite(rate[a:b]).any() else np.nan
        T = r.T[a:b] if r.has("T") else np.array([np.nan])
        if np.isnan(rr):
            dT = T[-1] - T[0] if T.size > 1 else 0
            typ = "heat" if dT > 2 else ("cool" if dT < -2 else "iso")
        else:
            typ = "heat" if rr > iso_tol else ("cool" if rr < -iso_tol else "iso")
        counts[typ] += 1
        seg = {"n": len(out) + 1, "i0": a, "i1": b, "type": typ, "label": f"{TYPE_PT[typ]} {counts[typ]}",
               "key": f"{typ}{counts[typ]}", "T0": float(T[0]), "T1": float(T[-1]),
               "rate": float(round(rr, 2)) if np.isfinite(rr) else None,
               "t0": float(r.time[a]) if r.has("time") else None, "t1": float(r.time[b - 1]) if r.has("time") else None}
        out.append(seg)
    # fundir segmentos vizinhos do mesmo tipo (sobra da detecção)
    merged = []
    for s in out:
        if merged and merged[-1]["type"] == s["type"] and not r.has("segment"):
            m = merged[-1]
            m.update(i1=s["i1"], T1=s["T1"], t1=s["t1"])
        else:
            merged.append(s)
    counts = {"heat": 0, "cool": 0, "iso": 0}
    for i, s in enumerate(merged):
        counts[s["type"]] += 1
        s.update(n=i + 1, label=f"{TYPE_PT[s['type']]} {counts[s['type']]}", key=f"{s['type']}{counts[s['type']]}")
        rr = rate[s["i0"]:s["i1"]]
        s["rate"] = float(round(np.nanmedian(rr), 2)) if np.isfinite(rr).any() else s["rate"]
    return merged


def select(r: Run, which=None) -> Run:
    """which: None (tudo), int (nº do segmento, 1-based), 'heat1', 'cool1', 'heat2', 'iso1',
    'last_heat', 'main_heat' (o aquecimento de maior faixa de T)."""
    if which in (None, "all", ""):
        return r
    segs = segments(r)
    s = None
    if isinstance(which, (int, np.integer)) or (isinstance(which, str) and which.isdigit()):
        k = int(which)
        s = next((x for x in segs if x["n"] == k), None)
    elif which == "last_heat":
        hs = [x for x in segs if x["type"] == "heat"]
        s = hs[-1] if hs else None
    elif which == "main_heat":
        # aquecimento com maior variação de massa (se houver TG); senão o de maior faixa de T
        hs = [x for x in segs if x["type"] == "heat" and abs(x["T1"] - x["T0"]) > 20]
        if hs and r.has("mass_pct"):
            def dm(x):
                mm = r.mass_pct[x["i0"]:x["i1"]]
                mm = mm[np.isfinite(mm)]
                return abs(mm[0] - mm[-1]) if mm.size else 0
            s = max(hs, key=dm)
        else:
            s = max(hs, key=lambda x: abs(x["T1"] - x["T0"])) if hs else None
    else:
        s = next((x for x in segs if x["key"] == which), None)
    if s is None:
        raise ValueError(f"segmento '{which}' não encontrado. Disponíveis: " +
                         ", ".join(f"{x['n']}={x['key']} ({x['T0']:.0f}→{x['T1']:.0f} °C)" for x in segs))
    m = np.zeros(r.n, bool)
    m[s["i0"]:s["i1"]] = True
    endo = endo_direction(r) if r.has("hf") else (None, None)   # decidir com o ensaio inteiro
    out = r.subset(m)
    out.meta["segment_info"] = s
    if endo[0]:
        out.meta["endo"], out.meta["endo_source"] = endo
    out.title = f"{r.title} – {s['label']}"
    return out


def heating_rate(r: Run) -> float | None:
    """Taxa média (K/min) do maior segmento de aquecimento (medida, não nominal)."""
    hs = [s for s in segments(r) if s["type"] == "heat" and s["rate"]]
    if hs:
        return max(hs, key=lambda s: s["i1"] - s["i0"])["rate"]
    nominal = r.meta.get("heating_rate")
    if isinstance(nominal, list) and nominal:
        return float(max(nominal))
    return float(nominal) if isinstance(nominal, (int, float)) else None


# ============================================================================ DSC: unidades e direção
def hf_specific(r: Run, mass_mg: float | None = None):
    """Fluxo de calor em W/g (= mW/mg) se possível. Retorna (array, unidade)."""
    if not r.has("hf"):
        return None, ""
    u = r.units.get("hf", "")
    y = r.hf.astype(float)
    m = mass_mg or r.sample_mass
    if u == "W/g":
        return y, "W/g"
    if u == "mW" and m:
        return y / m, "W/g"
    if u == "uV/mg":
        return y, "µV/mg"
    if u == "uV" and m:
        return y / m, "µV/mg"
    return y, {"uV": "µV", "K": "K"}.get(u, u)


def _noise(y):
    d2 = np.diff(y[np.isfinite(y)], 2)
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6)) + 1e-15 if d2.size > 5 else 1e-15


def endo_direction(r: Run, override: str | None = None):
    """Direção endotérmica do sinal ('up' = endo positivo). Prioridade:
    1) override do usuário; 2) rótulo explícito no arquivo; 3) dados: com aquecimento E
    resfriamento, a linha de base do aquecimento desloca-se no sentido endo (a amostra absorve
    calor ∝ Cp·β) e a do resfriamento no sentido exo; 4) evidência do arquivo (NETZSCH: área de
    fusão dos padrões de calibração); 5) convenção típica do fabricante; 6) desconhecida -> 'down'
    (exo para cima, o mais comum) com aviso."""
    if override in ("up", "down"):
        return override, "informado pelo usuário"
    src = r.meta.get("endo_source", "")
    declared = r.meta.get("endo")
    if declared and "rótulo" in src:
        return declared, src
    data_dir = None
    if r.has("hf") and r.has("T"):
        segs = segments(r)
        hs = [s for s in segs if s["type"] == "heat" and s["i1"] - s["i0"] > 20]
        cs = [s for s in segs if s["type"] == "cool" and s["i1"] - s["i0"] > 20]
        if hs and cs:
            y, _ = hf_specific(r)
            span = lambda s: abs(s["T1"] - s["T0"])
            h, c = max(hs, key=span), max(cs, key=span)
            Th, yh = r.T[h["i0"]:h["i1"]], y[h["i0"]:h["i1"]]
            Tc, yc = r.T[c["i0"]:c["i1"]], y[c["i0"]:c["i1"]]
            lo, hi = max(np.nanmin(Th), np.nanmin(Tc)), min(np.nanmax(Th), np.nanmax(Tc))
            if hi - lo > 20:
                mh = (Th > lo) & (Th < hi)
                mc = (Tc > lo) & (Tc < hi)
                diff = np.nanmedian(yh[mh]) - np.nanmedian(yc[mc])
                nz = max(_noise(yh), _noise(yc))
                if abs(diff) > 20 * nz:
                    data_dir = "up" if diff > 0 else "down"
    if data_dir:
        if declared and declared != data_dir:
            r.meta.setdefault("warnings", []).append(
                f"Direção endo pelos dados ({data_dir}) contradiz '{src}' ({declared}); usando os dados.")
        return data_dir, "dados (deslocamento da linha de base aquecimento × resfriamento)"
    if declared:
        return declared, src or "metadados"
    return "down", "DESCONHECIDA - assumido endo para baixo (exo para cima); confirme com um evento conhecido"


def endo_signal(r: Run, override: str | None = None):
    """(T, t, y_endo_pos, unidade, direção, fonte): sinal orientado com endo POSITIVO (uso interno)."""
    y, u = hf_specific(r)
    d, src = endo_direction(r, override)
    s = 1.0 if d == "up" else -1.0
    t = r.time if r.has("time") else None
    return r.T, t, s * y, u, d, src


# ============================================================================ TG: normalização e DTG
def mass_pct(r: Run, basis: str | float | None = None):
    """Massa em %. basis=None: como no arquivo (100% = massa inicial); basis=150 (°C) -> base seca
    (100% = massa a 150 °C no 1º aquecimento); basis='mg' -> massa absoluta."""
    if basis == "mg":
        return r.mass, "mg"
    m = r.mass_pct
    if m is None:
        return None, ""
    if basis is None:
        return m, "%"
    T = r.T
    idx = np.where(T >= float(basis))[0]
    if idx.size == 0:
        raise ValueError(f"temperatura de referência {basis} °C não atingida")
    return 100.0 * m / m[idx[0]], f"% (base: massa a {float(basis):g} °C)"


def dtg(r: Run, per: str = "T", window_K: float = 6.0, window_min: float | None = None, y=None):
    """DTG = d(massa %)/dT (%/°C, per='T') ou d(massa %)/dt (%/min, per='time'), alinhado com r.
    Calculada no tempo (sempre monótono) com Savitzky-Golay e dividida por dT/dt -> funciona em
    isotermas (per='time') e não amplifica ruído onde T 'anda para trás'. Negativa = perda.
    Janela: window_K/β (em min) para per='T'; window_min explícita se dada."""
    m = r.mass_pct if y is None else y
    if m is None:
        raise ValueError("sem canal de massa")
    if not r.has("time"):
        # sem tempo: derivar direto em T (só em trecho monótono)
        T = r.T
        o = np.argsort(T)
        w = max(7, int(round(window_K / max(np.median(np.abs(np.diff(T[o]))), 1e-6))) | 1)
        d = savgol_filter(m[o], min(w, (len(T) - 2) | 1), 2, deriv=1, delta=float(np.median(np.diff(T[o]))))
        out = np.empty_like(d)
        out[o] = d
        return out
    t = r.time.astype(float)
    ok = np.isfinite(t) & np.isfinite(m)
    tt, mm = t[ok], m[ok]
    o = np.argsort(tt, kind="stable")
    tt, mm = tt[o], mm[o]
    u = np.concatenate([[True], np.diff(tt) > 0])
    tt, mm = tt[u], mm[u]
    dt = float(np.median(np.diff(tt)))
    grid = np.arange(tt[0], tt[-1] + dt / 2, dt)
    mg = np.interp(grid, tt, mm)
    beta = heating_rate(r) or 10.0
    wmin = window_min if window_min else window_K / max(abs(beta), 0.1)
    w = max(7, int(round(wmin / dt)) | 1)
    w = min(w, (grid.size - 2) | 1)
    dmdt = savgol_filter(mg, w, 2, deriv=1, delta=dt)
    out_t = np.full(r.n, np.nan)
    out_t[ok] = np.interp(t[ok], grid, dmdt)
    if per == "time":
        return out_t
    rate = _rate(r, smooth_min=wmin)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.where(np.abs(rate) > 0.2, out_t / rate, np.nan)
    return out


# ============================================================================ branco / correção
def subtract_blank(r: Run, blank: Run, channel: str = "auto", on: str = "T") -> Run:
    """Subtrai uma curva de branco (cadinho vazio, mesmo programa): empuxo/deriva em TG e linha de
    base em DSC. Interpola o branco no eixo `on` ('T' só para trecho de aquecimento monótono;
    'time' para programas com isotermas/resfriamentos - exige mesmo programa e mesma origem)."""
    out = r.copy()
    chans = ["mass", "hf"] if channel == "auto" else [channel]
    for ch in chans:
        if ch == "mass":
            if not (r.has("mass") and (blank.has("mass") or "mass_change" in blank.extra)):
                continue
            bm = blank.mass if blank.has("mass") else blank.extra["mass_change"]
            bm = bm - bm[np.isfinite(bm)][0]
        elif ch == "hf":
            if not (r.has("hf") and blank.has("hf")):
                continue
            if r.units.get("hf") != blank.units.get("hf") and not (r.units.get("hf") in ("mW",) and blank.units.get("hf") in ("mW",)):
                raise ValueError(f"unidades diferentes: amostra {r.units.get('hf')} x branco {blank.units.get('hf')}")
            bm = blank.hf
        else:
            continue
        xs = blank.T if on == "T" else blank.time
        xr = r.T if on == "T" else r.time
        o = np.argsort(xs)
        b_on_r = np.interp(xr, xs[o], bm[o], left=np.nan, right=np.nan)
        out.data[ch] = r.data[ch] - b_on_r
        out.meta.setdefault("processing", []).append(f"branco subtraído de {ch} ({blank.title}, eixo {on})")
    if out.has("mass") and out.sample_mass:
        out.data["mass_pct"] = 100 * out.mass / out.mass[np.isfinite(out.mass)][0]
    return out


def smooth(y, window_pts: int = 11, poly: int = 2):
    w = max(5, int(window_pts) | 1)
    ok = np.isfinite(y)
    out = np.array(y, float)
    if ok.sum() > w:
        out[ok] = savgol_filter(y[ok], w, poly)
    return out
