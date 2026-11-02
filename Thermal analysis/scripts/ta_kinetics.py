"""
ta_kinetics.py - cinética não isotérmica isoconversional (recomendações ICTAC, Vyazovkin et al.,
Thermochim. Acta 520 (2011) 1-19) a partir de corridas em várias taxas de aquecimento

    from ta_kinetics import conversion, isoconversional, kissinger
    runs = [load(f) for f in arquivos]            # TG (ou DSC) em ≥3 taxas (ideal 4-5)
    res = isoconversional(runs, window=(250, 450)) # Friedman, OFW, KAS, Starink, Vyazovkin
    res["table"]                                   # Ea(α) por método ± IC95%
    kissinger(runs, window=(250, 450))             # Ea pelo pico da DTG (ou do DSC)

α para TG: α = (m0 − m)/(m0 − m∞) com m0 e m∞ nas bordas da janela de temperatura da etapa.
α para DSC: área parcial do pico / área total (janela do pico, linha de base linear).
β é MEDIDO (dT/dt médio na janela), não o nominal.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize_scalar
from scipy.signal import savgol_filter

from ta_io import Run
from ta_process import endo_signal, mass_pct, segments, _rate

R = 8.314462618


def conversion(r: Run, window, signal: str = "auto", endo: str | None = None):
    """α(T) no trecho de aquecimento dentro da janela (T1, T2). Retorna dict com T (K), t (s),
    alpha, dadt (1/s), beta (K/min medido), Tpeak (K, máximo de dα/dt)."""
    lo, hi = min(window), max(window)
    segs = [s for s in segments(r) if s["type"] == "heat"]
    if not segs:
        raise ValueError(f"{r.title}: sem segmento de aquecimento")
    s = max(segs, key=lambda x: x["i1"] - x["i0"])
    idx = np.arange(s["i0"], s["i1"])
    T = r.T[idx]
    m = (T >= lo) & (T <= hi)
    idx = idx[m]
    if idx.size < 30:
        raise ValueError(f"{r.title}: poucos pontos em {window}")
    T = r.T[idx].astype(float)
    t = r.time[idx].astype(float) * 60.0 if r.has("time") else None
    use = signal
    if use == "auto":
        use = "mass" if (r.has("mass_pct")) else "hf"
    if use == "mass":
        y, _ = mass_pct(r)
        y = y[idx]
        k = max(3, len(y) // 100)
        m0, mf = np.median(y[:k]), np.median(y[-k:])
        alpha = (m0 - y) / (m0 - mf)
    else:
        _, tt, q, u, d, src = endo_signal(r, endo)
        q = q[idx]
        if t is None:
            raise ValueError("DSC exige canal de tempo para α")
        b = np.interp(t, [t[0], t[-1]], [np.mean(q[:5]), np.mean(q[-5:])])
        net = q - b
        sgn = 1 if np.nanmax(net) >= -np.nanmin(net) else -1
        net = np.clip(sgn * net, 0, None)
        cum = np.concatenate([[0], np.cumsum(0.5 * (net[1:] + net[:-1]) * np.diff(t))])
        alpha = cum / cum[-1]
    if t is None:
        beta_meas = (r.meta.get("heating_rate") or [10.0])
        beta_meas = beta_meas[0] if isinstance(beta_meas, list) else beta_meas
        t = (T - T[0]) / beta_meas * 60.0
    beta = float(np.polyfit(t / 60.0, T, 1)[0])
    # dα/dt suavizado no tempo
    dt = np.median(np.diff(t))
    w = max(7, int(round(len(t) / 60)) | 1)
    w = min(w, (len(t) - 2) | 1)
    a_s = savgol_filter(alpha, w, 2)
    dadt = savgol_filter(alpha, w, 2, deriv=1, delta=dt)
    a_mono = np.maximum.accumulate(np.clip(a_s, -0.05, 1.05))
    jp = int(np.argmax(dadt))
    return dict(T=T + 273.15, t=t, alpha=a_mono, dadt=dadt, beta=beta, Tpeak=float(T[jp] + 273.15),
                alpha_peak=float(a_mono[jp]), title=r.title)


def _at_alpha(c, a):
    """T_α, t_α e (dα/dt)_α por interpolação no primeiro cruzamento."""
    al = c["alpha"]
    i = np.searchsorted(al, a)
    if i <= 0 or i >= len(al):
        return np.nan, np.nan, np.nan
    f = (a - al[i - 1]) / (al[i] - al[i - 1]) if al[i] != al[i - 1] else 0.0
    T = c["T"][i - 1] + f * (c["T"][i] - c["T"][i - 1])
    t = c["t"][i - 1] + f * (c["t"][i] - c["t"][i - 1])
    d = c["dadt"][i - 1] + f * (c["dadt"][i] - c["dadt"][i - 1])
    return T, t, d


def _p_senum_yang(x):
    """Integral de temperatura p(x), x = E/RT (aproximação de Senum-Yang de 4ª ordem)."""
    return np.exp(-x) / x * (x ** 3 + 18 * x ** 2 + 88 * x + 96) / (x ** 4 + 20 * x ** 3 + 120 * x ** 2 + 240 * x + 120)


def _linfit(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return np.nan, np.nan, np.nan, np.nan
    res = stats.linregress(x, y)
    tcrit = stats.t.ppf(0.975, len(x) - 2)
    return res.slope, res.stderr * tcrit, res.intercept, res.rvalue ** 2


def isoconversional(runs, window, alphas=None, signal: str = "auto", endo: str | None = None,
                    methods=("friedman", "ofw", "kas", "starink", "vyazovkin")) -> dict:
    """Ea(α) em kJ/mol por método, com IC95% (regressão). Retorna dict com 'table' (DataFrame),
    'summary' (média e faixa em 0.1 ≤ α ≤ 0.9), 'conv' (curvas α(T) por corrida) e 'warnings'."""
    if alphas is None:
        alphas = np.round(np.arange(0.05, 0.951, 0.05), 3)
    conv = [conversion(r, window, signal, endo) for r in runs]
    conv.sort(key=lambda c: c["beta"])
    betas = np.array([c["beta"] for c in conv])
    warn = []
    if len(conv) < 3:
        warn.append("Menos de 3 taxas de aquecimento: ICTAC recomenda 3-5 (Friedman e Vyazovkin mais robustos com 5).")
    if betas.max() / betas.min() < 3:
        warn.append("Faixa de β estreita (< 3×): Ea com incerteza alta.")
    rows = []
    for a in alphas:
        Ts, ts, ds = zip(*[_at_alpha(c, a) for c in conv])
        Ts, ts, ds = np.array(Ts), np.array(ts), np.array(ds)
        if np.isnan(Ts).any():
            continue
        x = 1.0 / Ts
        row = {"alpha": a}
        for meth in methods:
            if meth == "friedman":
                s, ci, b, r2 = _linfit(x, np.log(np.where(ds > 0, ds, np.nan)))
                E, dE = -s * R, ci * R
            elif meth == "ofw":
                s, ci, b, r2 = _linfit(x, np.log(betas))
                E, dE = -s * R / 1.052, ci * R / 1.052
            elif meth == "kas":
                s, ci, b, r2 = _linfit(x, np.log(betas / Ts ** 2))
                E, dE = -s * R, ci * R
            elif meth == "starink":
                s, ci, b, r2 = _linfit(x, np.log(betas / Ts ** 1.92))
                E, dE = -s * R / 1.0008, ci * R / 1.0008
            elif meth == "vyazovkin":
                bs = betas / 60.0
                def phi(Ek):
                    I = _p_senum_yang(Ek * 1e3 / (R * Ts))           # I(E,T) = (E/R)·p(x); E/R se cancela
                    M = np.outer(I / bs, bs / I)                       # M_ij = I_i β_j / (I_j β_i)
                    return float(np.sum(M) - len(bs))
                E0 = row.get("kas", 150.0) if np.isfinite(row.get("kas", np.nan)) else 150.0
                res = minimize_scalar(phi, bounds=(max(5.0, E0 * 0.3), E0 * 3 + 50), method="bounded")
                E, dE, r2 = res.x * 1e3, np.nan, np.nan
            row[meth] = E / 1e3
            row[f"{meth}_ic95"] = dE / 1e3 if np.isfinite(dE) else np.nan
            row[f"{meth}_r2"] = r2
        for c, T in zip(conv, Ts):
            row[f"T@{c['beta']:.3g}K/min"] = T - 273.15
        rows.append(row)
    df = pd.DataFrame(rows)
    summ = {}
    mid = df[(df.alpha >= 0.1) & (df.alpha <= 0.9)] if len(df) else df
    for meth in methods:
        if meth in mid:
            v = mid[meth].dropna()
            if len(v):
                summ[meth] = dict(media=float(v.mean()), dp=float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                                  min=float(v.min()), max=float(v.max()))
    ref = summ.get("friedman") or next(iter(summ.values()), None)
    if ref and ref["media"] > 0 and (ref["max"] - ref["min"]) / ref["media"] > 0.2:
        warn.append(f"Ea varia {100 * (ref['max'] - ref['min']) / ref['media']:.0f}% com α (0.1–0.9): processo "
                    "multietapa/complexo - não use um único Ea nem um modelo de etapa única para todo o processo.")
    for meth in methods:
        col = f"{meth}_r2"
        if col in mid and (mid[col] < 0.98).mean() > 0.3:
            warn.append(f"{meth}: R² < 0.98 em boa parte de α - confira sobreposição de etapas, janela e linha de base.")
    return dict(table=df, summary=summ, conv=conv, betas=betas, warnings=warn, window=window)


def kissinger(runs, window, signal: str = "auto", endo: str | None = None) -> dict:
    """Kissinger: ln(β/Tp²) = ln(AR/E) − E/(R·Tp), Tp = T no máximo de dα/dt (DTG ou DSC).
    Assume (rigorosamente) reação de 1ª ordem e α_p constante; ICTAC: use só como estimativa e
    prefira os isoconversionais. Retorna Ea (kJ/mol) ± IC95%, R², ln A (A em 1/s, 1ª ordem)."""
    conv = sorted([conversion(r, window, signal, endo) for r in runs], key=lambda c: c["beta"])
    b = np.array([c["beta"] for c in conv])
    Tp = np.array([c["Tpeak"] for c in conv])
    s, ci, icpt, r2 = _linfit(1 / Tp, np.log(b / Tp ** 2))
    E = -s * R
    lnA = icpt + np.log(E / R) - np.log(60.0)          # β em K/min -> A em 1/s
    ap = [c["alpha_peak"] for c in conv]
    warn = []
    if np.ptp(ap) > 0.1:
        warn.append(f"α no pico varia ({min(ap):.2f}–{max(ap):.2f}): Kissinger fica impreciso.")
    return dict(Ea=E / 1e3, ic95=ci * R / 1e3, r2=r2, lnA=lnA, betas=b, Tp_C=Tp - 273.15, alpha_peak=ap,
                warnings=warn)
