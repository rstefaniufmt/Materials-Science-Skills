"""
uvvis_peaks.py - λmax, ombros, largura de banda, métricas e pontos isosbésticos

    from uvvis_peaks import find_bands, band_metrics, isosbestic
    tab = find_bands(s)                     # DataFrame: λ (nm), E (eV), A, FWHM (nm e eV)...
    m = band_metrics(s, 664)                # A em λ fixo (opcional: linha de base local)
    m = band_metrics(s, (600, 720))         # máximo numa janela + área
    iso = isosbestic([s1, s2, s3])          # λ onde todos os espectros se cruzam

Posições: uma casa decimal se o intervalo de dados < 1 nm, senão inteiro.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, savgol_filter

from uvvis_io import HC_EV_NM, Spectrum
from uvvis_process import noise_level, resample

COLS = ["wavelength", "energy_eV", "kind", "value", "prominence", "rel_intensity", "fwhm_nm",
        "fwhm_eV", "snr", "flags"]


def _asc(s):
    o = np.argsort(s.x)
    return s.x[o].astype(float), s.y[o].astype(float)


def _signal(s: Spectrum) -> Spectrum:
    """Ordenada onde bandas apontam para cima."""
    if s.ykind in ("T", "%T"):
        return s.to_absorbance()
    if s.ykind in ("R", "%R"):
        return s.to_km()
    return s


def _rnd(v, step):
    return round(float(v), 1) if step < 1.0 else int(round(float(v)))


def _parabolic(x, y, i):
    if 0 < i < len(y) - 1:
        den = y[i - 1] - 2 * y[i] + y[i + 1]
        if den != 0:
            d = 0.5 * (y[i - 1] - y[i + 1]) / den
            if abs(d) <= 1:
                return x[i] + d * (x[i + 1] - x[i])
    return x[i]


def _local_noise(y, step, win_nm=40.0):
    w = max(5, int(round(3 / step)) | 1)
    if w >= y.size:
        return np.full(y.size, noise_level_arr(y))
    resid = y - savgol_filter(y, w, 2)
    k = max(11, int(round(win_nm / step)))
    mad = pd.Series(np.abs(resid)).rolling(k, center=True, min_periods=5).median().bfill().ffill().to_numpy()
    return 1.4826 * mad + 1e-12


def noise_level_arr(y):
    d2 = np.diff(y, 2)
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6)) + 1e-12


def find_bands(s: Spectrum, region=None, prominence=None, rel_prominence: float = 0.02,
               min_sep_nm: float = 5.0, shoulders: bool = True, top_n: int | None = None) -> pd.DataFrame:
    """Detecta máximos (λmax) e ombros. T/%T vira A e R/%R vira F(R) antes.
    prominence: limiar absoluto; padrão max(8×ruído, rel_prominence×amplitude).
    Flags: 'saturada?' (A>3), 'borda' (perto do limite do espectro: pode ser só subida de
    absorção, não máximo real), 'ruído?' (SNR<25), 'troca lâmpada/detector?' (λ em 340-360 ou
    790-870 com banda estreita)."""
    a = _signal(s)
    x, y = _asc(a)
    if region is not None:
        m = (x >= min(region)) & (x <= max(region))
        x, y = x[m], y[m]
    if x.size < 10:
        return pd.DataFrame(columns=COLS)
    step = float(np.median(np.diff(x)))
    noise = noise_level_arr(y)
    ys = savgol_filter(y, 5, 2) if y.size > 7 else y
    span = float(np.nanmax(ys) - np.nanpercentile(y, 1))
    if prominence is None:
        prominence = max(8 * noise, rel_prominence * span)
    dist = max(1, int(round(min_sep_nm / step)))
    idx, pr = find_peaks(y, prominence=prominence, distance=dist)
    ln = _local_noise(y, step)
    keep = pr["prominences"] >= 8 * ln[idx]
    idx = idx[keep]
    pr = {k: v[keep] for k, v in pr.items()}
    rows = []
    if idx.size:
        w = peak_widths(y, idx, rel_height=0.5, prominence_data=(pr["prominences"], pr["left_bases"], pr["right_bases"]))
        xi = np.arange(x.size)
        for k, i in enumerate(idx):
            xl = np.interp(w[2][k], xi, x)
            xr = np.interp(w[3][k], xi, x)
            rows.append(dict(wavelength=_parabolic(x, y, i), kind="peak", value=float(y[i]),
                             prominence=float(pr["prominences"][k]), fwhm_nm=float(xr - xl),
                             fwhm_eV=float(HC_EV_NM / xl - HC_EV_NM / xr),
                             snr=float(pr["prominences"][k] / ln[i])))
    if shoulders and x.size > 30:
        # janela e zona de exclusão proporcionais à largura das bandas principais: bandas largas
        # (SPR, d-d, corantes) geram mínimos espúrios de 2ª derivada nos flancos
        fw = np.array([r["fwhm_nm"] for r in rows]) if rows else np.array([20.0])
        win_nm = max(15.0, 0.3 * float(np.median(fw)))
        win = max(5, int(round(win_nm / step)) | 1)
        win = min(win, (x.size - 2) | 1)
        d2 = -savgol_filter(y, win, 3, deriv=2)
        j, _ = find_peaks(d2, prominence=0.12 * np.nanmax(d2), distance=max(dist, win // 2))
        main = np.array([r["wavelength"] for r in rows]) if rows else np.array([])
        mfw = np.array([r["fwhm_nm"] for r in rows]) if rows else np.array([])
        for i in j:
            if main.size:
                k = int(np.argmin(np.abs(main - x[i])))
                if abs(main[k] - x[i]) < max(10.0, 2 * min_sep_nm, 0.6 * mfw[k]):
                    continue
            lo_i, hi_i = max(0, i - 4 * win), min(x.size, i + 4 * win)
            h = y[i] - np.min(y[lo_i:hi_i])
            if h < 2 * prominence or h < 10 * ln[i]:
                continue
            rows.append(dict(wavelength=float(x[i]), kind="shoulder", value=float(y[i]), prominence=float(h),
                             fwhm_nm=np.nan, fwhm_eV=np.nan, snr=float(h / ln[i])))
    if not rows:
        return pd.DataFrame(columns=COLS)
    df = pd.DataFrame(rows)
    pmax = df.loc[df.kind == "peak", "prominence"].max() if (df.kind == "peak").any() else df.prominence.max()
    df["rel_intensity"] = (df.prominence / pmax).clip(upper=1).round(3)
    df["energy_eV"] = (HC_EV_NM / df.wavelength).round(3)
    flags = []
    for _, r in df.iterrows():
        f = []
        if a.ykind == "A" and r.value > 3.0:
            f.append("saturada?")
        width = r.fwhm_nm if np.isfinite(r.fwhm_nm) else 10
        if r.wavelength - x.min() < max(width / 2, 3 * step) or x.max() - r.wavelength < max(width / 2, 3 * step):
            f.append("borda")
        if r.snr < 25:
            f.append("ruído?")
        if (340 <= r.wavelength <= 360 or 790 <= r.wavelength <= 870) and (not np.isfinite(r.fwhm_nm) or r.fwhm_nm < 6):
            f.append("troca lâmpada/detector?")
        flags.append(",".join(f))
    df["flags"] = flags
    df["wavelength"] = [_rnd(v, step) for v in df.wavelength]
    df["fwhm_nm"] = df.fwhm_nm.round(1)
    df["fwhm_eV"] = df.fwhm_eV.round(3)
    df["value"] = df.value.round(4)
    df["prominence"] = df.prominence.round(4)
    df["snr"] = df.snr.round(0)
    df = df[COLS].sort_values("wavelength").reset_index(drop=True)
    if top_n:
        keep = df.sort_values("prominence", ascending=False).head(top_n).index
        df = df.loc[sorted(keep)].reset_index(drop=True)
    df.attrs.update(noise=noise, prominence_threshold=prominence, step=step, ykind=a.ykind,
                    title=s.title, source=s.source)
    return df


def band_metrics(s: Spectrum, where, base=None, half_window: float = 0.0) -> dict:
    """where = λ (nm) -> valor em λ (média em ±half_window);
       where = (lo, hi) -> máximo na janela (λmax, valor) e área.
    base = λ ou (λ1, λ2): subtrai linha de base local (constante em λ, ou reta entre λ1 e λ2).
    Use os MESMOS λ/janelas para todas as amostras de uma série."""
    a = _signal(s)
    x, y = _asc(a)
    step = float(np.median(np.diff(x)))

    def bl(xx):
        if base is None:
            return np.zeros_like(np.atleast_1d(xx), float)
        if np.isscalar(base):
            return np.full_like(np.atleast_1d(xx), np.interp(base, x, y), dtype=float)
        b1, b2 = min(base), max(base)
        y1, y2 = np.interp(b1, x, y), np.interp(b2, x, y)
        return y1 + (y2 - y1) * (np.atleast_1d(xx) - b1) / (b2 - b1)

    if np.isscalar(where):
        if half_window > 0:
            m = (x >= where - half_window) & (x <= where + half_window)
            v = float(np.mean(y[m] - bl(x[m])))
        else:
            v = float(np.interp(where, x, y) - bl(where)[0])
        return {"wavelength": where, "value": v, "ykind": a.ykind, "baseline": base}
    lo, hi = min(where), max(where)
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m] - bl(x[m])
    if xs.size < 3:
        raise ValueError(f"janela {where} fora do espectro")
    i = int(np.argmax(ys))
    return {"window": (lo, hi), "wavelength": _rnd(_parabolic(xs, ys, i), step), "value": float(ys[i]),
            "area": float(np.trapezoid(ys, xs)), "ykind": a.ykind, "baseline": base,
            "at_edge": i in (0, xs.size - 1)}


def isosbestic(spectra, region=None, step: float = 0.5, min_variation: float = 5.0):
    """Pontos isosbésticos: λ onde todos os espectros de uma série têm o MESMO valor enquanto
    variam em volta. Critérios: (1) a dispersão entre espectros (desvio-padrão) tem mínimo local
    pequeno comparado com a dispersão típica; (2) a maioria dos espectros troca de sinal em
    relação à média nesse λ (cruzamento real, não só região sem absorção).
    Retorna DataFrame: wavelength, value, spread (desvio entre espectros), quality (0-1)."""
    S = [_signal(s) for s in spectra]
    lo = max(s.xrange[0] for s in S)
    hi = min(s.xrange[1] for s in S)
    if region is not None:
        lo, hi = max(lo, min(region)), min(hi, max(region))
    R = [resample(s, step, lo, hi) for s in S]
    x = R[0].x
    Y = np.vstack([r.y for r in R])
    Y = savgol_filter(Y, max(5, int(3 / step) | 1), 2, axis=1)
    sd = Y.std(0)
    typical = np.percentile(sd, 75)
    noise = np.median([noise_level_arr(r.y) for r in R])
    if typical < min_variation * noise:
        return pd.DataFrame(columns=["wavelength", "value", "spread", "quality"])
    D = Y - Y.mean(0)
    idx, _ = find_peaks(-sd, distance=max(1, int(5 / step)))
    rows = []
    for i in idx:
        if sd[i] > 0.15 * typical:
            continue
        j0, j1 = max(0, i - int(6 / step)), min(x.size - 1, i + int(6 / step))
        flips = np.mean(np.sign(D[:, j0]) != np.sign(D[:, j1]))
        # precisa haver variação dos dois lados (senão é só região sem absorção)
        if flips < 0.6 or min(sd[j0], sd[j1]) < 2 * sd[i]:
            continue
        q = float(np.clip(1 - sd[i] / (0.15 * typical), 0, 1) * flips)
        rows.append(dict(wavelength=round(float(x[i]), 1), value=round(float(Y[:, i].mean()), 4),
                         spread=round(float(sd[i]), 5), quality=round(q, 2)))
    return pd.DataFrame(rows, columns=["wavelength", "value", "spread", "quality"])


# ------------------------------------------------------------------ verificação de comprimento de onda
# Valores nominais típicos (nm). Use SEMPRE os valores do certificado do seu filtro/padrão quando
# disponíveis: dependem da largura de banda espectral (SBW) e do lote.
STANDARDS = {
    "holmium_glass": [241.5, 279.4, 287.5, 333.7, 360.9, 385.9, 418.7, 453.2, 460.0, 484.5, 536.2, 637.5],
    "holmium_solution": [241.1, 249.9, 278.1, 287.2, 333.4, 345.5, 361.3, 385.7, 416.3, 451.3, 467.8,
                         485.3, 536.6, 640.5],   # Ho2O3 4% em HClO4 10%, SBW 1 nm (NIST SRM 2034)
    "didymium_glass": [329.8, 472.8, 513.5, 529.8, 573.0, 585.7, 684.8, 740.0, 807.7, 878.7],
}


def wavelength_check(s: Spectrum, standard="holmium_glass", tol: float = 3.0, certified=None) -> pd.DataFrame:
    """Compara λmax medidos com um padrão de calibração (filtro de hólmio, didímio...).
    Retorna tabela referência × medido × desvio; df.attrs['mean_offset'], ['sd'].
    Desvio sistemático (todas as bandas deslocadas igual) = erro de calibração do equipamento;
    desvios aleatórios grandes = resolução/ruído/SBW diferente do certificado."""
    ref = certified or STANDARDS[standard]
    t = find_bands(s, shoulders=False, rel_prominence=0.01)
    pk = t.wavelength.to_numpy(float)
    rows = []
    lo, hi = s.xrange
    for r in ref:
        if not (lo + 2 <= r <= hi - 2):
            continue
        if pk.size and np.min(np.abs(pk - r)) <= tol:
            m = pk[np.argmin(np.abs(pk - r))]
            rows.append(dict(reference=r, measured=m, deviation=round(m - r, 2)))
        else:
            rows.append(dict(reference=r, measured=np.nan, deviation=np.nan))
    df = pd.DataFrame(rows)
    d = df.deviation.dropna()
    df.attrs.update(mean_offset=float(d.mean()) if len(d) else np.nan,
                    sd=float(d.std(ddof=1)) if len(d) > 1 else np.nan, n_found=int(len(d)), standard=standard)
    return df
