"""
ftir_peaks.py - detecção de bandas, ombros e métricas de banda (altura/área/razões)

    from ftir_peaks import find_bands, band_metrics, band_ratio
    tab = find_bands(s)                       # DataFrame de bandas (sempre em absorbância)
    m = band_metrics(s, 1760, 1680)           # posição, altura e área com linha de base local
    r = band_ratio(s, (1760, 1680), (1490, 1440))   # ex.: índice de carbonila

Posições são reportadas como inteiros (cm-1); só com intervalo < 0.25 cm-1 (alta resolução,
fase gasosa) usa uma casa decimal. Mais casas seriam falsa precisão (a resolução instrumental típica é 4 cm-1).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, savgol_filter

from ftir_io import Spectrum
from ftir_process import atmospheric_flag, noise_level


def _asc(s):
    if s.x[0] > s.x[-1]:
        return s.x[::-1].copy(), s.y[::-1].copy()
    return s.x.copy(), s.y.copy()


def _round_pos(v, step):
    return int(round(v)) if step >= 0.25 else round(float(v), 1)


def _parabolic(x, y, i):
    if 0 < i < len(y) - 1:
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            d = 0.5 * (y0 - y2) / den
            if abs(d) <= 1:
                return x[i] + d * (x[i + 1] - x[i])
    return x[i]


def _local_noise(x, y, step, win_cm: float = 120.0):
    """Desvio robusto (MAD) do resíduo de alta frequência numa janela móvel de ~win_cm."""
    n = len(y)
    w = max(7, int(round(9 / step)) | 1)
    if w >= n:
        return np.full(n, np.std(y) * 0 + 1e-12)
    resid = y - savgol_filter(y, w, 2)
    k = max(15, int(round(win_cm / step)))
    s = pd.Series(np.abs(resid))
    mad = s.rolling(k, center=True, min_periods=5).median().bfill().ffill().to_numpy()
    return 1.4826 * mad


def intensity_class(rel: float) -> str:
    """Classe qualitativa pela intensidade relativa à banda mais forte da região."""
    if rel >= 0.75:
        return "vs"
    if rel >= 0.40:
        return "s"
    if rel >= 0.15:
        return "m"
    if rel >= 0.05:
        return "w"
    return "vw"


def find_bands(s: Spectrum, region=None, prominence=None, rel_prominence: float = 0.02,
               min_sep_cm: float = 6.0, shoulders: bool = True, top_n: int | None = None,
               broad_cm: float = 60.0) -> pd.DataFrame:
    """Detecta bandas em absorbância (T/%T é convertido antes).

    prominence      limiar absoluto (A). Padrão: max(5*ruído, rel_prominence*amplitude).
    min_sep_cm      separação mínima entre máximos (≈ resolução; 4-8 cm-1 para sólidos).
    shoulders       também procura ombros pela 2ª derivada (kind='shoulder').
    top_n           mantém só as N mais proeminentes (útil para rotular figuras).
    Colunas: position, kind, absorbance, prominence, rel_intensity, intensity, fwhm, shape, flags
    """
    a = s.to_absorbance() if s.ykind in ("T", "%T", "R", "%R") else s
    x, y = _asc(a)
    if region is not None:
        hi, lo = max(region), min(region)
        m = (x <= hi) & (x >= lo)
        x, y = x[m], y[m]
    if x.size < 10:
        return pd.DataFrame(columns=["position", "kind", "absorbance", "prominence", "rel_intensity",
                                     "intensity", "fwhm", "shape", "flags"])
    step = float(np.median(np.diff(x)))
    noise = noise_level(a)
    ys = savgol_filter(y, 5, 2) if y.size > 7 else y
    span = float(np.nanmax(ys) - np.nanpercentile(y, 1))
    if prominence is None:
        prominence = max(5 * noise, rel_prominence * span)
    dist = max(1, int(round(min_sep_cm / step)))
    idx, props = find_peaks(y, prominence=prominence, distance=dist)
    # ruído LOCAL: bordas do detector (<600 cm-1 no DTGS/MCT, >4000) e regiões de cristal
    # são muito mais ruidosas que o meio; exige proeminência >= 6x o ruído local
    lnoise = _local_noise(x, y, step)
    keep = props["prominences"] >= 8 * lnoise[idx]
    idx = idx[keep]
    props = {k: v[keep] for k, v in props.items()}
    rows = []
    if idx.size:
        w = peak_widths(y, idx, rel_height=0.5, prominence_data=(props["prominences"],
                                                                 props["left_bases"], props["right_bases"]))
        fwhm = w[0] * step
        for k, i in enumerate(idx):
            rows.append(dict(position=_parabolic(x, y, i), kind="peak", absorbance=float(y[i]),
                             prominence=float(props["prominences"][k]), fwhm=float(fwhm[k]),
                             snr=float(props["prominences"][k] / max(lnoise[i], 1e-12))))
    # ---- ombros pela 2ª derivada
    if shoulders and x.size > 30:
        win = max(5, int(round(13 / step)) | 1)
        if win >= x.size:
            win = (x.size - 2) | 1
        d2 = savgol_filter(y, win, 3, deriv=2)
        neg = -d2
        j, pj = find_peaks(neg, prominence=0.12 * np.nanmax(neg), distance=dist)
        main = np.array([r["position"] for r in rows]) if rows else np.array([])
        for i in j:
            pos = x[i]
            if main.size and np.min(np.abs(main - pos)) < max(10.0, 2 * min_sep_cm):
                continue
            # ombro precisa ter absorbância real acima do vale/linha de base local
            lo_i, hi_i = max(0, i - 4 * win), min(x.size, i + 4 * win)
            local_base = np.min(y[lo_i:hi_i])
            if y[i] - local_base < 2 * prominence or (y[i] - local_base) < 8 * lnoise[i]:
                continue
            rows.append(dict(position=float(pos), kind="shoulder", absorbance=float(y[i]),
                             prominence=float(y[i] - local_base), fwhm=np.nan,
                             snr=float((y[i] - local_base) / max(lnoise[i], 1e-12))))
    if not rows:
        return pd.DataFrame(columns=["position", "kind", "absorbance", "prominence", "rel_intensity",
                                     "intensity", "fwhm", "shape", "flags"])
    df = pd.DataFrame(rows)
    pmax = df.loc[df.kind == "peak", "prominence"].max() if (df.kind == "peak").any() else df.prominence.max()
    df["rel_intensity"] = (df["prominence"] / pmax).clip(upper=1.0).round(3)
    df["intensity"] = df["rel_intensity"].apply(intensity_class)
    shape = []
    for _, r in df.iterrows():
        if r.kind == "shoulder":
            shape.append("sh")
        elif r.fwhm > broad_cm:
            shape.append("br")
        else:
            shape.append("")
    df["shape"] = shape
    flags = []
    xmin, xmax = x.min(), x.max()
    for _, r in df.iterrows():
        f = []
        af = atmospheric_flag(r.position)
        # vapor d'água = linhas finas e fracas; banda larga/forte na mesma região é da amostra
        if af == "CO2?" or (af and np.isfinite(r.fwhm) and r.fwhm <= 8 and r.rel_intensity < 0.15):
            f.append(af)
        if r.absorbance > 2.5:
            f.append("saturada?")
        if r.snr < 25:
            f.append("ruído?")
        width = r.fwhm if np.isfinite(r.fwhm) else 10
        if r.position - xmin < width or xmax - r.position < width:
            f.append("borda")
        flags.append(",".join(f))
    df["flags"] = flags
    n_vapor = int(df["flags"].str.contains("H2O").sum())
    if n_vapor >= 6:
        df.attrs["vapor_warning"] = (f"{n_vapor} linhas finas em regiões de vapor d'água: provável "
                                     "compensação atmosférica ruim (purgar/novo background).")
    df["position"] = [_round_pos(v, step) for v in df["position"]]
    df["fwhm"] = df["fwhm"].round(0)
    df["snr"] = df["snr"].round(0)
    df["absorbance"] = df["absorbance"].round(4)
    df["prominence"] = df["prominence"].round(4)
    df = df.sort_values("position", ascending=False).reset_index(drop=True)
    if top_n:
        keep = df.sort_values("prominence", ascending=False).head(top_n).index
        df = df.loc[sorted(keep)].reset_index(drop=True)
    df.attrs.update(noise=noise, prominence_threshold=prominence, step=step, source=s.source,
                    title=s.title, converted_from=a.meta.get("converted_from"))
    return df


def band_metrics(s: Spectrum, hi: float, lo: float, base=None) -> dict:
    """Banda entre hi e lo (cm-1): posição do máximo, altura acima de uma linha de base reta
    ligando as bordas da janela (ou os pontos `base=(x1, x2)`), e área integrada acima dela.
    Altura NEGATIVA = nenhuma banda resolvida na janela (curvatura para o outro lado):
    isso é um resultado válido (ausência), não um erro - reporte o sinal.
    Use janelas IDÊNTICAS ao comparar uma série de amostras."""
    a = s.to_absorbance() if s.ykind in ("T", "%T", "R", "%R") else s
    x, y = _asc(a)
    hi, lo = max(hi, lo), min(hi, lo)
    m = (x >= lo) & (x <= hi)
    xs, ys = x[m], y[m]
    if xs.size < 3:
        raise ValueError(f"janela {hi}-{lo} fora do espectro")
    b1, b2 = (lo, hi) if base is None else (min(base), max(base))
    yb1, yb2 = np.interp(b1, x, y), np.interp(b2, x, y)
    bl = yb1 + (yb2 - yb1) * (xs - b1) / (b2 - b1)
    net = ys - bl
    i = int(np.argmax(net))
    area = np.trapezoid(net, xs) if hasattr(np, "trapezoid") else np.trapz(net, xs)
    step = float(np.median(np.diff(x)))
    return {"window": (hi, lo), "position": _round_pos(_parabolic(xs, net, i), step),
            "height": float(net[i]), "absorbance_at_max": float(ys[i]), "area": float(area)}


def band_ratio(s: Spectrum, num, den, use: str = "height") -> float:
    """Razão entre duas bandas (altura ou área com linha de base local).
    Ex.: índice de carbonila em poliolefinas = A(1760-1680)/A(1490-1440)."""
    n = band_metrics(s, *num)[use]
    d = band_metrics(s, *den)[use]
    return float(n / d) if d else float("nan")


def compare_series(spectra, windows: dict, use: str = "height") -> pd.DataFrame:
    """Tabela amostra × banda (mesmas janelas para todas). windows={'C=O': (1760, 1680), ...}"""
    rows = []
    for s in spectra:
        r = {"sample": s.title}
        for name, (hi, lo) in windows.items():
            met = band_metrics(s, hi, lo)
            r[f"{name} pos"] = met["position"]
            r[f"{name} {use}"] = round(met[use], 5)
        rows.append(r)
    return pd.DataFrame(rows)
