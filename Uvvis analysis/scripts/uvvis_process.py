"""
uvvis_process.py - pré-processamento de espectros UV-Vis

Todas as funções recebem/devolvem `Spectrum` (uvvis_io) sem alterar o original e trabalham
com x em nm. Correções de linha de base/espalhamento fazem sentido em ABSORBÂNCIA (ou F(R)).

    resample(s, step=1, lo=None, hi=None)            grade comum (subtrair, médias, séries)
    offset_correct(s, window=(750, 800))             subtrai a média de uma região sem absorção
    scatter_correct(s, fit=(600, 800), model='power') remove espalhamento A = a·λ^-n (+c)
    baseline(s, 'als'|'rubberband'|'poly'|'points', ...) -> (corrigido, linha de base)
    detect_steps(s) / splice_correct(s, at=860)      degrau na troca de detector (PMT/InGaAs...)
    smooth(s, window_nm=5) / derivative(s, order=2)  Savitzky-Golay
    normalize(s, 'max'|'at'|'area'|'minmax', wl=..)
    subtract(a, b, factor=None, region=None)         branco/solvente/matriz (fator ajustável)
    average(spectra)                                  média ± desvio em grade comum
    noise_level(s)
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve
from scipy.spatial import ConvexHull

from uvvis_io import Spectrum


def _asc(s: Spectrum):
    o = np.argsort(s.x)
    inv = np.argsort(o)
    return s.x[o], s.y[o], (lambda arr: np.asarray(arr)[inv])


def _pts(s: Spectrum, width_nm: float, minimum=5) -> int:
    n = max(int(round(width_nm / max(s.step, 1e-9))), minimum)
    return n + 1 if n % 2 == 0 else n


def _to_abs(s: Spectrum) -> Spectrum:
    return s.to_absorbance() if s.ykind in ("T", "%T") else s


def resample(s: Spectrum, step: float = 1.0, lo: float | None = None, hi: float | None = None) -> Spectrum:
    x, y, _ = _asc(s)
    lo = max(lo if lo is not None else x[0], x[0])
    hi = min(hi if hi is not None else x[-1], x[-1])
    nx = np.arange(np.ceil(lo / step) * step, hi + step / 2, step)
    out = s.copy(x=nx, y=np.interp(nx, x, y))
    out.meta["resampled_step"] = step
    return out


def offset_correct(s: Spectrum, window=(750, 800)) -> Spectrum:
    """Subtrai a média de y numa região onde a amostra não absorve (desloca a linha de base para 0).
    É a correção mais comum e mais segura para soluções; escolha uma janela realmente sem banda."""
    a = _to_abs(s)
    lo, hi = min(window), max(window)
    m = (a.x >= lo) & (a.x <= hi)
    if m.sum() < 2:
        raise ValueError(f"janela {window} fora do espectro ({a.xrange})")
    off = float(np.mean(a.y[m]))
    out = a.copy(y=a.y - off)
    out.meta["offset_correction"] = {"window": (lo, hi), "offset": off}
    return out


def scatter_correct(s: Spectrum, fit=(600, 800), model: str = "power", const: bool = False):
    """Remove espalhamento (turbidez, nanopartículas, suspensões) ajustando A = a·λ^-n (+ c)
    numa região SEM absorção e subtraindo a curva no espectro todo. n≈4 Rayleigh (partículas
    << λ); n menor (1-2) para partículas maiores. Retorna (corrigido, info)."""
    a = _to_abs(s)
    x, y, back = _asc(a)
    lo, hi = min(fit), max(fit)
    m = (x >= lo) & (x <= hi) & (y > 0)
    if m.sum() < 5:
        raise ValueError("poucos pontos (com A>0) na região de ajuste do espalhamento")
    if model == "power":
        if const:
            f = lambda xx, aa, nn, cc: aa * (xx / 500.0) ** (-nn) + cc
            p0 = [np.median(y[m]), 2.0, 0.0]
        else:
            f = lambda xx, aa, nn: aa * (xx / 500.0) ** (-nn)
            p0 = [np.median(y[m]), 2.0]
        p, _ = curve_fit(f, x[m], y[m], p0=p0, maxfev=20000)
        sc = f(x, *p)
        info = {"model": "a·(λ/500)^-n" + (" + c" if const else ""), "a": float(p[0]), "n": float(p[1]),
                "c": float(p[2]) if const else 0.0, "fit": (lo, hi)}
    else:
        raise ValueError(model)
    out = a.copy(y=back(y - sc))
    out.meta["scatter_correction"] = info
    return out, info


# ---------------------------------------------------------------- linha de base geral
def _als(y, lam=1e5, p=0.01, niter=15):
    L = len(y)
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(L, L - 2))
    D = lam * D.dot(D.transpose())
    w = np.ones(L)
    W = sparse.spdiags(w, 0, L, L)
    z = y
    for _ in range(niter):
        W.setdiag(w)
        z = spsolve(sparse.csc_matrix(W + D), w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def _rubberband(x, y):
    try:
        v = ConvexHull(np.column_stack([x, y])).vertices
    except Exception:
        return np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
    v = np.roll(v, -v.argmin())
    v = v[: v.argmax() + 1]
    v = np.unique(np.concatenate([[0], v, [len(x) - 1]]))
    return np.interp(x, x[v], y[v])


def baseline(s: Spectrum, method: str = "points", **kw):
    """Retorna (corrigido, linha_de_base). Em UV-Vis a linha de base raramente é 'curva sob as
    bandas' como no IR: prefira offset_correct ou scatter_correct. Use estas quando há deriva:
      'points'     anchors=[nm, ...] pontos sem banda, interpolação linear (padrão, controlável)
      'als'        lam (1e4-1e7), p (0.001-0.05)
      'rubberband' envoltória inferior - cuidado: 'come' bandas largas em espectros de sólidos
      'poly'       deg=2 ajustado só nos anchors"""
    a = _to_abs(s)
    x, y, back = _asc(a)
    if method == "points":
        anchors = np.sort(np.asarray(kw["anchors"], float))
        half = kw.get("half_nm", 2.0)
        ya = [np.mean(y[(x >= p - half) & (x <= p + half)]) if ((x >= p - half) & (x <= p + half)).any()
              else np.interp(p, x, y) for p in anchors]
        b = np.interp(x, anchors, ya)
    elif method == "poly":
        anchors = np.sort(np.asarray(kw["anchors"], float))
        ya = np.interp(anchors, x, y)
        b = np.polyval(np.polyfit(anchors, ya, kw.get("deg", 2)), x)
    elif method == "als":
        b = _als(y, kw.get("lam", 1e5), kw.get("p", 0.01))
    elif method == "rubberband":
        b = _rubberband(x, y)
    else:
        raise ValueError(method)
    out = a.copy(y=back(y - b))
    out.meta["baseline"] = {"method": method, **{k: v for k, v in kw.items()}}
    return out, a.copy(y=back(b), title=a.title + " (baseline)")


# ---------------------------------------------------------------- troca de detector
COMMON_SWITCHES = [(800, 870, "PMT→InGaAs/PbS (Cary ~800, Lambda ~860)"),
                   (1600, 1700, "InGaAs→PbS (~1650)"),
                   (310, 380, "troca de lâmpada D2→W (~320-350): ruído/degrau pequeno")]


def detect_steps(s: Spectrum, k: float = 8.0, side_nm: float = 15.0):
    """Procura degraus (descontinuidades) típicos de troca de detector/lâmpada. Para o maior
    salto ponto-a-ponto em cada faixa típica, ajusta retas nos dois lados (±side_nm) e compara a
    diferença extrapolada com o desvio residual. Retorna [(λ, salto, razão salto/ruído, descrição)]."""
    x, y, _ = _asc(s)
    d = np.diff(y)
    found = []
    for lo, hi, desc in COMMON_SWITCHES:
        m = np.where((x[:-1] >= lo) & (x[:-1] <= hi))[0]
        if m.size < 3:
            continue
        i = m[np.argmax(np.abs(d[m] - np.median(d[m])))]
        xm = (x[i] + x[i + 1]) / 2
        L = (x <= x[i]) & (x >= x[i] - side_nm)
        R = (x >= x[i + 1]) & (x <= x[i + 1] + side_nm)
        if L.sum() < 5 or R.sum() < 5:
            continue
        pl, pr = np.polyfit(x[L], y[L], 1), np.polyfit(x[R], y[R], 1)
        res = np.concatenate([y[L] - np.polyval(pl, x[L]), y[R] - np.polyval(pr, x[R])])
        sd = float(np.std(res, ddof=2)) + 1e-12
        jump = float(np.polyval(pr, xm) - np.polyval(pl, xm))
        if abs(jump) > k * sd:
            found.append((float(xm), jump, abs(jump) / sd, desc))
    return found


def splice_correct(s: Spectrum, at: float, width: float = 10.0, side: str = "high") -> Spectrum:
    """Corrige degrau em λ=at: extrapola linearmente os dois lados (janelas de `width` nm) e
    desloca o segmento `side` ('high' = λ > at, padrão; 'low') para encontrar o outro."""
    a = s.copy()
    x, y, back = _asc(a)
    L = (x >= at - width) & (x < at)
    R = (x > at) & (x <= at + width)
    if L.sum() < 2 or R.sum() < 2:
        raise ValueError("pontos insuficientes em volta do degrau")
    yl = np.polyval(np.polyfit(x[L], y[L], 1), at)
    yr = np.polyval(np.polyfit(x[R], y[R], 1), at)
    y2 = y.copy()
    if side == "high":
        y2[x > at] -= (yr - yl)
    else:
        y2[x < at] += (yr - yl)
    out = a.copy(y=back(y2))
    out.meta.setdefault("splice", []).append({"at": at, "shift": float(yr - yl), "side": side})
    return out


# ---------------------------------------------------------------- suavização / derivadas / normalização
def smooth(s: Spectrum, window_nm: float = 5.0, poly: int = 2) -> Spectrum:
    """Savitzky-Golay. Janela menor que a largura das bandas mais estreitas (bandas de solução
    tipicamente 20-100 nm de largura; linhas de Ho/Nd/lantanídeos são finas!)."""
    n = _pts(s, window_nm)
    x, y, back = _asc(s)
    out = s.copy(y=back(savgol_filter(y, n, poly)))
    out.meta["smoothing"] = f"SG {n} pts ({window_nm} nm), poly {poly}"
    return out


def derivative(s: Spectrum, order: int = 1, window_nm: float = 9.0, poly: int = 3) -> Spectrum:
    """Derivada SG dy/dλ. 1ª: borda de absorção (máximo) e ponto de inflexão; 2ª: separa bandas
    sobrepostas (bandas viram mínimos)."""
    n = _pts(s, window_nm)
    x, y, back = _asc(s)
    d = savgol_filter(y, n, poly, deriv=order, delta=float(np.mean(np.diff(x))))
    out = s.copy(y=back(d), ykind=f"d{order}")
    out.meta["derivative"] = f"SG d{order}, {n} pts, poly {poly}"
    return out


def normalize(s: Spectrum, method: str = "max", wl: float | None = None, region=None) -> Spectrum:
    """'max' (na região), 'at' (valor em λ=wl vira 1), 'minmax', 'area'."""
    y = s.y.astype(float)
    m = np.ones_like(y, bool)
    if region is not None:
        m = (s.x >= min(region)) & (s.x <= max(region))
    if method == "max":
        y = y / np.nanmax(y[m])
    elif method == "at":
        y = y / s.at(wl)
    elif method == "minmax":
        mn, mx = np.nanmin(y[m]), np.nanmax(y[m])
        y = (y - mn) / (mx - mn)
    elif method == "area":
        x, yy, _ = _asc(s)
        y = y / np.trapezoid(np.abs(yy), x)
    else:
        raise ValueError(method)
    out = s.copy(y=y)
    out.meta["normalization"] = method + (f"@{wl}" if wl else "")
    return out


def subtract(a: Spectrum, b: Spectrum, factor: float | None = 1.0, region=None, step: float | None = None):
    """a − f·b em absorbância numa grade comum. factor=None ajusta f por mínimos quadrados em
    `region` (onde só b absorve). Retorna (diferença, f)."""
    A, B = _to_abs(a), _to_abs(b)
    st = step or max(A.step, B.step)
    lo, hi = max(A.xrange[0], B.xrange[0]), min(A.xrange[1], B.xrange[1])
    A, B = resample(A, st, lo, hi), resample(B, st, lo, hi)
    if factor is None:
        m = np.ones_like(A.x, bool) if region is None else (A.x >= min(region)) & (A.x <= max(region))
        factor = float(np.dot(A.y[m], B.y[m]) / np.dot(B.y[m], B.y[m]))
    out = A.copy(y=A.y - factor * B.y, title=f"{a.title} − {factor:.3g}×{b.title}")
    out.meta["subtraction"] = {"factor": factor, "region": region}
    return out, factor


def average(spectra, step: float | None = None):
    """Média e desvio-padrão de réplicas numa grade comum. Retorna (média, desvio) como Spectrum."""
    st = step or max(s.step for s in spectra)
    lo = max(s.xrange[0] for s in spectra)
    hi = min(s.xrange[1] for s in spectra)
    R = [resample(_to_abs(s), st, lo, hi) for s in spectra]
    Y = np.vstack([r.y for r in R])
    m = R[0].copy(y=Y.mean(0), title=f"média de {len(R)}")
    sd = R[0].copy(y=Y.std(0, ddof=1) if len(R) > 1 else np.zeros_like(Y[0]), title="desvio-padrão")
    m.meta["n_replicates"] = len(R)
    return m, sd


def noise_level(s: Spectrum) -> float:
    """Ruído RMS robusto pela 2ª diferença (independe de bandas largas)."""
    _, y, _ = _asc(s)
    d2 = np.diff(y, 2)
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6))
