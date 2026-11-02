"""
ftir_process.py - pré-processamento de espectros FT-IR

Todas as funções recebem e devolvem `Spectrum` (de ftir_io) sem alterar o original.
Linha de base, derivadas e normalização devem ser feitas em ABSORBÂNCIA
(é a grandeza aditiva/proporcional à concentração); converta antes com s.to_absorbance().

Funções:
    resample(s, step=1.0, hi=None, lo=None)      grade uniforme (para comparar/subtrair)
    baseline(s, method='als'|'rubberband'|'poly'|'points', **kw) -> (corrigido, linha_de_base)
    smooth(s, window_cm=9, poly=2)                Savitzky-Golay
    derivative(s, order=2, window_cm=13, poly=3)  derivadas SG (2ª derivada revela ombros)
    normalize(s, method='max'|'area'|'band'|'minmax'|'snv'|'vector', band=(hi, lo))
    atr_correction(s, ref=1000.0)                 correção ATR simples (profundidade ∝ λ)
    subtract(a, b, factor=None, region=None)     a - f*b (f por mínimos quadrados na região)
    atmospheric_mask(x)                           pontos em regiões de CO2/H2O vapor
    co2_flatten(s)                                interpola sobre 2400-2280 cm-1
    noise_level(s)                                ruído RMS estimado (região 2200-1900 ou MAD da 2ª dif.)
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.signal import savgol_filter
from scipy.sparse.linalg import spsolve
from scipy.spatial import ConvexHull

from ftir_io import Spectrum

# regiões típicas de interferência atmosférica (cm-1)
CO2_REGIONS = [(2400, 2280), (680, 655)]
H2O_REGIONS = [(3950, 3500), (1950, 1300)]


def _asc(s: Spectrum):
    """x crescente (para algoritmos), e função para voltar à ordem original."""
    if s.x[0] > s.x[-1]:
        return s.x[::-1], s.y[::-1], (lambda arr: arr[::-1])
    return s.x, s.y, (lambda arr: arr)


def _pts(s: Spectrum, width_cm: float, odd=True, minimum=5) -> int:
    n = max(int(round(width_cm / max(s.step, 1e-9))), minimum)
    if odd and n % 2 == 0:
        n += 1
    return n


def resample(s: Spectrum, step: float = 1.0, hi: float | None = None, lo: float | None = None) -> Spectrum:
    x, y, _ = _asc(s)
    hi = min(hi if hi is not None else x[-1], x[-1])
    lo = max(lo if lo is not None else x[0], x[0])
    nx = np.arange(np.ceil(lo / step) * step, hi + step / 2, step)
    ny = np.interp(nx, x, y)
    out = s.copy(x=nx[::-1], y=ny[::-1])
    out.meta["resampled_step"] = step
    return out


# ---------------------------------------------------------------- linha de base
def _als(y, lam=1e6, p=0.01, niter=15):
    """Asymmetric Least Squares (Eilers & Boelens 2005)."""
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
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
    """Envoltória convexa inferior (rubberband), como no OPUS/Spectrum."""
    pts = np.column_stack([x, y])
    try:
        hull = ConvexHull(pts)
        v = hull.vertices
    except Exception:
        return np.interp(x, [x[0], x[-1]], [y[0], y[-1]])
    v = np.roll(v, -v.argmin())
    # segmento inferior: do menor x ao maior x
    v = v[: v.argmax() + 1]
    v = np.unique(np.concatenate([[0], v, [len(x) - 1]]))
    return np.interp(x, x[v], y[v])


def _modpoly(x, y, deg=3, niter=100, tol=1e-4):
    """Polinômio modificado (Lieber & Mahadevan-Jansen 2003): ajusta e corta por cima."""
    yy = y.copy()
    xn = (x - x.mean()) / (x.std() or 1)
    base = yy
    for _ in range(niter):
        c = np.polyfit(xn, yy, deg)
        base = np.polyval(c, xn)
        new = np.minimum(yy, base)
        if np.linalg.norm(new - yy) / (np.linalg.norm(yy) + 1e-12) < tol:
            yy = new
            break
        yy = new
    return base


def baseline(s: Spectrum, method: str = "als", **kw):
    """Retorna (espectro_corrigido, linha_de_base_como_Spectrum).
    method:
      'als'        lam (1e5-1e8, maior = mais rígida), p (0.001-0.05)  -> bom padrão geral
      'rubberband' sem parâmetros; ótimo para ATR/KBr com linha de base curva e bandas positivas
      'poly'       deg=3; polinômio modificado iterativo
      'points'     anchors=[cm-1, ...] pontos sem banda, interpolação linear (controle total;
                   o que você usaria para reproduzir uma correção feita 'à mão')
    Trabalha em absorbância. Se receber T/%T, converte antes (e avisa no meta)."""
    if s.ykind in ("T", "%T", "R", "%R"):
        s = s.to_absorbance()
    x, y, back = _asc(s)
    if method == "als":
        b = _als(y, lam=kw.get("lam", 1e6), p=kw.get("p", 0.01), niter=kw.get("niter", 15))
    elif method == "rubberband":
        b = _rubberband(x, y)
    elif method == "poly":
        b = _modpoly(x, y, deg=kw.get("deg", 3))
    elif method == "points":
        anchors = np.sort(np.asarray(kw["anchors"], float))
        idx = [int(np.argmin(np.abs(x - a))) for a in anchors]
        # média em ±half_cm para robustez ao ruído
        half = kw.get("half_cm", 3.0)
        ya = [np.mean(y[(x >= x[i] - half) & (x <= x[i] + half)]) for i in idx]
        b = np.interp(x, x[idx], ya)
    else:
        raise ValueError(f"método de linha de base desconhecido: {method}")
    corr = s.copy(y=back(y - b))
    corr.meta["baseline"] = {"method": method, **{k: v for k, v in kw.items() if k != "anchors"}}
    return corr, s.copy(y=back(b), title=s.title + " (baseline)")


# ---------------------------------------------------------------- suavização / derivadas
def smooth(s: Spectrum, window_cm: float = 9.0, poly: int = 2) -> Spectrum:
    """Savitzky-Golay. Mantenha a janela MENOR que a largura das bandas mais estreitas
    (bandas de sólidos ~10-20 cm-1 FWHM; gases muito mais estreitas)."""
    n = _pts(s, window_cm)
    out = s.copy(y=savgol_filter(s.y, n, poly))
    out.meta["smoothing"] = f"SG {n} pts ({window_cm} cm-1), poly {poly}"
    return out


def derivative(s: Spectrum, order: int = 2, window_cm: float = 13.0, poly: int = 3) -> Spectrum:
    """Derivada Savitzky-Golay em relação a x (cm-1). Na 2ª derivada, bandas viram mínimos
    negativos - útil para separar ombros (ex.: amida I, C=O sobrepostas)."""
    n = _pts(s, window_cm)
    d = savgol_filter(s.y, n, poly, deriv=order, delta=float(np.mean(np.diff(s.x))))
    out = s.copy(y=d, ykind=f"d{order}")
    out.meta["derivative"] = f"SG d{order}, {n} pts, poly {poly}"
    return out


# ---------------------------------------------------------------- normalização
def normalize(s: Spectrum, method: str = "max", band=None, region=None) -> Spectrum:
    """'max' -> maior valor = 1 (na região, se dada); 'minmax' -> 0..1;
    'area' -> área = 1; 'band' -> altura em `band=(hi, lo)` = 1 (banda de referência interna,
    ex.: 1450 cm-1 CH2 em polímeros); 'snv' -> média 0 desvio 1; 'vector' -> norma 2 = 1."""
    y = s.y.astype(float)
    m = np.ones_like(y, bool)
    if region is not None:
        hi, lo = max(region), min(region)
        m = (s.x <= hi) & (s.x >= lo)
    if method == "max":
        f = np.nanmax(y[m])
        y = y / f
    elif method == "minmax":
        mn, mx = np.nanmin(y[m]), np.nanmax(y[m])
        y = (y - mn) / (mx - mn)
    elif method == "area":
        x, yy, _ = _asc(s.copy(y=y))
        f = np.trapezoid(np.abs(yy), x) if hasattr(np, "trapezoid") else np.trapz(np.abs(yy), x)
        y = y / f
    elif method == "band":
        hi, lo = max(band), min(band)
        mb = (s.x <= hi) & (s.x >= lo)
        y = y / np.nanmax(y[mb])
    elif method == "snv":
        y = (y - np.nanmean(y[m])) / np.nanstd(y[m])
    elif method == "vector":
        y = y / np.linalg.norm(y[m])
    else:
        raise ValueError(method)
    out = s.copy(y=y)
    out.meta["normalization"] = method + (f" {band}" if band else "")
    return out


# ---------------------------------------------------------------- ATR, subtração, atmosfera
def atr_correction(s: Spectrum, ref: float = 1000.0) -> Spectrum:
    """Correção ATR simples: A_corr = A * (x / ref). Compensa a profundidade de penetração
    ∝ comprimento de onda (bandas de baixo número de onda ficam artificialmente fortes no ATR).
    NÃO corrige deslocamentos de banda por dispersão anômala (ATR 'avançado' precisa de n do
    cristal e da amostra). Use só para comparar com espectros de transmissão/biblioteca."""
    a = s.to_absorbance()
    out = a.copy(y=a.y * (a.x / ref))
    out.meta["atr_correction"] = f"simple, ref={ref} cm-1"
    return out


def subtract(a: Spectrum, b: Spectrum, factor: float | None = None, region=None, step: float | None = None):
    """a - f*b em absorbância, numa grade comum. Se factor=None, f é ajustado por mínimos
    quadrados em `region=(hi, lo)` (escolha uma região onde só `b` absorve, ex.: banda do
    solvente/matriz). Retorna (diferença, f)."""
    A, B = a.to_absorbance(), b.to_absorbance()
    st = step or max(A.step, B.step)
    hi = min(A.xrange[0], B.xrange[0])
    lo = max(A.xrange[1], B.xrange[1])
    A, B = resample(A, st, hi, lo), resample(B, st, hi, lo)
    if factor is None:
        m = np.ones_like(A.x, bool)
        if region is not None:
            m = (A.x <= max(region)) & (A.x >= min(region))
        factor = float(np.dot(A.y[m], B.y[m]) / np.dot(B.y[m], B.y[m]))
    out = A.copy(y=A.y - factor * B.y, title=f"{a.title} - {factor:.3g}×{b.title}")
    out.meta["subtraction"] = {"factor": factor, "region": region}
    return out, factor


def atmospheric_mask(x) -> np.ndarray:
    """True onde x cai em região de CO2 ou vapor d'água."""
    x = np.asarray(x)
    m = np.zeros(x.shape, bool)
    for hi, lo in CO2_REGIONS + H2O_REGIONS:
        m |= (x <= hi) & (x >= lo)
    return m


def atmospheric_flag(pos: float) -> str:
    for hi, lo in CO2_REGIONS:
        if lo <= pos <= hi:
            return "CO2?"
    for hi, lo in H2O_REGIONS:
        if lo <= pos <= hi:
            return "H2O(v)?"
    return ""


def co2_flatten(s: Spectrum, hi: float = 2400, lo: float = 2280) -> Spectrum:
    """Substitui 2400-2280 cm-1 por interpolação linear (remove o dubleto do CO2 só para
    visualização - informe isso na legenda/figura)."""
    y = s.y.copy()
    m = (s.x <= hi) & (s.x >= lo)
    if m.sum() < 2:
        return s.copy()
    x_out = s.x[~m]
    order = np.argsort(x_out)
    y[m] = np.interp(s.x[m], x_out[order], y[~m][order])
    out = s.copy(y=y)
    out.meta["co2_flattened"] = (hi, lo)
    return out


def noise_level(s: Spectrum) -> float:
    """Ruído RMS: menor desvio (após ajuste quadrático) entre janelas que costumam não ter
    bandas (2700-2400, 2200-1900, 3900-3750, 1900-1800); se nenhuma existir, estimativa
    robusta pela 2ª diferença. Tomar o mínimo evita superestimar quando uma janela tem banda
    (C≡N, C≡C) ou o ruído do cristal de diamante do ATR (~2300-1900)."""
    vals = []
    for hi, lo in ((2700, 2400), (2200, 1900), (3900, 3750), (1900, 1800)):
        m = (s.x <= hi) & (s.x >= lo)
        if m.sum() > 30:
            xs, seg = s.x[m], s.y[m]
            c = np.polyfit(xs - xs.mean(), seg, 2)
            vals.append(float(np.std(seg - np.polyval(c, xs - xs.mean()))))
    if vals:
        return min(vals)
    d2 = np.diff(s.y, 2)
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6))
