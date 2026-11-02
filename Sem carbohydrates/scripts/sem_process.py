"""
sem_process.py - pré-processamento e qualidade de imagens de MEV

    from sem_process import quality, flatten, denoise, enhance, um2px
    q = quality(im)                     # ruído, foco, saturação, carregamento, gradiente de iluminação
    f = flatten(im.data, sigma_px=200)   # remove gradiente lento de iluminação (não usar p/ medir brilho)
    d = denoise(im.data, 'nlm')          # 'gauss' | 'median' | 'nlm' (non-local means) | 'bilateral'
    e = enhance(im.data)                 # CLAHE - só para visualização
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi


def um2px(im, value_um, default_px=None):
    """Converte µm em pixels com a escala da imagem (ou devolve default_px se não calibrada)."""
    if value_um is None:
        return default_px
    if im.px_um:
        return float(value_um) / im.px_um
    return default_px if default_px is not None else float(value_um)


def noise_sigma(a: np.ndarray) -> float:
    """Estimativa robusta do desvio do ruído (MAD do Laplaciano, Immerkær)."""
    k = np.array([[1, -2, 1], [-2, 4, -2], [1, -2, 1]], float)
    c = ndi.convolve(a.astype(float), k)
    return float(np.sqrt(np.pi / 2) * np.mean(np.abs(c[1:-1, 1:-1])) / 6.0)


def focus_measure(a: np.ndarray) -> float:
    """Nitidez: variância do Laplaciano normalizada pelo ruído (comparativo entre imagens similares)."""
    lap = ndi.laplace(ndi.gaussian_filter(a.astype(float), 1.0))
    return float(lap.var() / (a.var() + 1e-9))


def charging_index(a: np.ndarray):
    """Indicadores de carregamento: fração saturada, linhas horizontais anômalas (riscas de varredura)
    e regiões brancas grandes. Retorna dict."""
    sat = a >= 0.98
    row_mean = a.mean(axis=1)
    d = np.abs(np.diff(row_mean))
    mad = np.median(np.abs(d - np.median(d))) + 1e-6
    streak_rows = int((d > np.median(d) + 12 * mad).sum())
    lab, n = ndi.label(sat)
    big = 0
    if n:
        sizes = ndi.sum(sat, lab, index=np.arange(1, n + 1))
        big = float(sizes.max() / a.size)
    return dict(saturado_pct=100 * float(sat.mean()), riscas=streak_rows, maior_regiao_saturada_pct=100 * big)


def illumination_gradient(a: np.ndarray) -> float:
    """Diferença (%) entre o brilho médio das metades opostas após forte suavização (máx. de esq/dir,
    cima/baixo). > ~15% sugere iluminação/detector desigual (sombra do detector ETD, carregamento)."""
    s = ndi.gaussian_filter(a.astype(float), max(a.shape) / 8)
    H, W = s.shape
    lr = abs(s[:, : W // 2].mean() - s[:, W // 2:].mean())
    tb = abs(s[: H // 2].mean() - s[H // 2:].mean())
    return float(100 * max(lr, tb) / (a.mean() + 1e-9))


def quality(im) -> dict:
    a = im.data
    q = dict(ruido_sigma=noise_sigma(a), snr=float(a.std() / (noise_sigma(a) + 1e-9)), foco=focus_measure(a),
             gradiente_iluminacao_pct=illumination_gradient(a), contraste_p1_p99=(float(np.percentile(a, 1)),
                                                                                    float(np.percentile(a, 99))))
    q.update(charging_index(a))
    notes = []
    if q["snr"] < 3:
        notes.append("imagem ruidosa (SNR baixo): aumente tempo de varredura/integração; suavize antes de segmentar")
    if q["saturado_pct"] > 2 or q["maior_regiao_saturada_pct"] > 0.5:
        notes.append("regiões saturadas - possível carregamento (amostra não condutora): recobrimento, kV menor, "
                     "baixo vácuo; não quantifique dentro delas")
    if q["riscas"] > 3:
        notes.append(f"{q['riscas']} saltos de brilho entre linhas (riscas de varredura: carregamento ou instabilidade)")
    if q["gradiente_iluminacao_pct"] > 15:
        notes.append(f"gradiente de iluminação {q['gradiente_iluminacao_pct']:.0f}% (sombra do detector/carregamento): "
                     "use limiar local ou --flatten")
    lo, hi = q["contraste_p1_p99"]
    if hi - lo < 0.25:
        notes.append("contraste baixo (faixa p1-p99 < 25% da escala)")
    q["notas"] = notes
    return q


def flatten(a: np.ndarray, sigma_px: float | None = None) -> np.ndarray:
    """Remove gradiente lento de iluminação dividindo por uma versão muito suavizada (mantém a média).
    sigma padrão = 1/6 do maior lado. Use para segmentar, não para comparar brilho absoluto."""
    a = a.astype(np.float32)
    s = sigma_px or max(a.shape) / 6
    bg = ndi.gaussian_filter(a, s) + 1e-3
    out = a / bg * float(a.mean())
    return np.clip(out, 0, 1).astype(np.float32)


def denoise(a: np.ndarray, method: str = "gauss", strength: float = 1.0) -> np.ndarray:
    a = a.astype(np.float32)
    if method in (None, "none"):
        return a
    if method == "gauss":
        return ndi.gaussian_filter(a, strength)
    if method == "median":
        return ndi.median_filter(a, size=max(3, int(2 * round(strength) + 1)))
    if method == "nlm":
        from skimage.restoration import denoise_nl_means, estimate_sigma
        sig = float(estimate_sigma(a))
        return denoise_nl_means(a, h=0.8 * sig * strength, sigma=sig, patch_size=5, patch_distance=6,
                                fast_mode=True).astype(np.float32)
    if method == "bilateral":
        from skimage.restoration import denoise_bilateral
        return denoise_bilateral(a, sigma_color=0.08 * strength, sigma_spatial=2).astype(np.float32)
    raise ValueError(f"método de suavização desconhecido: {method}")


def enhance(a: np.ndarray, clip: float = 0.01) -> np.ndarray:
    """CLAHE para visualização (não meça brilho depois disso)."""
    from skimage import exposure
    return exposure.equalize_adapthist(np.clip(a, 0, 1), clip_limit=clip).astype(np.float32)


def autolevel(a: np.ndarray, p_lo: float = 0.5, p_hi: float = 99.5) -> np.ndarray:
    lo, hi = np.percentile(a, [p_lo, p_hi])
    return np.clip((a - lo) / (hi - lo + 1e-9), 0, 1).astype(np.float32)


def threshold(a: np.ndarray, method="otsu", block_px: int | None = None):
    """Limiar global ('otsu', 'li', 'triangle', 'yen', 'mean' ou número 0..1) ou local ('local', bloco
    em px). Retorna (valor ou mapa de limiar, descrição)."""
    from skimage import filters
    if isinstance(method, (int, float)) or (isinstance(method, str) and method.replace(".", "", 1).isdigit()):
        return float(method), f"manual {float(method):.3f}"
    m = str(method).lower()
    if m == "local":
        b = int(block_px or max(31, (min(a.shape) // 8) | 1))
        b = b | 1
        return filters.threshold_local(a, b, method="gaussian", offset=-0.01), f"local (bloco {b} px)"
    f = {"otsu": filters.threshold_otsu, "li": filters.threshold_li, "triangle": filters.threshold_triangle,
         "yen": filters.threshold_yen, "mean": filters.threshold_mean}.get(m)
    if f is None:
        raise ValueError(f"limiar desconhecido: {method}")
    return float(f(a)), f"{m} {float(f(a)):.3f}"


def gray_close(a: np.ndarray, radius: float) -> np.ndarray:
    """Fechamento em tons de cinza com disco (OpenCV, rápido)."""
    import cv2
    r = max(1, int(round(radius)))
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    return cv2.morphologyEx(a.astype(np.float32), cv2.MORPH_CLOSE, k)


def gray_open(a: np.ndarray, radius: float) -> np.ndarray:
    import cv2
    r = max(1, int(round(radius)))
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    return cv2.morphologyEx(a.astype(np.float32), cv2.MORPH_OPEN, k)
