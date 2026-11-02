"""
sem_film.py - morfologia de filmes (amido, quitosana, alginato, pectina, CMC, celulose, blendas e
compósitos) por MEV: superfície e corte transversal (fratura)

    from sem_io import load
    from sem_film import surface_texture, surface_defects, cross_section
    im = load("filme_sup.tif")
    tx = surface_texture(im)            # CV de brilho, GLCM, comprimento de correlação, D fractal (PSD)
    df = surface_defects(im)            # poros, partículas/grânulos residuais, trincas (máscaras + tabelas)
    cs = cross_section(load("corte.tif"))   # espessura (média ± DP ao longo do filme), porosidade interna

IMPORTANTE: o brilho de MEV (SE) não é altura. Textura/“rugosidade” daqui é um índice COMPARATIVO
entre imagens adquiridas nas mesmas condições (kV, detector, aumento, contraste/brilho); rugosidade
real (Sa, Sq) exige AFM, perfilometria ou estereo-MEV.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from sem_process import flatten, gray_close, gray_open, um2px


# ----------------------------------------------------------------------------- textura

def _glcm(a, distances_px, levels=64):
    from skimage.feature import graycomatrix, graycoprops
    q = np.clip((a - np.percentile(a, 0.5)) / (np.percentile(a, 99.5) - np.percentile(a, 0.5) + 1e-9), 0, 1)
    q = np.round(q * (levels - 1)).astype(np.uint8)
    d = sorted(set(max(1, int(round(x))) for x in distances_px))
    P = graycomatrix(q, d, [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4], levels=levels, symmetric=True, normed=True)
    out = {}
    for prop in ("contrast", "homogeneity", "energy", "correlation", "dissimilarity"):
        out[prop] = graycoprops(P, prop).mean(axis=1)          # média nas 4 direções
    ent = []
    for i in range(len(d)):
        p = P[:, :, i, :].mean(axis=-1)
        p = p[p > 0]
        ent.append(float(-(p * np.log2(p)).sum()))
    out["entropy"] = np.array(ent)
    return d, out


def _radial_psd(a):
    a = a - a.mean()
    H, W = a.shape
    win = np.outer(np.hanning(H), np.hanning(W))
    F = np.fft.fftshift(np.abs(np.fft.fft2(a * win)) ** 2)
    yy, xx = np.indices(F.shape)
    r = np.hypot(yy - H / 2, xx - W / 2).astype(int)
    psd = np.bincount(r.ravel(), F.ravel()) / np.maximum(np.bincount(r.ravel()), 1)
    k = np.arange(psd.size) / max(H, W)            # ciclos/px
    return k[1:min(H, W) // 2], psd[1:min(H, W) // 2]


def _corr_length(a):
    """Comprimento de correlação (px): distância em que a autocorrelação radial cai para 1/e."""
    a = a - a.mean()
    F = np.fft.fft2(a)
    ac = np.fft.fftshift(np.real(np.fft.ifft2(np.abs(F) ** 2)))
    ac /= ac.max()
    H, W = a.shape
    yy, xx = np.indices(ac.shape)
    r = np.hypot(yy - H // 2, xx - W // 2).astype(int)
    prof = np.bincount(r.ravel(), ac.ravel()) / np.maximum(np.bincount(r.ravel()), 1)
    below = np.where(prof < 1 / np.e)[0]
    if not below.size:
        return np.nan
    i = below[0]
    return float(i - 1 + (prof[i - 1] - 1 / np.e) / (prof[i - 1] - prof[i] + 1e-12))


def surface_texture(im, distances_um=(0.5, 1, 2, 5), do_flatten: bool = True, mask=None) -> dict:
    """Índices de textura da superfície (comparativos). mask: região a usar (True)."""
    a = im.data.astype(np.float32)
    if do_flatten:
        a = flatten(a)
    sub = a if mask is None else np.where(mask, a, np.nan)
    v = sub[np.isfinite(sub)]
    out = dict(media=float(v.mean()), dp=float(v.std()), cv_pct=float(100 * v.std() / (v.mean() + 1e-9)),
               p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)))
    dpx = [um2px(im, d, d) for d in distances_um]
    dd, g = _glcm(a, dpx)
    out["glcm_dist"] = [round(d * (im.px_um or 1), 4) for d in dd]
    for k, val in g.items():
        out[f"glcm_{k}"] = [round(float(x), 4) for x in val]
    lc = _corr_length(a)
    out["comprimento_correlacao"] = lc * (im.px_um or 1) if np.isfinite(lc) else np.nan
    k, psd = _radial_psd(a)
    # inclinação da PSD na faixa intermediária (evita ruído de pixel e janela)
    sel = (k > 4 / max(a.shape)) & (k < 0.15)
    if sel.sum() > 5:
        beta = -np.polyfit(np.log(k[sel]), np.log(psd[sel] + 1e-30), 1)[0]
        out["psd_beta"] = float(beta)
        out["D_fractal_psd"] = float(np.clip((8 - beta) / 2, 2, 3)) if 2 < beta < 4 else np.nan
    out["unidade"] = "µm" if im.px_um else "px"
    return out


# ----------------------------------------------------------------------------- defeitos

def _blobs(a, polarity, rmax_px, k, min_d_px, mask=None):
    r = max(2, int(round(rmax_px)))
    if polarity == "dark":
        th = gray_close(a, r) - a
    else:
        th = a - gray_open(a, r)
    valid = np.ones_like(a, bool) if mask is None else mask
    v = th[valid]
    med = np.median(v)
    mad = np.median(np.abs(v - med)) * 1.4826 + 1e-6
    m = (th > med + k * mad) & valid
    lab, n = ndi.label(m)
    if n:
        sizes = np.bincount(lab.ravel())
        keep = sizes >= np.pi / 4 * min_d_px ** 2
        keep[0] = False
        m = keep[lab]
    return m, th


def _shape_filter(m, im, px, H, W):
    """Mantém só manchas compactas (AR ≤ 3, solidez ≥ 0.6) e devolve (máscara, tabela)."""
    from skimage import measure as skm
    lab, n = ndi.label(m)
    rows = []
    keep = np.zeros(n + 1, bool)
    for p in skm.regionprops(lab, intensity_image=im.data):
        el = p.axis_major_length / max(p.axis_minor_length, 1.0)
        if el > 3.0 or p.solidity < 0.6:
            continue
        keep[p.label] = True
        r0, c0, r1, c1 = p.bbox
        edge = r0 == 0 or c0 == 0 or r1 == H or c1 == W
        rows.append(dict(id=len(rows) + 1, ecd=np.sqrt(4 * p.area / np.pi) * px, area=p.area * px ** 2,
                         AR=el, circularidade=min(1, 4 * np.pi * p.area / max(p.perimeter, 1) ** 2),
                         borda=edge, x_px=p.centroid[1], y_px=p.centroid[0]))
    return keep[lab], pd.DataFrame(rows)


def surface_defects(im, pore_max_um=None, particle_max_um=None, crack_min_um=None, k: float = 5.0,
                    min_d_um=None, do_flatten: bool = True) -> dict:
    """Defeitos de superfície:
      poros/furos   - manchas escuras compactas (bolhas, poros, pites)
      partículas    - protuberâncias claras compactas (grânulos não gelatinizados, cristais, aglomerados
                      de carga, sal/plastificante aflorado)
      trincas       - estruturas escuras finas e longas (filtro de cristas de Sato)
    Tamanhos padrão relativos ao campo (2% da largura para poros, 5% para partículas, trinca ≥ 3% da
    largura). Retorna dict com tabelas (pandas), máscaras e resumo (densidades por mm²)."""
    from skimage import measure as skm
    from skimage.filters import sato
    from skimage.morphology import skeletonize
    a = im.data.astype(np.float32)
    if do_flatten:
        a = flatten(a)
    a = ndi.gaussian_filter(a, 0.8)
    H, W = a.shape
    px = im.px_um or 1.0
    pore_r = um2px(im, pore_max_um, 0.02 * W) / 2 if pore_max_um else 0.01 * W
    part_r = um2px(im, particle_max_um, 0.05 * W) / 2 if particle_max_um else 0.025 * W
    crack_min = um2px(im, crack_min_um, 0.05 * W) if crack_min_um else 0.05 * W
    min_d = um2px(im, min_d_um, 4) if min_d_um else 4
    # partículas claras primeiro: a sombra/anel escuro ao redor delas não é poro nem trinca
    qm_raw, _ = _blobs(a, "bright", part_r, k, min_d)
    qm0 = ndi.binary_fill_holes(qm_raw)
    # anel claro com interior escuro = borda de um poro, não partícula
    holes = qm0 & ~qm_raw
    lab_q, nq = ndi.label(qm0)
    if nq:
        idx = np.arange(1, nq + 1)
        h_area = np.asarray(ndi.sum(holes, lab_q, idx))
        h_mean = np.asarray(ndi.sum(np.where(holes, a, 0), lab_q, idx)) / np.maximum(h_area, 1)
        tot = np.asarray(ndi.sum(qm0, lab_q, idx))
        dark_lim = np.percentile(a, 5)
        rim = (h_area > 0.15 * tot) & (h_mean < dark_lim)
        keepq = np.ones(nq + 1, bool)
        keepq[0] = False
        keepq[1:][rim] = False
        qm0 = keepq[lab_q]
    qm0, ptab = _shape_filter(qm0, im, px, H, W)
    halo = ndi.binary_dilation(qm0, iterations=max(3, int(0.3 * part_r)))
    # trincas (antes dos poros, para não contá-las como poros)
    sig = [1, 1.5, 2, 3]
    rid = sato(a, sigmas=sig, black_ridges=True)
    rv = rid.ravel()
    med = np.median(rv)
    mad = np.median(np.abs(rv - med)) * 1.4826 + 1e-9
    cm = rid > med + 2 * k * mad
    # trinca tem de ser escura de verdade em relação à vizinhança
    loc = ndi.uniform_filter(a, int(max(9, 6 * max(sig))))
    cm &= a < loc - 0.5 * np.std(a)
    clab, nc = ndi.label(cm, structure=np.ones((3, 3)))
    cracks = []
    crack_mask = np.zeros_like(cm)
    total_len = 0.0
    for p in skm.regionprops(clab):
        sk = skeletonize(p.image)
        L = float(sk.sum())
        el = p.axis_major_length / max(p.axis_minor_length, 1.0)
        width = p.area / max(L, 1)
        r0_, c0_, r1_, c1_ = p.bbox
        in_halo = (halo[r0_:r1_, c0_:c1_] & p.image).sum() / max(p.area, 1)
        if in_halo > 0.3:
            continue
        if L >= crack_min and L / max(width, 1) >= 8 and (el >= 4 or p.solidity < 0.4):
            r0, c0, r1, c1 = p.bbox
            crack_mask[r0:r1, c0:c1] |= p.image
            width = p.area / max(L, 1)
            cracks.append(dict(id=len(cracks) + 1, comprimento=L * px, largura_media=width * px,
                               x_px=p.centroid[1], y_px=p.centroid[0]))
            total_len += L
    # poros escuros e partículas claras
    pm, _ = _blobs(a, "dark", pore_r, k, min_d, mask=~(ndi.binary_dilation(crack_mask, iterations=3) | halo))
    qm = qm0
    tabs = {"particulas": ptab}
    pm, tabs["poros"] = _shape_filter(pm, im, px, H, W)
    tabs["trincas"] = pd.DataFrame(cracks)
    area_mm2 = (H * W * px ** 2) * 1e-6 if im.px_um else np.nan
    summ = dict(area_analisada_mm2=area_mm2,
                n_poros=len(tabs["poros"]), n_particulas=len(tabs["particulas"]), n_trincas=len(cracks),
                poros_area_pct=float(100 * pm.mean()), particulas_area_pct=float(100 * qm.mean()),
                trincas_comprimento_total=total_len * px)
    if im.px_um:
        summ["poros_por_mm2"] = summ["n_poros"] / area_mm2
        summ["particulas_por_mm2"] = summ["n_particulas"] / area_mm2
        summ["densidade_trincas_mm_por_mm2"] = (total_len * px * 1e-3) / area_mm2
    for name in ("poros", "particulas"):
        df = tabs[name]
        if len(df):
            summ[f"{name}_ecd_media"] = float(df["ecd"].mean())
            summ[f"{name}_ecd_dp"] = float(df["ecd"].std(ddof=1)) if len(df) > 1 else 0.0
    return dict(resumo=summ, tabelas=tabs, mascaras=dict(poros=pm, particulas=qm, trincas=crack_mask),
                unidade="µm" if im.px_um else "px")


# ----------------------------------------------------------------------------- corte transversal

def cross_section(im, orientation: str = "auto", film_brighter: bool | None = None, smooth_um=None,
                  min_cover: float = 0.6, pore_k: float = 6.0) -> dict:
    """Espessura de filme em imagem de fratura/corte. O filme é a faixa contínua (mais clara e/ou mais
    texturizada) que atravessa a imagem; mede a distância entre as bordas superior e inferior em cada
    coluna (perpendicular à faixa, corrigindo a inclinação). orientation: 'auto' | 'h' (faixa
    horizontal) | 'v'. Retorna dict com espessura média ± DP, perfil, bordas e porosidade interna."""
    a = im.data.astype(np.float32)
    rotated = False
    if orientation == "auto":
        # faixa horizontal => perfil de média por linha varia muito mais que por coluna
        rv = np.var(ndi.gaussian_filter1d(a.mean(axis=1), 3))
        cv = np.var(ndi.gaussian_filter1d(a.mean(axis=0), 3))
        orientation = "h" if rv >= cv else "v"
    if orientation == "v":
        a = a.T
        rotated = True
    H, W = a.shape
    px = im.px_um or 1.0
    s = um2px(im, smooth_um, 3) if smooth_um else max(2.0, 0.004 * max(H, W))
    m = ndi.gaussian_filter(a, s)
    tex = np.sqrt(np.maximum(ndi.gaussian_filter(a ** 2, s * 2) - ndi.gaussian_filter(a, s * 2) ** 2, 0))
    from skimage.filters import threshold_otsu
    feats = []
    for f in (m, tex):
        f = (f - f.min()) / (f.max() - f.min() + 1e-9)
        feats.append(f)
    # combina brilho e textura (o filme costuma ser mais claro E mais texturizado que o vazio/fita)
    comb = 0.6 * feats[0] + 0.4 * feats[1]
    if film_brighter is False:
        comb = 1 - feats[0] * 0.6 + 0.4 * feats[1]
    th = threshold_otsu(comb)
    fm = comb > th
    fm = ndi.binary_closing(fm, structure=np.ones((int(3 * s) | 1, int(3 * s) | 1)))
    fm = ndi.binary_fill_holes(fm)
    lab, n = ndi.label(fm)
    if n == 0:
        raise ValueError("não encontrei a faixa do filme - recorte a região (--crop) ou informe a orientação")
    cover = np.array([np.unique(np.where(lab == i)[1]).size / W for i in range(1, n + 1)])
    areas = np.bincount(lab.ravel())[1:]
    score = cover * areas
    best = int(np.argmax(score)) + 1
    band = lab == best
    if cover[best - 1] < min_cover:
        warn = [f"a faixa detectada cobre só {100 * cover[best - 1]:.0f}% da largura - confira o recorte"]
    else:
        warn = []
    top = np.full(W, np.nan)
    bot = np.full(W, np.nan)
    for x in range(W):
        ys = np.where(band[:, x])[0]
        if ys.size:
            top[x], bot[x] = ys[0], ys[-1] + 1
    ok = np.isfinite(top) & (top > 1) & (bot < H - 1)
    if ok.sum() < 0.3 * W:
        warn.append("a faixa toca a borda da imagem em boa parte das colunas: espessura pode estar truncada "
                    "(use aumento menor)")
    # refina bordas com o gradiente local (sub-pixel ~ 1 px)
    g = ndi.gaussian_filter(a, 1.5)
    gy = ndi.sobel(g, axis=0)
    for arr, sign in ((top, 1), (bot, -1)):
        for x in np.where(ok)[0]:
            y0 = int(arr[x])
            lo, hi = max(1, y0 - int(2 * s)), min(H - 1, y0 + int(2 * s) + 1)
            seg = sign * gy[lo:hi, x] * (1 if np.mean(a[band]) > np.mean(a[~band]) else -1)
            if seg.size:
                arr[x] = lo + int(np.argmax(seg))
    # suaviza bordas ao longo de x e corrige inclinação (espessura perpendicular)
    xs = np.arange(W)
    tt = top.copy()
    bb = bot.copy()
    for arr in (tt, bb):
        good = np.isfinite(arr) & ok
        if good.sum() > 10:
            arr[:] = np.interp(xs, xs[good], ndi.median_filter(arr[good], size=max(7, int(W / 60) | 1)))
    mid = 0.5 * (tt + bb)
    slope = np.gradient(ndi.gaussian_filter1d(mid, max(5, W / 50)))
    vert = (bb - tt)
    thick = vert * np.cos(np.arctan(slope))
    t_ok = thick[ok]
    # porosidade interna: manchas escuras compactas dentro da faixa (longe das bordas)
    inner = band & ~ndi.binary_dilation(~band, iterations=int(max(4, 3 * s)))
    pores = np.zeros_like(band)
    pore_pct = np.nan
    pore_tab = pd.DataFrame()
    if inner.sum() > 100:
        r_max = max(3.0, 0.15 * float(np.nanmedian(vert[ok])) if ok.any() else 10.0)
        pm, _ = _blobs(ndi.gaussian_filter(a, 0.8), "dark", r_max, pore_k, 4, mask=inner)
        pores, pore_tab = _shape_filter(pm & inner, im, px, H, W)
        pore_pct = float(100 * pores.sum() / max(inner.sum(), 1))
    res = dict(n_colunas=int(ok.sum()), espessura_media=float(np.mean(t_ok)) * px if t_ok.size else np.nan,
               espessura_dp=float(np.std(t_ok, ddof=1)) * px if t_ok.size > 1 else np.nan,
               espessura_min=float(np.min(t_ok)) * px if t_ok.size else np.nan,
               espessura_max=float(np.max(t_ok)) * px if t_ok.size else np.nan,
               espessura_p5_p95=(float(np.percentile(t_ok, 5) * px), float(np.percentile(t_ok, 95) * px))
               if t_ok.size else (np.nan, np.nan),
               inclinacao_media_graus=float(np.degrees(np.arctan(np.median(slope[ok])))) if ok.any() else np.nan,
               porosidade_area_pct=pore_pct, n_poros=int(len(pore_tab)),
               poros_ecd_media=float(pore_tab["ecd"].mean()) if len(pore_tab) else np.nan,
               orientacao=orientation, unidade="µm" if im.px_um else "px",
               top=tt, bot=bb, ok=ok, perfil=thick * px, band=band, pores=pores, rotated=rotated, avisos=warn)
    if res["espessura_media"] and res["espessura_dp"] / res["espessura_media"] > 0.15:
        warn.append("espessura varia > 15% ao longo do campo (filme irregular, fratura inclinada ou faixa mal "
                    "detectada) - confira a figura")
    return res


def measure_line(im, x1, y1, x2, y2) -> float:
    """Distância entre dois pontos (px) em µm."""
    return float(np.hypot(x2 - x1, y2 - y1) * (im.px_um or 1.0))
