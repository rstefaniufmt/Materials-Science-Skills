"""
sem_granules.py - segmentação e morfometria de grânulos/partículas (amido, microcápsulas,
partículas de spray-drying, quitosana, celulose microcristalina...)

    from sem_io import load
    from sem_granules import segment, measure, size_stats, bimodal, shape_summary, pits
    im = load("milho.tif")
    seg = segment(im, min_d_um=2)                 # dict: labels, mask, info
    df = measure(im, seg["labels"])               # tabela por grânulo (µm)
    st = size_stats(df)                           # n, média, DP, D10/D50/D90 (número e volume), D[4,3]...
    bm = bimodal(df)                              # duas populações? (amido de trigo/cevada: A e B)
    pits(im, seg["labels"], df)                   # poros/pites superficiais (hidrólise enzimática)

Definições (ISO 9276-6 / ISO 13322-1):
  ECD (d_eq)      diâmetro do círculo de mesma área projetada
  Feret máx/mín   maior/menor distância entre retas paralelas tangentes ao contorno
  AR              razão de aspecto = Feret máx / Feret mín
  circularidade   4πA/P²  (1 = círculo; sensível à rugosidade do contorno e à resolução)
  arredondamento  4A/(π·Feret_máx²)
  solidez         A / A_envoltória_convexa (baixa = reentrâncias, fragmentos, agregados)
Objetos que tocam a borda são excluídos das estatísticas (tamanho truncado), mas contados.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import ndimage as ndi

from sem_process import denoise, flatten, gray_close, threshold, um2px


# ----------------------------------------------------------------------------- segmentação

def segment(im, method="otsu", dark_objects: bool = False, min_d_um=None, max_d_um=None, smooth: float = 1.0,
            do_flatten: bool = False, split: bool = True, split_h: float = 0.25, fill_holes: bool = True,
            block_um=None, open_px: int = 1) -> dict:
    """Segmenta objetos claros sobre fundo escuro (SE típico: grânulos sobre fita de carbono).
    dark_objects=True para o inverso. split: separa grânulos encostados por watershed na transformada
    de distância (marcadores = máximos com proeminência split_h × raio). Retorna dict(labels, mask, info)."""
    from skimage import morphology, segmentation
    a = im.data.astype(np.float32)
    if do_flatten:
        a = flatten(a)
    a = denoise(a, "gauss", smooth) if smooth else a
    if dark_objects:
        a = 1.0 - a
    min_px = um2px(im, min_d_um, None)
    if min_px is None:
        min_px = 6.0
    max_px = um2px(im, max_d_um, None)
    block_px = um2px(im, block_um, None)
    th, th_desc = threshold(a, method, int(block_px) if block_px else None)
    fg = a > th
    if open_px:
        fg = ndi.binary_opening(fg, structure=morphology.disk(open_px))
    bg_level = float(np.median(a[~fg])) if (~fg).any() else 0.0
    fg_level = float(np.median(a[fg])) if fg.any() else 1.0
    holes = np.zeros_like(fg)
    if fill_holes:
        filled = ndi.binary_fill_holes(fg)
        holes = filled & ~fg
        fg = filled
    fg = _remove_small(fg, max(4, int(np.pi / 4 * min_px ** 2 * 0.5)))
    if split:
        dist = ndi.distance_transform_edt(fg)
        dist_s = ndi.gaussian_filter(dist, 1.0)
        r_typ = max(2.0, min_px / 2)
        h = max(1.5, split_h * r_typ)
        lab0, n0 = ndi.label(fg)
        # proeminência proporcional ao raio de cada objeto (evita supersegmentar objetos grandes)
        rmax = ndi.maximum(dist_s, lab0, index=np.arange(1, n0 + 1)) if n0 else np.array([])
        hmap = np.zeros(n0 + 1)
        hmap[1:] = np.maximum(h, split_h * rmax)
        hloc = hmap[lab0]
        markers = _hmax_markers(dist_s, fg, hloc)
        labels = segmentation.watershed(-dist_s, markers, mask=fg)
    else:
        labels, _ = ndi.label(fg)
    if fill_holes and holes.any():
        # buraco cercado por 2+ grânulos (vazio entre grânulos encostados) volta a ser fundo;
        # buraco dentro de um único grânulo (pite, centro escuro, hilo) fica preenchido
        lab_h, nh = ndi.label(holes)
        ring = ndi.grey_dilation(lab_h, size=3) * (lab_h == 0)
        multi = np.zeros(nh + 1, bool)
        for h_id in range(1, nh + 1):
            neigh = np.unique(labels[(ring == h_id) & (labels > 0)])
            if neigh.size >= 2:
                multi[h_id] = True
        if multi.any():
            labels = labels.copy()
            labels[multi[lab_h]] = 0
    labels = _relabel_sizes(labels, min_px, max_px)
    info = dict(threshold=th_desc, bg_level=bg_level, fg_level=fg_level, min_d_px=min_px,
                n=int(labels.max()), dark_objects=dark_objects, split=split)
    return dict(labels=labels, mask=labels > 0, info=info)


def _remove_small(mask, min_area):
    lab, n = ndi.label(mask)
    if n == 0:
        return mask
    sizes = np.bincount(lab.ravel())
    keep = sizes >= min_area
    keep[0] = False
    return keep[lab]


def _hmax_markers(dist, mask, hloc):
    """Marcadores = componentes de máximos regionais cuja proeminência local excede hloc (h-maxima
    com h variável por objeto)."""
    from skimage.morphology import reconstruction
    seed = dist - hloc
    seed = np.where(mask, seed, 0)
    rec = reconstruction(np.minimum(seed, dist), dist, method="dilation")
    peaks = (dist - rec) >= (hloc - 1e-6) * 0.999
    peaks &= mask
    markers, n = ndi.label(peaks, structure=np.ones((3, 3)))
    # objetos sem marcador recebem o máximo global como marcador
    lab0, n0 = ndi.label(mask)
    have = np.unique(lab0[markers > 0])
    miss = np.setdiff1d(np.arange(1, n0 + 1), have)
    if miss.size:
        pos = ndi.maximum_position(dist, lab0, index=miss)
        for k, p in enumerate(pos, start=1):
            markers[p] = n + k
    return markers


def _relabel_sizes(labels, min_px, max_px):
    from skimage.segmentation import relabel_sequential
    if labels.max() == 0:
        return labels
    areas = np.bincount(labels.ravel())
    ecd = np.sqrt(4 * areas / np.pi)
    bad = ecd < min_px
    if max_px:
        bad |= ecd > max_px
    bad[0] = False
    labels = labels.copy()
    labels[bad[labels]] = 0
    return relabel_sequential(labels)[0]


# ----------------------------------------------------------------------------- medidas

def _min_feret(coords):
    """Feret mínimo (calibre mínimo) pelo método das pinças rotativas na envoltória convexa."""
    from scipy.spatial import ConvexHull
    if len(coords) < 4:
        return np.nan, np.nan
    try:
        hull = coords[ConvexHull(coords).vertices]
    except Exception:
        return np.nan, np.nan
    best = (np.inf, 0.0)
    n = len(hull)
    for i in range(n):
        p, q = hull[i], hull[(i + 1) % n]
        e = q - p
        L = np.hypot(*e)
        if L == 0:
            continue
        nrm = np.array([-e[1], e[0]]) / L
        w = (hull - p) @ nrm
        width = w.max() - w.min()
        if width < best[0]:
            best = (width, np.degrees(np.arctan2(e[0], e[1])))
    return best


def measure(im, labels: np.ndarray, exclude_edge: bool = True) -> pd.DataFrame:
    """Tabela de medidas por objeto, em µm (ou px se não calibrado)."""
    from skimage import measure as skm
    import cv2
    px = im.px_um or 1.0
    H, W = labels.shape
    rows = []
    for p in skm.regionprops(labels, intensity_image=im.data):
        minr, minc, maxr, maxc = p.bbox
        edge = minr == 0 or minc == 0 or maxr == H or maxc == W
        A = p.area
        P = p.perimeter_crofton if p.perimeter_crofton > 0 else p.perimeter
        # pontos do contorno (cantos de pixel) para Feret
        img = p.image.astype(np.uint8)
        cnts, _ = cv2.findContours(img, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        cnt = max(cnts, key=cv2.contourArea) if cnts else None
        if cnt is not None and len(cnt) >= 3:
            pts = cnt[:, 0, :].astype(float)
            pts = np.vstack([pts, pts + [1, 0], pts + [0, 1], pts + [1, 1]])
            fmin, _ = _min_feret(pts)
            hull = cv2.convexHull(cnt)
            hp = cv2.arcLength(hull, True)
            approx = cv2.approxPolyDP(hull, 0.04 * hp, True)
            nvert = len(approx)
        else:
            fmin, hp, nvert = np.nan, np.nan, 0
        fmax = p.feret_diameter_max + 1.0 if A > 1 else 1.0
        ecd = np.sqrt(4 * A / np.pi)
        circ = min(1.0, 4 * np.pi * A / P ** 2) if P > 0 else np.nan
        rows.append(dict(
            id=p.label, borda=edge, x_px=p.centroid[1], y_px=p.centroid[0],
            area=A * px ** 2, ecd=ecd * px, feret_max=fmax * px, feret_min=fmin * px if np.isfinite(fmin) else np.nan,
            eixo_maior=p.axis_major_length * px, eixo_menor=p.axis_minor_length * px,
            perimetro=P * px, AR=fmax / fmin if fmin and np.isfinite(fmin) and fmin > 0 else np.nan,
            circularidade=circ, arredondamento=4 * A / (np.pi * fmax ** 2),
            solidez=p.solidity, convexidade=(hp / P) if P > 0 and np.isfinite(hp) else np.nan,
            vertices=nvert, orientacao=np.degrees(p.orientation), brilho=p.intensity_mean))
    df = pd.DataFrame(rows)
    if len(df):
        # fragmentos de watershed de um grânulo cortado pela borda: contato longo com vizinho que toca a borda
        border_ids = set(df.loc[df.borda, "id"])
        if border_ids:
            from skimage.segmentation import find_boundaries
            edge_lab = np.isin(labels, list(border_ids))
            near = ndi.binary_dilation(edge_lab, iterations=2)
            sl = ndi.find_objects(labels)
            for i, r in df.iterrows():
                if r["borda"]:
                    continue
                ys, xs = sl[int(r["id"]) - 1]
                ys = slice(max(0, ys.start - 3), ys.stop + 3)
                xs = slice(max(0, xs.start - 3), xs.stop + 3)
                if not near[ys, xs].any():
                    continue
                m = labels[ys, xs] == r["id"]
                contact = (find_boundaries(m, mode="inner") & near[ys, xs]).sum()
                if contact > 0.35 * r["ecd"] / px:
                    df.at[i, "borda"] = True
        df["classe"] = [classify(r) for r in df.itertuples()]
    df.attrs["unit"] = "µm" if im.px_um else "px"
    df.attrs["exclude_edge"] = exclude_edge
    return df


def classify(r) -> str:
    """Classe de forma (descritiva) para grânulos de amido/partículas."""
    ar, sol, circ, nv = r.AR, r.solidez, r.circularidade, r.vertices
    if not np.isfinite(ar):
        return "indefinido"
    if sol < 0.80:
        return "agregado/sobreposto"
    if sol < 0.90:
        return "irregular/fragmentado"
    if ar >= 2.0:
        return "alongado/lenticular (de perfil)"
    if 3 <= nv <= 7 and circ < 0.92 and ar < 1.5:
        return "poliédrico/angular"
    if ar < 1.2:
        return "esférico/arredondado"
    return "oval/elíptico"


def valid(df: pd.DataFrame) -> pd.DataFrame:
    """Objetos válidos para estatística (exclui os que tocam a borda)."""
    if not len(df):
        return df
    return df[~df["borda"]] if df.attrs.get("exclude_edge", True) else df


# ----------------------------------------------------------------------------- estatística de tamanho

def _wq(x, w, q):
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = np.cumsum(w) - 0.5 * w
    c /= w.sum()
    return float(np.interp(q, c, x))


def size_stats(df: pd.DataFrame, col: str = "ecd") -> dict:
    """Estatística de tamanho: em número e em volume (esfera equivalente, comparável à difração laser)."""
    d = valid(df)[col].dropna().to_numpy(float)
    n = d.size
    out = dict(n=int(n), n_borda=int(df["borda"].sum()) if len(df) else 0, unidade=df.attrs.get("unit", "µm"))
    if n == 0:
        return out
    v = d ** 3
    out.update(media=float(d.mean()), dp=float(d.std(ddof=1)) if n > 1 else 0.0,
               cv_pct=float(100 * d.std(ddof=1) / d.mean()) if n > 1 else 0.0,
               mediana=float(np.median(d)), min=float(d.min()), max=float(d.max()),
               D10=float(np.percentile(d, 10)), D50=float(np.percentile(d, 50)), D90=float(np.percentile(d, 90)),
               Dv10=_wq(d, v, 0.1), Dv50=_wq(d, v, 0.5), Dv90=_wq(d, v, 0.9),
               D43=float((d ** 4).sum() / v.sum()), D32=float(v.sum() / (d ** 2).sum()))
    out["span_n"] = (out["D90"] - out["D10"]) / out["D50"]
    out["span_v"] = (out["Dv90"] - out["Dv10"]) / out["Dv50"]
    if n > 1:
        tcrit = 1.96 if n > 30 else {2: 12.7, 3: 4.30, 5: 2.78, 10: 2.26, 20: 2.09, 30: 2.04}.get(
            min([k for k in (2, 3, 5, 10, 20, 30) if k >= n], default=30), 2.0)
        out["ic95_media"] = float(tcrit * out["dp"] / np.sqrt(n))
    warn = []
    if n < 100:
        warn.append(f"só {n} objetos medidos: use várias imagens/campos (ideal ≥ 300-500 por amostra)")
    elif n < 300:
        warn.append(f"{n} objetos: aceitável para média, pouco para D10/D90 e populações (ideal ≥ 300-500)")
    if out["n_borda"] > 0.3 * (n + out["n_borda"]):
        warn.append("muitos objetos cortados pela borda: aumento alto demais para o tamanho dos grânulos")
    out["avisos"] = warn
    return out


def bimodal(df: pd.DataFrame, col: str = "ecd", min_n: int = 40) -> dict:
    """Testa uma vs duas populações em log(ECD) por mistura gaussiana (BIC). Para trigo, cevada,
    centeio, triticale (grânulos A lenticulares e B esféricos) e misturas."""
    d = valid(df)[col].dropna().to_numpy(float)
    if d.size < min_n:
        return dict(bimodal=False, motivo=f"poucos objetos ({d.size} < {min_n})")
    from sklearn.mixture import GaussianMixture
    x = np.log(d).reshape(-1, 1)
    g1 = GaussianMixture(1, random_state=0).fit(x)
    g2 = GaussianMixture(2, random_state=0, n_init=5).fit(x)
    b1, b2 = g1.bic(x), g2.bic(x)
    mu = g2.means_.ravel()
    sd = np.sqrt(g2.covariances_.ravel())
    o = np.argsort(mu)
    mu, sd, w = mu[o], sd[o], g2.weights_[o]
    sep = abs(mu[1] - mu[0]) / (2 * np.sqrt(sd[0] * sd[1]))   # Ashman D/2
    is_bi = (b1 - b2) > 10 and sep > 1.0 and w.min() > 0.03
    out = dict(bimodal=bool(is_bi), dBIC=float(b1 - b2), separacao_ashman=float(2 * sep))
    if is_bi:
        lab = g2.predict(x)
        lab = np.where(lab == o[0], 0, 1)
        grid = np.linspace(mu[0], mu[1], 400)
        dens = [w[k] * np.exp(-0.5 * ((grid - mu[k]) / sd[k]) ** 2) / sd[k] for k in (0, 1)]
        cut = float(np.exp(grid[np.argmin(np.abs(dens[0] - dens[1]))]))
        pops = []
        for k, name in ((0, "pequenos (tipo B)"), (1, "grandes (tipo A)")):
            dk = d[lab == k]
            pops.append(dict(populacao=name, n=int(dk.size), frac_numero=float(dk.size / d.size),
                             frac_volume=float((dk ** 3).sum() / (d ** 3).sum()), media=float(dk.mean()),
                             dp=float(dk.std(ddof=1)) if dk.size > 1 else 0.0, mediana=float(np.median(dk))))
        out.update(corte=cut, populacoes=pops, media_log=[float(np.exp(m)) for m in mu])
    return out


def shape_summary(df: pd.DataFrame) -> dict:
    v = valid(df)
    if not len(v):
        return {}
    cls = v["classe"].value_counts(normalize=True).mul(100).round(1).to_dict()
    agg = {c: (float(v[c].mean()), float(v[c].std(ddof=1)) if len(v) > 1 else 0.0)
           for c in ("AR", "circularidade", "arredondamento", "solidez")}
    return dict(classes_pct=cls, medias=agg)


# ----------------------------------------------------------------------------- pites / poros superficiais

def pits(im, labels: np.ndarray, df: pd.DataFrame | None = None, pit_max_um: float = 1.5, k: float = 5.0,
         min_d_px: float = 5.0, max_rel: float = 0.35, min_depth: float = 0.35) -> dict:
    """Detecta pites/poros escuros na superfície dos grânulos (amido hidrolisado por amilases, poroso,
    danificado). Top-hat preto com elemento do tamanho máximo do pite, dentro dos grânulos (erodidos
    para evitar a borda). Acrescenta colunas n_pites e area_pites_pct ao df. Retorna resumo e máscara."""
    from skimage import morphology
    r = max(2, int(round(um2px(im, pit_max_um, 8) / 2)))
    a = ndi.gaussian_filter(im.data.astype(np.float32), 0.8)
    inner = ndi.binary_erosion(labels > 0, iterations=max(2, r // 2))
    th = gray_close(a, r) - a
    vals = th[inner]
    if vals.size == 0:
        return dict(n_pites=0)
    med = np.median(vals)
    mad = np.median(np.abs(vals - med)) * 1.4826 + 1e-6
    pm = (th > med + k * mad) & inner
    pm = _remove_small(pm, int(np.pi / 4 * min_d_px ** 2))
    host_ecd = np.sqrt(4 * np.bincount(labels.ravel(), minlength=labels.max() + 1) / np.pi)
    # só pites compactos (descarta sulcos/arestas de facetas e fissuras alongadas)
    from skimage import measure as skm
    plab, n = ndi.label(pm)
    keep = np.zeros(n + 1, bool)
    dmax = 2 * r + 2
    for p in skm.regionprops(plab):
        el = p.axis_major_length / max(p.axis_minor_length, 1.0)
        host = labels[int(p.centroid[0]), int(p.centroid[1])]
        d_eq = np.sqrt(4 * p.area / np.pi)
        if (el <= 2.5 and p.solidity >= 0.75 and p.axis_major_length <= 1.5 * dmax
                and host > 0 and d_eq <= max_rel * host_ecd[host]):
            keep[p.label] = True
    pm = keep[plab]
    plab, n = ndi.label(pm)
    # contraste de borda: pite real é um degrau nítido (buraco escuro), não um sombreamento suave
    if n:
        ring = ndi.binary_dilation(pm, iterations=2) & ~pm & (labels > 0)
        rlab = np.where(ring, ndi.grey_dilation(plab, size=5), 0)
        idx = np.arange(1, n + 1)
        m_in = np.asarray(ndi.minimum(a, plab, index=idx))
        m_ring = np.asarray(ndi.mean(a, rlab, index=idx))
        span = max(1e-3, float(np.percentile(a[labels > 0], 90) - np.percentile(a[labels == 0], 50))) \
            if (labels == 0).any() else 0.5
        depth = (m_ring - m_in) / span
        keep = np.zeros(n + 1, bool)
        keep[1:] = depth >= min_depth
        pm = keep[plab]
        plab, n = ndi.label(pm)
    per = np.zeros(labels.max() + 1)
    area = np.zeros(labels.max() + 1)
    if n:
        owner = ndi.maximum(labels, plab, index=np.arange(1, n + 1)).astype(int)
        sizes = ndi.sum(pm, plab, index=np.arange(1, n + 1))
        np.add.at(per, owner, 1)
        np.add.at(area, owner, sizes)
    garea = np.bincount(labels.ravel(), minlength=labels.max() + 1).astype(float)
    pct = np.where(garea > 0, 100 * area / np.maximum(garea, 1), 0)
    if df is not None and len(df):
        df["n_pites"] = per[df["id"].to_numpy()].astype(int)
        df["area_pites_pct"] = pct[df["id"].to_numpy()]
    v = valid(df) if df is not None else None
    frac_with = float((v["n_pites"] > 0).mean()) if v is not None and len(v) else np.nan
    pit_d = 2 * np.sqrt(ndi.sum(pm, plab, index=np.arange(1, n + 1)) / np.pi) * (im.px_um or 1) if n else np.array([])
    return dict(n_pites=int(n), frac_graos_com_pites=frac_with,
                pites_por_grao=float(v["n_pites"].mean()) if v is not None and len(v) else np.nan,
                area_pites_pct=float(v["area_pites_pct"].mean()) if v is not None and len(v) else np.nan,
                d_pite_medio=float(pit_d.mean()) if n else np.nan, mask=pm)


# ----------------------------------------------------------------------------- comparação entre amostras

def compare(samples: dict, col: str = "ecd") -> dict:
    """Compara distribuições entre amostras {nome: df}. Kruskal-Wallis + Mann-Whitney par a par com
    correção de Holm (não paramétrico: distribuições de tamanho são assimétricas)."""
    from scipy import stats
    names = list(samples)
    data = [valid(samples[k])[col].dropna().to_numpy() for k in names]
    res = dict(variavel=col, resumo=[dict(amostra=k, n=len(x), media=float(np.mean(x)), dp=float(np.std(x, ddof=1)),
                                         mediana=float(np.median(x))) for k, x in zip(names, data) if len(x) > 1])
    if len(names) < 2:
        return res
    if len(names) > 2:
        H, p = stats.kruskal(*data)
        res["kruskal"] = dict(H=float(H), p=float(p))
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            U, p = stats.mannwhitneyu(data[i], data[j], alternative="two-sided")
            r_rb = 1 - 2 * U / (len(data[i]) * len(data[j]))        # correlação rank-biserial
            pairs.append(dict(a=names[i], b=names[j], U=float(U), p=float(p), efeito_rb=float(r_rb)))
    ps = np.array([q["p"] for q in pairs])
    o = np.argsort(ps)
    m = len(ps)
    adj = np.empty(m)
    run = 0.0
    for rank, k in enumerate(o):
        run = max(run, (m - rank) * ps[k])
        adj[k] = min(1.0, run)
    for q, pa in zip(pairs, adj):
        q["p_holm"] = float(pa)
    res["pares"] = pairs
    return res
