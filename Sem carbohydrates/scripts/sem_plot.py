"""
sem_plot.py - figuras de MEV: imagem com barra de escala, sobreposições de segmentação,
distribuições de tamanho/forma, defeitos de filme, espessura, painéis para artigo e comparações

    from sem_plot import show, plot_granules, plot_film, plot_section, panel, plot_compare
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from scipy import ndimage as ndi

OKABE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
CLASS_COLORS = {"esférico/arredondado": "#56B4E9", "oval/elíptico": "#009E73", "poliédrico/angular": "#E69F00",
                "alongado/lenticular (de perfil)": "#CC79A7", "irregular/fragmentado": "#D55E00",
                "agregado/sobreposto": "#F0E442", "indefinido": "#999999"}
plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300, "figure.dpi": 100})
T = {"pt": dict(n="n", ecd="Diâmetro equivalente (µm)", ecd_px="Diâmetro equivalente (px)", freq="Frequência (%)",
                cum="Acumulado (%)", numero="número", volume="volume", ar="Razão de aspecto (Feret máx/mín)",
                circ="Circularidade (4πA/P²)", thick="Espessura (µm)", pos="Posição ao longo do filme (µm)",
                poros="poros", part="partículas/protuberâncias", trincas="trincas", borda="tocando a borda (excluído)"),
     "en": dict(n="n", ecd="Equivalent diameter (µm)", ecd_px="Equivalent diameter (px)", freq="Frequency (%)",
                cum="Cumulative (%)", numero="number", volume="volume", ar="Aspect ratio (Feret max/min)",
                circ="Circularity (4πA/P²)", thick="Thickness (µm)", pos="Position along film (µm)",
                poros="pores", part="particles/protrusions", trincas="cracks", borda="touching border (excluded)")}


def nice_length(target_um: float) -> float:
    """Comprimento 'redondo' (1, 2, 5 × 10ⁿ) ≤ alvo."""
    if target_um <= 0:
        return 1.0
    e = np.floor(np.log10(target_um))
    for m in (5, 2, 1):
        v = m * 10 ** e
        if v <= target_um:
            return float(v)
    return float(10 ** e)


def _fmt_um(v):
    if v >= 1000:
        return f"{v / 1000:g} mm"
    if v >= 1:
        return f"{v:g} µm"
    return f"{v * 1000:g} nm"


def scalebar(ax, px_um, shape, length_um=None, loc="lower right", color="white", frac=0.2, fontsize=9,
             box=True, lw=None):
    """Desenha barra de escala em µm nos eixos de uma imagem (coordenadas em px)."""
    if not px_um:
        return
    H, W = shape[:2]
    L = length_um or nice_length(frac * W * px_um)
    Lpx = L / px_um
    m = 0.04 * min(H, W)
    x1 = W - m if "right" in loc else m + Lpx
    x0 = x1 - Lpx
    y = H - m if "lower" in loc else m + 0.05 * H
    th = lw or max(2.5, 0.008 * H * 72 / 100)
    other = "black" if color == "white" else "white"
    if box:
        pad = 0.012 * min(H, W)
        ax.add_patch(plt.Rectangle((x0 - pad, y - 0.07 * H - pad), Lpx + 2 * pad, 0.07 * H + 2 * pad,
                                   facecolor=other, alpha=0.55, edgecolor="none", zorder=4))
    ax.plot([x0, x1], [y, y], color=color, lw=th, solid_capstyle="butt", zorder=5)
    ax.text((x0 + x1) / 2, y - 0.018 * H, _fmt_um(L), color=color, ha="center", va="bottom",
            fontsize=fontsize, zorder=6, fontweight="bold")


def show(im, ax=None, title=None, bar=True, bar_um=None, enhance=False, letter=None, lang="pt", cmap="gray"):
    from sem_process import autolevel, enhance as clahe
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.4, 6.4 * im.data.shape[0] / im.data.shape[1]))
    a = clahe(im.data) if enhance else autolevel(im.data, 0.1, 99.9)
    ax.imshow(a, cmap=cmap, vmin=0, vmax=1, interpolation="nearest" if max(a.shape) < 800 else "antialiased")
    ax.set_axis_off()
    if bar and im.px_um:
        scalebar(ax, im.px_um, a.shape, bar_um)
    if letter:
        ax.text(0.025, 0.975, letter, transform=ax.transAxes, ha="left", va="top", fontsize=13, fontweight="bold",
                color="black", bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.9))
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def _boundaries(labels):
    from skimage.segmentation import find_boundaries
    return find_boundaries(labels, mode="inner")


def overlay_granules(ax, im, labels, df, color_by="classe", numbers=False, lang="pt"):
    show(im, ax=ax, bar=True, lang=lang)
    b = _boundaries(labels)
    H, W = labels.shape
    rgba = np.zeros((H, W, 4))
    lut = {}
    if len(df):
        for r in df.itertuples():
            if r.borda:
                c = (0.6, 0.6, 0.6)
            elif color_by == "pop" and hasattr(r, "pop"):
                c = matplotlib.colors.to_rgb(OKABE[0] if r.pop == 0 else OKABE[1])
            else:
                c = matplotlib.colors.to_rgb(CLASS_COLORS.get(r.classe, "#56B4E9"))
            lut[r.id] = c
    lab_b = np.where(b, labels, 0)
    lw = max(2, int(round(max(H, W) / 450)))
    lab_b = ndi.grey_dilation(lab_b, size=lw)
    ids = np.unique(lab_b)
    for i in ids[ids > 0]:
        m = lab_b == i
        rgba[m, :3] = lut.get(i, (0.34, 0.71, 0.91))
        rgba[m, 3] = 1.0
    ax.imshow(rgba, interpolation="nearest")
    if numbers and len(df):
        for r in df.itertuples():
            if not r.borda:
                ax.text(r.x_px, r.y_px, str(r.id), color="yellow", fontsize=5, ha="center", va="center",
                        path_effects=[pe.withStroke(linewidth=1.2, foreground="black")])
    if len(df):
        present = [c for c in CLASS_COLORS if c in set(df.loc[~df.borda, "classe"])] if color_by == "classe" else []
        handles = [Line2D([], [], color=CLASS_COLORS[c], lw=2, label=c) for c in present]
        if df.borda.any():
            handles.append(Line2D([], [], color="0.6", lw=2, label=T[lang]["borda"]))
        if handles:
            ax.legend(handles=handles, loc="upper right", fontsize=6, framealpha=0.85, handlelength=1.2)


def _plain_log(ax):
    from matplotlib.ticker import FuncFormatter, LogLocator
    ax.xaxis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(1, 10)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(FuncFormatter(lambda v, _: ""))


def plot_granules(im, labels, df, stats=None, bm=None, pits_mask=None, lang="pt", out=None, title=None,
                  numbers=False):
    """Figura 2×2: (a) segmentação, (b) histograma em número (+ populações), (c) acumuladas em número e
    volume com D10/D50/D90, (d) forma (AR × circularidade)."""
    from sem_granules import valid
    t = T[lang]
    v = valid(df)
    unit = df.attrs.get("unit", "µm")
    xl = t["ecd"] if unit == "µm" else t["ecd_px"]
    fig = plt.figure(figsize=(10, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.28, wspace=0.22)
    ax0 = fig.add_subplot(gs[0, 0])
    overlay_granules(ax0, im, labels, df, numbers=numbers, lang=lang)
    if pits_mask is not None and pits_mask.any():
        rg = np.zeros(pits_mask.shape + (4,))
        rg[pits_mask] = (1, 0.1, 0.1, 0.9)
        ax0.imshow(rg, interpolation="nearest")
    ax0.set_title(title or im.title, fontsize=9)
    ax0.text(0.01, 0.01, f"n = {len(v)} (+{int(df.borda.sum()) if len(df) else 0} na borda)", transform=ax0.transAxes,
             fontsize=7, color="white", va="bottom", path_effects=[pe.withStroke(linewidth=1.5, foreground="black")])
    d = v["ecd"].to_numpy() if len(v) else np.array([])
    ax1 = fig.add_subplot(gs[0, 1])
    if d.size:
        logb = d.max() / max(d.min(), 1e-9) > 6
        bins = np.geomspace(d.min() * 0.95, d.max() * 1.05, 25) if logb else np.linspace(0, d.max() * 1.05, 25)
        w = np.ones_like(d) * 100 / d.size
        ax1.hist(d, bins=bins, weights=w, color="#9ecae1", edgecolor="#08519c", lw=0.6)
        if logb:
            ax1.set_xscale("log")
            _plain_log(ax1)
        if bm and bm.get("bimodal"):
            ax1.axvline(bm["corte"], color="k", ls="--", lw=0.8)
            for p, c in zip(bm["populacoes"], (OKABE[0], OKABE[5])):
                ax1.text(p["media"], ax1.get_ylim()[1] * 0.6, f"{p['populacao']}\n{p['media']:.1f} ({100 * p['frac_numero']:.0f}% n)",
                         color=c, ha="center", va="center", fontsize=7, fontweight="bold",
                         bbox=dict(fc="white", ec="none", alpha=0.8))
        if stats:
            ax1.text(0.98, 0.97, f"n = {stats['n']}\nmédia = {stats['media']:.2f} ± {stats['dp']:.2f}\n"
                                 f"D50 = {stats['D50']:.2f}\nCV = {stats['cv_pct']:.0f}%",
                     transform=ax1.transAxes, ha="right", va="top", fontsize=7,
                     bbox=dict(fc="white", ec="0.8", alpha=0.9))
    ax1.set_xlabel(xl)
    ax1.set_ylabel(t["freq"])
    ax2 = fig.add_subplot(gs[1, 0])
    if d.size:
        o = np.sort(d)
        cn = (np.arange(1, o.size + 1) - 0.5) / o.size * 100
        w3 = o ** 3
        cv = (np.cumsum(w3) - 0.5 * w3) / np.sum(w3) * 100
        ax2.plot(o, cn, color=OKABE[0], lw=1.6, label=t["numero"])
        ax2.plot(o, cv, color=OKABE[5], lw=1.6, label=t["volume"])
        if stats:
            for key, c, yy in (("Dv10", OKABE[5], 10), ("Dv50", OKABE[5], 50), ("Dv90", OKABE[5], 90)):
                ax2.plot(stats[key], yy, "o", ms=3.5, color=c)
            for q_, yy in ((10, 10), (50, 50), (90, 90)):
                ax2.plot(np.interp(yy, cn, o), yy, "o", ms=3.5, color=OKABE[0])
            ax2.text(0.98, 0.03, f"D10/50/90 (n) = {stats['D10']:.1f} / {stats['D50']:.1f} / {stats['D90']:.1f}\n"
                                 f"Dv10/50/90 = {stats['Dv10']:.1f} / {stats['Dv50']:.1f} / {stats['Dv90']:.1f}\n"
                                 f"D[4,3] = {stats['D43']:.1f}   D[3,2] = {stats['D32']:.1f}",
                     transform=ax2.transAxes, ha="right", va="bottom", fontsize=7,
                     bbox=dict(fc="white", ec="0.8", alpha=0.9))
        if d.max() / max(d.min(), 1e-9) > 6:
            ax2.set_xscale("log")
            _plain_log(ax2)
        ax2.legend(loc="upper left", fontsize=7, frameon=False)
    ax2.set_xlabel(xl)
    ax2.set_ylabel(t["cum"])
    ax2.set_ylim(0, 102)
    ax3 = fig.add_subplot(gs[1, 1])
    if len(v):
        for c, g in v.groupby("classe"):
            ax3.scatter(g["circularidade"], g["AR"], s=10, color=CLASS_COLORS.get(c, "0.5"), label=f"{c} ({len(g)})",
                        alpha=0.8, edgecolor="none")
        ax3.legend(fontsize=6, frameon=False, loc="upper left")
        ax3.set_ylim(0.95, max(2.2, float(v["AR"].quantile(0.99)) * 1.1))
    ax3.set_xlabel(t["circ"])
    ax3.set_ylabel(t["ar"])
    for a_, L in zip((ax0, ax1, ax2, ax3), "abcd"):
        a_.text(-0.02 if a_ is ax0 else -0.14, 1.03, f"({L})", transform=a_.transAxes, fontsize=10, fontweight="bold")
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_film(im, res, lang="pt", out=None, title=None):
    t = T[lang]
    fig, ax = plt.subplots(figsize=(7.5, 7.5 * im.data.shape[0] / im.data.shape[1] + 0.6))
    show(im, ax=ax, bar=True)
    cols = {"poros": "#00E5FF", "particulas": "#FFD400", "trincas": "#FF2D2D"}
    handles = []
    for k, lab in (("poros", t["poros"]), ("particulas", t["part"]), ("trincas", t["trincas"])):
        m = res["mascaras"][k]
        if m.any():
            mm = ndi.binary_dilation(m, iterations=1) if k == "trincas" else m
            ax.contour(mm, [0.5], colors=cols[k], linewidths=0.8)
        n = res["resumo"].get(f"n_{k}", 0)
        handles.append(Line2D([], [], color=cols[k], lw=2, label=f"{lab} ({n})"))
    ax.legend(handles=handles, loc="upper right", fontsize=7, framealpha=0.85)
    ax.set_title(title or im.title, fontsize=9)
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_section(im, cs, lang="pt", out=None, title=None):
    t = T[lang]
    a = im.data.T if cs["rotated"] else im.data
    from types import SimpleNamespace
    im2 = SimpleNamespace(data=a, px_um=im.px_um, title=im.title)
    H, W = a.shape
    fig = plt.figure(figsize=(7.5, 7.5 * H / W + 2.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[H / W * 3, 1], hspace=0.12)
    ax = fig.add_subplot(gs[0])
    show(im2, ax=ax, bar=True)
    x = np.arange(W)
    ok = cs["ok"]
    for arr in (cs["top"], cs["bot"]):
        ax.plot(np.where(ok, x, np.nan), np.where(ok, arr, np.nan), color="#FFD400", lw=1.0)
    if cs["pores"].any():
        ax.contour(cs["pores"], [0.5], colors="#00E5FF", linewidths=0.6)
    for xi in np.linspace(0.1 * W, 0.9 * W, 5).astype(int):
        if ok[xi]:
            ax.annotate("", xy=(xi, cs["bot"][xi]), xytext=(xi, cs["top"][xi]),
                        arrowprops=dict(arrowstyle="<->", color="#FFD400", lw=0.9))
            ax.text(xi + 0.01 * W, 0.5 * (cs["top"][xi] + cs["bot"][xi]), f"{cs['perfil'][xi]:.1f}",
                    color="#FFD400", fontsize=7, va="center",
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="black")])
    ax.set_title(title or im.title, fontsize=9)
    ax2 = fig.add_subplot(gs[1])
    xs = x * (im.px_um or 1)
    ax2.plot(xs[ok], cs["perfil"][ok], color=OKABE[0], lw=1)
    m, s = cs["espessura_media"], cs["espessura_dp"]
    ax2.axhline(m, color="k", lw=0.8)
    ax2.axhspan(m - s, m + s, color=OKABE[0], alpha=0.12)
    ax2.set_xlabel(t["pos"] if im.px_um else "px")
    ax2.set_ylabel(t["thick"] if im.px_um else "px")
    ax2.set_xlim(xs[0], xs[-1])
    ax2.text(0.99, 0.95, f"{m:.1f} ± {s:.1f} {cs['unidade']} (n = {cs['n_colunas']} colunas)",
             transform=ax2.transAxes, ha="right", va="top", fontsize=8)
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def panel(images, letters=True, titles=None, ncols=2, bar_um=None, out=None, width_in=7.0, same_bar=True,
          enhance=False):
    """Painel de figuras para artigo: imagens lado a lado com barra de escala e letras (a), (b)..."""
    n = len(images)
    ncols = min(ncols, n)
    nrows = int(np.ceil(n / ncols))
    asp = max(im.data.shape[0] / im.data.shape[1] for im in images)
    fig, axes = plt.subplots(nrows, ncols, figsize=(width_in, width_in / ncols * asp * nrows + 0.3 * nrows),
                             squeeze=False)
    if bar_um is None and same_bar:
        cal = [im for im in images if im.px_um]
        if cal and len({round(im.px_um, 6) for im in cal}) == 1:
            bar_um = nice_length(0.2 * cal[0].data.shape[1] * cal[0].px_um)
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.set_axis_off()
            continue
        im = images[k]
        show(im, ax=ax, bar=True, bar_um=bar_um, enhance=enhance,
             letter=f"({chr(97 + k)})" if letters else None, title=titles[k] if titles and k < len(titles) else None)
    fig.subplots_adjust(wspace=0.03, hspace=0.08 if not titles else 0.18, left=0.01, right=0.99, top=0.97, bottom=0.01)
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_compare(samples: dict, col="ecd", res=None, lang="pt", out=None, ylabel=None):
    from sem_granules import valid
    names = list(samples)
    data = [valid(samples[k])[col].dropna().to_numpy() for k in names]
    fig, ax = plt.subplots(figsize=(max(3.5, 1.1 * len(names) + 1.5), 3.6))
    parts = ax.violinplot(data, showextrema=False, widths=0.8)
    for b, c in zip(parts["bodies"], OKABE * 3):
        b.set_facecolor(c)
        b.set_alpha(0.35)
    ax.boxplot(data, widths=0.18, showfliers=False, medianprops=dict(color="k"))
    rng = np.random.default_rng(0)
    for i, (x, c) in enumerate(zip(data, OKABE * 3), start=1):
        ax.scatter(i + rng.uniform(-0.07, 0.07, x.size), x, s=4, color=c, alpha=0.5, edgecolor="none")
        ax.text(i, ax.get_ylim()[0], f"n={x.size}", ha="center", va="bottom", fontsize=7)
    ax.set_xticks(range(1, len(names) + 1), names, rotation=20 if len(names) > 3 else 0)
    unit = samples[names[0]].attrs.get("unit", "µm")
    lab = {"ecd": T[lang]["ecd"], "AR": T[lang]["ar"], "circularidade": T[lang]["circ"]}.get(col, col)
    ax.set_ylabel(ylabel or (lab if unit == "µm" or col != "ecd" else T[lang]["ecd_px"]))
    if res and res.get("kruskal"):
        ax.set_title(f"Kruskal-Wallis p = {res['kruskal']['p']:.2g}", fontsize=8)
    elif res and res.get("pares"):
        ax.set_title(f"Mann-Whitney p = {res['pares'][0]['p']:.2g}", fontsize=8)
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig
