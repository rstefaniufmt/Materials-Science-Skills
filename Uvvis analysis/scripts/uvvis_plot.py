"""
uvvis_plot.py - figuras de UV-Vis: espectros, séries, Tauc, calibração e cinética

    plot_spectra(specs, mode="A", peaks=[tab], series_values=[0, 5, 10], series_label="t (min)",
                 ev_axis=True, lang="pt", out="fig.png")
    plot_tauc(res, out="tauc.png")            # res = uvvis_bandgap.tauc(...)
    plot_calibration(cal, out="curva.png")    # cal = uvvis_quant.calibrate(...)
    plot_kinetics([kin1, kin2], names=[...], out="cinetica.png")

Convenções: eixo x em nm CRESCENTE (padrão de UV-Vis; diferente do IR), eixo superior opcional
em eV, séries (tempo/concentração) com gradiente de cor e barra de cores, rótulos de λmax
verticais sem sobreposição.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402
from matplotlib import cm, colors as mcolors  # noqa: E402

from uvvis_io import HC_EV_NM, Spectrum  # noqa: E402

TXT = {
    "en": {"x": "Wavelength (nm)", "ev": "Energy (eV)", "A": "Absorbance", "%T": "Transmittance (%)",
           "T": "Transmittance", "%R": "Reflectance (%)", "R": "Reflectance", "KM": "F(R)",
           "eps": "ε (L mol$^{-1}$ cm$^{-1}$)", "logeps": "log ε", "alpha": "α (cm$^{-1}$)",
           "counts": "Counts", "arb": "Intensity", "d1": "dA/dλ", "d2": "d²A/dλ²", "norm": "normalized",
           "offset": "offset", "au": "a.u.", "conc": "Concentration", "time": "Time",
           "hv": "Photon energy, hν (eV)", "resid": "Residual"},
    "pt": {"x": "Comprimento de onda (nm)", "ev": "Energia (eV)", "A": "Absorbância",
           "%T": "Transmitância (%)", "T": "Transmitância", "%R": "Reflectância (%)", "R": "Reflectância",
           "KM": "F(R)", "eps": "ε (L mol$^{-1}$ cm$^{-1}$)", "logeps": "log ε", "alpha": "α (cm$^{-1}$)",
           "counts": "Contagens", "arb": "Intensidade", "d1": "dA/dλ", "d2": "d²A/dλ²", "norm": "normalizada",
           "offset": "deslocado", "au": "u.a.", "conc": "Concentração", "time": "Tempo",
           "hv": "Energia do fóton, hν (eV)", "resid": "Resíduo"},
}
PALETTE = ["#1f1f1f", "#c0392b", "#1f5fa8", "#2e8b57", "#8e44ad", "#d35400", "#16a085", "#7f6000"]


def _style():
    plt.rcParams.update({
        "font.size": 10, "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": False, "ytick.right": True, "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.minor.visible": True, "ytick.minor.visible": True, "xtick.minor.size": 2, "ytick.minor.size": 2,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.dpi": 300,
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def _prep(s: Spectrum, mode: str) -> Spectrum:
    if mode == "A" and s.ykind in ("T", "%T"):
        return s.to_absorbance()
    if mode == "A" and s.ykind in ("R", "%R"):
        return s.to_absorbance()       # log(1/R)
    if mode == "%T" and s.ykind in ("A", "T"):
        return s.to_transmittance(True)
    if mode == "KM" and s.ykind in ("R", "%R"):
        return s.to_km()
    if mode == "%R" and s.ykind == "R":
        return s.copy(y=s.y * 100, ykind="%R")
    return s


def _spread(xs, width, lo, hi, iters=300):
    order = np.argsort(xs)
    p = np.asarray(xs, float)[order].copy()
    for _ in range(iters):
        moved = False
        for i in range(len(p) - 1):
            gap = p[i + 1] - p[i]
            if gap < width:
                d = (width - gap) / 2 + 1e-9
                p[i] -= d
                p[i + 1] += d
                moved = True
        p = np.clip(p, lo + width / 2, hi - width / 2)
        if not moved:
            break
    out = np.empty_like(p)
    out[order] = p
    return out


def add_ev_axis(ax, lang="en"):
    """Eixo superior em eV (E = 1239.84/λ)."""
    f = lambda x: HC_EV_NM / np.clip(x, 1e-6, None)
    sec = ax.secondary_xaxis("top", functions=(f, f))
    sec.set_xlabel(TXT[lang]["ev"])
    sec.tick_params(direction="in")
    return sec


def plot_spectra(spectra, mode: str = "auto", offset=None, peaks=None, label_fmt: str = "{wavelength}",
                 label_extra: str | None = None, xlim=None, ylim=None, labels=None, legend: str = "auto",
                 series_values=None, series_label: str | None = None, cmap: str = "viridis",
                 highlight=None, ev_axis: bool = False, lang: str = "en", colors=None, lw: float = 1.2,
                 fontsize: float = 7.5, figsize=(6.4, 4.2), title=None, normalized: bool = False,
                 vlines=None, out: str | None = None, ax=None):
    """mode: 'A' | 'T' (=%T) | 'R' (=%R) | 'KM' | 'auto' (tipo do 1º espectro) | 'raw'
    peaks: DataFrame (find_bands) ou lista - rotula λmax acima dos máximos (A/KM) ou abaixo (%T).
    series_values: um número por espectro (tempo, concentração, pH...) -> cores em gradiente
    contínuo (cmap) + barra de cores com `series_label`, no lugar da legenda.
    offset: None | 'auto' | número.  highlight: [(λ1, λ2, 'texto'), ...].
    vlines: [(λ, 'texto'), ...] linhas verticais tracejadas (ex.: pontos isosbésticos)."""
    _style()
    if isinstance(spectra, Spectrum):
        spectra = [spectra]
    n = len(spectra)
    if peaks is not None and not isinstance(peaks, (list, tuple)):
        peaks = [peaks]
    peaks = list(peaks) if peaks is not None else [None] * n
    peaks += [None] * (n - len(peaks))
    if mode == "auto":
        k0 = spectra[0].ykind
        mode = {"T": "%T", "%T": "%T", "R": "%R", "%R": "%R"}.get(k0, "A" if k0 in ("A",) else "raw")
    mode = {"T": "%T", "R": "%R"}.get(mode, mode)
    data = [_prep(s, mode) if mode != "raw" else s for s in spectra]
    T = TXT.get(lang, TXT["en"])
    labels = labels or [s.title for s in data]
    if xlim is None:
        xlim = (min(s.xrange[0] for s in data), max(s.xrange[1] for s in data))
    xlo, xhi = min(xlim), max(xlim)
    if series_values is not None:
        norm = mcolors.Normalize(vmin=min(series_values), vmax=max(series_values))
        cmap_ = matplotlib.colormaps[cmap]
        cols = [cmap_(0.05 + 0.85 * norm(v)) for v in series_values]
    else:
        cols = [(colors or PALETTE)[i % len(colors or PALETTE)] for i in range(n)]

    def vis(s):
        o = np.argsort(s.x)
        x, y = s.x[o], s.y[o]
        m = (x >= xlo) & (x <= xhi)
        return x[m], y[m]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    step = 0.0
    if offset is not None and n > 1:
        step = (0.7 * max(np.ptp(vis(s)[1]) for s in data)) if offset == "auto" else float(offset)
    offs = np.array([(n - 1 - i) * step for i in range(n)]) if step else np.zeros(n)
    ax.set_xlim(xlo, xhi)
    if highlight:
        tr = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for h in highlight:
            a_, b_ = min(h[0], h[1]), max(h[0], h[1])
            ax.axvspan(a_, b_, color="#f1c40f", alpha=0.18, lw=0, zorder=0)
            if len(h) > 2 and h[2]:
                ax.text((a_ + b_) / 2, 0.985, h[2], transform=tr, ha="center", va="top", fontsize=fontsize,
                        color="#7a5c00")
    for i, s in enumerate(data):
        x, y = vis(s)
        ax.plot(x, y + offs[i], color=cols[i], lw=lw, label=labels[i], zorder=3)
    if vlines:
        tr2 = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for xv, lab in vlines:
            ax.axvline(xv, color="0.35", lw=0.8, ls="--", zorder=2)
            if lab:
                ax.text(xv, 0.985, f" {lab}", transform=tr2, ha="left", va="top", fontsize=fontsize, color="0.25")
    ally = np.concatenate([vis(s)[1] + offs[i] for i, s in enumerate(data)])
    ymin, ymax = np.nanmin(ally), np.nanmax(ally)
    yr = (ymax - ymin) or 1.0
    has_lab = any(p is not None and len(p) for p in peaks)
    ax.set_ylim(ymin - 0.05 * yr * (3 if mode == "%T" and has_lab else 1),
                ymax + 0.05 * yr * (3 if mode != "%T" and has_lab else 1))
    # ---- rótulos de λmax
    texts = []
    if has_lab:
        fig.canvas.draw()
        ppx = abs(ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0])
        lw_nm = (fontsize * 1.35 / 72 * fig.dpi) / ppx
        pad = 0.02 * yr
        for i, (s, p) in enumerate(zip(data, peaks)):
            if p is None or not len(p):
                continue
            x, y = vis(s)
            y = y + offs[i]
            p = p[(p.wavelength >= xlo) & (p.wavelength <= xhi)]
            xa = p.wavelength.to_numpy(float)
            xl = _spread(xa, lw_nm, xlo, xhi)
            for (_, row), x0, x1 in zip(p.iterrows(), xa, xl):
                y0 = float(np.interp(x0, x, y))
                win = (x >= min(x0, x1) - lw_nm / 2) & (x <= max(x0, x1) + lw_nm / 2)
                txt = label_fmt.format(**row.to_dict())
                if row.get("kind", "") == "shoulder":
                    txt += " sh"
                if label_extra and label_extra in row and str(row[label_extra]) not in ("", "nan"):
                    txt += "  " + str(row[label_extra])
                if mode == "%T":
                    yl, va, sg = (np.nanmin(y[win]) if win.any() else y0) - pad, "top", -1
                else:
                    yl, va, sg = (np.nanmax(y[win]) if win.any() else y0) + pad, "bottom", 1
                if abs(x1 - x0) > 0.25 * lw_nm or abs(yl - y0) > 2 * pad:
                    ax.plot([x0, x0, x1], [y0 + sg * 0.3 * pad, y0 + sg * 0.3 * pad + (yl - y0) / 2, yl],
                            lw=0.5, color="0.45", zorder=4)
                texts.append(ax.text(x1, yl + sg * 0.15 * pad, txt, rotation=90, ha="center", va=va,
                                     fontsize=fontsize, color=cols[i] if n > 1 and series_values is None else "0.1",
                                     zorder=6))
        for _ in range(6):
            fig.canvas.draw()
            inv = ax.transData.inverted()
            ys = []
            for t in texts:
                b = t.get_window_extent()
                ys += [inv.transform((b.x0, b.y0))[1], inv.transform((b.x0, b.y1))[1]]
            lo_, hi_ = ax.get_ylim()
            cur = hi_ - lo_
            nl, nh = min(lo_, min(ys) - 0.02 * cur), max(hi_, max(ys) + 0.02 * cur)
            if nl >= lo_ and nh <= hi_:
                break
            ax.set_ylim(nl, nh)
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_xlabel(T["x"])
    ylab = T.get(mode if mode != "raw" else data[0].ykind, T["arb"])
    if mode == "A" and spectra[0].ykind in ("R", "%R"):
        ylab = "log(1/R)"
    if normalized:
        ylab += f" ({T['norm']})"
    if step:
        ylab = f"{ylab.split(' (')[0]} ({T['au']}, {T['offset']})"
        ax.set_yticklabels([])
        ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_ylabel(ylab)
    if ev_axis:
        add_ev_axis(ax, lang)
    if title:
        ax.set_title(title, fontsize=10, pad=24 if ev_axis else 6)
    if series_values is not None:
        sm = cm.ScalarMappable(norm=norm, cmap=matplotlib.colormaps[cmap])
        cb = fig.colorbar(sm, ax=ax, pad=0.02, fraction=0.05)
        cb.set_label(series_label or "")
        cb.ax.tick_params(direction="in")
    elif legend == "auto" and n > 1:
        if step:
            for i, s in enumerate(data):
                x, y = vis(s)
                ax.text(xhi - 0.01 * (xhi - xlo), y[-max(1, len(y) // 15):].max() + offs[i] + 0.01 * yr,
                        labels[i], ha="right", va="bottom", fontsize=8, color=cols[i])
        else:
            ax.legend(fontsize=8, loc="best")
    elif legend == "box":
        ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax


def plot_tauc(res: dict, out: str | None = None, lang: str = "en", xlim=None, show_baseline: bool = True,
              title: str | None = None, ax=None, color="#1f5fa8", name: str | None = None):
    """Gráfico de Tauc com a reta ajustada estendida até o eixo, região de ajuste sombreada e Eg.
    Para sobrepor várias amostras, passe o mesmo `ax` e cores/nomes diferentes."""
    _style()
    T = TXT[lang]
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
    else:
        fig = ax.figure
    E, Y = res["E"], res["Y"]
    e1, e2 = res["fit_range"]
    if xlim is None:
        xlim = (max(E.min(), res["Eg"] - 1.0), min(E.max(), e2 + 0.8))
    m = (E >= xlim[0]) & (E <= xlim[1])
    ax.plot(E[m], Y[m], color=color, lw=1.3, label=name)
    ax.axvspan(e1, e2, color=color, alpha=0.10, lw=0)
    xx = np.linspace(res["Eg"], e2 + 0.15 * (e2 - e1) + 0.05, 50)
    ax.plot(xx, res["slope"] * xx + res["intercept"], "--", color="#c0392b", lw=1.1)
    if show_baseline and res.get("base_line") and np.isfinite(res.get("Eg_baseline", np.nan)):
        mb, bb, (b1, b2) = res["base_line"]
        xb = np.linspace(b1, res["Eg_baseline"] + 0.05, 30)
        ax.plot(xb, mb * xb + bb, ":", color="0.4", lw=1.0)
    ymax = np.nanmax(Y[m])
    ax.set_ylim(bottom=min(0, np.nanmin(Y[m])), top=ymax * 1.08)
    ax.set_xlim(*xlim)
    ax.axhline(0, color="0.6", lw=0.6)
    lab = f"$E_g$ = {res['Eg']:.2f} eV"
    if np.isfinite(res.get("Eg_spread", np.nan)) and res["Eg_spread"] >= 0.01:
        lab += f" (±{res['Eg_spread'] / 2:.2f})"
    ax.annotate(lab, xy=(res["Eg"], 0), xytext=(res["Eg"] - 0.05 * (xlim[1] - xlim[0]), 0.55 * ymax),
                ha="right", fontsize=9, color="#c0392b",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=0.8))
    ax.set_xlabel(T["hv"])
    sym = "α" if "α" in res["source"] else ("F(R)" if "Kubelka" in res["source"] else "A")
    pw = {2.0: "2", 0.5: "1/2"}.get(res["power"], f"{res['power']:.2g}")
    ax.set_ylabel(f"$({sym}\\cdot h\\nu)^{{{pw}}}$ ({T['au']})")
    if title:
        ax.set_title(title, fontsize=10)
    if name:
        ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax


def plot_calibration(cal: dict, out: str | None = None, lang: str = "en", unknowns=None,
                     xlabel: str | None = None, ylabel: str | None = None, title=None):
    """Curva analítica com faixa de confiança 95% da reta, equação, R², LOD/LOQ e resíduos.
    unknowns: lista de dicts de predict() para marcar amostras na curva."""
    from scipy import stats
    _style()
    T = TXT[lang]
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(5.2, 5.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    x, y = cal["conc"], cal["A"]
    xx = np.linspace(0 if cal["through_origin"] or x.min() < 0.2 * x.max() else x.min(), x.max() * 1.05, 100)
    yy = cal["slope"] * xx + cal["intercept"]
    t = stats.t.ppf(0.975, cal["dof"])
    if cal["through_origin"]:
        se = cal["syx"] * np.sqrt(xx ** 2 / np.sum(x ** 2))
    else:
        se = cal["syx"] * np.sqrt(1 / cal["n"] + (xx - cal["x_mean"]) ** 2 / cal["sxx"])
    ax.fill_between(xx, yy - t * se, yy + t * se, color="#1f5fa8", alpha=0.15, lw=0)
    ax.plot(xx, yy, color="#1f5fa8", lw=1.2)
    ax.plot(x, y, "o", color="0.1", ms=5, zorder=5)
    if unknowns:
        for u in unknowns:
            c = u["conc"] / u.get("dilution", 1.0)
            ax.errorbar(c, u["A"], xerr=(u["ci95"][1] - u["ci95"][0]) / 2 / u.get("dilution", 1.0), fmt="s",
                        color="#c0392b", ms=5, capsize=3, zorder=6)
    sign = "+" if cal["intercept"] >= 0 else "−"
    eq = f"A = {cal['slope']:.4g}·c" + ("" if cal["through_origin"] else f" {sign} {abs(cal['intercept']):.3g}")
    txt = f"{eq}\nR² = {cal['r2']:.5f}\nLOD = {cal['lod']:.3g}   LOQ = {cal['loq']:.3g} {cal.get('conc_unit', '')}"
    if "epsilon" in cal:
        txt += f"\nε = {cal['epsilon']:.4g} {cal['epsilon_unit']}"
    ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=8.5)
    ax.set_ylabel(ylabel or T["A"])
    axr.axhline(0, color="0.5", lw=0.7)
    axr.plot(x, cal["residuals"], "o", color="0.1", ms=4)
    axr.set_ylabel(T["resid"])
    cu = cal.get("conc_unit", "")
    axr.set_xlabel(xlabel or (f"{T['conc']} ({cu})" if cu else T["conc"]))
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, (ax, axr)


def plot_kinetics(kins, names=None, out: str | None = None, lang: str = "en", time_unit: str = "min",
                  colors=None):
    """Dois painéis: C/C₀ vs t (com fase escura, se houver) e forma linearizada com a reta.
    kins: um dict de kinetics() ou lista (várias amostras/catalisadores)."""
    _style()
    T = TXT[lang]
    if isinstance(kins, dict):
        kins = [kins]
    names = names or [f"#{i + 1}" for i in range(len(kins))]
    cols = colors or PALETTE
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.4, 3.6))
    for i, k in enumerate(kins):
        c = cols[i % len(cols)]
        t, C = k["t"], k["C_rel"]
        a1.plot(t, C, "o-", color=c, ms=4, lw=1, label=names[i])
        tt = np.linspace(0, t.max(), 100)
        a2.plot(t, k["linear_y"], "o", color=c, ms=4)
        a2.plot(tt, k["k"] * tt + k.get("intercept", 0.0), "-", color=c, lw=1,
                label=f"{names[i]}: k = {k['k']:.3g} {time_unit}$^{{-1}}$, R² = {k['r2']:.3f}"
                if k["model"] == "pfo" else f"{names[i]}: k = {k['k']:.3g}, R² = {k['r2']:.3f}")
    a1.set_xlabel(f"{T['time']} ({time_unit})")
    a1.set_ylabel("C/C$_0$")
    a1.set_ylim(bottom=0)
    a2.set_xlabel(f"{T['time']} ({time_unit})")
    a2.set_ylabel(kins[0]["linear_label"])
    a1.legend(fontsize=8)
    a2.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, (a1, a2)
