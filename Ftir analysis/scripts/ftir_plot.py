"""
ftir_plot.py - figuras de FT-IR com bandas rotuladas, sem sobreposição de rótulos

    from ftir_plot import plot_spectra
    fig, ax = plot_spectra([s1, s2], mode="A", offset="auto",
                           peaks=[tab1, tab2],            # DataFrames de find_bands (ou None)
                           highlight=[(1760, 1680, "C=O")],
                           lang="pt", out="figura.png")

Convenções que a função já segue (e que a figura final deve manter):
  * eixo x invertido (número de onda decrescente, 4000 -> 400 cm-1);
  * modo 'A' rotula acima dos máximos; modo 'T' (%T) rotula abaixo dos mínimos;
  * com offset, os valores absolutos do eixo y perdem significado -> ticks de y são ocultados
    e o rótulo do eixo diz 'deslocado'/'offset';
  * rótulos verticais, afastados entre si automaticamente, com linha-guia quando deslocados.
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms as mtransforms  # noqa: E402

from ftir_io import Spectrum  # noqa: E402

TXT = {
    "en": {"x": "Wavenumber (cm$^{-1}$)", "A": "Absorbance", "%T": "Transmittance (%)",
           "T": "Transmittance", "KM": "Kubelka–Munk", "%R": "Reflectance (%)", "R": "Reflectance",
           "norm": "normalized", "offset": "offset", "au": "a.u.", "eps": "Absorptivity",
           "d2": "2nd derivative", "d1": "1st derivative", "arb": "Intensity"},
    "pt": {"x": "Número de onda (cm$^{-1}$)", "A": "Absorbância", "%T": "Transmitância (%)",
           "T": "Transmitância", "KM": "Kubelka–Munk", "%R": "Reflectância (%)", "R": "Reflectância",
           "norm": "normalizada", "offset": "deslocado", "au": "u.a.", "eps": "Absortividade",
           "d2": "2ª derivada", "d1": "1ª derivada", "arb": "Intensidade"},
}

# paleta qualitativa sóbria e distinguível em P&B (ordem pensada para pilhas)
PALETTE = ["#1f1f1f", "#c0392b", "#1f5fa8", "#2e8b57", "#8e44ad", "#d35400", "#16a085", "#7f6000"]


def _style():
    plt.rcParams.update({
        "font.size": 10, "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "xtick.major.size": 4, "ytick.major.size": 4,
        "xtick.minor.visible": True, "xtick.minor.size": 2, "legend.frameon": False,
        "savefig.bbox": "tight", "savefig.dpi": 300, "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def _prep(s: Spectrum, mode: str) -> Spectrum:
    if mode == "A" and s.ykind in ("T", "%T", "R", "%R"):
        return s.to_absorbance()
    if mode in ("T", "%T") and s.ykind == "A":
        return s.to_transmittance(percent=True)
    if mode in ("T", "%T") and s.ykind == "T":
        return s.to_transmittance(percent=True)
    return s


def _spread(xs, width, lo, hi, iters=200):
    """Afasta posições 1D para que fiquem a >= width umas das outras, mantendo a ordem."""
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


def _label_text(row, fmt, assign_col):
    txt = fmt.format(position=row["position"], **{k: row[k] for k in row.index if k != "position"})
    if row.get("kind", "") == "shoulder" and "sh" not in txt:
        txt += " sh"
    if assign_col and assign_col in row and str(row[assign_col]) not in ("", "nan"):
        a = str(row[assign_col])
        if assign_col == "assignment":          # nome da tabela: corta o comentário entre parênteses
            a = a.split(" (")[0].split(";")[0].strip()
        txt += "  " + (a if len(a) <= 30 else a[:29] + "…")
    return txt


def plot_spectra(spectra, mode: str = "auto", offset=None, peaks=None, label_fmt: str = "{position}",
                 label_assign: str | None = None, highlight=None, xlim=(4000, 400), ylim=None,
                 labels=None, legend: str = "auto", lang: str = "en", colors=None, lw: float = 1.0,
                 fontsize: float = 7.5, figsize=(7.0, 4.2), title: str | None = None,
                 normalized: bool = False, out: str | None = None, ax=None):
    """Plota 1+ espectros.

    mode        'A' | 'T' (=%T) | 'auto' (usa o tipo do 1º espectro) | 'raw' (não converte)
    offset      None (sobrepostos) | 'auto' | número (deslocamento vertical entre traços)
    peaks       DataFrame (1 espectro) ou lista de DataFrames/None (um por espectro)
    label_fmt   formato do rótulo, ex. '{position}' ou '{position} ({intensity})'
    label_assign nome de coluna para anexar ao rótulo (ex. 'assignment' ou uma coluna sua)
    highlight   [(hi, lo, 'texto'), ...] faixas sombreadas
    labels      nomes dos traços (padrão: s.title)
    legend      'auto' | 'inline' (nome ao lado do traço) | 'box' | None
    """
    _style()
    if isinstance(spectra, Spectrum):
        spectra = [spectra]
    if peaks is not None and not isinstance(peaks, (list, tuple)):
        peaks = [peaks]
    peaks = list(peaks) if peaks is not None else [None] * len(spectra)
    peaks += [None] * (len(spectra) - len(peaks))
    if mode == "auto":
        k0 = spectra[0].ykind
        mode = "A" if k0 in ("A", "eps", "KM") else ("T" if k0 in ("T", "%T") else "raw")
    if mode == "T":
        mode = "%T"
    data = [_prep(s, mode) if mode != "raw" else s for s in spectra]
    T = TXT.get(lang, TXT["en"])
    colors = colors or PALETTE
    labels = labels or [s.title for s in data]
    xhi, xlo = max(xlim), min(xlim)

    def vis(s):
        m = (s.x <= xhi) & (s.x >= xlo)
        return s.x[m], s.y[m]

    n = len(data)
    any_labels = any(p is not None and len(p) for p in peaks)
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure
    step = 0.0
    if offset is not None and n > 1:
        if offset == "auto":
            spans = [np.nanmax(vis(s)[1]) - np.nanmin(vis(s)[1]) for s in data]
            step = 0.6 * max(spans)
        else:
            step = float(offset)
    # com offset 'auto' e rótulos, aumenta o deslocamento até os rótulos de um traço
    # não invadirem o traço vizinho (acima no modo A, abaixo no modo %T)
    for _ in range(8):
        ax.cla()
        offs = np.array([(n - 1 - i) * step for i in range(n)])
        placed = _draw(ax, fig, data, peaks, offs, mode, xhi, xlo, vis, colors, lw, labels, highlight,
                       label_fmt, label_assign, fontsize, n, any_labels)
        if not (offset == "auto" and n > 1 and any_labels):
            break
        need = _overlap_need(ax, placed, data, offs, mode, vis)
        if need <= 0:
            break
        step += need * 1.1
    if ylim:
        ax.set_ylim(*ylim)

    # ---- eixos e legenda
    ax.set_xlim(xhi, xlo)
    ax.set_xlabel(T["x"])
    yl_txt = T.get(mode if mode != "raw" else data[0].ykind, T["arb"])
    if normalized:
        yl_txt += f" ({T['norm']})"
    if offset is not None and n > 1:
        yl_txt = f"{yl_txt.split(' (')[0]} ({T['au']}, {T['offset']})"
        ax.set_yticklabels([])
        ax.tick_params(axis="y", which="both", left=False, right=False)
    ax.set_ylabel(yl_txt)
    if title:
        ax.set_title(title, fontsize=10)
    if legend == "auto":
        legend = "inline" if (offset is not None and n > 1) else ("box" if n > 1 else None)
    if legend == "box":
        ax.legend(fontsize=8, loc="best")
    elif legend == "inline":
        yr = np.diff(ax.get_ylim())[0]
        for i, s in enumerate(data):
            x, y = vis(s)
            if not len(x):
                continue
            # nome sobre o trecho de alto número de onda (normalmente plano, sem rótulos)
            m = x >= xhi - 0.09 * (xhi - xlo)
            seg = y[m] if m.any() else y[:max(1, len(y) // 12)]
            if mode == "%T":
                ax.text(xhi - 0.01 * (xhi - xlo), np.nanmin(seg) + offs[i] - 0.015 * yr, labels[i], ha="left",
                        va="top", fontsize=8, color=colors[i % len(colors)])
            else:
                ax.text(xhi - 0.01 * (xhi - xlo), np.nanmax(seg) + offs[i] + 0.015 * yr, labels[i], ha="left",
                        va="bottom", fontsize=8, color=colors[i % len(colors)])
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax


def _draw(ax, fig, data, peaks, offs, mode, xhi, xlo, vis, colors, lw, labels, highlight,
          label_fmt, label_assign, fontsize, n, any_labels):
    """Desenha traços, faixas e rótulos; devolve [(índice_do_traço, Text)]."""
    ax.set_xlim(xhi, xlo)
    if highlight:
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        for h in highlight:
            hi_, lo_ = max(h[0], h[1]), min(h[0], h[1])
            ax.axvspan(lo_, hi_, color="#f1c40f", alpha=0.18, lw=0, zorder=0)
            if len(h) > 2 and h[2]:
                ax.text((hi_ + lo_) / 2, 0.985, h[2], transform=trans, ha="center", va="top",
                        fontsize=fontsize, color="#7a5c00", zorder=5)
    for i, s in enumerate(data):
        x, y = vis(s)
        ax.plot(x, y + offs[i], color=colors[i % len(colors)], lw=lw, label=labels[i], zorder=3)
    ally = np.concatenate([vis(s)[1] + offs[i] for i, s in enumerate(data)])
    ymin, ymax = np.nanmin(ally), np.nanmax(ally)
    yr = ymax - ymin if ymax > ymin else 1.0
    head = 0.08 * yr
    ax.set_ylim(ymin - head * (2.5 if mode == "%T" and any_labels else 1),
                ymax + head * (2.5 if mode != "%T" and any_labels else 1))
    fig.canvas.draw()
    px_per_x = abs(ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0])
    label_w = (fontsize * 1.35 / 72 * fig.dpi) / px_per_x  # largura do texto vertical em cm-1
    pad = 0.02 * yr
    placed = []
    for i, (s, p) in enumerate(zip(data, peaks)):
        if p is None or len(p) == 0:
            continue
        x, y = s.x, s.y + offs[i]
        p = p[(p.position <= xhi) & (p.position >= xlo)]
        if not len(p):
            continue
        xa = p.position.to_numpy(float)
        xl = _spread(xa, label_w, xlo, xhi)
        xs_, ys_ = (x[::-1], y[::-1]) if x[0] > x[-1] else (x, y)
        for (_, row), x0, x1 in zip(p.iterrows(), xa, xl):
            y0 = float(np.interp(x0, xs_, ys_))
            a_, b_ = min(x0, x1) - label_w / 2, max(x0, x1) + label_w / 2
            win = (x >= a_) & (x <= b_)
            txt = _label_text(row, label_fmt, label_assign)
            if mode == "%T":
                yl = (np.nanmin(y[win]) if win.any() else y0) - pad
                va, sgn = "top", -1
            else:
                yl = (np.nanmax(y[win]) if win.any() else y0) + pad
                va, sgn = "bottom", 1
            if abs(x1 - x0) > 0.25 * label_w or abs(yl - y0) > 2 * pad:
                ax.plot([x0, x0, x1], [y0 + sgn * pad * 0.3, y0 + sgn * pad * 0.3 + (yl - y0) * 0.5, yl],
                        lw=0.5, color="0.45", zorder=4)
            t = ax.text(x1, yl + sgn * pad * 0.15, txt, rotation=90, ha="center", va=va,
                        fontsize=fontsize, color=colors[i % len(colors)] if n > 1 else "0.1", zorder=6)
            placed.append((i, t))
    # o texto tem altura fixa em pontos: ao ampliar ylim ele "encolhe" em unidades de dados,
    # então itera até todos os rótulos caberem dentro dos eixos
    for _ in range(6):
        if not placed:
            break
        fig.canvas.draw()
        inv = ax.transData.inverted()
        bb = [t.get_window_extent() for _, t in placed]
        ys = [inv.transform((b.x0, b.y1))[1] for b in bb] + [inv.transform((b.x0, b.y0))[1] for b in bb]
        lo_, hi_ = ax.get_ylim()
        cur = hi_ - lo_
        new_lo, new_hi = min(lo_, min(ys) - 0.02 * cur), max(hi_, max(ys) + 0.02 * cur)
        if new_lo >= lo_ and new_hi <= hi_:
            break
        ax.set_ylim(new_lo, new_hi)
    return placed


def _overlap_need(ax, placed, data, offs, mode, vis):
    """Quanto (em unidades de y) falta para os rótulos não tocarem o traço vizinho."""
    if not placed:
        return 0.0
    ax.figure.canvas.draw()
    inv = ax.transData.inverted()
    yr = np.diff(ax.get_ylim())[0]
    margin = 0.01 * yr
    need = 0.0
    for i, t in placed:
        j = i - 1 if mode != "%T" else i + 1          # vizinho acima (A) ou abaixo (%T)
        if j < 0 or j >= len(data):
            continue
        b = t.get_window_extent()
        (xa, ya), (xb, yb) = inv.transform((b.x0, b.y0)), inv.transform((b.x1, b.y1))
        lo_x, hi_x = min(xa, xb), max(xa, xb)
        x, y = vis(data[j])
        m = (x >= lo_x) & (x <= hi_x)
        if not m.any():
            continue
        if mode != "%T":
            need = max(need, max(ya, yb) - (np.nanmin(y[m]) + offs[j]) + margin)
        else:
            need = max(need, (np.nanmax(y[m]) + offs[j]) - min(ya, yb) + margin)
    return need


def plot_with_derivative(s: Spectrum, d2: Spectrum, peaks=None, out=None, lang="en", xlim=(4000, 400), **kw):
    """Painel duplo: espectro (A) em cima, 2ª derivada embaixo (ombros = mínimos)."""
    _style()
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 5.6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    plot_spectra(s, mode="A", peaks=peaks, ax=a1, lang=lang, xlim=xlim, **kw)
    a1.set_xlabel("")
    m = (d2.x <= max(xlim)) & (d2.x >= min(xlim))
    a2.plot(d2.x[m], d2.y[m], color="0.2", lw=0.9)
    a2.axhline(0, color="0.6", lw=0.5)
    a2.set_ylabel(TXT[lang]["d2"])
    a2.set_yticklabels([])
    a2.set_xlabel(TXT[lang]["x"])
    a2.set_xlim(max(xlim), min(xlim))
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, (a1, a2)
