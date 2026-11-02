"""
nmr_plot.py - figuras de RMN: ¹H com integrais e rótulos de multipletos, zoom dos multipletos com
árvore de J, ¹³C com multiplicidade DEPT, comparação experimental × predito das candidatas

    from nmr_plot import plot_1h, plot_multiplets, plot_13c, plot_candidate, plot_stack
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe

plt.rcParams.update({"font.size": 9, "axes.linewidth": 0.8, "savefig.dpi": 300})
MCOL = {"CH3": "#0072B2", "CH2": "#D55E00", "CH": "#009E73", "C": "#000000", "CH/CH3": "#56B4E9", "C/CH2": "#E69F00",
        "": "#333333", "?": "#999999"}
T = {"pt": dict(ppm="δ (ppm)", sol="solvente", imp="impureza"), "en": dict(ppm="δ (ppm)", sol="solvent", imp="impurity")}


def _xlim(s, res_tab, pad=0.4, nucleus="1H"):
    if res_tab is not None and len(res_tab):
        hi = float(res_tab.ppm.max()) + pad
        lo = float(res_tab.ppm.min()) - pad
        lo = min(lo, 0.0 if nucleus == "1H" else lo)
        return max(hi, lo + 1), lo
    return float(s.ppm.max()), float(s.ppm.min())


def plot_1h(s, res, out=None, xlim=None, lang="pt", title=None, figsize=(10, 4.2)):
    t = res["table"]
    fig, ax = plt.subplots(figsize=figsize)
    hi, lo = xlim or _xlim(s, t)
    m = (s.ppm <= hi) & (s.ppm >= lo)
    ax.plot(s.ppm[m], s.y[m], color="k", lw=0.6)
    ymax = np.percentile(s.y[m], 99.95) if m.any() else 1
    samp = t[t.amostra] if len(t) else t
    big = samp.apply(lambda r: s.y[r.i0:r.i1 + 1].max(), axis=1).max() if len(samp) else ymax
    top = 1.18 * big
    ax.set_ylim(-0.2 * top, top * 1.18)
    # integrais (só sinais da amostra)
    nmax = max(1.0, float(t.nH.max())) if len(t) else 1.0
    placed = []
    for r in t.itertuples():
        if not r.amostra or r.nH < 0.3:
            continue
        seg = s.y[r.i0:r.i1 + 1]
        x = s.ppm[r.i0:r.i1 + 1]
        cum = np.cumsum(seg)
        if cum[-1] <= 0:
            continue
        col = "#0072B2"
        base = top * 0.95
        ax.plot(x, base + 0.15 * top * cum / cum[-1] * min(1, max(0.25, r.nH / nmax)), color=col, lw=0.8)
        xm = float(np.mean(x))
        row = 0
        while any(abs(xm - px) < 0.035 * (hi - lo) and pr == row for px, pr in placed):
            row += 1
        placed.append((xm, row))
        ax.text(xm, -0.04 * top - row * 0.05 * top, f"{r.nH:.2f}", ha="center", va="top", fontsize=6.5, color=col)
    # rótulos
    lab_y = []
    for r in t.itertuples():
        if r.ppm > hi or r.ppm < lo:
            continue
        h = s.y[r.i0:r.i1 + 1].max()
        if r.amostra and r.nH >= 0.5:
            txt = f"{r.ppm:.2f} {r.mult}" + (f" ({', '.join(f'{j:.1f}' for j in r.J)})" if r.J else "")
            col = "#0072B2"
        elif r.rotulo:
            txt = r.rotulo.replace("possível ", "") if "satélite" not in r.rotulo else ""
            col = "#999999"
        else:
            continue
        if not txt:
            continue
        y = min(h + 0.03 * top, 0.86 * top)
        ax.text(r.ppm, y, txt, rotation=90, ha="center", va="bottom", fontsize=6.5, color=col,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_xlim(hi, lo)
    ax.set_xlabel(T[lang]["ppm"])
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    from nmr_data import SOLVENTS
    sol = SOLVENTS.get(s.solvent or "", {}).get("name", s.solvent or "?")
    ax.set_title(title or f"{s.name}: ¹H ({s.sf:.0f} MHz, {sol})", fontsize=9, loc="left")
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_multiplets(s, res, out=None, max_n=24, ncols=6):
    """Zoom de cada multipleto (eixo em Hz relativo ao centro) com a árvore de acoplamentos."""
    t = res["table"]
    t = t[t.amostra & (t.nH >= 0.5)]
    t = t.iloc[:max_n]
    n = len(t)
    if n == 0:
        return None
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(1.9 * ncols, 1.7 * nrows), squeeze=False)
    for ax, r in zip(axs.ravel(), t.itertuples()):
        x = (s.ppm[r.i0:r.i1 + 1] - r.ppm) * s.sf
        y = s.y[r.i0:r.i1 + 1]
        ax.plot(x, y, color="k", lw=0.7)
        ax.set_xlim(x.max(), x.min())
        ax.set_yticks([])
        ax.tick_params(labelsize=6)
        ax.set_title(f"{r.ppm:.2f} {r.mult}" + (f"\nJ = {', '.join(f'{j:.1f}' for j in r.J)} Hz" if r.J else "")
                     + f"  ({r.nH_int:.0f}H)", fontsize=6.5)
        # árvore
        if r.J and r.mult not in ("m", "AB"):
            from nmr_data import MULT_NAMES
            inv = {v: k for k, v in MULT_NAMES.items()}
            letters = []
            k = 0
            mult = r.mult
            while k < len(mult):
                for L in sorted(inv, key=len, reverse=True):
                    if mult.startswith(L, k):
                        letters.append(inv[L])
                        k += len(L)
                        break
                else:
                    k += 1
            ytop = y.max()
            L = sum(letters)
            step = 0.12 * ytop
            lev = ytop * 1.08 + L * step          # raiz no alto, ramos descendo até as linhas
            ax.set_ylim(-0.05 * ytop, lev + 0.1 * ytop)
            pos = [0.0]
            for J, cnt in zip(r.J, letters):
                for _ in range(cnt):
                    new = []
                    for p_ in pos:
                        for q in (p_ - J / 2, p_ + J / 2):
                            ax.plot([p_, q], [lev, lev - step], color="#0072B2", lw=0.5)
                            new.append(q)
                    pos = new
                    lev -= step
            # inverte a árvore (de cima para baixo)
    for ax in axs.ravel()[n:]:
        ax.set_axis_off()
    fig.supxlabel("Hz (relativo ao centro)", fontsize=7)
    fig.tight_layout()
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_13c(s, res, dept135=None, dept90=None, out=None, lang="pt", title=None, figsize=(10, 4.2)):
    t = res["table"]
    n_tr = 1 + (dept135 is not None) + (dept90 is not None)
    fig, axs = plt.subplots(n_tr, 1, figsize=(figsize[0], figsize[1] + 1.2 * (n_tr - 1)), sharex=True,
                            gridspec_kw=dict(height_ratios=[3] + [1] * (n_tr - 1), hspace=0.05))
    axs = np.atleast_1d(axs)
    samp = t[t.rotulo.eq("")] if len(t) else t
    hi = (samp.ppm.max() + 8) if len(samp) else s.ppm.max()
    lo = (samp.ppm.min() - 8) if len(samp) else s.ppm.min()
    hi, lo = min(hi, s.ppm.max()), max(lo, s.ppm.min())
    ax = axs[0]
    m = (s.ppm <= hi) & (s.ppm >= lo)
    ax.plot(s.ppm[m], s.y[m], color="k", lw=0.6)
    top = samp.altura.max() * 1.15 if len(samp) else s.y[m].max()
    ax.set_ylim(-0.05 * top, top * 1.35)
    for r in t.itertuples():
        if r.ppm > hi or r.ppm < lo:
            continue
        if r.rotulo == "":
            col = MCOL.get(r.mult, "#333333")
            txt = f"{r.ppm:.1f}" + (f" {r.mult}" if r.mult else "")
        elif r.rotulo == "solvente" or r.rotulo.startswith("pequeno"):
            continue
        else:
            col, txt = "#999999", f"{r.ppm:.1f} {r.rotulo.replace('possível ', '')}"[:28]
        ax.text(r.ppm, min(r.altura, top) + 0.03 * top, txt, rotation=90, ha="center", va="bottom", fontsize=6.5,
                color=col, path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    ax.set_yticks([])
    from nmr_data import SOLVENTS
    sol = SOLVENTS.get(s.solvent or "", {}).get("name", s.solvent or "?")
    ax.set_title(title or f"{s.name}: ¹³C ({s.sf:.0f} MHz, {sol})", fontsize=9, loc="left")
    k = 1
    for d, lab in ((dept135, "DEPT-135"), (dept90, "DEPT-90")):
        if d is None:
            continue
        a = axs[k]
        mm = (d.ppm <= hi) & (d.ppm >= lo)
        a.plot(d.ppm[mm], d.y[mm], color="#555555", lw=0.6)
        lim = np.abs(d.y[mm]).max() * 1.1
        a.set_ylim(-lim, lim)
        a.set_yticks([])
        a.text(0.005, 0.85, lab, transform=a.transAxes, fontsize=7)
        k += 1
    for a in axs:
        for sp in ("left", "right", "top"):
            a.spines[sp].set_visible(False)
    axs[-1].set_xlim(hi, lo)
    axs[-1].set_xlabel(T[lang]["ppm"])
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_candidate(name, det, out=None, title=None):
    """Experimental (acima) × predito (abaixo), com linhas da atribuição, para ¹³C e ¹H."""
    have_c = bool(det.get("c13"))
    have_h = bool(det.get("h1"))
    n = have_c + have_h
    if n == 0:
        return None
    fig, axs = plt.subplots(n, 1, figsize=(10, 2.6 * n), squeeze=False)
    k = 0
    for nuc, key in (("¹³C", "c13"), ("¹H", "h1")):
        lst = det.get(key) or []
        if not lst:
            continue
        ax = axs[k, 0]
        k += 1
        exp_all = sorted({a["exp"] for a in lst if np.isfinite(a.get("exp", np.nan))} |
                         set(det.get("extra_c" if key == "c13" else "extra_h", [])))
        for e in exp_all:
            ax.plot([e, e], [0.15, 1.0], color="k", lw=1.0)
        for e in det.get("extra_c" if key == "c13" else "extra_h", []):
            ax.plot([e, e], [0.15, 1.0], color="#D55E00", lw=1.2)
            ax.text(e, 1.03, "?", color="#D55E00", ha="center", fontsize=7)
        for a in lst:
            if a.get("exch"):
                continue
            p = a["pred"]
            if not np.isfinite(p):
                continue
            ok = np.isfinite(a.get("exp", np.nan))
            col = "#0072B2" if ok and abs(a["err"]) <= (3 if key == "c13" else 0.3) else "#D55E00"
            ax.plot([p, p], [-1.0, -0.15], color=col, lw=1.0)
            if ok:
                ax.plot([p, a["exp"]], [-0.15, 0.15], color=col, lw=0.6, alpha=0.8)
            ax.text(p, -1.05, a["label"], rotation=90, ha="center", va="top", fontsize=5.5, color=col)
        allx = [a["pred"] for a in lst if np.isfinite(a["pred"])] + exp_all
        if allx:
            pad = 5 if key == "c13" else 0.4
            ax.set_xlim(max(allx) + pad, min(allx) - pad)
        ax.set_ylim(-1.6, 1.15)
        ax.set_yticks([0.6, -0.6], ["exp.", "prev."])
        ax.axhline(0, color="0.85", lw=0.5)
        ax.set_xlabel(f"δ {nuc} (ppm)")
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)
    axs[0, 0].set_title(title or f"Atribuição: {name}  (azul: |erro| ≤ 3 ppm ¹³C / 0,3 ppm ¹H; laranja: maior ou sem sinal)",
                        fontsize=8.5, loc="left")
    fig.tight_layout()
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig


def plot_stack(specs, names=None, out=None, xlim=None, offset=1.1):
    fig, ax = plt.subplots(figsize=(10, 1.4 + 1.1 * len(specs)))
    for k, s in enumerate(specs):
        m = slice(None)
        y = s.y / (np.percentile(s.y, 99.9) or 1)
        ax.plot(s.ppm, y + k * offset, lw=0.6, color=plt.cm.viridis(k / max(1, len(specs) - 1)))
        ax.text(1.0, k * offset + 0.3, (names[k] if names else s.name), transform=ax.get_yaxis_transform(), fontsize=7,
                ha="right")
    hi = xlim[0] if xlim else max(s.ppm.max() for s in specs)
    lo = xlim[1] if xlim else min(s.ppm.min() for s in specs)
    ax.set_xlim(hi, lo)
    ax.set_yticks([])
    ax.set_xlabel("δ (ppm)")
    if out:
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    return fig
