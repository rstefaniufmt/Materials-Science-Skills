"""
ta_plot.py - figuras de análise térmica

    plot_tga(runs, steps=[df], lang="pt", out="tg.png")        # TG + DTG, etapas anotadas
    plot_dsc(r, segments=["heat1","cool1","heat2"], results=..., display="exo_up", out="dsc.png")
    plot_sta(r, steps=df, events=..., out="sta.png")           # TG/DTG em cima, DSC embaixo
    plot_kinetics(res, out="cinetica.png")                     # α(T) por β + Ea(α) por método

Convenções: DSC com seta 'exo ↑' (ou 'endo ↑') no canto - obrigatória em figura de DSC;
segmentos (aquecimento/resfriamento) empilhados com deslocamento e rotulados; DTG em %/°C no
eixo direito (picos para baixo = perda de massa).
"""
from __future__ import annotations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ta_io import Run  # noqa: E402
from ta_process import dtg, endo_signal, select, segments, mass_pct, _rate  # noqa: E402

PALETTE = ["#1f1f1f", "#c0392b", "#1f5fa8", "#2e8b57", "#8e44ad", "#d35400", "#16a085", "#7f6000"]
TXT = {
    "pt": {"T": "Temperatura (°C)", "t": "Tempo (min)", "mass": "Massa (%)", "dtg": "DTG (%/°C)",
           "hf": "Fluxo de calor", "exo": "exo", "endo": "endo", "res": "Resíduo", "Tprog": "Temperatura (°C)",
           "alpha": "Conversão, α", "Ea": "$E_a$ (kJ/mol)", "heat": "aquecimento", "cool": "resfriamento",
           "dta": "DTA", "au": "u.a."},
    "en": {"T": "Temperature (°C)", "t": "Time (min)", "mass": "Mass (%)", "dtg": "DTG (%/°C)",
           "hf": "Heat flow", "exo": "exo", "endo": "endo", "res": "Residue", "Tprog": "Temperature (°C)",
           "alpha": "Conversion, α", "Ea": "$E_a$ (kJ/mol)", "heat": "heating", "cool": "cooling",
           "dta": "DTA", "au": "a.u."},
}
METHOD_NAMES = {"friedman": "Friedman", "ofw": "Ozawa-Flynn-Wall", "kas": "KAS", "starink": "Starink",
                "vyazovkin": "Vyazovkin"}


def _style():
    plt.rcParams.update({
        "font.size": 10, "axes.linewidth": 0.9, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": False, "xtick.minor.visible": True, "ytick.minor.visible": True,
        "legend.frameon": False, "savefig.bbox": "tight", "savefig.dpi": 300, "svg.fonttype": "none",
        "pdf.fonttype": 42})


def _hf_label(u, lang, kind="DSC"):
    T = TXT[lang]
    name = T["dta"] if kind == "DTA" else T["hf"]
    return f"{name} ({u})" if u else name


def _arrow(ax, display, lang, loc="upper left"):
    T = TXT[lang]
    word = T["exo"] if display == "exo_up" else T["endo"]
    x = 0.015 if "left" in loc else 0.985
    ax.annotate(f"{word} ↑", xy=(x, 0.97), xycoords="axes fraction", ha="left" if "left" in loc else "right",
                va="top", fontsize=9, fontweight="bold")


# ============================================================================ TG / DTG
def plot_tga(runs, steps=None, labels=None, x: str = "auto", show_dtg: bool = True, annotate: bool = True,
             basis=None, lang: str = "pt", out: str | None = None, figsize=(6.6, 4.4), title=None,
             window_K: float = 6.0, colors=None, xlim=None):
    """TG (massa %) e DTG. x='T' | 'time' | 'auto' (tempo se houver resfriamento ou isotermas
    longas). Com uma curva e `steps`, anota Δm de cada etapa e o pico da DTG."""
    _style()
    if isinstance(runs, Run):
        runs = [runs]
    steps = steps if isinstance(steps, list) else [steps] + [None] * (len(runs) - 1)
    labels = labels or [r.title for r in runs]
    cols = colors or PALETTE
    T = TXT[lang]
    if x == "auto":
        segs = segments(runs[0])
        iso_t = sum((s["t1"] or 0) - (s["t0"] or 0) for s in segs if s["type"] == "iso")
        tot = (runs[0].time[-1] - runs[0].time[0]) if runs[0].has("time") else 1
        x = "time" if (any(s["type"] == "cool" and abs(s["T1"] - s["T0"]) > 20 for s in segs)
                       or sum(s["type"] == "heat" for s in segs) > 1 and iso_t > 0.3 * tot) else "T"
        if not runs[0].has("time"):
            x = "T"
    fig, ax = plt.subplots(figsize=figsize)
    ax2 = ax.twinx() if (show_dtg or x == "time") else None
    for i, r in enumerate(runs):
        m, unit = mass_pct(r, basis)
        xv = r.T if x == "T" else r.time
        ax.plot(xv, m, color=cols[i % len(cols)], lw=1.4, label=labels[i])
        if x == "T" and show_dtg:
            d = dtg(r, "T", window_K=window_K, y=m)
            heat = np.zeros(r.n, bool)
            rr = _rate(r) if r.has("time") else np.full(r.n, np.nan)
            for s in segments(r):
                if s["type"] == "heat" and abs(s["T1"] - s["T0"]) > 20:
                    heat[s["i0"]:s["i1"]] = True
                    if s["rate"]:
                        heat[s["i0"]:s["i1"]] &= ~(np.abs(rr[s["i0"]:s["i1"]]) < 0.8 * abs(s["rate"]))
            ax2.plot(np.where(heat, xv, np.nan), np.where(heat, d, np.nan), color=cols[i % len(cols)], lw=1.0,
                     ls="--", alpha=0.85)
        if x == "time" and i == 0 and r.has("T"):
            ax2.plot(r.time, r.T, color="0.55", lw=0.9, ls=":")
    ax.set_xlabel(T["T"] if x == "T" else T["t"])
    ax.set_ylabel(T["mass"] if basis is None else f"{T['mass'][:-4]} ({unit})")
    if ax2 is not None:
        ax2.set_ylabel(T["dtg"] if x == "T" else T["Tprog"], color="0.35")
        ax2.tick_params(axis="y", colors="0.35", direction="in")
    if xlim:
        ax.set_xlim(*xlim)
    # anotações (uma curva)
    if annotate and len(runs) == 1 and steps[0] is not None and len(steps[0]) and x == "T":
        r, st = runs[0], steps[0]
        m, _ = mass_pct(r, basis)
        ylo, yhi = ax.get_ylim()
        span = yhi - ylo
        for _, s in st.iterrows():
            if not str(s["tipo"]).startswith(("perda", "ganho")):
                continue
            xm = s["T_pico"]
            ya, yb = s["massa_inicio"], s["massa_fim"]
            xb = s["T_fim"] if np.isfinite(s["T_fim"]) else xm
            ax.annotate("", xy=(xb, yb), xytext=(xb, ya),
                        arrowprops=dict(arrowstyle="<->", color="0.25", lw=0.8, shrinkA=0, shrinkB=0))
            ax.text(xb + 0.01 * (ax.get_xlim()[1] - ax.get_xlim()[0]), (ya + yb) / 2,
                    f"{'−' if s['perda_pct'] > 0 else '+'}{abs(s['perda_pct']):.1f}%", fontsize=8.5, va="center")
            if ax2 is not None and np.isfinite(xm):
                dmin = s["DTG_max_pct_C"]
                ax2.annotate(f"{xm:.0f} °C", xy=(xm, dmin), xytext=(0, -10), textcoords="offset points",
                             ha="center", va="top", fontsize=7.5, color="0.3")
        fin = m[np.isfinite(m)][-1]
        ax.text(0.985, 0.97, f"{T['res']} ({r.T[np.isfinite(m)][-1]:.0f} °C): {fin:.1f}%", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5)
        if ax2 is not None:
            lo2, hi2 = ax2.get_ylim()
            ax2.set_ylim(lo2 - 0.12 * (hi2 - lo2), hi2 + 0.6 * (hi2 - lo2))   # DTG na parte de baixo
    elif ax2 is not None and x == "T":
        lo2, hi2 = ax2.get_ylim()
        ax2.set_ylim(lo2 - 0.05 * (hi2 - lo2), hi2 + 0.6 * (hi2 - lo2))
    if len(runs) > 1:
        ax.legend(fontsize=8, loc="lower left")
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax


# ============================================================================ DSC / DTA
def _disp(y_endo_pos, display):
    return y_endo_pos if display == "endo_up" else -y_endo_pos


def plot_dsc(r: Run, segs=None, results=None, display: str = "exo_up", offset="auto", endo=None,
             lang: str = "pt", out: str | None = None, figsize=(6.6, 4.6), title=None, xlim=None, labels=None):
    """Curvas DSC/DTA por segmento (empilhadas) com eventos anotados.
    segs: lista de chaves ('heat1','cool1','heat2'...) ou None (todos os aquecimentos/resfriamentos).
    results: dict {chave_segmento: [resultado de integrate_peak / glass_transition, ...]}.
    display: 'exo_up' (padrão TA/Mettler/NETZSCH) ou 'endo_up'."""
    _style()
    T = TXT[lang]
    allsegs = [s for s in segments(r) if s["type"] in ("heat", "cool") and abs(s["T1"] - s["T0"]) > 10]
    if segs is None:
        keys = [s["key"] for s in allsegs] or [None]
    else:
        keys = segs
    parts = []
    for k in keys:
        sr = select(r, k) if k else r
        Tt, t, y, u, d, src = endo_signal(sr, endo)
        yd = _disp(y, display)
        info = sr.meta.get("segment_info") or {}
        if info.get("rate"):
            # oculta o transiente inicial/final (taxa ainda fora da programada): vira traço vertical
            rr = _rate(sr)
            yd = np.where(np.abs(rr) >= 0.8 * abs(info["rate"]), yd, np.nan)
        parts.append((k, sr, Tt, yd, u))
    fig, ax = plt.subplots(figsize=figsize)
    step = 0.0
    if len(parts) > 1 and offset is not None:
        # empilhar: cada curva relativa à própria mediana (aquecimento e resfriamento têm linhas de
        # base em lados opostos), espaçamento = maior excursão
        parts = [(k, sr, Tt, yd - np.nanmedian(yd), u) for k, sr, Tt, yd, u in parts]
        ups = [np.nanmax(p[3]) for p in parts]
        dns = [np.nanmin(p[3]) for p in parts]
        step = (1.1 * max(a - b for a, b in zip(ups, dns))) if offset == "auto" else float(offset)
    spans = [np.nanpercentile(p[3], 99.5) - np.nanpercentile(p[3], 0.5) for p in parts]
    order = list(range(len(parts)))
    for j, (k, sr, Tt, yd, u) in enumerate(parts):
        off = (len(parts) - 1 - j) * step
        col = PALETTE[j % len(PALETTE)]
        ax.plot(Tt, yd + off, color=col, lw=1.3)
        info = sr.meta.get("segment_info", {})
        lab = (labels[j] if labels else None) or (info.get("label") if info else sr.title)
        if info.get("rate"):
            lab += f" ({abs(info['rate']):.0f} K/min)"
        okp = np.isfinite(yd)
        xl_ = np.nanmin(Tt[okp]) if xlim is None else min(xlim)
        startzone = okp & (Tt <= xl_ + 0.25 * (np.nanmax(Tt) - np.nanmin(Tt)))
        ystart = np.nanmin(yd[startzone]) if startzone.any() else np.nanmedian(yd)
        ax.text(xl_, ystart + off - 0.05 * max(spans), lab, fontsize=8, color=col, va="top", ha="left")
        for res in (results or {}).get(k, []):
            if "dH" in res:
                p = res["_plot"]
                xs = p["x"]
                # curva no referencial plotado; linha de base = curva − sinal líquido (endo positivo)
                yy = np.interp(xs, Tt[np.argsort(Tt)], yd[np.argsort(Tt)])
                net_e = (p["y"] - p["b"]) * (1 if res["endo"] == "up" else -1)
                base_plot = yy - _disp(net_e, display)
                ax.fill_between(xs, yy + off, base_plot + off, color=col, alpha=0.15, lw=0)
                ax.plot(xs, base_plot + off, color=col, lw=0.6, ls="--")
                ip = int(np.argmax(np.abs(net_e)))
                lab2 = f"{res['T_pico']:.1f} °C\n{res['dH']:.1f} {res['unidade_dH']}"
                up = (yy[ip] - base_plot[ip]) > 0
                ax.annotate(lab2, xy=(xs[ip], yy[ip] + off), xytext=(6, 8 if up else -8), textcoords="offset points",
                            fontsize=7.5, va="bottom" if up else "top", color=col)
            elif "T_mid" in res:
                Tm = res["T_mid"] if np.isfinite(res["T_mid"]) else res["T_inflexao"]
                yy = float(np.interp(Tm, Tt[np.argsort(Tt)], yd[np.argsort(Tt)]))
                ax.plot([Tm], [yy + off], "o", ms=4, color=col)
                txt = f"$T_g$ = {Tm:.1f} °C" + (f"\nΔ$C_p$ = {res['dCp']:.2f} J/(g·K)" if np.isfinite(res.get("dCp", np.nan)) else "")
                ax.annotate(txt, xy=(Tm, yy + off), xytext=(-8, 12), textcoords="offset points", fontsize=7.5,
                            ha="right", color=col, arrowprops=dict(arrowstyle="-", color=col, lw=0.6))
    ax.margins(y=0.14)
    u = parts[0][4]
    kind = r.kind if r.kind in ("DSC", "DTA") else "DSC"
    ax.set_ylabel(_hf_label(u, lang, kind) + (f" — {T['au']}, deslocado" if step and lang == "pt" else
                                               (f" — {T['au']}, offset" if step else "")))
    if step:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", which="both", left=False)
    ax.set_xlabel(T["T"])
    if xlim:
        ax.set_xlim(*xlim)
    _arrow(ax, display, lang, "upper right")
    if title:
        ax.set_title(title, fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax


# ============================================================================ STA
def plot_sta(r: Run, steps=None, results=None, display: str = "exo_up", endo=None, lang: str = "pt",
             out: str | None = None, figsize=(6.6, 6.4), title=None, window_K: float = 6.0):
    """Dois painéis com o mesmo eixo de temperatura: TG (+DTG) em cima, DSC/DTA embaixo
    (só o maior aquecimento)."""
    _style()
    T = TXT[lang]
    h = select(r, "main_heat")
    fig, (a1, a2) = plt.subplots(2, 1, figsize=figsize, sharex=True, gridspec_kw={"height_ratios": [1.1, 1]})
    m, _ = mass_pct(h)
    a1.plot(h.T, m, color=PALETTE[0], lw=1.4)
    a1b = a1.twinx()
    d = dtg(h, "T", window_K=window_K, y=m)
    info = h.meta.get("segment_info") or {}
    if info.get("rate") and h.has("time"):
        rr = _rate(h)
        d = np.where(np.abs(rr) >= 0.8 * abs(info["rate"]), d, np.nan)
    a1b.plot(h.T, d, color=PALETTE[2], lw=1.0, ls="--")
    a1.set_ylabel(T["mass"])
    a1b.set_ylabel(T["dtg"], color=PALETTE[2])
    a1b.tick_params(axis="y", colors=PALETTE[2], direction="in")
    if steps is not None and len(steps):
        for _, s in steps.iterrows():
            if str(s["tipo"]) in ("perda", "ganho") and np.isfinite(s["T_pico"]):
                a1.axvline(s["T_pico"], color="0.75", lw=0.6, ls=":")
                a2.axvline(s["T_pico"], color="0.75", lw=0.6, ls=":")
                a1.text(s["T_pico"], 0.97, f" −{s['perda_pct']:.1f}%", fontsize=7.5, va="top", ha="right", rotation=90,
                        transform=a1.get_xaxis_transform())
    Tt, t, y, u, dd, src = endo_signal(h, endo)
    yd = _disp(y, display)
    a2.plot(Tt, yd, color=PALETTE[1], lw=1.3)
    for res in (results or []):
        if "dH" in res:
            p = res["_plot"]
            net_e = (p["y"] - p["b"]) * (1 if res["endo"] == "up" else -1)
            yy = np.interp(p["x"], Tt[np.argsort(Tt)], yd[np.argsort(Tt)])
            base_plot = yy - _disp(net_e, display)
            a2.fill_between(p["x"], yy, base_plot, color=PALETTE[1], alpha=0.15, lw=0)
            ip = int(np.argmax(np.abs(net_e)))
            a2.annotate(f"{res['T_pico']:.0f} °C ({res['tipo']})\n{res['dH']:.0f} {res['unidade_dH']}",
                        xy=(p["x"][ip], yy[ip]), xytext=(5, 6), textcoords="offset points", fontsize=7.5)
    a2.margins(y=0.12)
    a2.set_ylabel(_hf_label(u, lang, "DTA" if r.kind == "DTA" else "DSC"))
    a2.set_xlabel(T["T"])
    _arrow(a2, display, lang, "upper right")
    if title:
        a1.set_title(title, fontsize=10)
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, (a1, a2)


# ============================================================================ cinética
def plot_kinetics(res: dict, methods=("friedman", "ofw", "kas", "vyazovkin"), lang: str = "pt",
                  out: str | None = None, figsize=(9.0, 3.8)):
    _style()
    T = TXT[lang]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=figsize)
    for i, c in enumerate(res["conv"]):
        a1.plot(c["T"] - 273.15, c["alpha"], color=PALETTE[i % len(PALETTE)], lw=1.3, label=f"{c['beta']:.1f} K/min")
    a1.set_xlabel(T["T"])
    a1.set_ylabel(T["alpha"])
    a1.set_ylim(-0.02, 1.02)
    a1.legend(fontsize=8, loc="upper left", title="β", title_fontsize=8)
    tb = res["table"]
    mk = {"friedman": "o", "ofw": "s", "kas": "^", "starink": "v", "vyazovkin": "D"}
    for j, m in enumerate(methods):
        if m not in tb:
            continue
        err = tb.get(f"{m}_ic95") if j == 0 else None          # barras só no 1º método (legibilidade)
        a2.errorbar(tb["alpha"], tb[m], yerr=err if err is not None and np.isfinite(err).any() else None,
                    marker=mk.get(m, "o"), ms=4, lw=1, capsize=2, color=PALETTE[j % len(PALETTE)],
                    label=METHOD_NAMES.get(m, m))
    a2.set_xlabel(T["alpha"])
    a2.set_ylabel(T["Ea"])
    a2.set_xlim(0, 1)
    a2.legend(fontsize=8)
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, (a1, a2)


def plot_oit(r: Run, res: dict, endo=None, lang: str = "pt", out: str | None = None, figsize=(6.4, 4.0)):
    _style()
    Tt, t, y, u, d, src = endo_signal(r, endo)
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(t, -y, color=PALETTE[0], lw=1.2)
    ax.axvline(res["t_switch"], color=PALETTE[2], ls="--", lw=0.8)
    ax.axvline(res["t_onset"], color=PALETTE[1], ls="--", lw=0.8)
    ax.annotate("", xy=(res["t_onset"], 0.9), xytext=(res["t_switch"], 0.9), xycoords=("data", "axes fraction"),
                arrowprops=dict(arrowstyle="<->", color="0.2"))
    ax.text((res["t_onset"] + res["t_switch"]) / 2, 0.92, f"OIT = {res['OIT']:.1f} min", transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=9)
    ax.set_xlabel(TXT[lang]["t"])
    ax.set_ylabel(_hf_label(u, lang))
    _arrow(ax, "exo_up", lang, "upper left")
    fig.tight_layout()
    if out:
        fig.savefig(out)
    return fig, ax
