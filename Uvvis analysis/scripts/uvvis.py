#!/usr/bin/env python3
"""
uvvis.py - linha de comando da skill de UV-Vis

  info      arquivos...                          metadados, faixa, tipo de ordenada, avisos
  plot      arquivos... -o fig.png [--mode A|T|R|KM] [--label] [--top 5] [--series 0 10 20
            --series-label "t (min)"] [--offset auto] [--ev-axis] [--xlim 300 800]
            [--offset-correct 750 800] [--scatter 600 800] [--smooth 5] [--normalize max]
  peaks     arquivo [-o bandas.csv] [--region 300 800] [--top 10]
  calib     (--table curva.csv | padroes... --conc 1 2 5 10) --wl 664 [--base 800]
            [--path 1 --unit mol/L] [--origin] [--unknown amostra.csv ... | --unknown-A 0.42 0.43]
            [--dilution 10] -o curva.png
  kinetics  (--table dados.csv | espectros... --times -30 0 10 20 ...) [--wl 664 | --window 600 700]
            [--model pfo|pso|zero|all] -o cinetica.png
  tauc      arquivo [--transition direct|indirect|both] [--window 2.5 4] [--fit 3.2 3.5]
            [--thickness 2e-5] [--ykind %R] -o tauc.png
  isosbestic espectros... [--region 300 700]
  wlcheck   arquivo_holmio [--standard holmium_glass|holmium_solution|didymium_glass]
  convert   arquivo -o saida.csv [--to A|T|KM] [--decimal ,] [--ev]
  report    arquivo -d pasta [--lang pt]       figura com λmax + tabela + resumo .md

Formatos: .csv/.txt/.prn/.asc (inclusive Cary pareado, Shimadzu, OceanView, vírgula decimal),
.xlsx, .jdx/.dx, .spc (Galactic), .sp (PerkinElmer). Curingas funcionam ("dados/*.csv").
--ykind força o tipo de ordenada quando o arquivo não diz (ex.: %R de reflectância difusa).
"""
from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from uvvis_io import HC_EV_NM, load_all, save_csv  # noqa: E402
from uvvis_process import normalize, offset_correct, scatter_correct, smooth, detect_steps  # noqa: E402
from uvvis_peaks import band_metrics, find_bands, isosbestic, wavelength_check  # noqa: E402
from uvvis_quant import calibrate, compare_kinetic_models, kinetics, predict  # noqa: E402
from uvvis_bandgap import derivative_edge, onset_wavelength, tauc, urbach  # noqa: E402
from uvvis_plot import plot_calibration, plot_kinetics, plot_spectra, plot_tauc  # noqa: E402


def _expand(paths):
    out = []
    for p in paths:
        g = sorted(glob.glob(p))
        out += g if g else [p]
    return out


def _load(paths, a):
    specs = []
    for p in _expand(paths):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            got = load_all(p)
        for ww in w:
            print(f"[aviso] {os.path.basename(p)}: {ww.message}", file=sys.stderr)
        if getattr(a, "index", None) is not None:
            got = [got[a.index]]
        for s in got:
            if getattr(a, "ykind", None):
                s.ykind = a.ykind
                s.meta["ykind_guessed"] = False
        specs += got
    return specs


def _stem(s):
    """Nome curto: arquivo para espectro único; título (coluna/amostra) quando o arquivo tem vários."""
    base = os.path.splitext(os.path.basename(s.source))[0]
    fmt = s.meta.get("format", "")
    multi = (s.title and s.title != base and (s.meta.get("column") is not None or "pareado" in fmt
             or "largo" in fmt or "Excel" in fmt)) or bool(s.meta.get("subfile"))
    return s.title if multi else base


def _pre(s, a):
    if getattr(a, "offset_correct", None):
        s = offset_correct(s, tuple(a.offset_correct))
    if getattr(a, "scatter", None):
        s, info = scatter_correct(s, tuple(a.scatter))
        print(f"[espalhamento] {s.title}: n = {info['n']:.2f}", file=sys.stderr)
    if getattr(a, "smooth", None):
        s = smooth(s, a.smooth)
    if getattr(a, "normalize", None):
        s = normalize(s, a.normalize, wl=getattr(a, "norm_wl", None))
    return s


def cmd_info(a):
    for s in _load(a.files, a):
        print(s.summary())
        st = detect_steps(s) if s.ykind not in ("counts",) else []
        for lam, j, r, d in st:
            print(f"DEGRAU      : {lam:.1f} nm (salto {j:+.4g}, {r:.0f}× ruído) - {d}")
        print("-" * 60)


def cmd_plot(a):
    specs = [_pre(s, a) for s in _load(a.files, a)]
    peaks = None
    if a.label:
        peaks = []
        for s in specs:
            t = find_bands(s, region=tuple(a.xlim) if a.xlim else None, shoulders=a.shoulders)
            t = t[~t["flags"].str.contains("ruído|borda|troca", regex=True)]
            if a.top:
                t = t.sort_values("prominence", ascending=False).head(a.top).sort_values("wavelength")
            peaks.append(t)
        if len(specs) > 3 and not a.label_all:     # séries: rótulos só nos extremos, senão empilham
            peaks = [p if i in (0, len(specs) - 1) else None for i, p in enumerate(peaks)]
    hl = [(float(h.split(":")[0]), float(h.split(":")[1]), h.split(":")[2] if h.count(":") > 1 else "")
          for h in (a.highlight or [])]
    offset = None if a.offset is None else ("auto" if a.offset == "auto" else float(a.offset))
    vl = None
    if a.isosbestic:
        iso = isosbestic(specs, region=tuple(a.xlim) if a.xlim else None)
        vl = [(r.wavelength, f"{r.wavelength:.0f} nm") for _, r in iso.iterrows()]
        print("pontos isosbésticos:", ", ".join(f"{r.wavelength} nm" for _, r in iso.iterrows()) or "nenhum")
    plot_spectra(specs, vlines=vl, mode=a.mode, offset=offset, peaks=peaks, xlim=tuple(a.xlim) if a.xlim else None,
                 labels=a.names or [_stem(s) for s in specs], series_values=a.series, series_label=a.series_label,
                 highlight=hl, ev_axis=a.ev_axis, lang=a.lang, normalized=bool(a.normalize), title=a.title,
                 figsize=tuple(a.size), out=a.output, cmap=a.cmap)
    print(f"gravado: {a.output}")


def cmd_peaks(a):
    pd.set_option("display.width", 200)
    for s in _load([a.file], a):
        s = _pre(s, a)
        t = find_bands(s, region=tuple(a.region) if a.region else None, top_n=a.top)
        print(f"# {_stem(s)}  ({t.attrs.get('ykind', '')}; ruído ~{t.attrs.get('noise', 0):.2g})")
        print(t.to_string(index=False))
        if a.output:
            t.to_csv(a.output, index=False, sep=";" if a.decimal == "," else ",", decimal=a.decimal)
            print(f"gravado: {a.output}")


def _measure(specs, a):
    where = a.wl if a.wl is not None else tuple(a.window)
    base = None
    if a.base:
        base = a.base[0] if len(a.base) == 1 else tuple(a.base)
    vals = []
    for s in specs:
        m = band_metrics(s, where, base=base)
        vals.append(m["value"])
    return np.array(vals), where


def _read_table(path):
    """Tabela de 2 colunas (x, y) com cabeçalho opcional; separador e vírgula decimal detectados."""
    from uvvis_io import _decode_bytes, _has_decimal_comma, _try_parse
    with open(path, "rb") as fh:
        lines = _decode_bytes(fh.read()).splitlines()
    best = None
    for delim in (";", "\t", ",", "ws"):
        for dec in (".", ","):
            if delim == "," and dec == ",":
                continue
            if dec == "," and not _has_decimal_comma(lines, delim):
                continue
            rows = [r for r in _try_parse(lines, delim, dec) if r is not None and len(r) >= 2]
            if best is None or len(rows) > len(best[0]):
                best = (rows, delim)
    if not best or len(best[0]) < 3:
        raise ValueError(f"não consegui ler 2 colunas numéricas de {path}")
    arr = np.array([r[:2] for r in best[0]], float)
    first = lines[0].split("\t" if best[1] == "\t" else (best[1] if best[1] != "ws" else None))
    cols = [c.strip().strip('"') for c in first[:2]] if not _try_parse([lines[0]], best[1], ".")[0] else ["x", "y"]
    return arr[:, 0], arr[:, 1], cols


def cmd_calib(a):
    if a.table:
        conc, A, cols = _read_table(a.table)
        where = cols[1]
    else:
        specs = [_pre(s, a) for s in _load(a.files, a)]
        if len(a.conc) != len(specs):
            sys.exit(f"--conc tem {len(a.conc)} valores mas há {len(specs)} espectros")
        conc = np.array(a.conc, float)
        A, where = _measure(specs, a)
        if a.wl is None:
            print("Aviso: usando o MÁXIMO na janela; para quantificação prefira um λ fixo (--wl).")
    cal = calibrate(conc, A, through_origin=a.origin, path_cm=a.path, conc_unit=a.unit or "")
    print(f"Medida: {where}" + (f" (linha de base {a.base})" if a.base else ""))
    print(pd.DataFrame({"conc": conc, "A": np.round(A, 5), "residuo": np.round(cal['residuals'], 5)}).to_string(index=False))
    sign = "+" if cal["intercept"] >= 0 else "-"
    print(f"\nA = {cal['slope']:.5g}·c {sign} {abs(cal['intercept']):.4g}   (b = {cal['slope']:.5g} ± {cal['se_slope']:.2g})")
    print(f"R² = {cal['r2']:.5f}  R²aj = {cal['r2_adj']:.5f}  s(y/x) = {cal['syx']:.3g}  n = {cal['n']}")
    print(f"LOD = {cal['lod']:.4g}  LOQ = {cal['loq']:.4g} {cal['conc_unit']}")
    if "epsilon" in cal:
        print(f"ε = {cal['epsilon']:.5g} {cal['epsilon_unit']}")
    if "mandel_p" in cal:
        print(f"Teste de Mandel p = {cal['mandel_p']:.3f}   intercepto p = {cal.get('intercept_p', float('nan')):.3f}")
    for w in cal["warnings"]:
        print("AVISO:", w)
    unk = []
    if a.unknown_A:
        unk.append(("A informada", predict(cal, a.unknown_A, dilution=a.dilution)))
    if a.unknown:
        for s in [_pre(s, a) for s in _load(a.unknown, a)]:
            v, _ = _measure([s], a)
            unk.append((_stem(s), predict(cal, v, dilution=a.dilution)))
    for name, u in unk:
        f = ("  [" + "; ".join(u["flags"]) + "]") if u["flags"] else ""
        print(f"Amostra {name}: A = {u['A']:.4f} -> c = {u['conc']:.4g} ± {u['s_conc']:.2g} "
              f"(IC95% {u['ci95'][0]:.4g} - {u['ci95'][1]:.4g}) {cal['conc_unit']}{f}")
    if a.output:
        plot_calibration(cal, out=a.output, lang=a.lang, unknowns=[u for _, u in unk] or None)
        print(f"gravado: {a.output}")


def cmd_kinetics(a):
    if a.table:
        t, A, _ = _read_table(a.table)
        specs = None
    else:
        specs = [_pre(s, a) for s in _load(a.files, a)]
        if len(a.times) != len(specs):
            sys.exit(f"--times tem {len(a.times)} valores mas há {len(specs)} espectros")
        t = np.array(a.times, float)
        if a.wl is None and a.window is None:
            i0 = int(np.argmin(np.abs(t)))
            tb = find_bands(specs[i0], region=(400, 800))
            a.wl = float(tb.sort_values("prominence").iloc[-1].wavelength)
            print(f"λ de acompanhamento = {a.wl} nm (λmax mais proeminente em t=0, 400-800 nm)")
        A, where = _measure(specs, a)
    models = ["zero", "pfo", "pso"] if a.model == "all" else [a.model]
    if a.model == "all":
        print("Comparação de modelos (R²):")
        for m, r2, k in compare_kinetic_models(t, A):
            print(f"  {m:5s} R² = {r2:.4f}  k = {k:.4g}" if np.isfinite(r2) else f"  {m}: {k}")
    res = kinetics(t, A, model="pfo" if a.model == "all" else a.model)
    print(pd.DataFrame({"t": t, "A": np.round(A, 5)}).to_string(index=False))
    print(f"\nModelo {res['model']}: k = {res['k']:.4g} {a.time_unit}^-1 (pfo)  R² = {res['r2']:.4f}  "
          f"t½ = {res['half_life']:.3g} {a.time_unit}  degradação final = {res['degradation_pct']:.1f}%")
    if "dark_adsorption_pct" in res:
        print(f"Adsorção no escuro (até t=0): {res['dark_adsorption_pct']:.1f}%")
    if a.output:
        plot_kinetics(res, names=[a.name or "amostra"], out=a.output, lang=a.lang, time_unit=a.time_unit)
        print(f"gravado: {a.output}")
        if specs is not None:
            base, ext = os.path.splitext(a.output)
            plot_spectra(specs, mode="A", series_values=list(t), series_label=f"t ({a.time_unit})", lang=a.lang,
                         out=base + "_espectros" + ext)
            print(f"gravado: {base}_espectros{ext}")


def cmd_tauc(a):
    s = _pre(_load([a.file], a)[0], a)
    trans = ["direct", "indirect"] if a.transition == "both" else [a.transition]
    for tr in trans:
        r = tauc(s, tr, thickness_cm=a.thickness, energy_window=tuple(a.window) if a.window else None,
                 fit_window=tuple(a.fit) if a.fit else None, r2_min=a.r2)
        print(f"[{tr}] F = {r['source']}")
        print(f"  Eg (interseção y=0)        = {r['Eg']:.3f} eV  ({r['lambda_edge_nm']:.0f} nm)")
        print(f"  Eg (interseção c/ sub-gap) = {r['Eg_baseline']:.3f} eV")
        print(f"  ajuste {r['fit_range'][0]:.3f}-{r['fit_range'][1]:.3f} eV, {r['n_fit']} pts, R² = {r['r2']:.4f}, "
              f"variação com a janela ±{r['Eg_spread'] / 2:.3f} eV ({'auto' if r['auto'] else 'manual'})")
        for w in r["warnings"]:
            print("  AVISO:", w)
        try:
            u = urbach(s, Eg=r["Eg"], thickness_cm=a.thickness)
            note = ""
            if r.get("other_edges_eV") or u["Eu_meV"] > 300:
                note = "  <- duvidoso: absorção sub-gap não exponencial (dopante/defeito/espalhamento)"
            print(f"  Urbach E_U = {u['Eu_meV']:.0f} meV (R² = {u['r2']:.3f}, "
                  f"{u['fit_range'][0]:.2f}-{u['fit_range'][1]:.2f} eV){note}")
        except ValueError as e:
            print(f"  Urbach: {e}")
        d = derivative_edge(s, thickness_cm=a.thickness, Eg=r["Eg"])
        print(f"  conferência: inflexão dF/dE = {d['E_inflection']:.3f} eV")
        if a.output:
            base, ext = os.path.splitext(a.output)
            outp = a.output if len(trans) == 1 else f"{base}_{tr}{ext}"
            plot_tauc(r, out=outp, lang=a.lang, title=a.title)
            print(f"  gravado: {outp}")
    try:
        o = onset_wavelength(s, Eg=r["Eg"])
        print(f"Onset (tangente em A vs λ): λ = {o['lambda_onset_nm']:.0f} nm -> {o['Eg_onset']:.3f} eV")
    except Exception as e:  # noqa
        print(f"Onset: {e}")


def cmd_isosbestic(a):
    specs = [_pre(s, a) for s in _load(a.files, a)]
    iso = isosbestic(specs, region=tuple(a.region) if a.region else None)
    print(iso.to_string(index=False) if len(iso) else "Nenhum ponto isosbéstico claro.")


def cmd_wlcheck(a):
    s = _load([a.file], a)[0]
    df = wavelength_check(s, a.standard)
    print(df.to_string(index=False))
    print(f"\nDesvio médio = {df.attrs['mean_offset']:+.2f} nm  (s = {df.attrs['sd']:.2f}; {df.attrs['n_found']} bandas)")
    print("Valores de referência são nominais; use os do certificado do seu filtro.")


def cmd_convert(a):
    s = _load([a.file], a)[0]
    if a.to == "A":
        s = s.to_absorbance()
    elif a.to == "T":
        s = s.to_transmittance(True)
    elif a.to == "KM":
        s = s.to_km()
    save_csv(s, a.output, decimal=a.decimal, energy=a.ev)
    print(f"gravado: {a.output} ({s.npoints} pontos, {s.ykind})")


def cmd_report(a):
    os.makedirs(a.dir, exist_ok=True)
    for s0 in _load([a.file], a):
        s = _pre(s0, a)
        name = "".join(c if c.isalnum() or c in "-_." else "_" for c in _stem(s))[:60] or "espectro"
        base = os.path.join(a.dir, name)
        t = find_bands(s)
        lab = t[~t["flags"].str.contains("ruído|borda|troca", regex=True)]
        lab = lab.sort_values("prominence", ascending=False).head(a.top).sort_values("wavelength")
        plot_spectra(s, peaks=lab, labels=[_stem(s)], ev_axis=True, lang=a.lang, out=base + "_lambda_max.png")
        t.to_csv(base + "_bandas.csv", index=False)
        steps = detect_steps(s) if s.ykind != "counts" else []
        with open(base + "_relatorio.md", "w", encoding="utf-8") as fh:
            fh.write(f"# UV-Vis: {_stem(s)}\n\n```\n{s0.summary()}\n```\n\n")
            if steps:
                fh.write("**Degraus instrumentais:** " + "; ".join(f"{l:.0f} nm ({d})" for l, _, _, d in steps) + "\n\n")
            fh.write("| λ (nm) | E (eV) | tipo | valor | FWHM (nm) | flags |\n|---|---|---|---|---|---|\n")
            for _, r in t.iterrows():
                fh.write(f"| {r.wavelength} | {r.energy_eV} | {r.kind} | {r.value} | {r.fwhm_nm} | {r['flags']} |\n")
        print(f"gravados: {base}_lambda_max.png, {base}_bandas.csv, {base}_relatorio.md")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--ykind", choices=["A", "T", "%T", "R", "%R", "KM", "eps", "alpha", "counts"])
        sp.add_argument("--index", type=int, help="espectro a usar quando o arquivo tem vários")
        sp.add_argument("--offset-correct", nargs=2, type=float, metavar=("NM1", "NM2"),
                        help="subtrai a média nessa região sem absorção")
        sp.add_argument("--scatter", nargs=2, type=float, metavar=("NM1", "NM2"),
                        help="remove espalhamento λ^-n ajustado nessa região sem absorção")
        sp.add_argument("--smooth", type=float, help="janela Savitzky-Golay em nm")
        sp.add_argument("--normalize", choices=["max", "at", "minmax", "area"])
        sp.add_argument("--norm-wl", type=float)
        sp.add_argument("--lang", default="en", choices=["en", "pt"])

    def meas(sp):
        sp.add_argument("--wl", type=float, help="λ fixo de medida (nm)")
        sp.add_argument("--window", nargs=2, type=float, help="máximo dentro da janela (nm)")
        sp.add_argument("--base", nargs="+", type=float, help="λ de linha de base (1 valor) ou reta entre 2 λ")

    s = sub.add_parser("info"); s.add_argument("files", nargs="+"); common(s); s.set_defaults(f=cmd_info)
    s = sub.add_parser("plot"); s.add_argument("files", nargs="+"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--mode", default="auto", choices=["auto", "A", "T", "R", "KM", "raw"])
    s.add_argument("--offset"); s.add_argument("--label", action="store_true"); s.add_argument("--top", type=int)
    s.add_argument("--shoulders", action="store_true")
    s.add_argument("--series", nargs="+", type=float); s.add_argument("--series-label")
    s.add_argument("--cmap", default="viridis"); s.add_argument("--names", nargs="+")
    s.add_argument("--highlight", nargs="*"); s.add_argument("--xlim", nargs=2, type=float)
    s.add_argument("--ev-axis", action="store_true"); s.add_argument("--title")
    s.add_argument("--isosbestic", action="store_true", help="detecta e marca pontos isosbésticos")
    s.add_argument("--label-all", action="store_true", help="rotula todos os espectros de uma série (padrão: só 1º e último)")
    s.add_argument("--size", nargs=2, type=float, default=[6.4, 4.2]); common(s); s.set_defaults(f=cmd_plot)
    s = sub.add_parser("peaks"); s.add_argument("file"); s.add_argument("-o", "--output")
    s.add_argument("--region", nargs=2, type=float); s.add_argument("--top", type=int)
    s.add_argument("--decimal", default="."); common(s); s.set_defaults(f=cmd_peaks)
    s = sub.add_parser("calib"); s.add_argument("files", nargs="*"); s.add_argument("--table")
    s.add_argument("--conc", nargs="+", type=float); meas(s)
    s.add_argument("--origin", action="store_true"); s.add_argument("--path", type=float)
    s.add_argument("--unit", default=""); s.add_argument("--unknown", nargs="+")
    s.add_argument("--unknown-A", nargs="+", type=float); s.add_argument("--dilution", type=float, default=1.0)
    s.add_argument("-o", "--output"); common(s); s.set_defaults(f=cmd_calib)
    s = sub.add_parser("kinetics"); s.add_argument("files", nargs="*"); s.add_argument("--table")
    s.add_argument("--times", nargs="+", type=float); meas(s)
    s.add_argument("--model", default="pfo", choices=["pfo", "pso", "zero", "all"])
    s.add_argument("--time-unit", default="min"); s.add_argument("--name"); s.add_argument("-o", "--output")
    common(s); s.set_defaults(f=cmd_kinetics)
    s = sub.add_parser("tauc"); s.add_argument("file")
    s.add_argument("--transition", default="direct",
                   choices=["direct", "indirect", "direct_forbidden", "indirect_forbidden", "both"])
    s.add_argument("--window", nargs=2, type=float, help="faixa de energia (eV) onde procurar a borda")
    s.add_argument("--fit", nargs=2, type=float, help="força a região de ajuste (eV)")
    s.add_argument("--thickness", type=float, help="espessura do filme em cm (α = 2.303A/d)")
    s.add_argument("--r2", type=float, default=0.995); s.add_argument("--title"); s.add_argument("-o", "--output")
    common(s); s.set_defaults(f=cmd_tauc)
    s = sub.add_parser("isosbestic"); s.add_argument("files", nargs="+"); s.add_argument("--region", nargs=2, type=float)
    common(s); s.set_defaults(f=cmd_isosbestic)
    s = sub.add_parser("wlcheck"); s.add_argument("file")
    s.add_argument("--standard", default="holmium_glass", choices=["holmium_glass", "holmium_solution", "didymium_glass"])
    common(s); s.set_defaults(f=cmd_wlcheck)
    s = sub.add_parser("convert"); s.add_argument("file"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--to", choices=["A", "T", "KM"]); s.add_argument("--decimal", default=".")
    s.add_argument("--ev", action="store_true"); common(s); s.set_defaults(f=cmd_convert)
    s = sub.add_parser("report"); s.add_argument("file"); s.add_argument("-d", "--dir", default="uvvis_report")
    s.add_argument("--top", type=int, default=8); common(s); s.set_defaults(f=cmd_report)
    a = p.parse_args(argv)
    try:
        a.f(a)
    except (ValueError, FileNotFoundError) as e:
        sys.exit(f"ERRO: {e}")


if __name__ == "__main__":
    main()
