#!/usr/bin/env python3
"""
ftir.py - linha de comando para a skill de FT-IR

  python ftir.py info     arquivos...                         metadados, faixa, tipo de ordenada
  python ftir.py convert  arquivo -o saida.csv [--to A|T] [--decimal ,] [--jcamp]
  python ftir.py peaks    arquivo [-o bandas.csv] [--region 4000 400] [--prom 0.01]
                          [--min-sep 6] [--top 25] [--assign] [--baseline als]
  python ftir.py plot     arquivos... -o fig.png [--mode A|T] [--offset auto|0.5]
                          [--label] [--top 15] [--assign-labels] [--highlight 1760:1680:C=O]
                          [--xlim 4000 400] [--baseline rubberband] [--normalize max]
                          [--smooth 9] [--lang pt] [--no-co2]
  python ftir.py report   arquivo -d pasta_saida [--lang pt] [--top 20] [...]
                          (figura rotulada + tabela de bandas com candidatos + resumo .md)
  python ftir.py compare  arquivos... --band C=O:1760:1680 --band ref:1490:1440 [-o tabela.csv]
                          (mesmas janelas para todas as amostras; altura/área com linha de base local)

Formatos de entrada: .sp (PerkinElmer), .jdx/.dx/.jcamp, .csv/.txt/.dpt/.prn/.xy/.asc.
Pode passar curingas (ex.: "dados/*.sp"). Execute de qualquer pasta.
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

from ftir_io import load_all, save_csv, save_jcamp  # noqa: E402
from ftir_process import baseline, co2_flatten, normalize, smooth  # noqa: E402
from ftir_peaks import compare_series, find_bands  # noqa: E402
from ftir_assign import assign, screen_groups  # noqa: E402
from ftir_plot import plot_spectra  # noqa: E402


def _expand(paths):
    out = []
    for p in paths:
        g = sorted(glob.glob(p))
        out += g if g else [p]
    return out


def _load(paths, index=None):
    specs = []
    for p in _expand(paths):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            got = load_all(p)
        for ww in w:
            print(f"[aviso] {os.path.basename(p)}: {ww.message}", file=sys.stderr)
        specs += got if index is None else [got[index]]
    return specs


def _pre(s, a):
    """Pipeline comum: absorbância -> CO2 -> linha de base -> suavização -> normalização."""
    if getattr(a, "baseline", None) or getattr(a, "normalize", None) or getattr(a, "smooth", None):
        if s.ykind in ("T", "%T", "R", "%R"):
            s = s.to_absorbance()
    if getattr(a, "no_co2", False):
        s = co2_flatten(s)
    if getattr(a, "baseline", None):
        s, _ = baseline(s, a.baseline)
    if getattr(a, "smooth", None):
        s = smooth(s, a.smooth)
    if getattr(a, "normalize", None):
        s = normalize(s, a.normalize, band=tuple(a.norm_band) if a.norm_band else None)
    return s


def _stem(s):
    """Nome curto para legenda: nome do arquivo (+ coluna/bloco quando o arquivo tem vários)."""
    base = os.path.splitext(os.path.basename(s.source))[0]
    if s.meta.get("column") is not None and s.meta.get("y_label"):
        return s.meta["y_label"]
    if "(largo" in s.meta.get("format", ""):
        return s.title
    return base


def cmd_info(a):
    for s in _load(a.files):
        print(s.summary())
        print("-" * 60)


def cmd_convert(a):
    s = _load([a.file], a.index or 0)[0]
    if a.to == "A":
        s = s.to_absorbance()
    elif a.to == "T":
        s = s.to_transmittance(percent=True)
    if a.jcamp or a.output.lower().endswith((".jdx", ".dx")):
        save_jcamp(s, a.output)
    else:
        save_csv(s, a.output, decimal=a.decimal)
    print(f"gravado: {a.output}  ({s.npoints} pontos, ykind={s.ykind})")


def _bands(s, a):
    region = tuple(a.region) if a.region else None
    tab = find_bands(s, region=region, prominence=a.prom, min_sep_cm=a.min_sep,
                     shoulders=not a.no_shoulders)
    if a.assign:
        tab = assign(tab)   # sobre a lista completa: a checagem de bandas parceiras precisa dela
    if a.top:
        attrs = dict(tab.attrs)
        tab = tab.sort_values("prominence", ascending=False).head(a.top).sort_values("position", ascending=False)
        tab.attrs.update(attrs)
    return tab


def cmd_peaks(a):
    for s in _load([a.file], a.index):
        s = _pre(s, a)
        tab = _bands(s, a)
        pd.set_option("display.width", 200, "display.max_colwidth", 70, "display.max_rows", 500)
        print(f"# {s.title}  (ruído~{tab.attrs.get('noise', float('nan')):.2g}, limiar de proeminência "
              f"{tab.attrs.get('prominence_threshold', float('nan')):.3g})")
        if tab.attrs.get("vapor_warning"):
            print("AVISO: " + tab.attrs["vapor_warning"])
        cols = [c for c in ["position", "kind", "absorbance", "rel_intensity", "intensity", "shape", "fwhm", "snr",
                            "flags", "assignment"] if c in tab.columns]
        print(tab[cols].to_string(index=False))
        if a.assign:
            print("\nTriagem:\n" + screen_groups(tab))
        if a.output:
            tab.to_csv(a.output, index=False, sep=";" if a.decimal == "," else ",", decimal=a.decimal)
            print(f"\ngravado: {a.output}")


def cmd_plot(a):
    specs = [_pre(s, a) for s in _load(a.files, a.index)]
    peaks = None
    if a.label:
        peaks = []
        for s in specs:
            t = find_bands(s, region=tuple(a.region) if a.region else tuple(a.xlim), prominence=a.prom,
                           min_sep_cm=a.min_sep, shoulders=a.label_shoulders)
            if a.assign_labels:
                t = assign(t)
            t = t[~t["flags"].str.contains("ruído|CO2|H2O", regex=True)]
            if a.top:
                t = t.sort_values("prominence", ascending=False).head(a.top).sort_values("position", ascending=False)
            peaks.append(t)
    hl = []
    for h in a.highlight or []:
        parts = h.split(":")
        hl.append((float(parts[0]), float(parts[1]), parts[2] if len(parts) > 2 else ""))
    offset = None
    if a.offset is not None:
        offset = "auto" if a.offset == "auto" else float(a.offset)
    label_col = "assignment" if a.assign_labels else None
    if a.label_map:
        # rótulos revisados: '1715=ν(C=O)' -> banda detectada mais próxima (±10 cm-1) recebe o texto;
        # sem banda detectada perto, o rótulo é colocado na posição indicada
        mp = []
        for item in a.label_map:
            k, _, v = item.partition("=")
            mp.append((float(k), v))
        peaks = peaks or [None] * len(specs)
        new = []
        for s, t in zip(specs, peaks):
            full = find_bands(s, region=tuple(a.xlim), prominence=a.prom, min_sep_cm=a.min_sep, shoulders=True)
            rows = []
            for pos, txt in mp:
                near = full[(full.position - pos).abs() <= 10] if len(full) else full
                if len(near):
                    r = near.iloc[(near.position - pos).abs().argmin()].copy()
                else:
                    r = pd.Series({"position": int(round(pos)), "kind": "peak"})
                r["user_label"] = txt
                rows.append(r)
            new.append(pd.DataFrame(rows))
        peaks, label_col = new, "user_label"
    names = a.names or [_stem(s) for s in specs]
    plot_spectra(specs, mode=a.mode, offset=offset, peaks=peaks, highlight=hl, xlim=tuple(a.xlim),
                 lang=a.lang, labels=names, label_assign=label_col,
                 normalized=bool(a.normalize), title=a.title, figsize=tuple(a.size), out=a.output)
    print(f"gravado: {a.output}")


def cmd_report(a):
    os.makedirs(a.dir, exist_ok=True)
    for s0 in _load([a.file], a.index):
        s = _pre(s0, a)
        base = os.path.join(a.dir, "".join(c if c.isalnum() or c in "-_." else "_" for c in _stem(s))[:60] or "spectrum")
        full = assign(find_bands(s, prominence=a.prom, min_sep_cm=a.min_sep))
        pk = full[(full.kind == "peak") & ~full["flags"].str.contains("ruído|CO2|H2O", regex=True)]
        lab = pk.sort_values("prominence", ascending=False).head(a.top).sort_values("position", ascending=False)
        plot_spectra(s, mode=a.mode, peaks=lab, labels=[_stem(s)], xlim=tuple(a.xlim), lang=a.lang, out=base + "_labeled.png",
                     normalized=bool(a.normalize))
        full.to_csv(base + "_bands.csv", index=False)
        with open(base + "_report.md", "w", encoding="utf-8") as fh:
            fh.write(f"# FT-IR: {s.title}\n\n```\n{s0.summary()}\n```\n\n")
            fh.write(f"Processamento: {', '.join(f'{k}={s.meta[k]}' for k in ('converted_from', 'baseline', 'smoothing', 'normalization', 'co2_flattened') if k in s.meta) or 'nenhum'}\n\n")
            fh.write(f"Ruído estimado {full.attrs.get('noise', 0):.2g}; limiar de proeminência "
                     f"{full.attrs.get('prominence_threshold', 0):.3g}; intervalo {full.attrs.get('step', 0):.3g} cm-1.\n\n")
            if full.attrs.get("vapor_warning"):
                fh.write("**Aviso:** " + full.attrs["vapor_warning"] + "\n\n")
            fh.write("## Triagem de grupos (automática - conferir no espectro)\n\n```\n" + screen_groups(full) + "\n```\n\n")
            fh.write("## Bandas e candidatos de atribuição\n\n| cm⁻¹ | int. | forma | flags | candidatos |\n|---|---|---|---|---|\n")
            for _, r in full.iterrows():
                fh.write(f"| {r.position} | {r.intensity} | {r['shape']} | {r['flags']} | {r.candidates} |\n")
        print(f"gravados: {base}_labeled.png, {base}_bands.csv, {base}_report.md")


def cmd_compare(a):
    specs = [_pre(s, a) for s in _load(a.files, a.index)]
    win = {}
    for b in a.band:
        name, hi, lo = b.split(":")
        win[name] = (float(hi), float(lo))
    df = compare_series(specs, win, use=a.use)
    print(df.to_string(index=False))
    if a.output:
        df.to_csv(a.output, index=False)
        print(f"gravado: {a.output}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common_proc(sp):
        sp.add_argument("--baseline", choices=["als", "rubberband", "poly"])
        sp.add_argument("--smooth", type=float, help="janela Savitzky-Golay em cm-1")
        sp.add_argument("--normalize", choices=["max", "minmax", "area", "band", "snv", "vector"])
        sp.add_argument("--norm-band", nargs=2, type=float, help="faixa hi lo para --normalize band")
        sp.add_argument("--no-co2", action="store_true", help="interpola sobre 2400-2280 cm-1")
        sp.add_argument("--index", type=int, help="bloco/coluna a usar quando o arquivo tem vários espectros")

    def common_peaks(sp):
        sp.add_argument("--prom", type=float, help="proeminência mínima (absorbância)")
        sp.add_argument("--min-sep", type=float, default=6.0)
        sp.add_argument("--region", nargs=2, type=float)
        sp.add_argument("--top", type=int)
        sp.add_argument("--no-shoulders", action="store_true")

    s = sub.add_parser("info"); s.add_argument("files", nargs="+"); s.set_defaults(f=cmd_info)
    s = sub.add_parser("convert"); s.add_argument("file"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--to", choices=["A", "T"]); s.add_argument("--decimal", default=".")
    s.add_argument("--jcamp", action="store_true"); s.add_argument("--index", type=int); s.set_defaults(f=cmd_convert)
    s = sub.add_parser("peaks"); s.add_argument("file"); s.add_argument("-o", "--output")
    s.add_argument("--assign", action="store_true"); s.add_argument("--decimal", default=".")
    common_peaks(s); common_proc(s); s.set_defaults(f=cmd_peaks)
    s = sub.add_parser("plot"); s.add_argument("files", nargs="+"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--mode", default="auto", choices=["auto", "A", "T", "raw"])
    s.add_argument("--offset"); s.add_argument("--label", action="store_true")
    s.add_argument("--assign-labels", action="store_true"); s.add_argument("--highlight", nargs="*")
    s.add_argument("--xlim", nargs=2, type=float, default=[4000, 400]); s.add_argument("--lang", default="en")
    s.add_argument("--title"); s.add_argument("--size", nargs=2, type=float, default=[7.0, 4.2])
    s.add_argument("--names", nargs="+", help="nomes dos traços na legenda (padrão: nome do arquivo)")
    s.add_argument("--label-shoulders", action="store_true", help="também rotula ombros (2ª derivada)")
    s.add_argument("--label-map", nargs="+", help="rótulos revisados, ex.: 1715='ν(C=O)' 1240='ν(C–O)' (só essas bandas)")
    common_peaks(s); common_proc(s); s.set_defaults(f=cmd_plot)
    s = sub.add_parser("report"); s.add_argument("file"); s.add_argument("-d", "--dir", default="ftir_report")
    s.add_argument("--mode", default="auto", choices=["auto", "A", "T"]); s.add_argument("--lang", default="en")
    s.add_argument("--xlim", nargs=2, type=float, default=[4000, 400])
    s.add_argument("--prom", type=float); s.add_argument("--min-sep", type=float, default=6.0)
    s.add_argument("--top", type=int, default=20); common_proc(s); s.set_defaults(f=cmd_report)
    s = sub.add_parser("compare"); s.add_argument("files", nargs="+"); s.add_argument("--band", action="append", required=True)
    s.add_argument("--use", default="height", choices=["height", "area"]); s.add_argument("-o", "--output")
    common_proc(s); s.set_defaults(f=cmd_compare)
    a = p.parse_args(argv)
    a.f(a)


if __name__ == "__main__":
    main()
