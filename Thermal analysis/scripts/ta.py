#!/usr/bin/env python3
"""
ta.py - linha de comando da skill de análise térmica (TG/DTG, DSC, DTA, STA)

  info      arquivos...                         formato, tipo, segmentos do programa, massa, direção endo
  tga       arquivo(s) [-o fig.png] [--csv etapas.csv] [--basis 150] [--min-loss 0.5]
            [--theory "CaC2O4·H2O:H2O,CO,CO2"] [--blank branco.txt] [--x T|time] [--lang pt]
  dsc       arquivo [--segment heat2 ...] [--endo up|down] [--tg 60:105] [--peak 215:280 ...]
            [--polymer PET | --dh0 140] [--w 0.7] [--baseline linear|sigmoidal] [--display exo_up|endo_up]
            [-o fig.png] [--csv eventos.csv]       (sem --tg/--peak: detecção e integração automáticas)
  sta       arquivo [-o fig.png] [--endo ...]   TG/DTG + DSC/DTA; caráter endo/exo de cada etapa
  kinetics  arquivos... --window 250:450 [--signal mass|hf] [-o fig.png] [--csv ea.csv]
  oit       arquivo [--t-switch 13.8] [-o fig.png]
  cp        --sample a --blank b --sapphire c [--m-sample mg] [--m-sapphire mg] [--range 50:300] [-o cp.csv]
  stoich    "CuSO4·5H2O" --release 2H2O 2H2O H2O      perdas teóricas por etapa
  hydrate   --loss 14.4 --anhydrous CuSO4             nº de águas a partir da perda medida
  plot      arquivos... -o fig.png [--kind tga|dsc] [--segment heat2] [--names ...]   sobreposição
  convert   arquivo -o saida.csv [--decimal ,]
  report    arquivo -d pasta [--lang pt]             análise completa conforme o tipo (TGA/DSC/STA)

Leitura: NETZSCH .ngb-* (binário; pip install pyngb) e ASCII; TA Universal Analysis/TRIOS (txt,
xls); Mettler STARe; PerkinElmer Pyris; Setaram; CSV/TXT/Excel genéricos. Arquivo sem cabeçalho:
--columns time:min,T:°C,mass:mg. NETZSCH com correção separada: --correction arquivo.ngb-bs3.
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

from ta_io import load_all, save_csv  # noqa: E402
from ta_process import endo_direction, heating_rate, segments, select, subtract_blank  # noqa: E402
from ta_tga import (char_temps, compare_with_theory, find_steps, hydration_number, molar_mass,  # noqa: E402
                    residues, theoretical_losses)
from ta_dsc import (DH0, cp_sapphire, crystallinity, find_events, glass_transition, integrate_peak,  # noqa: E402
                    oit)
from ta_kinetics import isoconversional, kissinger  # noqa: E402
from ta_plot import plot_dsc, plot_kinetics, plot_oit, plot_sta, plot_tga  # noqa: E402

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 40)


def _expand(paths):
    out = []
    for p in paths:
        g = sorted(glob.glob(p))
        out += g if g else [p]
    return out


def _load(paths, a, one=False):
    runs = []
    for p in _expand(paths):
        kw = {}
        if getattr(a, "columns", None):
            kw["columns"] = a.columns
        if getattr(a, "correction", None):
            kw["baseline"] = a.correction
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            got = load_all(p, **kw)
        for ww in w:
            print(f"[aviso] {os.path.basename(p)}: {ww.message}", file=sys.stderr)
        if getattr(a, "index", None) is not None:
            got = [got[a.index]]
        for r in got:
            if getattr(a, "mass", None):
                r.meta["sample_mass_mg"] = a.mass
                from ta_io import _finish
                r.data.pop("mass_pct", None) if r.has("mass") else None
                _finish(r)
        runs += got
    if one:
        if len(runs) > 1:
            print(f"[aviso] {len(runs)} curvas; usando a primeira (--index N para outra).", file=sys.stderr)
        return runs[0]
    return runs


def _win(s):
    a, b = s.replace(",", ".").split(":")
    return float(a), float(b)


def _stem(r):
    return os.path.splitext(os.path.basename(r.source))[0]


def _print_segments(r):
    segs = segments(r)
    print("Segmentos do programa:")
    for s in segs:
        dur = (s["t1"] - s["t0"]) if s["t0"] is not None else float("nan")
        print(f"  {s['n']:>2}  {s['key']:<7} {s['T0']:7.1f} → {s['T1']:7.1f} °C   "
              f"{(s['rate'] if s['rate'] is not None else float('nan')):7.2f} K/min   {dur:6.1f} min")


# ============================================================================ comandos
def cmd_info(a):
    for r in _load(a.files, a):
        print(r.summary())
        _print_segments(r)
        if r.has("hf"):
            d, src = endo_direction(r, a.endo)
            print(f"Direção endo: {'para cima' if d == 'up' else 'para baixo'}  ({src})")
        print("-" * 70)


def _tga_one(r, a, out_fig=None, out_csv=None, quiet=False):
    if getattr(a, "blank", None):
        b = _load([a.blank], a, one=True)
        r = subtract_blank(r, b, "mass", on="T" if not a.blank_time else "time")
    basis = a.basis
    st = find_steps(r, min_loss=a.min_loss, basis=basis)
    tc = char_temps(r, basis=basis)
    rs = residues(r, basis=basis)
    if not quiet:
        beta = heating_rate(r)
        print(f"# {r.title}  ({r.kind}; β ≈ {beta} K/min; massa inicial {r.sample_mass or float('nan'):.4g} mg"
              + (f"; base: massa a {basis} °C" if basis else "") + ")")
        for w in r.meta.get("warnings", []):
            print("AVISO:", w)
        cols = [c for c in ("etapa", "tipo", "segmento", "T_inicio", "T_onset", "T_pico", "T_endset", "T_fim",
                            "perda_pct", "perda_mg", "DTG_max_pct_C", "DTG_max_pct_min", "sobreposta") if c in st]
        print(st[cols].to_string(index=False) if len(st) else "(nenhuma etapa ≥ limite)")
        print("Temperaturas características:", ", ".join(f"{k} = {v}" for k, v in tc.items()))
        print("Massa residual:", ", ".join(f"{k} = {v}" for k, v in rs.items()))
        if len(st) and st["sobreposta"].any():
            print("AVISO: etapas sobrepostas (DTG não volta à linha de base entre picos) - Δm por etapa é aproximado.")
    if getattr(a, "theory", None):
        reac, _, rel = a.theory.partition(":")
        th = theoretical_losses(reac.strip(), [x.strip() for x in rel.split(",") if x.strip()])
        print(f"\nTeórico para {reac.strip()} (M = {th.attrs['M']:.2f} g/mol):")
        print(compare_with_theory(st, th).to_string(index=False))
    if out_csv:
        st.to_csv(out_csv, index=False, sep=";" if a.decimal == "," else ",", decimal=a.decimal)
        print(f"gravado: {out_csv}")
    if out_fig:
        plot_tga(r, steps=st, x=a.x, basis=basis, lang=a.lang, out=out_fig)
        print(f"gravado: {out_fig}")
    return r, st, tc, rs


def cmd_tga(a):
    runs = _load(a.files, a)
    if len(runs) > 1:
        rows = []
        for r in runs:
            _, st, tc, rs = _tga_one(r, a, quiet=False)
            rows.append({"amostra": r.title, **tc, **{k: v for k, v in rs.items()}})
            print()
        print(pd.DataFrame(rows).to_string(index=False))
        if a.output:
            plot_tga(runs, x=a.x, basis=a.basis, lang=a.lang, out=a.output, labels=a.names)
            print(f"gravado: {a.output}")
        return
    _tga_one(runs[0], a, a.output, a.csv)


def _auto_dsc(seg, a):
    """Eventos automáticos -> Tg e picos integrados. Janelas recortadas para não invadir eventos
    vizinhos; picos < 8% do maior evento do segmento são só listados (não quantificados)."""
    ev = find_events(seg, endo=a.endo)
    res = []
    if not len(ev):
        return ev, res
    hmax = float(ev.loc[ev["evento"] == "pico", "altura"].max()) if (ev["evento"] == "pico").any() else float(ev["altura"].max())
    Tmin, Tmax = float(np.nanmin(seg.T)) + 2, float(np.nanmax(seg.T)) - 2
    is_tg_ = lambda e: str(e["evento"]).startswith("degrau") or e["evento"] == "pico+degrau"
    evs = [e for e in ev.to_dict("records") if (is_tg_(e) or e["altura"] >= 0.08 * hmax)
           and not str(e["evento"]).startswith("borda")]
    for i, e in enumerate(evs):
        # limites: ponto médio até os centros dos eventos vizinhos quantificados
        prev_end = (evs[i - 1]["T_pico"] + e["T_pico"]) / 2 if i > 0 else Tmin
        next_ini = (evs[i + 1]["T_pico"] + e["T_pico"]) / 2 if i + 1 < len(evs) else Tmax
        is_tg = is_tg_(e)
        try:
            if is_tg:
                lo, hi = max(e["T_pico"] - 25, prev_end, Tmin), min(e["T_pico"] + 25, next_ini, Tmax)
                res.append(glass_transition(seg, (lo, hi), endo=a.endo))
            else:
                width = e["T_fim"] - e["T_ini"]
                lo = max(e["T_ini"] - 0.6 * width - 5, prev_end, Tmin)
                hi = min(e["T_fim"] + 0.6 * width + 5, next_ini, Tmax)
                res.append(integrate_peak(seg, (lo, hi), endo=a.endo, baseline=a.baseline))
        except ValueError as ex:
            print(f"  (evento em {e['T_pico']} °C não quantificado: {ex})")
    return ev, res


def _fmt_res(res):
    if "dH" in res:
        return (f"pico {res['tipo']:<4}  janela {res['janela'][0]:.0f}–{res['janela'][1]:.0f} °C  "
                f"T_onset {res['T_onset']}  T_pico {res['T_pico']}  T_endset {res['T_endset']}  "
                f"ΔH = {res['dH']:.2f} {res['unidade_dH']}  (linha de base {res['linha_base']})")
    return (f"Tg           janela {res['janela'][0]:.0f}–{res['janela'][1]:.0f} °C  onset {res['T_onset']}  "
            f"meia altura {res['T_mid']}  inflexão {res['T_inflexao']}  endset {res['T_endset']}  "
            f"ΔCp = {res['dCp']} J/(g·K)" + ("  [relaxação entálpica]" if res.get("relaxacao") else ""))


def cmd_dsc(a):
    r = _load([a.file], a, one=True)
    d, src = endo_direction(r, a.endo)
    a.endo = d
    print(f"# {r.title}  ({r.kind}; massa {r.sample_mass or float('nan'):.4g} mg)")
    print(f"Direção endo: {'↑' if d == 'up' else '↓'}  ({src})")
    for w in r.meta.get("warnings", []):
        print("AVISO:", w)
    _print_segments(r)
    allsegs = [s for s in segments(r) if s["type"] in ("heat", "cool") and abs(s["T1"] - s["T0"]) > 10]
    keys = a.segment or [s["key"] for s in allsegs]
    results, rows = {}, []
    manual = bool(a.tg or a.peak)
    for k in keys:
        seg = select(r, k)
        print(f"\n== {seg.meta['segment_info']['label']}  ({seg.meta['segment_info']['T0']:.0f} → "
              f"{seg.meta['segment_info']['T1']:.0f} °C, {seg.meta['segment_info']['rate']} K/min)")
        res = []
        if manual:
            for w in a.tg or []:
                res.append(glass_transition(seg, _win(w), endo=d))
            for w in a.peak or []:
                lo, hi = _win(w)
                if lo < np.nanmax(seg.T) and hi > np.nanmin(seg.T):
                    res.append(integrate_peak(seg, (lo, hi), endo=d, baseline=a.baseline))
        else:
            ev, res = _auto_dsc(seg, a)
            if len(ev):
                print("Eventos detectados (candidatos):")
                print(ev.drop(columns=["unidade"]).to_string(index=False))
        for x in res:
            print("  " + _fmt_res(x))
            for w in x.get("avisos", []):
                print("    AVISO:", w)
            rows.append({"segmento": k, **{kk: (str(v) if isinstance(v, tuple) else v) for kk, v in x.items()
                                              if not kk.startswith("_") and kk not in ("avisos", "endo_source")}})
        results[k] = res
        # cristalinidade
        if a.polymer or a.dh0:
            pk = [x for x in res if "dH" in x]
            melt = [x for x in pk if x["tipo"] == "endo"]
            cold = [x for x in pk if x["tipo"] == "exo" and (seg.meta["segment_info"]["type"] == "heat")]
            if melt and seg.meta["segment_info"]["type"] == "heat":
                mm = max(melt, key=lambda x: x["dH"])
                cc = sum(x["dH"] for x in cold if x["T_pico"] < mm["T_pico"])
                xc = crystallinity(mm["dH"], a.dh0 or a.polymer, dHcc=cc, w=a.w)
                print(f"  Cristalinidade: Xc = {xc['Xc']} %  (ΔHm {xc['dHm']:.2f} − ΔHcc {xc['dHcc']:.2f}) / "
                      f"(ΔH0 {xc['dH0']} × w {xc['w']})")
            elif seg.meta["segment_info"]["type"] == "cool":
                cr = [x for x in pk if x["tipo"] == "exo"]
                if cr:
                    xc = crystallinity(max(cr, key=lambda x: x["dH"])["dH"], a.dh0 or a.polymer, w=a.w)
                    print(f"  Cristalinidade desenvolvida no resfriamento: Xc = {xc['Xc']} % (ΔHc/ΔH0)")
    if not manual:
        print("\nJanelas AUTOMÁTICAS: confira no gráfico e refaça com --tg T1:T2 / --peak T1:T2 se necessário.")
    if a.csv and rows:
        pd.DataFrame(rows).to_csv(a.csv, index=False)
        print(f"gravado: {a.csv}")
    if a.output:
        plot_dsc(r, keys, results=results, display=a.display, endo=d, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")


def cmd_sta(a):
    r = _load([a.file], a, one=True)
    if r.kind != "STA":
        sys.exit(f"ERRO: {r.title} é {r.kind}, não STA (TG + DSC/DTA simultâneos).")
    d, src = endo_direction(r, a.endo)
    a.endo = d
    r2, st, tc, rs = _tga_one(r, a)
    print(f"\nDireção endo: {'↑' if d == 'up' else '↓'}  ({src})")
    h = select(r, "main_heat")
    res = []
    print("\nCaráter térmico de cada etapa de massa (sinal DSC/DTA integrado na faixa da etapa):")
    for _, s in st.iterrows():
        if str(s["tipo"]) not in ("perda", "ganho"):
            continue
        hw = s.get("largura_DTG_K", np.nan)
        hw = hw if np.isfinite(hw) else (s["T_fim"] - s["T_inicio"]) / 3
        lo = max(s["T_inicio"], s["T_pico"] - max(1.2 * hw, 10))
        hi = min(s["T_fim"], s["T_pico"] + max(1.2 * hw, 10))
        try:
            p = integrate_peak(h, (lo, hi), endo=d)
            yy = h.hf[np.isfinite(h.hf)]
            if abs(p["altura"]) >= 0.05 * (np.nanpercentile(yy, 99) - np.nanpercentile(yy, 1)) * (1 / h.sample_mass if h.units.get("hf") == "mW" and h.sample_mass else 1):
                res.append(p)                  # só anota no gráfico eventos visíveis
            print(f"  etapa {s['etapa']}: {s['T_inicio']:.0f}–{s['T_fim']:.0f} °C, Δm = {s['perda_pct']:.2f}%  →  "
                  f"{p['tipo'].upper()} (DSC integrado em {lo:.0f}–{hi:.0f} °C: pico {p['T_pico']} °C, "
                  f"área {p['dH']:.3g} {p['unidade_dH']})")
        except ValueError as e:
            print(f"  etapa {s['etapa']}: {e}")
    ev = find_events(h, endo=d)
    no_mass = []
    for _, e in ev.iterrows():
        if not any(s["T_inicio"] <= e["T_pico"] <= s["T_fim"] for _, s in st.iterrows()):
            no_mass.append(e)
    if no_mass:
        print("\nEventos DSC/DTA SEM perda de massa (fusão, transição de fase, Tg, cristalização...):")
        for e in no_mass:
            print(f"  {e['tipo']} {e['evento']} em {e['T_pico']} °C ({e['T_ini']}–{e['T_fim']} °C): {e['sugestao']}")
    print("\nNota: com DSC de STA, áreas em mW·s/mg só são entalpias se o sinal estiver calibrado; "
          "com perda de massa simultânea, normalize pela massa inicial e trate como semiquantitativo.")
    if a.output:
        plot_sta(r, steps=st, results=res, display=a.display, endo=d, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")


def cmd_kinetics(a):
    runs = _load(a.files, a)
    w = _win(a.window)
    res = isoconversional(runs, w, signal=a.signal, endo=a.endo)
    tb = res["table"]
    meths = [m for m in ("friedman", "ofw", "kas", "starink", "vyazovkin") if m in tb]
    show = ["alpha"] + [c for m in meths for c in (m, f"{m}_ic95") if c in tb and not tb[c].isna().all()]
    print(f"Taxas medidas (K/min): {', '.join(f'{b:.2f}' for b in res['betas'])}   janela {w[0]:g}–{w[1]:g} °C")
    print(tb[show].round(1).to_string(index=False))
    print("\nEa médio (0.1 ≤ α ≤ 0.9), kJ/mol:")
    for m, v in res["summary"].items():
        print(f"  {m:<10} {v['media']:7.1f} ± {v['dp']:.1f} (faixa {v['min']:.0f}–{v['max']:.0f})")
    try:
        k = kissinger(runs, w, signal=a.signal, endo=a.endo)
        print(f"  Kissinger  {k['Ea']:7.1f} ± {k['ic95']:.1f} (IC95%), R² = {k['r2']:.4f}, ln A = {k['lnA']:.2f} (A em s⁻¹, 1ª ordem)"
              f"; Tp = {', '.join(f'{x:.1f}' for x in k['Tp_C'])} °C")
        res["warnings"] += k["warnings"]
    except ValueError as e:
        print("  Kissinger:", e)
    for x in res["warnings"]:
        print("AVISO:", x)
    if a.csv:
        tb.to_csv(a.csv, index=False)
        print(f"gravado: {a.csv}")
    if a.output:
        plot_kinetics(res, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")


def cmd_oit(a):
    r = _load([a.file], a, one=True)
    res = oit(r, t_switch=a.t_switch, endo=a.endo)
    for k, v in res.items():
        print(f"{k:<12}: {v}")
    if a.output and np.isfinite(res.get("OIT", np.nan)):
        plot_oit(r, res, endo=a.endo, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")


def cmd_cp(a):
    s = _load([a.sample], a, one=True)
    b = _load([a.blank], a, one=True)
    z = _load([a.sapphire], a, one=True)
    df = cp_sapphire(s, b, z, a.m_sample, a.m_sapphire, _win(a.range) if a.range else None)
    print(df.iloc[::max(1, len(df) // 25)].to_string(index=False))
    if a.output:
        df.to_csv(a.output, index=False)
        print(f"gravado: {a.output}")


def cmd_stoich(a):
    th = theoretical_losses(a.formula, a.release)
    print(f"{a.formula}: M = {th.attrs['M']:.3f} g/mol")
    print(th.to_string(index=False))
    if a.compare:
        r = _load([a.compare], a, one=True)
        print(compare_with_theory(find_steps(r, min_loss=a.min_loss), th).to_string(index=False))


def cmd_hydrate(a):
    n = hydration_number(a.loss, a.anhydrous)
    print(f"{a.anhydrous}·nH2O com perda de {a.loss}% -> n = {n:.2f}  (M anidro = {molar_mass(a.anhydrous):.2f} g/mol)")


def cmd_plot(a):
    runs = _load(a.files, a)
    kind = a.kind or ("tga" if all(r.has("mass_pct") for r in runs) else "dsc")
    if kind == "tga":
        plot_tga(runs, x=a.x, lang=a.lang, out=a.output, labels=a.names, annotate=False)
    else:
        import matplotlib.pyplot as plt
        from ta_plot import _style, _disp, PALETTE, _arrow, _hf_label, TXT
        from ta_process import endo_signal
        _style()
        fig, ax = plt.subplots(figsize=(6.6, 4.6))
        step = None
        curves = []
        for i, r in enumerate(runs):
            sr = select(r, a.segment[0]) if a.segment else r
            Tt, t, y, u, d, src = endo_signal(sr, a.endo)
            yd = _disp(y, a.display)
            curves.append((Tt, yd - np.nanmedian(yd), u, (a.names[i] if a.names else r.title)))
        step = 1.1 * max(np.nanmax(c[1]) - np.nanmin(c[1]) for c in curves) if a.offset else 0
        for i, (Tt, yd, u, lab) in enumerate(curves):
            off = (len(curves) - 1 - i) * step
            ax.plot(Tt, yd + off, color=PALETTE[i % len(PALETTE)], lw=1.3, label=lab)
        ax.set_xlabel(TXT[a.lang]["T"])
        ax.set_ylabel(_hf_label(curves[0][2], a.lang) + (" — deslocado" if step else ""))
        if step:
            ax.set_yticklabels([])
        ax.legend(fontsize=8, loc="lower left")
        _arrow(ax, a.display, a.lang, "upper right")
        fig.tight_layout()
        fig.savefig(a.output)
    print(f"gravado: {a.output}")


def cmd_convert(a):
    r = _load([a.file], a, one=True)
    save_csv(r, a.output, decimal=a.decimal, extra=a.extra)
    print(f"gravado: {a.output}  ({r.n} pontos; colunas canônicas: tempo min, T °C, massa mg/%, fluxo {r.units.get('hf', '-')})")


def cmd_report(a):
    os.makedirs(a.dir, exist_ok=True)
    r = _load([a.file], a, one=True)
    base = os.path.join(a.dir, "".join(c if c.isalnum() or c in "-_." else "_" for c in _stem(r))[:60])
    import contextlib
    import io
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print(r.summary())
        _print_segments(r)
        if r.kind in ("TGA", "STA"):
            print()
            _tga_one(r, a, base + "_TG.png", base + "_etapas.csv", quiet=(r.kind == "STA"))
        if r.kind in ("DSC", "DTA"):
            a.file, a.output, a.csv = a.file, base + "_DSC.png", base + "_eventos.csv"
            a.segment, a.tg, a.peak = None, None, None
            print()
            cmd_dsc(a)
        if r.kind == "STA":
            a.output = base + "_STA.png"
            print()
            cmd_sta(a)
    txt = buf.getvalue()
    with open(base + "_relatorio.md", "w", encoding="utf-8") as fh:
        fh.write(f"# Análise térmica: {r.title}\n\n```\n{txt}\n```\n")
    print(txt)
    print(f"relatório: {base}_relatorio.md")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, endo=True):
        sp.add_argument("--lang", default="pt", choices=["pt", "en"])
        sp.add_argument("--columns", help="papéis das colunas p/ arquivo sem cabeçalho: time:min,T:°C,mass:mg")
        sp.add_argument("--index", type=int, help="curva a usar quando o arquivo tem várias")
        sp.add_argument("--mass", type=float, help="massa da amostra (mg) se o arquivo não informa")
        sp.add_argument("--correction", help="NETZSCH: arquivo de correção .ngb-bs3 a subtrair")
        sp.add_argument("--decimal", default=".")
        if endo:
            sp.add_argument("--endo", choices=["up", "down"], help="direção endotérmica do sinal no arquivo")
            sp.add_argument("--display", default="exo_up", choices=["exo_up", "endo_up"])

    def tgaopts(sp):
        sp.add_argument("--basis", type=float, help="base seca: 100%% = massa a esta T (°C)")
        sp.add_argument("--min-loss", type=float, default=0.5)
        sp.add_argument("--theory", help='"Reagente:espécie1,espécie2,..." ex. "CaC2O4·H2O:H2O,CO,CO2"')
        sp.add_argument("--blank", help="curva de branco (cadinho vazio) para subtrair da massa")
        sp.add_argument("--blank-time", action="store_true", help="alinhar branco pelo tempo (programas com isotermas)")
        sp.add_argument("--x", default="auto", choices=["auto", "T", "time"])

    s = sub.add_parser("info"); s.add_argument("files", nargs="+"); common(s); s.set_defaults(f=cmd_info)
    s = sub.add_parser("tga"); s.add_argument("files", nargs="+"); s.add_argument("-o", "--output")
    s.add_argument("--csv"); s.add_argument("--names", nargs="+"); tgaopts(s); common(s, False); s.set_defaults(f=cmd_tga)
    s = sub.add_parser("dsc"); s.add_argument("file"); s.add_argument("--segment", nargs="+")
    s.add_argument("--tg", action="append"); s.add_argument("--peak", action="append")
    s.add_argument("--polymer"); s.add_argument("--dh0", type=float); s.add_argument("--w", type=float, default=1.0)
    s.add_argument("--baseline", default="linear", choices=["linear", "sigmoidal"])
    s.add_argument("-o", "--output"); s.add_argument("--csv"); common(s); s.set_defaults(f=cmd_dsc)
    s = sub.add_parser("sta"); s.add_argument("file"); s.add_argument("-o", "--output"); tgaopts(s)
    common(s); s.set_defaults(f=cmd_sta, csv=None)
    s = sub.add_parser("kinetics"); s.add_argument("files", nargs="+"); s.add_argument("--window", required=True)
    s.add_argument("--signal", default="auto", choices=["auto", "mass", "hf"]); s.add_argument("-o", "--output")
    s.add_argument("--csv"); common(s); s.set_defaults(f=cmd_kinetics)
    s = sub.add_parser("oit"); s.add_argument("file"); s.add_argument("--t-switch", type=float)
    s.add_argument("-o", "--output"); common(s); s.set_defaults(f=cmd_oit)
    s = sub.add_parser("cp"); s.add_argument("--sample", required=True); s.add_argument("--blank", required=True)
    s.add_argument("--sapphire", required=True); s.add_argument("--m-sample", type=float)
    s.add_argument("--m-sapphire", type=float); s.add_argument("--range"); s.add_argument("-o", "--output")
    common(s); s.set_defaults(f=cmd_cp)
    s = sub.add_parser("stoich"); s.add_argument("formula"); s.add_argument("--release", nargs="+", required=True)
    s.add_argument("--compare"); s.add_argument("--min-loss", type=float, default=0.5); common(s, False)
    s.set_defaults(f=cmd_stoich)
    s = sub.add_parser("hydrate"); s.add_argument("--loss", type=float, required=True)
    s.add_argument("--anhydrous", required=True); s.set_defaults(f=cmd_hydrate)
    s = sub.add_parser("plot"); s.add_argument("files", nargs="+"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--kind", choices=["tga", "dsc"]); s.add_argument("--segment", nargs="+")
    s.add_argument("--names", nargs="+"); s.add_argument("--offset", action="store_true")
    s.add_argument("--x", default="auto", choices=["auto", "T", "time"]); common(s); s.set_defaults(f=cmd_plot)
    s = sub.add_parser("convert"); s.add_argument("file"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--extra", action="store_true"); common(s, False); s.set_defaults(f=cmd_convert)
    s = sub.add_parser("report"); s.add_argument("file"); s.add_argument("-d", "--dir", default="ta_report")
    tgaopts(s); s.add_argument("--polymer"); s.add_argument("--dh0", type=float); s.add_argument("--w", type=float, default=1.0)
    s.add_argument("--baseline", default="linear", choices=["linear", "sigmoidal"])
    common(s); s.set_defaults(f=cmd_report, csv=None, output=None, segment=None, tg=None, peak=None)
    a = p.parse_args(argv)
    np.seterr(all="ignore")
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    try:
        a.f(a)
    except (ValueError, FileNotFoundError, ImportError) as e:
        sys.exit(f"ERRO: {e}")


if __name__ == "__main__":
    main()
