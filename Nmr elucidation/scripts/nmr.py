#!/usr/bin/env python3
"""
nmr.py - linha de comando da skill de RMN ¹H/¹³C (JCAMP-DX e outros)

  info        arquivos...                        núcleo, campo, solvente, eixo, referência, ruído, fase
  h1          1h.jdx [--formula C9H8O4] [--norm 2.29=3] [-o fig.png] [--multiplets mult.png] [--csv h.csv]
  c13         13c.jdx [--dept135 d135.jdx] [--dept90 d90.jdx] [--apt apt.jdx] [-o fig.png] [--csv c.csv]
  report      --h1 1h.jdx --c13 13c.jdx [--dept135 ...] [--dept90 ...] [--formula ...] -d pasta
  candidates  --h1 1h.jdx --c13 13c.jdx [--dept135 ...] --smiles "nome=SMILES" ... | --file cand.smi
              [--formula C11H14N2O] [-d pasta]       ranqueia estruturas candidatas (predição CSP5)
  predict     "SMILES" [--solvent cdcl3] [-o estrutura.png]   δ previstos por átomo + sinais esperados
  impurities  [--solvent cdcl3] [--nucleus 1H]       tabela de solventes/impurezas comuns
  process     fid.jdx -o espectro.jdx|.csv [--lb 0.3] [--zf 2]   FID -> espectro (fase automática)
  stack       arquivos... -o fig.png [--xlim 10:0]   sobreposição de espectros

Comuns: --solvent cdcl3 (se o arquivo não diz), --ref auto|file|tms|7.31=7.26, --sf MHz (texto),
--no-baseline, --lang pt|en, --decimal , (CSV/listas com vírgula)
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
np.seterr(all="ignore")

from nmr_io import load, save_csv, save_jcamp  # noqa: E402
from nmr_process import prepare, reference, noise  # noqa: E402


def _load(path, a, ref=None):
    s = load(path, nucleus=getattr(a, "nucleus", None), sf=a.sf, solvent=a.solvent)
    s, info = prepare(s, do_baseline=not a.no_baseline, ref=ref or a.ref)
    return s, info


def _ref_msgs(s, info):
    out = []
    if info.get("rephased"):
        out.append(info["rephased"])
    if info.get("ref"):
        out.append("Referência: " + info["ref"].get("note", ""))
    return out


def cmd_info(a):
    for f in a.files:
        s0 = load(f, sf=a.sf, solvent=a.solvent)
        print(s0.summary())
        if s0.kind == "spectrum":
            s, info = prepare(s0, do_baseline=not a.no_baseline, ref=a.ref)
            for m in _ref_msgs(s, info):
                print(m)
            ymax = np.max(s.y)
            print(f"Ruído σ = {info['noise']:.3g}; maior sinal S/N = {ymax / info['noise']:.0f}; "
                  f"resolução digital {s.hz_per_pt:.3f} Hz/ponto"
                  + ("  (grossa para J < 1 Hz: use zero-filling maior)" if s.hz_per_pt > 0.5 and s.nucleus == "1H" else ""))
            if s.nucleus == "1H" and s.sf < 100:
                print("  - campo baixo (bancada): multipletos de 2ª ordem frequentes; J e multiplicidades com cautela")
        elif s0.peaks is not None:
            print(pd.DataFrame(s0.peaks, columns=["ppm", "intensidade"]).head(30).to_string(index=False))
        print()


def _h1(a, path):
    from nmr_peaks import analyze_1h
    s, info = _load(path, a)
    if s.nucleus != "1H":
        print(f"AVISO: {path} é {s.nucleus}, não ¹H")
    r = analyze_1h(s, formula=a.formula, norm=a.norm, k=a.k)
    return s, info, r


def _print_h1(s, info, r, a):
    from nmr_peaks import acs_1h
    for m in _ref_msgs(s, info):
        print(m)
    t = r["table"]
    print(f"Integrais normalizadas: {r['norm_desc']}   (σ = {r['sigma']:.3g})")
    cols = ["ppm", "mult", "J", "nH", "nH_int", "ppm_max", "ppm_min", "rotulo", "nota", "parceiros_J"]
    tt = t[cols].copy()
    tt["J"] = tt.J.apply(lambda j: ", ".join(f"{x:.1f}" for x in j))
    tt["nH"] = tt.nH.round(2)
    tt["ppm"] = tt.ppm.round(3)
    tt["ppm_max"] = tt.ppm_max.round(2)
    tt["ppm_min"] = tt.ppm_min.round(2)
    print(tt.to_string(index=False))
    for w in r.get("warnings", []):
        print(f"AVISO: {w}")
    print()
    print(acs_1h(r, s, decimal=a.decimal))


def cmd_h1(a):
    from nmr_plot import plot_1h, plot_multiplets
    s, info, r = _h1(a, a.file)
    _print_h1(s, info, r, a)
    if a.csv:
        t = r["table"].drop(columns=["linhas_hz", "alturas", "i0", "i1"])
        t["J"] = t.J.apply(lambda j: "; ".join(f"{x:.1f}" for x in j))
        t.to_csv(a.csv, index=False, sep=";" if a.decimal == "," else ",", decimal=a.decimal, float_format="%.4f")
        print(f"gravado: {a.csv}")
    if a.output:
        plot_1h(s, r, out=a.output, lang=a.lang, xlim=_xl(a))
        print(f"gravado: {a.output}")
    if a.multiplets:
        plot_multiplets(s, r, out=a.multiplets)
        print(f"gravado: {a.multiplets}")
    return s, r


def _xl(a):
    if getattr(a, "xlim", None):
        hi, lo = [float(v) for v in a.xlim.split(":")]
        return max(hi, lo), min(hi, lo)
    return None


def _c13(a, path, dept135=None, dept90=None, apt=None):
    from nmr_peaks import analyze_13c
    s, info = _load(path, a)
    delta = info.get("ref", {}).get("delta", 0.0) if "re-referenciado" in info.get("ref", {}).get("note", "") else 0.0
    ds = {}
    for key, p in (("d135", dept135), ("d90", dept90), ("apt", apt)):
        if p:
            d, _ = _load(p, a, ref="file")
            ds[key] = reference(d, delta) if delta else d
    r = analyze_13c(s, k=a.k13, dept135=ds.get("d135"), dept90=ds.get("d90"), apt=ds.get("apt"))
    return s, info, r, ds


def _print_c13(s, info, r, a):
    from nmr_peaks import acs_13c
    for m in _ref_msgs(s, info):
        print(m)
    t = r["table"]
    if len(t):
        tt = t[["ppm", "snr", "mult", "rotulo", "ambientes"]].copy()
        tt["ppm"] = tt.ppm.round(2)
        tt["snr"] = tt.snr.round(0)
        tt["ambientes"] = tt.ambientes.str.slice(0, 70)
        print(tt.to_string(index=False))
    if r.get("counts"):
        print("Contagem: " + ", ".join(f"{k} {v}" for k, v in r["counts"].items()) +
              (f"; H em carbono (DEPT) = {r['nH_CH']}" if r.get("nH_CH") is not None else ""))
    if r.get("apt"):
        print("APT: " + r["apt"])
    print()
    print(acs_13c(r, s))


def cmd_c13(a):
    from nmr_plot import plot_13c
    s, info, r, ds = _c13(a, a.file, a.dept135, a.dept90, a.apt)
    _print_c13(s, info, r, a)
    if a.csv:
        r["table"].drop(columns=["idx"], errors="ignore").to_csv(a.csv, index=False, sep=";" if a.decimal == "," else ",",
                                                              decimal=a.decimal, float_format="%.3f")
        print(f"gravado: {a.csv}")
    if a.output:
        plot_13c(s, r, dept135=ds.get("d135"), dept90=ds.get("d90"), out=a.output, lang=a.lang)
        print(f"gravado: {a.output}")
    return s, r


def cmd_report(a):
    from nmr_peaks import acs_1h, acs_13c, inventory
    from nmr_plot import plot_1h, plot_13c, plot_multiplets
    import contextlib
    import io
    os.makedirs(a.dir, exist_ok=True)
    buf = io.StringIO()
    md = ["# Relatório de RMN\n"]
    r1 = r13 = s1 = s13 = None
    with contextlib.redirect_stdout(buf):
        if a.h1:
            s1, i1, r1 = _h1(a, a.h1)
            print(s1.summary())
            _print_h1(s1, i1, r1, a)
            plot_1h(s1, r1, out=os.path.join(a.dir, "h1.png"), lang=a.lang)
            plot_multiplets(s1, r1, out=os.path.join(a.dir, "h1_multipletos.png"))
            print()
        if a.c13:
            s13, i13, r13, ds = _c13(a, a.c13, a.dept135, a.dept90, a.apt)
            print(s13.summary())
            _print_c13(s13, i13, r13, a)
            plot_13c(s13, r13, dept135=ds.get("d135"), dept90=ds.get("d90"), out=os.path.join(a.dir, "c13.png"), lang=a.lang)
    txt = buf.getvalue()
    print(txt)
    inv = inventory(r1, r13, a.formula)
    print("Resumo estrutural (dos dados, sem propor estrutura):")
    for l in inv:
        print("  - " + l)
    md.append("## Dados para a parte experimental\n")
    if r1 is not None:
        md.append(acs_1h(r1, s1, decimal=a.decimal) + "\n")
    if r13 is not None:
        md.append(acs_13c(r13, s13) + "\n")
    md.append("\n## Resumo estrutural\n" + "\n".join("- " + l for l in inv) + "\n")
    md.append("\n## Saída completa\n```\n" + txt + "```\n")
    for f in ("h1.png", "h1_multipletos.png", "c13.png"):
        if os.path.exists(os.path.join(a.dir, f)):
            md.append(f"![{f}]({f})\n")
    p = os.path.join(a.dir, "relatorio.md")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("\n".join(md))
    print(f"gravado: {p} (+ figuras em {a.dir})")


def cmd_candidates(a):
    from nmr_struct import Candidate, rank, read_candidates, draw_assigned
    from nmr_plot import plot_candidate
    cands = []
    if a.file:
        cands += read_candidates(a.file)
    if a.smiles:
        cands += read_candidates(a.smiles)
    if not cands:
        print("ERRO: informe --smiles 'nome=SMILES' ... ou --file candidatas.smi")
        sys.exit(2)
    import contextlib
    import io
    r1 = r13 = None
    sol = a.solvent
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        if a.h1:
            s1, i1, r1 = _h1(a, a.h1)
            sol = sol or s1.solvent
        if a.c13:
            s13, i13, r13, _ = _c13(a, a.c13, a.dept135, a.dept90, a.apt)
            sol = sol or s13.solvent
    res = rank(cands, h1=r1, c13=r13, solvent=sol, formula=a.formula, num_conformers=a.conformers)
    t = res["table"]
    if not res["predicted"]:
        print("AVISO: CSP5 não instalado (pip install csp5) - só contagens/simetria, sem ranquear por deslocamento.")
    show = [c for c in ["candidata", "formula", "formula_ok", "prob_rel_pct", "MAE_13C", "max_13C", "MAE_1H", "max_1H",
                        "C_sem_sinal", "sinais_C_sobrando", "H_sem_sinal", "sinais_C_previstos", "sinais_C_obs",
                        "H_total", "H_obs"] if c in t]
    tt = t[show].copy()
    for c in ("prob_rel_pct", "MAE_13C", "max_13C"):
        if c in tt:
            tt[c] = tt[c].round(1)
    for c in ("MAE_1H", "max_1H"):
        if c in tt:
            tt[c] = tt[c].round(2)
    print(tt.to_string(index=False))
    for r in t.itertuples():
        if r.alertas:
            print(f"  {r.candidata}: {r.alertas}")
    print("\nprob_rel_pct = probabilidade relativa ENTRE as candidatas (estilo DP4: t de Student, σ = 2,0 ppm ¹³C e "
          "0,22 ppm ¹H). Não diz nada sobre estruturas fora da lista.")
    if a.dir and res["predicted"]:
        os.makedirs(a.dir, exist_ok=True)
        t.to_csv(os.path.join(a.dir, "ranking.csv"), index=False, float_format="%.4f")
        for k, r in enumerate(t.itertuples()):
            if k >= a.top:
                break
            det = res["details"][r.candidata]
            safe = "".join(ch if ch.isalnum() else "_" for ch in r.candidata)
            plot_candidate(r.candidata, det, out=os.path.join(a.dir, f"{k + 1}_{safe}_atribuicao.png"))
            try:
                draw_assigned(res["cands"][r.candidata], det, os.path.join(a.dir, f"{k + 1}_{safe}_estrutura.png"))
            except Exception as e:
                print(f"(estrutura não desenhada: {e})")
            rows = []
            for x in det["c13"]:
                rows.append(dict(nucleo="13C", atomo=x["label"], tipo=x["mult"], n=x["n"], previsto=round(x["pred"], 1),
                                 experimental=round(x["exp"], 2) if np.isfinite(x["exp"]) else None,
                                 erro=round(x["err"], 2) if np.isfinite(x["err"]) else None))
            for x in det["h1"]:
                rows.append(dict(nucleo="1H", atomo=x["label"], tipo="trocável" if x["exch"] else "", n=x["n"],
                                 previsto=round(x["pred"], 2), experimental=round(x["exp"], 3) if np.isfinite(x["exp"]) else None,
                                 erro=round(x["err"], 2) if np.isfinite(x["err"]) else None))
            pd.DataFrame(rows).to_csv(os.path.join(a.dir, f"{k + 1}_{safe}_atribuicao.csv"), index=False)
            if k == 0:
                print(f"\nAtribuição da 1ª colocada ({r.candidata}):")
                print(pd.DataFrame(rows).to_string(index=False))
        print(f"gravado: {a.dir}/ (ranking.csv, figuras e tabelas de atribuição das {min(a.top, len(t))} primeiras)")


def cmd_predict(a):
    from nmr_struct import Candidate, mol_info, predict, predicted_signals
    c = Candidate("molecula", a.smiles)
    inf = mol_info(c)
    ok = predict([c], solvent=a.solvent, num_conformers=a.conformers)
    print(f"Fórmula {inf['formula']}  MM {inf['mw']:.2f}  IDH {inf['dbe']:g}  sinais de ¹³C esperados: {inf['n_sig_c']} "
          f"({', '.join(f'{k} {v}' for k, v in inf['mult_counts'].items())}); H: {inf['n_h']} ({inf['n_h_exch']} trocáveis)")
    if not ok:
        print("CSP5 não instalado (pip install csp5): sem deslocamentos previstos.")
        return
    C, H = predicted_signals(c)
    print("\n¹³C previsto:")
    print(pd.DataFrame([dict(atomo=x["label"], tipo=x["mult"], n=x["n"], ppm=round(x["pred"], 1)) for x in
                        sorted(C, key=lambda d: -d["pred"])]).to_string(index=False))
    print("\n¹H previsto:")
    print(pd.DataFrame([dict(atomo=x["label"], nH=x["n"], ppm=round(x["pred"], 2), trocavel=x["exch"]) for x in
                        sorted(H, key=lambda d: -d["pred"])]).to_string(index=False))
    if a.output:
        from nmr_struct import draw_assigned
        det = dict(c13=[dict(x, exp=np.nan, err=np.nan) for x in C], h1=[dict(x, exp=np.nan, err=np.nan) for x in H])
        draw_assigned(c, det, a.output)
        print(f"gravado: {a.output} (δ previstos entre parênteses)")


def cmd_impurities(a):
    from nmr_data import IMPURITIES, SOLVENTS
    from nmr_io import norm_solvent
    sol = norm_solvent(a.solvent) or "cdcl3"
    s = SOLVENTS.get(sol, {})
    print(f"Solvente {s.get('name', sol)}: ¹H residual {', '.join(f'{d} ({m})' for d, m in s.get('h', []))}; "
          f"¹³C {', '.join(str(c) for c in s.get('c', []))}; água {s.get('water')}")
    rows = [dict(composto=k, H=", ".join(f"{x:.2f}" for x in v[0]), C=", ".join(f"{x:.2f}" for x in v[1]))
            for k, v in IMPURITIES.get(sol, {}).items()]
    print(pd.DataFrame(rows).to_string(index=False))
    print("Fonte: Gottlieb et al. 1997; Fulmer et al. 2010 (variam ±0,02-0,1 ppm com T e concentração).")


def cmd_process(a):
    s = load(a.file, sf=a.sf, solvent=a.solvent, lb=a.lb, zf=a.zf)
    s, info = prepare(s, do_baseline=not a.no_baseline, ref=a.ref)
    print(s.summary())
    for m in _ref_msgs(s, info):
        print(m)
    if a.output.lower().endswith((".csv", ".txt")):
        save_csv(s, a.output, decimal=a.decimal)
    else:
        save_jcamp(s, a.output)
    print(f"gravado: {a.output}")


def cmd_stack(a):
    from nmr_plot import plot_stack
    specs = [_load(f, a)[0] for f in a.files]
    plot_stack(specs, names=a.names, out=a.output, xlim=_xl(a))
    print(f"gravado: {a.output}")


def main(argv=None):
    p = argparse.ArgumentParser(prog="nmr.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = argparse.ArgumentParser(add_help=False)
    c.add_argument("--solvent")
    c.add_argument("--ref", default="auto", help="auto | file | tms | 7.31=7.26")
    c.add_argument("--sf", type=float, help="frequência (MHz) para arquivos de texto")
    c.add_argument("--no-baseline", action="store_true")
    c.add_argument("--lang", default="pt", choices=["pt", "en"])
    c.add_argument("--decimal", default=".")
    c.add_argument("-k", type=float, default=6.0, help="limiar de picos ¹H (× σ)")
    c.add_argument("--k13", type=float, default=5.0, help="limiar de picos ¹³C (× σ)")
    c.add_argument("--formula")
    c.add_argument("--norm", help="fixa integral: '2.29=3'")
    c.add_argument("--xlim")

    s = sub.add_parser("info", parents=[c]); s.add_argument("files", nargs="+")
    s = sub.add_parser("h1", parents=[c]); s.add_argument("file"); s.add_argument("-o", "--output")
    s.add_argument("--multiplets"); s.add_argument("--csv")
    s = sub.add_parser("c13", parents=[c]); s.add_argument("file"); s.add_argument("-o", "--output")
    s.add_argument("--csv"); s.add_argument("--dept135"); s.add_argument("--dept90"); s.add_argument("--apt")
    for name in ("report", "candidates"):
        s = sub.add_parser(name, parents=[c])
        s.add_argument("--h1"); s.add_argument("--c13"); s.add_argument("--dept135"); s.add_argument("--dept90")
        s.add_argument("--apt"); s.add_argument("-d", "--dir")
        if name == "candidates":
            s.add_argument("--smiles", nargs="+"); s.add_argument("--file"); s.add_argument("--top", type=int, default=3)
            s.add_argument("--conformers", type=int, default=1)
    s = sub.add_parser("predict", parents=[c]); s.add_argument("smiles"); s.add_argument("-o", "--output")
    s.add_argument("--conformers", type=int, default=1)
    s = sub.add_parser("impurities", parents=[c]); s.add_argument("--nucleus", default="1H")
    s = sub.add_parser("process", parents=[c]); s.add_argument("file"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--lb", type=float); s.add_argument("--zf", type=int, default=2)
    s = sub.add_parser("stack", parents=[c]); s.add_argument("files", nargs="+"); s.add_argument("-o", "--output", required=True)
    s.add_argument("--names", nargs="+")
    a = p.parse_args(argv)
    if a.cmd == "report" and not a.dir:
        p.error("report precisa de -d pasta")
    try:
        {"info": cmd_info, "h1": cmd_h1, "c13": cmd_c13, "report": cmd_report, "candidates": cmd_candidates,
         "predict": cmd_predict, "impurities": cmd_impurities, "process": cmd_process, "stack": cmd_stack}[a.cmd](a)
    except (ValueError, FileNotFoundError, ImportError) as e:
        print(f"ERRO: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
