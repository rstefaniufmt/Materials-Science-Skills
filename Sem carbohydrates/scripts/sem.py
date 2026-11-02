#!/usr/bin/env python3
"""
sem.py - linha de comando da skill de MEV para carboidratos (grânulos de amido, partículas e filmes)

  info      arquivos...                     escala (e de onde veio), kV, aumento, WD, detector, qualidade
  scale     arquivo [--bar-value 10um] [-o calibrado.tif]   detecta databar/barra/OCR; grava TIF calibrado (ImageJ)
  granules  arquivos... [-o fig.png] [--csv medidas.csv] [--min-d 2um] [--max-d 60um] [--thresh otsu|li|local|0.4]
            [--dark-objects] [--no-split] [--split-h 0.25] [--flatten] [--pits] [--pit-max 1.5um] [--numbers]
  film      arquivo [-o fig.png] [--csv defeitos.csv] [--pore-max 5um] [--particle-max 15um] [--crack-min 10um]
  section   arquivo [-o fig.png] [--csv perfil.csv] [--orientation auto|h|v] [--film-darker]
  panel     arquivos... -o painel.png [--titles A B ...] [--ncols 2] [--bar 10um] [--no-letters] [--enhance]
  compare   a.csv b.csv ... --names A B ... [--col ecd|AR|circularidade|solidez] [-o comp.png]
  measure   arquivo --points x1,y1,x2,y2      distância entre dois pontos (px da imagem sem databar)
  report    arquivo(s) --kind granules|film|section -d pasta    análise completa + relatorio.md

Escala (todas): --pixel-size 0.05um | --scale 10um=235px | --bar-value 10um (mede a barra desenhada)
Recorte: --databar N (linhas da base a remover; 0 = não remover) | --crop x0,y0,x1,y1  |  --lang pt|en
"""
from __future__ import annotations

import argparse
import os
import sys
import warnings

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
np.seterr(all="ignore")

from sem_io import fmt_len, load, parse_length, save_tif  # noqa: E402


def _len(s):
    return None if s is None else parse_length(s)


def _load(path, a):
    db = a.databar if a.databar is not None else "auto"
    crop = [int(v) for v in a.crop.split(",")] if a.crop else None
    return load(path, pixel_size=a.pixel_size, scale=a.scale, bar_value=a.bar_value, databar=db, crop=crop,
                ocr=not a.no_ocr)


def _u(im):
    return "µm" if im.px_um else "px"


def _r(v, nd=2):
    return f"{v:.{nd}f}" if isinstance(v, (float, np.floating)) and np.isfinite(v) else str(v)


# ----------------------------------------------------------------------------- comandos

def cmd_info(a):
    from sem_process import quality
    for f in a.files:
        im = _load(f, a)
        print(im.summary())
        q = quality(im)
        print(f"Qualidade   : SNR ≈ {q['snr']:.1f}, ruído σ = {q['ruido_sigma']:.3f}, saturado {q['saturado_pct']:.2f}%, "
              f"gradiente de iluminação {q['gradiente_iluminacao_pct']:.0f}%")
        for n in q["notas"]:
            print(f"  - {n}")
        if im.meta.get("hv_kV") and im.meta["hv_kV"] > 10:
            print(f"  - {im.meta['hv_kV']:g} kV é alto para carboidratos (dano por feixe/carregamento, perda de detalhe de "
                  "superfície): 1-5 kV costuma ser melhor")
        print()


def cmd_scale(a):
    im = _load(a.file, a)
    print(im.summary())
    if a.output:
        save_tif(im, a.output)
        print(f"gravado: {a.output}  (TIF 8 bits sem databar, calibrado para ImageJ/Fiji)")


def _seg_one(f, a):
    from sem_granules import measure, pits, segment
    im = _load(f, a)
    seg = segment(im, method=a.thresh, dark_objects=a.dark_objects, min_d_um=_len(a.min_d), max_d_um=_len(a.max_d),
                  smooth=a.smooth, do_flatten=a.flatten, split=not a.no_split, split_h=a.split_h)
    df = measure(im, seg["labels"], exclude_edge=not a.keep_edge)
    pr = None
    if a.pits and len(df):
        pr = pits(im, seg["labels"], df, pit_max_um=_len(a.pit_max) or 1.5)
    df.insert(0, "imagem", im.title)
    return im, seg, df, pr


def _print_gran(df, st, bm, sh, pr, info=None):
    u = st.get("unidade", "µm")
    print(f"Objetos: {st['n']} válidos (+{st['n_borda']} tocando a borda, excluídos)")
    if st["n"] == 0:
        print("Nenhum objeto válido - ajuste --thresh, --min-d, --dark-objects ou o recorte.")
        return
    print(f"ECD ({u}): média {st['media']:.2f} ± {st['dp']:.2f} (DP; IC95% ± {st.get('ic95_media', np.nan):.2f}), "
          f"CV {st['cv_pct']:.0f}%, mediana {st['mediana']:.2f}, faixa {st['min']:.2f}–{st['max']:.2f}")
    print(f"  em número: D10 {st['D10']:.2f}  D50 {st['D50']:.2f}  D90 {st['D90']:.2f}  span {st['span_n']:.2f}")
    print(f"  em volume: Dv10 {st['Dv10']:.2f}  Dv50 {st['Dv50']:.2f}  Dv90 {st['Dv90']:.2f}  span {st['span_v']:.2f}"
          f"   D[4,3] {st['D43']:.2f}  D[3,2] {st['D32']:.2f}")
    v = df[~df["borda"]]
    print(f"Forma (média ± DP): Feret máx {v.feret_max.mean():.2f} ± {v.feret_max.std():.2f} {u}; "
          f"AR {v.AR.mean():.2f} ± {v.AR.std():.2f}; circularidade {v.circularidade.mean():.3f} ± {v.circularidade.std():.3f}; "
          f"arredondamento {v.arredondamento.mean():.3f}; solidez {v.solidez.mean():.3f}")
    print("Classes de forma (%): " + ", ".join(f"{k} {p:.0f}" for k, p in sh["classes_pct"].items()))
    if bm.get("bimodal"):
        print(f"DISTRIBUIÇÃO BIMODAL (ΔBIC = {bm['dBIC']:.0f}; corte ≈ {bm['corte']:.2f} {u}):")
        for p in bm["populacoes"]:
            print(f"  {p['populacao']:<18} n = {p['n']:<5} {100 * p['frac_numero']:.0f}% em número, "
                  f"{100 * p['frac_volume']:.0f}% em volume; média {p['media']:.2f} ± {p['dp']:.2f} {u}")
    elif "motivo" in bm:
        print(f"Bimodalidade: não testada ({bm['motivo']})")
    else:
        print(f"Distribuição unimodal (ΔBIC = {bm['dBIC']:.0f})")
    if pr:
        if pr.get("n_pites", 0):
            print(f"Pites/poros superficiais: {pr['n_pites']} no total; {100 * pr['frac_graos_com_pites']:.0f}% dos grânulos "
                  f"com pites; {pr['pites_por_grao']:.2f} pites/grânulo; d médio {pr['d_pite_medio']:.2f} {u}; "
                  f"área {pr['area_pites_pct']:.2f}% da superfície projetada")
        else:
            print("Pites/poros superficiais: nenhum detectado")
        print("  (contagem por imagem 2D da face visível - compare com o amido nativo nas mesmas condições)")
    for w in st.get("avisos", []):
        print(f"AVISO: {w}")


def cmd_granules(a):
    from sem_granules import bimodal, shape_summary, size_stats
    from sem_plot import plot_granules
    ims, dfs, segs, prs = [], [], [], []
    for f in a.files:
        im, seg, df, pr = _seg_one(f, a)
        ims.append(im)
        segs.append(seg)
        dfs.append(df)
        prs.append(pr)
        print(f"# {im.title}: pixel {fmt_len(im.px_um) if im.px_um else 'NÃO CALIBRADO'} ({im.px_source}); "
              f"limiar {seg['info']['threshold']}; {int((~df.borda).sum()) if len(df) else 0} objetos válidos")
        for w in im.warnings:
            print(f"  AVISO: {w}")
    units = {_u(im) for im in ims}
    if len(units) > 1:
        print("ERRO: há imagens calibradas e não calibradas juntas - informe a escala de todas.")
        sys.exit(2)
    pxs = {round(im.px_um, 6) for im in ims if im.px_um}
    if len(pxs) > 1:
        print("Nota: imagens com aumentos diferentes - medidas em µm são combináveis, mas o limite inferior de "
              "detecção muda (objetos < ~10 px na imagem de menor aumento são perdidos).")
    df = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    df.attrs["unit"] = units.pop()
    df.attrs["exclude_edge"] = not a.keep_edge
    st = size_stats(df)
    bm = bimodal(df)
    sh = shape_summary(df)
    pr = None
    if a.pits:
        tot = [p for p in prs if p]
        if tot:
            v = df[~df.borda]
            pr = dict(n_pites=sum(p["n_pites"] for p in tot),
                      frac_graos_com_pites=float((v["n_pites"] > 0).mean()) if "n_pites" in v else np.nan,
                      pites_por_grao=float(v["n_pites"].mean()) if "n_pites" in v else np.nan,
                      area_pites_pct=float(v["area_pites_pct"].mean()) if "area_pites_pct" in v else np.nan,
                      d_pite_medio=float(np.nanmean([p["d_pite_medio"] for p in tot])))
    print(f"\n== Resultado ({len(ims)} imagem(ns))")
    _print_gran(df, st, bm, sh, pr)
    if bm.get("bimodal"):
        lab = (df["ecd"] >= bm["corte"]).astype(int)
        df["pop"] = lab
    if a.csv:
        df.drop(columns=[c for c in ("x_px", "y_px") if c in df and not a.keep_xy], errors="ignore") \
          .to_csv(a.csv, index=False, sep=a.sep, decimal=a.decimal, float_format="%.4g")
        print(f"gravado: {a.csv}")
    if a.output:
        d0 = df[df.imagem == ims[0].title]
        plot_granules(ims[0], segs[0]["labels"], d0, stats=st, bm=bm,
                      pits_mask=prs[0]["mask"] if prs[0] else None, lang=a.lang, out=a.output, numbers=a.numbers,
                      title=ims[0].title + (f"  (+{len(ims) - 1} imagens nas estatísticas)" if len(ims) > 1 else ""))
        print(f"gravado: {a.output}")
        if len(ims) > 1 and a.overlays:
            from sem_plot import overlay_granules
            import matplotlib.pyplot as plt
            os.makedirs(a.overlays, exist_ok=True)
            for im, seg, d in zip(ims, segs, dfs):
                fig, ax = plt.subplots(figsize=(7, 7 * im.data.shape[0] / im.data.shape[1]))
                overlay_granules(ax, im, seg["labels"], d, numbers=a.numbers, lang=a.lang)
                p = os.path.join(a.overlays, f"{im.title}_segmentacao.png")
                fig.savefig(p, bbox_inches="tight", dpi=200)
                plt.close(fig)
                print(f"gravado: {p}")
    return ims, df, st, bm, sh, pr


def cmd_film(a):
    from sem_film import surface_defects, surface_texture
    from sem_plot import plot_film
    im = _load(a.file, a)
    print(f"# {im.title}: pixel {fmt_len(im.px_um) if im.px_um else 'NÃO CALIBRADO'} ({im.px_source})")
    for w in im.warnings:
        print(f"  AVISO: {w}")
    tx = surface_texture(im)
    res = surface_defects(im, pore_max_um=_len(a.pore_max), particle_max_um=_len(a.particle_max),
                          crack_min_um=_len(a.crack_min), k=a.k)
    u = res["unidade"]
    r = res["resumo"]
    print(f"Área analisada: {r['area_analisada_mm2']:.4f} mm²" if im.px_um else "Área analisada: (não calibrada)")
    print(f"Poros/furos escuros : {r['n_poros']}" + (f"  ({r['poros_por_mm2']:.0f}/mm²)" if im.px_um else "") +
          (f", ECD {r['poros_ecd_media']:.2f} ± {r['poros_ecd_dp']:.2f} {u}" if r['n_poros'] else "") +
          f", {r['poros_area_pct']:.2f}% da área")
    print(f"Partículas/protuber.: {r['n_particulas']}" + (f"  ({r['particulas_por_mm2']:.0f}/mm²)" if im.px_um else "") +
          (f", ECD {r['particulas_ecd_media']:.2f} ± {r['particulas_ecd_dp']:.2f} {u}" if r['n_particulas'] else "") +
          f", {r['particulas_area_pct']:.2f}% da área")
    print(f"Trincas             : {r['n_trincas']}, comprimento total {r['trincas_comprimento_total']:.1f} {u}" +
          (f" ({r['densidade_trincas_mm_por_mm2']:.2f} mm/mm²)" if im.px_um else ""))
    print(f"Textura (comparativa): CV de brilho {tx['cv_pct']:.1f}%; comprimento de correlação "
          f"{_r(tx['comprimento_correlacao'])} {tx['unidade']}; GLCM (d = {tx['glcm_dist'][1]} {tx['unidade']}): "
          f"contraste {tx['glcm_contrast'][1]:.2f}, homogeneidade {tx['glcm_homogeneity'][1]:.3f}, "
          f"entropia {tx['glcm_entropy'][1]:.2f}" + (f"; D fractal (PSD) {tx['D_fractal_psd']:.2f}"
                                                     if np.isfinite(tx.get('D_fractal_psd', np.nan)) else ""))
    print("Nota: brilho de MEV não é altura - textura é índice comparativo (mesmas condições de aquisição); "
          "rugosidade real → AFM/perfilometria.")
    if a.csv:
        rows = []
        for k, df in res["tabelas"].items():
            if len(df):
                d = df.copy()
                d.insert(0, "tipo", k)
                rows.append(d)
        if rows:
            pd.concat(rows).to_csv(a.csv, index=False, sep=a.sep, decimal=a.decimal, float_format="%.4g")
            print(f"gravado: {a.csv}")
    if a.output:
        plot_film(im, res, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")
    return im, tx, res


def cmd_section(a):
    from sem_film import cross_section
    from sem_plot import plot_section
    im = _load(a.file, a)
    print(f"# {im.title}: pixel {fmt_len(im.px_um) if im.px_um else 'NÃO CALIBRADO'} ({im.px_source})")
    for w in im.warnings:
        print(f"  AVISO: {w}")
    cs = cross_section(im, orientation=a.orientation, film_brighter=False if a.film_darker else None)
    u = cs["unidade"]
    print(f"Espessura: {cs['espessura_media']:.2f} ± {cs['espessura_dp']:.2f} {u} (média ± DP ao longo de "
          f"{cs['n_colunas']} colunas); mín {cs['espessura_min']:.2f}, máx {cs['espessura_max']:.2f}; "
          f"p5–p95 {cs['espessura_p5_p95'][0]:.2f}–{cs['espessura_p5_p95'][1]:.2f} {u}")
    print(f"Orientação da faixa: {'horizontal' if cs['orientacao'] == 'h' else 'vertical'}; inclinação média "
          f"{cs['inclinacao_media_graus']:.1f}° (corrigida: espessura perpendicular à faixa)")
    print(f"Porosidade aparente no corte: {_r(cs['porosidade_area_pct'])}% da área do filme; {cs['n_poros']} poros"
          + (f", ECD médio {cs['poros_ecd_media']:.2f} {u}" if cs['n_poros'] else ""))
    if im.px_um:
        print(f"Resolução da medida: ±{2 * im.px_um:.2f} µm por borda (≈ 2 px); espessura = "
              f"{cs['espessura_media'] / im.px_um:.0f} px")
    for w in cs["avisos"]:
        print(f"AVISO: {w}")
    print("Nota: confirme com micrômetro (média de ≥ 5 pontos); MEV mede o filme seco/sob vácuo e só no campo visto.")
    if a.csv:
        x = np.arange(cs["perfil"].size) * (im.px_um or 1)
        pd.DataFrame(dict(posicao=x[cs["ok"]], espessura=cs["perfil"][cs["ok"]])).to_csv(
            a.csv, index=False, sep=a.sep, decimal=a.decimal, float_format="%.4g")
        print(f"gravado: {a.csv}")
    if a.output:
        plot_section(im, cs, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")
    return im, cs


def cmd_panel(a):
    from sem_plot import panel
    ims = [_load(f, a) for f in a.files]
    for im in ims:
        if not im.px_um:
            print(f"AVISO: {im.title} sem escala - painel sem barra para ela")
    panel(ims, letters=not a.no_letters, titles=a.titles, ncols=a.ncols, bar_um=_len(a.bar), out=a.output,
          enhance=a.enhance, width_in=a.width)
    print(f"gravado: {a.output}")


def cmd_compare(a):
    from sem_granules import compare
    from sem_plot import plot_compare
    names = a.names or [os.path.splitext(os.path.basename(f))[0] for f in a.files]
    samples = {}
    for n, f in zip(names, a.files):
        sep = ";" if open(f, encoding="utf-8").readline().count(";") > 0 else ","
        d = pd.read_csv(f, sep=sep, decimal="," if sep == ";" else ".")
        d.attrs["unit"] = "µm"
        if "borda" in d:
            d["borda"] = d["borda"].astype(str).str.lower().isin(["true", "1", "verdadeiro"])
        else:
            d["borda"] = False
        samples[n] = d
    res = compare(samples, a.col)
    print(f"Variável: {a.col}")
    for r in res["resumo"]:
        print(f"  {r['amostra']:<20} n = {r['n']:<5} média {r['media']:.3f} ± {r['dp']:.3f}  mediana {r['mediana']:.3f}")
    if "kruskal" in res:
        print(f"Kruskal-Wallis: H = {res['kruskal']['H']:.2f}, p = {res['kruskal']['p']:.3g}")
    for p in res.get("pares", []):
        print(f"  {p['a']} × {p['b']}: Mann-Whitney U = {p['U']:.0f}, p = {p['p']:.3g}, p(Holm) = {p['p_holm']:.3g}, "
              f"efeito (rank-biserial) = {p['efeito_rb']:.2f}")
    print("Nota: cada grânulo é tratado como observação; com várias amostras/lotes por tratamento, compare as médias "
          "por lâmina/réplica (n = réplicas), não os grânulos.")
    if a.output:
        plot_compare(samples, a.col, res, lang=a.lang, out=a.output)
        print(f"gravado: {a.output}")


def cmd_measure(a):
    from sem_film import measure_line
    im = _load(a.file, a)
    x1, y1, x2, y2 = [float(v) for v in a.points.split(",")]
    d = measure_line(im, x1, y1, x2, y2)
    print(f"Distância: {d:.3f} {_u(im)}  ({np.hypot(x2 - x1, y2 - y1):.1f} px; pixel {fmt_len(im.px_um) if im.px_um else '-'})")


def cmd_report(a):
    os.makedirs(a.dir, exist_ok=True)
    import contextlib
    import io
    base = os.path.splitext(os.path.basename(a.files[0]))[0]
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        if a.kind == "granules":
            a.output = os.path.join(a.dir, f"{base}_granulos.png")
            a.csv = os.path.join(a.dir, f"{base}_granulos.csv")
            cmd_info(a)
            cmd_granules(a)
        elif a.kind == "film":
            a.file = a.files[0]
            a.output = os.path.join(a.dir, f"{base}_superficie.png")
            a.csv = os.path.join(a.dir, f"{base}_defeitos.csv")
            cmd_info(argparse.Namespace(**{**vars(a), "files": a.files[:1]}))
            cmd_film(a)
        else:
            a.file = a.files[0]
            a.output = os.path.join(a.dir, f"{base}_espessura.png")
            a.csv = os.path.join(a.dir, f"{base}_espessura.csv")
            cmd_info(argparse.Namespace(**{**vars(a), "files": a.files[:1]}))
            cmd_section(a)
    txt = buf.getvalue()
    print(txt)
    md = os.path.join(a.dir, f"{base}_relatorio.md")
    with open(md, "w", encoding="utf-8") as f:
        f.write(f"# Relatório MEV - {base}\n\n```\n{txt}```\n\n![figura]({os.path.basename(a.output)})\n")
    print(f"gravado: {md}")


# ----------------------------------------------------------------------------- argumentos

def main(argv=None):
    p = argparse.ArgumentParser(prog="sem.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--pixel-size")
    common.add_argument("--scale", help="ex.: 10um=235px")
    common.add_argument("--bar-value", help="valor da barra de escala desenhada, ex.: 10um")
    common.add_argument("--databar", type=int, help="linhas da databar a remover (0 = nenhuma)")
    common.add_argument("--crop", help="x0,y0,x1,y1 (px, após remover a databar)")
    common.add_argument("--no-ocr", action="store_true")
    common.add_argument("--lang", default="pt", choices=["pt", "en"])
    common.add_argument("--sep", default=",")
    common.add_argument("--decimal", default=".")

    s = sub.add_parser("info", parents=[common])
    s.add_argument("files", nargs="+")
    s = sub.add_parser("scale", parents=[common])
    s.add_argument("file")
    s.add_argument("-o", "--output")

    def gran_opts(s):
        s.add_argument("--min-d", default=None, help="diâmetro mínimo (ex.: 2um); padrão 6 px")
        s.add_argument("--max-d", default=None)
        s.add_argument("--thresh", default="otsu", help="otsu | li | triangle | yen | local | valor 0-1")
        s.add_argument("--dark-objects", action="store_true", help="objetos mais escuros que o fundo")
        s.add_argument("--no-split", action="store_true", help="não separar grânulos encostados")
        s.add_argument("--split-h", type=float, default=0.25, help="sensibilidade do watershed (menor = separa mais)")
        s.add_argument("--smooth", type=float, default=1.0)
        s.add_argument("--flatten", action="store_true", help="corrigir gradiente de iluminação antes do limiar")
        s.add_argument("--keep-edge", action="store_true", help="incluir objetos que tocam a borda (não recomendado)")
        s.add_argument("--pits", action="store_true", help="contar pites/poros na superfície dos grânulos")
        s.add_argument("--pit-max", default="1.5um")
        s.add_argument("--numbers", action="store_true")
        s.add_argument("--keep-xy", action="store_true")
        s.add_argument("--overlays", help="pasta para salvar a segmentação de cada imagem")

    s = sub.add_parser("granules", parents=[common])
    s.add_argument("files", nargs="+")
    s.add_argument("-o", "--output")
    s.add_argument("--csv")
    gran_opts(s)

    def film_opts(s):
        s.add_argument("--pore-max", default=None)
        s.add_argument("--particle-max", default=None)
        s.add_argument("--crack-min", default=None)
        s.add_argument("-k", type=float, default=5.0, help="limiar de detecção (× desvio robusto)")
        s.add_argument("--orientation", default="auto", choices=["auto", "h", "v"])
        s.add_argument("--film-darker", action="store_true")

    s = sub.add_parser("film", parents=[common])
    s.add_argument("file")
    s.add_argument("-o", "--output")
    s.add_argument("--csv")
    film_opts(s)
    s = sub.add_parser("section", parents=[common])
    s.add_argument("file")
    s.add_argument("-o", "--output")
    s.add_argument("--csv")
    film_opts(s)
    s = sub.add_parser("panel", parents=[common])
    s.add_argument("files", nargs="+")
    s.add_argument("-o", "--output", required=True)
    s.add_argument("--titles", nargs="+")
    s.add_argument("--ncols", type=int, default=2)
    s.add_argument("--bar")
    s.add_argument("--no-letters", action="store_true")
    s.add_argument("--enhance", action="store_true", help="CLAHE (só visual)")
    s.add_argument("--width", type=float, default=7.0, help="largura da figura (polegadas)")
    s = sub.add_parser("compare", parents=[common])
    s.add_argument("files", nargs="+")
    s.add_argument("--names", nargs="+")
    s.add_argument("--col", default="ecd")
    s.add_argument("-o", "--output")
    s = sub.add_parser("measure", parents=[common])
    s.add_argument("file")
    s.add_argument("--points", required=True)
    s = sub.add_parser("report", parents=[common])
    s.add_argument("files", nargs="+")
    s.add_argument("--kind", required=True, choices=["granules", "film", "section"])
    s.add_argument("-d", "--dir", required=True)
    gran_opts(s)
    film_opts(s)

    a = p.parse_args(argv)
    try:
        {"info": cmd_info, "scale": cmd_scale, "granules": cmd_granules, "film": cmd_film, "section": cmd_section,
         "panel": cmd_panel, "compare": cmd_compare, "measure": cmd_measure, "report": cmd_report}[a.cmd](a)
    except (ValueError, FileNotFoundError, KeyError) as e:
        print(f"ERRO: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
