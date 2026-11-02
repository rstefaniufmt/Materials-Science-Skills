"""
nmr_struct.py - estruturas candidatas: simetria, sinais esperados, predição de deslocamentos (CSP5),
atribuição ótima aos sinais experimentais e ranqueamento (escore e probabilidade relativa estilo DP4)

    from nmr_struct import Candidate, predict, rank
    cands = [Candidate("aspirina", "CC(=O)Oc1ccccc1C(=O)O"), Candidate("isômero", "COC(=O)c1ccccc1O...")]
    res = rank(cands, h1=r1, c13=r13, solvent="cdcl3", formula="C9H8O4")
    res["table"]     # escores por candidata
    res["details"]   # atribuição átomo a átomo

Predição: CSP5 (pip install csp5; modelos inclusos, roda offline; MIT). Sem CSP5 a skill ainda calcula
simetria/sinais esperados/fórmula e confere contagens (DEPT, nH), mas não ranqueia por deslocamento.
IMPORTANTE: o ranqueamento só escolhe ENTRE as candidatas dadas; a estrutura certa precisa estar na
lista. Probabilidade "estilo DP4" é relativa ao conjunto e supõe erros de predição típicos.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy import stats

SIGMA = {"13C": 2.0, "1H": 0.22}      # erro típico (ppm) assumido para a predição
NU = 4                                 # graus de liberdade da t de Student (caudas para erros grandes)
MISS = {"13C": 8.0, "1H": 0.9}         # custo de não atribuir (ppm)


@dataclass
class Candidate:
    name: str
    smiles: str
    mol: object = None
    info: dict = field(default_factory=dict)
    pred: dict = field(default_factory=dict)       # {'13C': {atom_idx: δ}, '1H': {...}}


def _rdkit():
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        from rdkit.Chem import Descriptors, rdMolDescriptors
        return Chem, Descriptors, rdMolDescriptors
    except ImportError as e:
        raise ImportError("RDKit não instalado: pip install rdkit") from e


def mol_info(c: Candidate) -> dict:
    """Fórmula, IDH, classes de equivalência (simetria topológica) de C e H, multiplicidade dos C."""
    Chem, Descriptors, rdm = _rdkit()
    m = Chem.MolFromSmiles(c.smiles)
    if m is None:
        raise ValueError(f"SMILES inválido para {c.name}: {c.smiles}")
    mh = Chem.AddHs(m)
    c.mol = mh
    ranks = list(Chem.CanonicalRankAtoms(mh, breakTies=False))
    cls_c, cls_h = {}, {}
    for a in mh.GetAtoms():
        i = a.GetIdx()
        if a.GetSymbol() == "C":
            cls_c.setdefault(ranks[i], []).append(i)
        elif a.GetSymbol() == "H":
            heavy = a.GetNeighbors()[0]
            cls_h.setdefault(ranks[i], []).append(i)
    carbons = []
    for r, idx in cls_c.items():
        a = mh.GetAtomWithIdx(idx[0])
        nH = sum(1 for nb in a.GetNeighbors() if nb.GetSymbol() == "H")
        carbons.append(dict(cls=r, atoms=idx, n=len(idx), mult={0: "C", 1: "CH", 2: "CH2", 3: "CH3"}.get(nH, "CH4"),
                            label=f"C{min(idx) + 1}"))
    hyd = []
    chiral = len(Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)) > 0
    for r, idx in cls_h.items():
        heavy = mh.GetAtomWithIdx(idx[0]).GetNeighbors()[0]
        exch = heavy.GetSymbol() in ("O", "N", "S")
        lab = f"H-{heavy.GetSymbol()}{heavy.GetIdx() + 1}"
        # CH2 diastereotópico (anel ou molécula quiral): os dois H podem ter δ diferentes
        by_c = {}
        for i in idx:
            by_c.setdefault(mh.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx(), []).append(i)
        split = (not exch and heavy.GetSymbol() == "C" and heavy.GetHybridization() == Chem.HybridizationType.SP3
                 and all(len(v) == 2 for v in by_c.values()) and (chiral or heavy.IsInRing()))
        if split:
            for cidx, hs in by_c.items():
                for k, h in enumerate(hs):
                    hyd.append(dict(cls=r, atoms=[h], n=1, on=cidx, on_sym="C", exch=False, diastereo=True,
                                    label=f"H-C{cidx + 1}{'ab'[k]}"))
            continue
        hyd.append(dict(cls=r, atoms=idx, n=len(idx), on=heavy.GetIdx(), on_sym=heavy.GetSymbol(), exch=exch,
                        label=lab))
    form = rdm.CalcMolFormula(m)
    from nmr_data import dbe
    info = dict(formula=form, dbe=dbe(form), mw=Descriptors.MolWt(m), n_sig_c=len(carbons),
                n_h=sum(h["n"] for h in hyd), n_h_exch=sum(h["n"] for h in hyd if h["exch"]),
                carbons=sorted(carbons, key=lambda d: d["atoms"][0]), hydrogens=sorted(hyd, key=lambda d: d["atoms"][0]),
                mult_counts=pd.Series([cc["mult"] for cc in carbons]).value_counts().to_dict())
    c.info = info
    return info


def predict(cands, solvent=None, nuclei=("13C", "1H"), num_conformers=1):
    """Predição CSP5 para todas as candidatas (em lote). Preenche c.pred[nuc] = {atom_idx: δ}."""
    try:
        from csp5 import predict_smiles
    except ImportError:
        for c in cands:
            c.pred = {}
        return False
    smi = [c.smiles for c in cands]
    for nuc in nuclei:
        r = predict_smiles(smi, nucleus=nuc, solvent=solvent or "unknown", num_conformers=num_conformers,
                           random_seed=0)
        df = r.predictions
        for k, c in enumerate(cands):
            d = df[df.molecule_id == k]
            c.pred[nuc] = dict(zip(d.atom_index.astype(int), d.shift_ppm.astype(float)))
        for f in getattr(r, "failures", []) or []:
            pass
    return True


def predicted_signals(c: Candidate):
    """Sinais previstos por classe de simetria: ¹³C (δ, mult, n) e ¹H (δ, nH, trocável)."""
    pc = c.pred.get("13C", {})
    ph = c.pred.get("1H", {})
    C = []
    for cc in c.info["carbons"]:
        vals = [pc[i] for i in cc["atoms"] if i in pc]
        C.append(dict(label=cc["label"], atoms=cc["atoms"], n=cc["n"], mult=cc["mult"],
                      pred=float(np.mean(vals)) if vals else np.nan))
    H = []
    for hh in c.info["hydrogens"]:
        vals = [ph[i] for i in hh["atoms"] if i in ph]
        H.append(dict(label=hh["label"], atoms=hh["atoms"], n=hh["n"], exch=hh["exch"], on=hh["on"],
                      pred=float(np.mean(vals)) if vals else np.nan))
    return C, H


def _t_logpdf(err, sigma):
    return float(stats.t.logpdf(err / sigma, NU) - math.log(sigma))


def match_13c(C, exp_ppm, exp_mult=None, max_share=2):
    """Atribuição ótima classes previstas -> picos experimentais. Um pico pode receber até max_share
    classes (sobreposição). Retorna lista de dicts e erros."""
    rows = [c for c in C if np.isfinite(c["pred"])]
    E = np.asarray(exp_ppm, float)
    if not rows or E.size == 0:
        return [], {}
    cols = []
    for j in range(E.size):
        for k in range(max_share):
            cols.append((j, k))
    nR = len(rows)
    nC = len(cols) + nR
    cost = np.full((nR, nC), 1e6)
    for i, r in enumerate(rows):
        for c_i, (j, k) in enumerate(cols):
            d = abs(r["pred"] - E[j])
            pen = 0.0
            if k > 0:
                pen += 3.0
            if exp_mult is not None:
                em = exp_mult[j]
                if em and em not in ("?", "") and not _mult_ok(r["mult"], em):
                    pen += 15.0
            cost[i, c_i] = d + pen
        cost[i, len(cols) + i] = MISS["13C"]
    ri, ci = linear_sum_assignment(cost)
    out = []
    for i, cidx in zip(ri, ci):
        r = rows[i]
        if cidx < len(cols):
            j = cols[cidx][0]
            out.append(dict(r, exp=float(E[j]), err=float(E[j] - r["pred"]), exp_idx=int(j),
                            exp_mult=(exp_mult[j] if exp_mult is not None else "")))
        else:
            out.append(dict(r, exp=np.nan, err=np.nan, exp_idx=-1, exp_mult=""))
    used = {o["exp_idx"] for o in out if o["exp_idx"] >= 0}
    extra = [float(E[j]) for j in range(E.size) if j not in used]
    return out, dict(extra=extra)


def _mult_ok(pm, em):
    em = em.replace("₃", "3").replace("₂", "2")
    if em == "CH/CH3":
        return pm in ("CH", "CH3")
    if em == "C/CH2":
        return pm in ("C", "CH2")
    return pm == em


def match_1h(H, exp_table, fill_pen=0.6):
    """Atribuição de GRUPOS de H equivalentes (simetria) aos multipletos experimentais: cada grupo vai
    inteiro para um multipleto (H equivalentes têm um só δ); a integral de cada multipleto é a
    capacidade, com penalidade por H faltando/sobrando. Programação linear inteira (scipy.milp).
    Trocáveis (OH/NH) casam preferencialmente com sinais largos, ou ficam sem atribuição sem custo."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    t = exp_table.reset_index(drop=True)
    G = [g for g in H if np.isfinite(g["pred"])]
    if not G or not len(t):
        return [], {}
    nG, nJ = len(G), len(t)
    cap = np.maximum(np.round(t.nH_int.values.astype(float)), 0)
    C = np.zeros((nG, nJ))
    for i, g in enumerate(G):
        for j, r in t.iterrows():
            if g["exch"]:
                C[i, j] = 0.2 if r.mult == "sl" else 1.0
            elif r.mult == "m" and r.ppm_min <= g["pred"] <= r.ppm_max:
                C[i, j] = 0.0
            else:
                d = abs(g["pred"] - r.ppm)
                if r.mult == "m":
                    d = min(d, abs(g["pred"] - r.ppm_min), abs(g["pred"] - r.ppm_max))
                C[i, j] = d
    miss = np.array([0.3 if g["exch"] else MISS["1H"] for g in G])
    ng = np.array([g["n"] for g in G], float)
    # variáveis: x (nG*nJ), y_miss (nG), under (nJ), over (nJ)
    nx = nG * nJ
    nv = nx + nG + 2 * nJ
    c = np.concatenate([(C * ng[:, None] ** 0.5).ravel(), miss * ng ** 0.5, np.full(nJ, fill_pen), np.full(nJ, fill_pen)])
    A, lb, ub = [], [], []
    for i in range(nG):                           # cada grupo: um multipleto ou "sem atribuição"
        row = np.zeros(nv)
        row[i * nJ:(i + 1) * nJ] = 1
        row[nx + i] = 1
        A.append(row); lb.append(1); ub.append(1)
    for j in range(nJ):                           # Σ n·x − cap = over − under
        row = np.zeros(nv)
        for i in range(nG):
            if not G[i]["exch"] or t.mult[j] == "sl":
                row[i * nJ + j] = ng[i]
        row[nx + nG + j] = 1          # under
        row[nx + nG + nJ + j] = -1    # over
        A.append(row); lb.append(cap[j]); ub.append(cap[j])
    integ = np.concatenate([np.ones(nx + nG), np.zeros(2 * nJ)])
    bnd = Bounds(np.zeros(nv), np.concatenate([np.ones(nx + nG), np.full(2 * nJ, np.inf)]))
    res = milp(c, constraints=LinearConstraint(np.array(A), lb, ub), integrality=integ, bounds=bnd,
               options=dict(time_limit=20))
    if res.x is None:
        return [], {}
    X = np.round(res.x[:nx]).reshape(nG, nJ)
    out = []
    for i, g in enumerate(G):
        js = np.where(X[i] > 0.5)[0]
        if js.size:
            r = t.iloc[js[0]]
            inside = r.mult == "m" and r.ppm_min <= g["pred"] <= r.ppm_max
            err = 0.0 if inside else float(r.ppm - g["pred"])
            if r.mult == "m" and not inside:
                err = float((r.ppm_min if g["pred"] < r.ppm_min else r.ppm_max) - g["pred"])
            out.append(dict(g, exp=float(r.ppm), err=err, exp_rows=[int(js[0])], n_assigned=g["n"],
                            exp_mult=r.mult))
        else:
            out.append(dict(g, exp=np.nan, err=np.nan, exp_rows=[], n_assigned=0, exp_mult=""))
    under = res.x[nx + nG:nx + nG + nJ]
    extra = [float(t.ppm[j]) for j in range(nJ) if under[j] >= 0.5 and t.nH_int[j] >= 1]
    return out, dict(extra=extra, under=float(under.sum()), over=float(res.x[nx + nG + nJ:].sum()))


def score(match_c, match_h, extra_c=None, extra_h=None, h_fill=None):
    """Métricas e log-verossimilhança (t de Student) para uma candidata."""
    ec = np.array([m["err"] for m in match_c if np.isfinite(m["err"])])
    eh = np.array([m["err"] for m in match_h if np.isfinite(m["err"]) and not m["exch"]])
    miss_c = sum(1 for m in match_c if not np.isfinite(m["err"]))
    miss_h = sum(1 for m in match_h if not np.isfinite(m["err"]) and not m["exch"])
    ll = 0.0
    ll += sum(_t_logpdf(e, SIGMA["13C"]) for e in ec)
    ll += sum(_t_logpdf(e, SIGMA["1H"]) for e in eh)
    ll += miss_c * _t_logpdf(MISS["13C"], SIGMA["13C"]) + miss_h * _t_logpdf(MISS["1H"], SIGMA["1H"])
    n_extra_c = len(extra_c or [])
    if h_fill:
        bad = h_fill.get("under", 0) + h_fill.get("over", 0)
        ll += bad * (_t_logpdf(0.6, SIGMA["1H"]) - _t_logpdf(0.0, SIGMA["1H"]))
    ll += 0.5 * n_extra_c * _t_logpdf(MISS["13C"], SIGMA["13C"])        # sinal experimental sem dono
    return dict(MAE_13C=float(np.mean(np.abs(ec))) if ec.size else np.nan,
                max_13C=float(np.max(np.abs(ec))) if ec.size else np.nan,
                MAE_1H=float(np.mean(np.abs(eh))) if eh.size else np.nan,
                max_1H=float(np.max(np.abs(eh))) if eh.size else np.nan,
                C_sem_sinal=miss_c, H_sem_sinal=miss_h, sinais_C_sobrando=n_extra_c,
                sinais_H_sobrando=len(extra_h or []), logL=ll)


def rank(cands, h1=None, c13=None, solvent=None, formula=None, num_conformers=1):
    """Ranqueia candidatas. h1: resultado de analyze_1h; c13: resultado de analyze_13c."""
    for c in cands:
        mol_info(c)
    have_pred = predict(cands, solvent=solvent, num_conformers=num_conformers)
    rows, details = [], {}
    exp_c = exp_m = None
    if c13 is not None and len(c13["table"]):
        t = c13["table"]
        t = t[t.rotulo.eq("")]
        exp_c = t.ppm.values
        exp_m = t["mult"].values if t["mult"].astype(bool).any() else None
    exp_h = None
    if h1 is not None and len(h1["table"]):
        t = h1["table"]
        exp_h = t[t.amostra & (t.nH >= 0.5)].copy()
    for c in cands:
        inf = c.info
        row = dict(candidata=c.name, formula=inf["formula"], IDH=inf["dbe"], sinais_C_previstos=inf["n_sig_c"],
                   H_total=inf["n_h"], H_trocaveis=inf["n_h_exch"])
        if formula:
            from nmr_data import parse_formula
            ok = parse_formula(formula) == parse_formula(inf["formula"])
            row["formula_ok"] = ok
        flags = []
        if exp_c is not None:
            row["sinais_C_obs"] = len(exp_c)
            if len(exp_c) > inf["n_sig_c"]:
                flags.append(f"{len(exp_c) - inf['n_sig_c']} sinal(is) de ¹³C a mais que o previsto pela simetria")
            if exp_m is not None:
                obs = pd.Series(exp_m).value_counts().to_dict()
                pred = inf["mult_counts"]
                diff = {k: obs.get(k, 0) - pred.get(k, 0) for k in set(obs) | set(pred) if obs.get(k, 0) != pred.get(k, 0)}
                if diff:
                    flags.append("DEPT ≠ previsto: " + ", ".join(f"{k} {v:+d}" for k, v in diff.items()))
        if exp_h is not None:
            nexp = float(exp_h.nH_int.sum())
            row["H_obs"] = nexp
            if abs(nexp - inf["n_h"]) > 0.5 and abs(nexp - (inf["n_h"] - inf["n_h_exch"])) > 0.5:
                flags.append(f"integrais somam {nexp:.0f} H (candidata: {inf['n_h']} H, {inf['n_h_exch']} trocáveis)")
        if have_pred:
            C, H = predicted_signals(c)
            mc, xc = match_13c(C, exp_c, exp_m) if exp_c is not None else ([], {})
            mh, xh = match_1h(H, exp_h) if exp_h is not None else ([], {})
            sc = score(mc, mh, xc.get("extra"), xh.get("extra"), xh)
            sc["H_integral_desajuste"] = xh.get("under", 0) + xh.get("over", 0)
            row.update(sc)
            details[c.name] = dict(c13=mc, h1=mh, extra_c=xc.get("extra", []), extra_h=xh.get("extra", []),
                                   pred_C=C, pred_H=H)
        row["alertas"] = "; ".join(flags)
        rows.append(row)
    df = pd.DataFrame(rows)
    if have_pred and len(df) and "logL" in df:
        ll = df.logL.values.astype(float)
        if formula and "formula_ok" in df:
            ll = np.where(df.formula_ok.values, ll, -np.inf)
        if np.isfinite(ll).any():
            p = np.exp(ll - np.nanmax(ll[np.isfinite(ll)]))
            p[~np.isfinite(ll)] = 0
            df["prob_rel_pct"] = 100 * p / p.sum()
        df = df.sort_values("logL", ascending=False).reset_index(drop=True)
    return dict(table=df, details=details, predicted=have_pred, cands={c.name: c for c in cands})


def read_candidates(path_or_list):
    """Lê candidatas de .smi/.csv/.txt (nome e SMILES, em qualquer ordem) ou lista 'nome=SMILES'."""
    out = []
    Chem = _rdkit()[0]
    if isinstance(path_or_list, (list, tuple)):
        for k, item in enumerate(path_or_list, 1):
            n, sm = f"cand{k}", item
            if "=" in item:
                left, right = item.split("=", 1)
                if Chem.MolFromSmiles(right.strip()) is not None and (Chem.MolFromSmiles(item) is None or
                                                                       not any(ch in left for ch in "()[]")):
                    n, sm = left, right
            out.append(Candidate(n.strip(), sm.strip()))
        return out
    import re
    with open(path_or_list, encoding="utf-8") as f:
        for k, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = [p.strip() for p in re.split(r"[,;\t]| {2,}|\s(?=\S+$)", line) if p.strip()]
            smi, name = None, None
            for p in parts:
                if smi is None and Chem.MolFromSmiles(p, sanitize=True) is not None and not p.isalpha():
                    smi = p
                elif name is None:
                    name = p
            if smi is None and len(parts) == 1 and Chem.MolFromSmiles(parts[0]) is not None:
                smi = parts[0]
            if smi is None or smi.lower() in ("smiles",):
                continue
            out.append(Candidate(name or f"cand{len(out) + 1}", smi))
    return out


def draw_assigned(c: Candidate, det: dict, out: str, size=(700, 500)):
    """Estrutura com δC (azul) e δH (vermelho) experimentais atribuídos (ou previstos entre parênteses)."""
    Chem = _rdkit()[0]
    from rdkit.Chem import Draw, AllChem
    from rdkit.Chem.Draw import rdMolDraw2D
    m = Chem.RemoveHs(c.mol)
    notes = {}
    for a in det.get("c13", []):
        txt = f"{a['exp']:.1f}" if np.isfinite(a.get("exp", np.nan)) else f"({a['pred']:.1f})"
        for i in a["atoms"]:
            notes.setdefault(i, []).append(txt)
    hnotes = {}
    for g in det.get("h1", []):
        txt = f"{g['exp']:.2f}" if np.isfinite(g.get("exp", np.nan)) else f"({g['pred']:.2f})"
        hnotes.setdefault(g["on"], txt)
    for a in m.GetAtoms():
        i = a.GetIdx()
        parts = []
        if i in notes:
            parts.append(notes[i][0])
        if i in hnotes:
            parts.append("H " + hnotes[i])
        if parts:
            a.SetProp("atomNote", " / ".join(parts))
    AllChem.Compute2DCoords(m)
    d = rdMolDraw2D.MolDraw2DCairo(*size)
    op = d.drawOptions()
    op.annotationFontScale = 0.55
    op.addAtomIndices = False
    d.DrawMolecule(m)
    d.FinishDrawing()
    with open(out, "wb") as f:
        f.write(d.GetDrawingText())
