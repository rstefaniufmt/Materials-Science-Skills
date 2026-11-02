"""
nmr_io.py - leitura de espectros de RMN 1D (¹H, ¹³C, DEPT, APT) em JCAMP-DX e outros formatos

    from nmr_io import load
    s = load("1h.jdx")            # objeto NMRSpectrum (FID é transformado automaticamente)
    s.ppm, s.y                    # eixo em ppm (decrescente) e intensidade real
    s.nucleus, s.sf, s.solvent    # '1H', MHz, 'cdcl3'
    s.summary()

Formatos:
  JCAMP-DX 4.24/5/6 (.jdx, .dx, .jcamp): XYDATA (X++(Y..Y)) e NTUPLES (Bruker/MestReNova/JEOL/
    Varian/Nanalysis/Oxford/Magritek...), espectro processado OU FID (NMR FID), compressões ASDF,
    blocos LINK (usa o espectro; lê PEAK TABLE se houver)
  Pasta Bruker (fid/ser + acqus, ou pdata/1/1r + procs) e Varian/Agilent (fid + procpar) via nmrglue
  CSV/TXT de duas colunas (ppm, intensidade) ou lista de picos (ppm[, intensidade])

Eixo de deslocamento (ordem de prioridade):
  XUNITS=PPM -> direto
  Bruker ($OFFSET + $SF): ppm = OFFSET − (x₀ − x)/SF
  .SHIFT REFERENCE= ..., ponto, ppm  -> ppm(ponto) fixo
  Hz sem referência (JEOL, MestReNova): ppm = Hz / .OBSERVE FREQUENCY
Depois disso nmr_process.check_reference() confere o sinal residual do solvente.
"""
from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

warnings.filterwarnings("ignore", module="nmrglue")

SOLVENT_ALIASES = {
    "cdcl3": "cdcl3", "chloroform-d": "cdcl3", "chloroform": "cdcl3", "chcl3": "cdcl3",
    "dmso": "dmso-d6", "dmso-d6": "dmso-d6", "d6-dmso": "dmso-d6", "dmsod6": "dmso-d6", "(cd3)2so": "dmso-d6",
    "cd3od": "cd3od", "meod": "cd3od", "methanol-d4": "cd3od", "methanol": "cd3od", "meod-d4": "cd3od",
    "d2o": "d2o", "water": "d2o", "h2o+d2o": "d2o", "h2o": "d2o",
    "acetone": "acetone-d6", "acetone-d6": "acetone-d6", "(cd3)2co": "acetone-d6",
    "c6d6": "c6d6", "benzene-d6": "c6d6", "benzene": "c6d6",
    "cd3cn": "cd3cn", "acetonitrile": "cd3cn", "acetonitrile-d3": "cd3cn",
    "cd2cl2": "cd2cl2", "dichloromethane": "cd2cl2", "methylenechloride-d2": "cd2cl2",
    "thf": "thf-d8", "thf-d8": "thf-d8", "pyridine": "pyridine-d5", "pyridine-d5": "pyridine-d5",
    "toluene": "toluene-d8", "toluene-d8": "toluene-d8", "dmf": "dmf-d7", "dmf-d7": "dmf-d7",
    "tfa": "tfa-d", "neat": "neat",
}


def norm_solvent(s) -> Optional[str]:
    if not s:
        return None
    t = str(s).strip().strip("<>").lower().replace(" ", "")
    t = re.sub(r"\d+vol%$", "", t)
    if t in SOLVENT_ALIASES:
        return SOLVENT_ALIASES[t]
    for k, v in SOLVENT_ALIASES.items():
        if t.startswith(k):
            return v
    return t or None


@dataclass
class NMRSpectrum:
    ppm: np.ndarray
    y: np.ndarray
    nucleus: str = "1H"
    sf: float = 400.0                   # MHz (frequência de referência para Hz <-> ppm)
    solvent: Optional[str] = None
    title: str = ""
    path: str = ""
    kind: str = "spectrum"              # spectrum | peaks
    experiment: str = "1D"              # 1D | DEPT-135 | DEPT-90 | DEPT-45 | DEPTQ | APT | ...
    meta: dict = field(default_factory=dict)
    peaks: Optional[np.ndarray] = None  # lista de picos (ppm, int) se o arquivo trouxer
    yi: Optional[np.ndarray] = None
    warnings: list = field(default_factory=list)

    @property
    def name(self):
        if self.path:
            b = os.path.basename(os.path.normpath(self.path))
            return os.path.splitext(b)[0] if not os.path.isdir(self.path) else b
        return self.title or "espectro"

    @property
    def hz_per_pt(self):
        return float(abs(np.median(np.diff(self.ppm))) * self.sf) if self.ppm.size > 1 else np.nan

    def hz(self, ppm):
        return np.asarray(ppm) * self.sf

    def copy(self, **kw):
        d = dict(self.__dict__)
        d.update(kw)
        d["meta"] = dict(self.meta)
        d["warnings"] = list(self.warnings)
        return NMRSpectrum(**d)

    def summary(self) -> str:
        m = self.meta
        L = [f"Arquivo      : {self.path}",
             f"Título       : {self.title}",
             f"Núcleo       : {self.nucleus}   experimento: {self.experiment}   ({m.get('pulse', '?')})",
             f"Frequência   : {self.sf:.4f} MHz" + (f"  (campo ≈ {self.sf / {'1H': 1, '13C': 0.25145, '19F': 0.94077, '31P': 0.40481}.get(self.nucleus, 1):.0f} MHz ¹H)" if self.nucleus != '1H' else ""),
             f"Solvente     : {self.solvent or '? (informe --solvent)'}",
             f"Origem       : {m.get('origin', '?')}  ({m.get('format', '?')})"]
        if self.kind == "spectrum" and self.ppm.size:
            L.append(f"Faixa        : {self.ppm.max():.2f} a {self.ppm.min():.2f} ppm, {self.ppm.size} pontos, "
                     f"{self.hz_per_pt:.3f} Hz/ponto")
        if m.get("from_fid"):
            L.append(f"FID          : transformado aqui (LB = {m.get('lb')} Hz, zero-filling ×{m.get('zf')}, fase "
                     f"{m.get('phase_desc', '?')})")
        for k, lab in (("ns", "Varreduras"), ("temp", "Temperatura"), ("date", "Data"), ("axis", "Eixo ppm")):
            if m.get(k) not in (None, ""):
                L.append(f"{lab:<13}: {m[k]}")
        if self.peaks is not None:
            L.append(f"Lista de picos no arquivo: {len(self.peaks)} picos")
        for w in self.warnings:
            L.append(f"AVISO: {w}")
        return "\n".join(L)


# ----------------------------------------------------------------------------- utilidades

def _f(v, default=None):
    if v is None:
        return default
    if isinstance(v, (list, tuple)):
        v = v[0] if v else None
        if v is None:
            return default
    s = str(v).split("$$")[0].strip().strip("<>").strip()
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group(0)) if m else default


def _s(v, default=""):
    if v is None:
        return default
    if isinstance(v, (list, tuple)):
        v = v[0] if v else ""
    return str(v).split("$$")[0].strip().strip("<>").strip()


def _nucleus(dic):
    for k in (".OBSERVENUCLEUS", "$NUC1", "$AXNUC", ".NUCLEUS"):
        v = _s(dic.get(k))
        if v and v.lower() != "off":
            v = v.replace("^", "")
            m = re.match(r"(\d+)([A-Za-z]+)", v)
            if m:
                return f"{m.group(1)}{m.group(2).capitalize()}"
    t = _s(dic.get("TITLE")).upper()
    for n in ("13C", "1H", "19F", "31P", "15N"):
        if n in t:
            return n
    return "1H"


def _experiment(dic, title=""):
    p = (_s(dic.get(".PULSESEQUENCE")) or _s(dic.get("$PULPROG")) or "").lower()
    t = (title or "").lower()
    for key, lab in (("deptq", "DEPTQ"), ("dept135", "DEPT-135"), ("dept90", "DEPT-90"), ("dept45", "DEPT-45"),
                     ("jmod", "APT"), ("apt", "APT"), ("uday", "APT")):
        if key in p or key in t.replace(" ", "").replace("-", ""):
            return lab, p
    if "dept" in p or "dept" in t:
        return "DEPT", p
    return "1D", p


def _common_meta(dic, path):
    title = _s(dic.get("TITLE"))
    exp, pulse = _experiment(dic, title)
    meta = dict(origin=_s(dic.get("ORIGIN")), pulse=pulse or "?", ns=_s(dic.get(".AVERAGES")) or _s(dic.get("$NS")),
                temp=_s(dic.get("$TE")) or _s(dic.get("TEMPERATURE")), date=_s(dic.get("LONGDATE")) or _s(dic.get("DATE")))
    solvent = norm_solvent(_s(dic.get(".SOLVENTNAME")) or _s(dic.get("$SOLVENT")))
    return title, exp, meta, solvent


def _sf(dic):
    sf = _f(dic.get("$SF")) or _f(dic.get("$BF1")) or _f(dic.get(".OBSERVEFREQUENCY"))
    return sf


def _shift_reference(dic):
    v = _s(dic.get(".SHIFTREFERENCE"))
    if not v:
        return None
    parts = [p.strip() for p in v.strip("()").split(",")]
    try:
        return float(parts[-2]), float(parts[-1])
    except (ValueError, IndexError):
        return None


def _axis_from_hz(dic, x, n, warnings_):
    """Converte eixo em Hz para ppm conforme o dialeto."""
    sf = _sf(dic)
    if sf is None:
        raise ValueError("arquivo sem frequência do espectrômetro (.OBSERVE FREQUENCY/$SF)")
    off = _f(dic.get("$OFFSET"))
    if off is not None and ("$SF" in dic or "$BF1" in dic):
        return off - (x[0] - x) / sf, sf, "Bruker $OFFSET/$SF"
    ref = _shift_reference(dic)
    if ref is not None:
        pt, val = ref
        i = int(round(pt)) - 1 if pt >= 1 else int(round(pt))
        i = min(max(i, 0), n - 1)
        return val + (x - x[i]) / sf, sf, f".SHIFT REFERENCE (ponto {pt:g} = {val:g} ppm)"
    return x / sf, sf, "Hz / frequência observada"


# ----------------------------------------------------------------------------- JCAMP

def _peak_table_text(txt):
    """Extrai uma PEAK TABLE (XY..XY ou XYW..XYW) do texto bruto."""
    m = re.search(r"##PEAK\s*TABLE=\s*\((XYW?)\.\.(XYW?)\)\s*\n(.*?)(?=\n##)", txt, re.S | re.I)
    if not m:
        return None
    vals = [float(t) for t in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", m.group(3))]
    k = 3 if m.group(1).upper() == "XYW" else 2
    arr = np.array(vals[: len(vals) // k * k]).reshape(-1, k)
    return arr


def read_jcamp(path: str) -> NMRSpectrum:
    import nmrglue as ng
    with open(path, "rb") as fh:
        raw = fh.read()
    txt = raw.decode("utf-8", errors="replace")
    try:
        dic, data = ng.jcampdx.read(path)
    except Exception as e:
        dic, data = {}, None
        warn0 = f"nmrglue não leu ({e}); usando leitor simples"
    else:
        warn0 = None
    warn = [warn0] if warn0 else []
    # bloco LINK: o espectro principal vem em data; parâmetros no topo
    if not dic.get("TITLE") and "_datatype_NMRSPECTRUM" in dic:
        sub = dic["_datatype_NMRSPECTRUM"]
        dic = sub[0] if isinstance(sub, list) else sub
    title, exp, meta, solvent = _common_meta(dic, path)
    nuc = _nucleus(dic)
    dtype = _s(dic.get("DATATYPE")).upper()
    meta["format"] = "JCAMP-DX " + _s(dic.get("JCAMPDX")) + (f" {_s(dic.get('DATACLASS'))}" if dic.get("DATACLASS") else "")
    is_fid = "FID" in dtype
    peaks = _peak_table_text(txt)
    if data is None or (isinstance(data, np.ndarray) and data.size == 0):
        # PEAK TABLE densa (alguns simuladores gravam o espectro assim)
        if peaks is not None and len(peaks) > 2000:
            xu = _s(dic.get("XUNITS")) or (re.search(r"##XUNITS=\s*(.*)", txt).group(1) if "##XUNITS" in txt else "")
            sf = _f(dic.get(".OBSERVEFREQUENCY")) or _f(re.search(r"##\.OBSERVE FREQUENCY=\s*(\S+)", txt).group(1))
            x, y = peaks[:, 0], peaks[:, 1]
            if "HZ" in xu.upper():
                x = x / sf
            o = np.argsort(x)[::-1]
            meta["axis"] = "ppm (PEAK TABLE densa)"
            return NMRSpectrum(ppm=x[o], y=y[o].astype(float), nucleus=nuc, sf=sf or 400.0, solvent=solvent, title=title,
                               path=path, experiment=exp, meta=meta, warnings=warn)
        if peaks is not None:
            sf = _sf(dic) or 400.0
            x = peaks[:, 0]
            if "HZ" in _s(dic.get("XUNITS")).upper():
                x = x / sf
            return NMRSpectrum(ppm=np.array([]), y=np.array([]), nucleus=nuc, sf=sf, solvent=solvent, title=title,
                               path=path, kind="peaks", experiment=exp, meta=meta,
                               peaks=np.column_stack([x, peaks[:, 1]]), warnings=warn)
        raise ValueError("JCAMP sem dados de espectro reconhecíveis")
    if isinstance(data, list):
        re_, im_ = np.asarray(data[0], float), np.asarray(data[1], float) if len(data) > 1 else None
    else:
        re_, im_ = np.asarray(data, float), None
    n = re_.size
    if is_fid:
        from nmr_process import fid_to_spectrum
        fid = re_ + 1j * (im_ if im_ is not None else 0)
        return fid_to_spectrum(fid, dic, nucleus=nuc, title=title, path=path, solvent=solvent, experiment=exp,
                               meta=meta, warnings=warn)
    # eixo
    xunits = (_s(dic.get("XUNITS")) or "").upper()
    first, last = _f(dic.get("FIRSTX")), _f(dic.get("LASTX"))
    if first is None:                                          # NTUPLES
        fl = _s(dic.get("FIRST")).split(",")
        ll = _s(dic.get("LAST")).split(",")
        units = _s(dic.get("UNITS")).split(",")
        try:
            first, last = float(fl[0]), float(ll[0])
            xunits = units[0].strip().upper()
        except (ValueError, IndexError):
            first, last = None, None
    if first is None:
        raise ValueError("não encontrei FIRSTX/LASTX")
    x = np.linspace(first, last, n)
    sf = _sf(dic)
    if "PPM" in xunits:
        ppm, axis = x, "ppm no arquivo"
        sf = sf or 400.0
    else:
        ppm, sf, axis = _axis_from_hz(dic, x, n, warn)
    meta["axis"] = axis
    if ppm[0] < ppm[-1]:
        ppm, re_ = ppm[::-1], re_[::-1]
        im_ = im_[::-1] if im_ is not None else None
    sp = NMRSpectrum(ppm=ppm, y=re_, nucleus=nuc, sf=float(sf), solvent=solvent, title=title, path=path,
                     experiment=exp, meta=meta, yi=im_, warnings=warn)
    if peaks is not None and len(peaks) < 5000:
        px = peaks[:, 0]
        if "PPM" not in _s(dic.get("XUNITS")).upper() and px.max() > 300:
            px = px / sf
        sp.peaks = np.column_stack([px, peaks[:, 1]])
    return sp


# ----------------------------------------------------------------------------- Bruker / Varian / texto

def read_bruker(path: str) -> NMRSpectrum:
    import nmrglue as ng
    p = path
    pdata = None
    if os.path.basename(p) in ("1r", "1i") or os.path.exists(os.path.join(p, "1r")):
        pdata = p if os.path.isdir(p) else os.path.dirname(p)
    elif os.path.isdir(os.path.join(p, "pdata", "1")) and os.path.exists(os.path.join(p, "pdata", "1", "1r")):
        pdata = os.path.join(p, "pdata", "1")
    if pdata:
        dic, data = ng.bruker.read_pdata(pdata)
        udic = ng.bruker.guess_udic(dic, data)
        uc = ng.fileiobase.uc_from_udic(udic)
        ppm = uc.ppm_scale()
        acq = dic.get("acqus", {})
        nuc = re.sub(r"[<>]", "", str(acq.get("NUC1", "1H")))
        exp_dir = os.path.dirname(os.path.dirname(os.path.abspath(pdata)))
        d2 = {"$PULPROG": acq.get("PULPROG", ""), "TITLE": os.path.basename(exp_dir)}
        exp, pulse = _experiment(d2)
        sp = NMRSpectrum(ppm=ppm, y=np.asarray(data, float), nucleus=nuc, sf=float(dic["procs"]["SF"]),
                         solvent=norm_solvent(acq.get("SOLVENT")), title=os.path.basename(exp_dir), path=path,
                         experiment=exp, meta=dict(origin="Bruker", format="Bruker pdata (1r)", pulse=pulse,
                                                   ns=acq.get("NS"), temp=acq.get("TE"), axis="procs OFFSET/SF"))
        return sp
    dic, data = ng.bruker.read(p)
    acq = dic["acqus"]
    d2 = {"$SF": acq.get("BF1"), "$SFO1": acq.get("SFO1"), "$SW_h": acq.get("SW_h"), "$O1": acq.get("O1"),
          "$GRPDLY": acq.get("GRPDLY"), "$DECIM": acq.get("DECIM"), "$DSPFVS": acq.get("DSPFVS"),
          "$PULPROG": acq.get("PULPROG"), "TITLE": os.path.basename(os.path.normpath(p))}
    nuc = re.sub(r"[<>]", "", str(acq.get("NUC1", "1H")))
    exp, pulse = _experiment(d2)
    from nmr_process import fid_to_spectrum
    data = ng.bruker.remove_digital_filter(dic, data)
    d2["$GRPDLY"] = 0
    return fid_to_spectrum(np.asarray(data), d2, nucleus=nuc, title=d2["TITLE"], path=path,
                           solvent=norm_solvent(acq.get("SOLVENT")), experiment=exp,
                           meta=dict(origin="Bruker", format="Bruker fid", pulse=pulse, ns=acq.get("NS")), warnings=[],
                           filter_removed=True)


def read_varian(path: str) -> NMRSpectrum:
    import nmrglue as ng
    d = path if os.path.isdir(path) else os.path.dirname(path)
    dic, data = ng.varian.read(d)
    pp = dic["procpar"]

    def g(k, default=None):
        try:
            return float(pp[k]["values"][0])
        except Exception:
            return default
    sfrq = g("sfrq", 400.0)
    sw = g("sw", 5000.0)
    rfl, rfp = g("rfl", 0.0), g("rfp", 0.0)
    nuc = pp.get("tn", {}).get("values", ["1H"])[0]
    d2 = {"$SF": sfrq, "$SW_h": sw, "TITLE": os.path.basename(os.path.normpath(d)), "_varian": (rfl, rfp, sw, sfrq)}
    from nmr_process import fid_to_spectrum
    seq = pp.get("seqfil", {}).get("values", [""])[0]
    d2["$PULPROG"] = seq
    exp, pulse = _experiment(d2)
    return fid_to_spectrum(np.asarray(data).ravel(), d2, nucleus=nuc, title=d2["TITLE"], path=path,
                           solvent=norm_solvent(pp.get("solvent", {}).get("values", [""])[0]), experiment=exp,
                           meta=dict(origin="Varian/Agilent", format="VnmrJ fid", pulse=pulse), warnings=[],
                           filter_removed=True)


def read_text(path: str, nucleus=None, sf=None) -> NMRSpectrum:
    txt = open(path, encoding="utf-8", errors="replace").read()
    rows = []
    for line in txt.splitlines():
        t = line.strip().replace(";", " ").replace("\t", " ")
        if not t or t[0] in "#%":
            continue
        if re.match(r"^[-+\d.,eE ]+$", t):
            parts = t.split()
            if len(parts) == 1 and t.count(",") == 1:
                parts = t.split(",")
            elif "," in t and " " not in t.strip():
                parts = t.split(",")
            try:
                rows.append([float(p.replace(",", ".")) if p.count(",") == 1 and "." not in p else float(p) for p in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError("arquivo de texto sem números")
    w = min(len(r) for r in rows)
    arr = np.array([r[:w] for r in rows])
    nuc = nucleus or ("13C" if np.nanmax(arr[:, 0]) > 20 else "1H")
    if len(arr) < 400:     # lista de picos
        return NMRSpectrum(ppm=np.array([]), y=np.array([]), nucleus=nuc, sf=sf or 400.0, path=path, kind="peaks",
                           title=os.path.basename(path), meta=dict(format="lista de picos (texto)"),
                           peaks=np.column_stack([arr[:, 0], arr[:, 1] if w > 1 else np.ones(len(arr))]))
    x, y = arr[:, 0], arr[:, 1]
    o = np.argsort(x)[::-1]
    return NMRSpectrum(ppm=x[o], y=y[o], nucleus=nuc, sf=sf or 400.0, path=path, title=os.path.basename(path),
                       meta=dict(format="texto ppm,intensidade", axis="ppm no arquivo"),
                       warnings=[] if sf else ["frequência do espectrômetro desconhecida: informe --sf (MHz) para J em Hz"])


def load(path: str, nucleus=None, sf=None, solvent=None, **kw) -> NMRSpectrum:
    """Carrega um espectro 1D. kw vão para o processamento de FID (lb, zf, phase)."""
    p = path
    if os.path.isdir(p):
        if os.path.exists(os.path.join(p, "procpar")):
            s = read_varian(p)
        else:
            s = read_bruker(p)
    else:
        ext = os.path.splitext(p)[1].lower()
        base = os.path.basename(p)
        if base in ("fid", "ser") and os.path.exists(os.path.join(os.path.dirname(p), "acqus")):
            s = read_bruker(os.path.dirname(p))
        elif base == "fid" and os.path.exists(os.path.join(os.path.dirname(p), "procpar")):
            s = read_varian(p)
        elif base in ("1r",):
            s = read_bruker(p)
        elif ext in (".jdx", ".dx", ".jcamp", ".jcm") or open(p, "rb").read(200).lstrip().startswith(b"##"):
            s = read_jcamp(p)
        else:
            s = read_text(p, nucleus, sf)
    if nucleus:
        s.nucleus = nucleus
    if sf:
        s.sf = float(sf)
    if solvent:
        s.solvent = norm_solvent(solvent)
    return s


def save_csv(s: NMRSpectrum, path: str, decimal="."):
    import pandas as pd
    df = pd.DataFrame({"ppm": s.ppm, "intensidade": s.y})
    df.to_csv(path, index=False, sep=";" if decimal == "," else ",", decimal=decimal, float_format="%.6g")


def save_jcamp(s: NMRSpectrum, path: str):
    """Grava espectro processado em JCAMP-DX 5.01 (XYDATA, AFFN, eixo em ppm)."""
    y = np.asarray(s.y, float)
    yf = (np.abs(y).max() / 2e9) if np.abs(y).max() > 0 else 1.0
    yi = np.round(y / yf).astype(np.int64)
    L = [f"##TITLE={s.title or s.name}", "##JCAMP-DX=5.01", "##DATA TYPE=NMR SPECTRUM", "##DATA CLASS=XYDATA",
         "##ORIGIN=nmr-elucidation skill", f"##.OBSERVE FREQUENCY={s.sf:.6f}", f"##.OBSERVE NUCLEUS=^{s.nucleus}",
         f"##.SOLVENT NAME={s.solvent or ''}", "##XUNITS=PPM", "##YUNITS=ARBITRARY UNITS",
         f"##XFACTOR=1", f"##YFACTOR={yf:.10g}", f"##FIRSTX={s.ppm[0]:.8f}", f"##LASTX={s.ppm[-1]:.8f}",
         f"##NPOINTS={len(y)}", f"##FIRSTY={y[0]:.6g}", "##XYDATA=(X++(Y..Y))"]
    step = (s.ppm[-1] - s.ppm[0]) / (len(y) - 1)
    for i in range(0, len(y), 10):
        L.append(f"{s.ppm[0] + i * step:.8f} " + " ".join(str(v) for v in yi[i:i + 10]))
    L.append("##END=")
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
