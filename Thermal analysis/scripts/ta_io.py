"""
ta_io.py - leitura de dados de análise térmica (TG/TGA, DTG, DSC, DTA, STA)

    from ta_io import load, load_all
    r = load("amostra.txt")            # um Run
    r.T, r.time, r.mass, r.mass_pct, r.hf   # arrays (None se o canal não existe)

Formatos reconhecidos pelo conteúdo (não pela extensão):
    NETZSCH binário .ngb-ss3/.ngb-bs3/.ngb-ds3 (via pacote `pyngb`; pip install pyngb)
    NETZSCH Proteus ASCII (#cabeçalho, ##Temp./°C;Time/min;...)
    TA Instruments Universal Analysis (Nsig/SigN/StartOfData; UTF-16) e TRIOS (.txt [step]; Excel)
    Mettler Toledo STARe (Index t Ts Tr Value + linha de unidades [s] [°C] ...)
    PerkinElmer Pyris (cabeçalho de 2 linhas, blocos por etapa do método)
    Setaram Calisto/SETSOFT (Index Time(s) Sample Temperature(°C) HeatFlow(mW); UTF-16)
    Genérico CSV/TXT/Excel (nomes com unidade "(mg)", "[°C]", "/min" ou linha de unidades; vírgula decimal)

Canais canônicos (sempre nestas unidades):
    time  [min]   T [°C] (amostra)   T_ref [°C] (programa/referência/forno)
    mass  [mg]    mass_pct [%]       hf (DSC/DTA; unidade em r.units['hf']: 'W/g', 'mW', 'µV',
                                         'µV/mg', 'K' (ΔT de DTA))
    dtg   (DTG exportado pelo software, unidade em r.units['dtg'])   segment [int]
Outros canais (vazões de gás, pressão, potência...) ficam em r.extra.

Direção endotérmica (essencial para DSC/DTA): r.meta['endo'] = 'up' | 'down' | None e
r.meta['endo_source'] diz de onde veio. Determinação final em ta_process.endo_direction().
"""
from __future__ import annotations

import io
import json
import os
import re
import struct
import unicodedata
import warnings
import zipfile
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

__all__ = ["Run", "load", "load_all", "detect_format", "save_csv"]

CANON = ("time", "T", "T_ref", "mass", "mass_pct", "hf", "dtg", "segment")


# ============================================================================ modelo
@dataclass
class Run:
    data: Dict[str, np.ndarray]
    units: Dict[str, str] = field(default_factory=dict)
    extra: Dict[str, np.ndarray] = field(default_factory=dict)
    extra_units: Dict[str, str] = field(default_factory=dict)
    title: str = ""
    source: str = ""
    meta: dict = field(default_factory=dict)

    def __getattr__(self, k):
        if k in CANON:
            return self.data.get(k)
        raise AttributeError(k)

    def has(self, k) -> bool:
        return k in self.data and self.data[k] is not None and np.isfinite(self.data[k]).any()

    @property
    def n(self) -> int:
        for v in self.data.values():
            return len(v)
        return 0

    @property
    def kind(self) -> str:
        return self.meta.get("kind", "?")

    @property
    def sample_mass(self) -> Optional[float]:
        m = self.meta.get("sample_mass_mg")
        return float(m) if m not in (None, "") and np.isfinite(m) and m > 0 else None

    def copy(self, **kw) -> "Run":
        d = dict(data={k: v.copy() for k, v in self.data.items()}, units=dict(self.units),
                 extra={k: v.copy() for k, v in self.extra.items()}, extra_units=dict(self.extra_units),
                 title=self.title, source=self.source, meta=json.loads(json.dumps(self.meta, default=str)))
        d.update(kw)
        return Run(**d)

    def subset(self, mask) -> "Run":
        r = self.copy()
        r.data = {k: v[mask] for k, v in self.data.items()}
        r.extra = {k: v[mask] for k, v in self.extra.items()}
        return r

    def summary(self) -> str:
        L = [f"Título      : {self.title or '(sem título)'}",
             f"Arquivo     : {self.source}",
             f"Formato     : {self.meta.get('format', '?')}",
             f"Tipo        : {self.kind}   ({self.n} pontos)"]
        if self.has("T"):
            L.append(f"Temperatura : {np.nanmin(self.T):.1f} → {np.nanmax(self.T):.1f} °C")
        if self.has("time"):
            L.append(f"Tempo       : {np.nanmin(self.time):.2f} → {np.nanmax(self.time):.2f} min")
        if self.sample_mass:
            L.append(f"Massa amostra: {self.sample_mass:.4g} mg")
        if self.has("mass_pct"):
            L.append(f"Massa       : {self.mass_pct[np.isfinite(self.mass_pct)][0]:.2f} → "
                     f"{self.mass_pct[np.isfinite(self.mass_pct)][-1]:.2f} %")
        if self.has("hf"):
            L.append(f"Sinal DSC/DTA: {np.nanmin(self.hf):.4g} … {np.nanmax(self.hf):.4g} {self.units.get('hf', '')}"
                     + (f"   endo {'↑' if self.meta.get('endo') == 'up' else '↓'} ({self.meta.get('endo_source', '')})"
                        if self.meta.get("endo") else "   direção endo: DESCONHECIDA"))
        for k in ("instrument", "operator", "date", "atmosphere", "crucible", "program", "heating_rate"):
            v = self.meta.get(k)
            if v not in (None, "", [], {}):
                s = v if isinstance(v, str) else json.dumps(v, ensure_ascii=False, default=str)
                L.append(f"{k:<12}: {s[:150]}")
        if self.extra:
            L.append("Outros canais: " + ", ".join(f"{k} [{self.extra_units.get(k, '')}]" for k in self.extra))
        for w in self.meta.get("warnings", []):
            L.append(f"AVISO       : {w}")
        return "\n".join(L)


# ============================================================================ texto: utilitários
def _decode_bytes(raw: bytes) -> str:
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    if len(raw) > 40 and raw[1:40:2].count(b"\x00") > 15:      # UTF-16 LE sem BOM
        try:
            return raw.decode("utf-16-le")
        except UnicodeDecodeError:
            pass
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


_NUM = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?$")
_NUMC = re.compile(r"^[+-]?(?:\d+,?\d*|,\d+)(?:[Ee][+-]?\d+)?$")


def _split(line: str, delim: str):
    if delim == "ws":
        return line.split()
    parts = [p.strip().strip('"').strip("'") for p in line.split(delim)]
    return parts


def _tofloat(s: str, dec: str):
    s = s.strip().strip('"')
    if s == "":
        return np.nan
    if dec == ",":
        if not _NUMC.match(s):
            return None
        return float(s.replace(",", "."))
    if not _NUM.match(s):
        if s.lower() in ("nan", "-nan", "inf", "-inf", "na", "n/a", "---"):
            return np.nan
        return None
    return float(s)


def _row(line, delim, dec, strip_prefix=None):
    if strip_prefix and line.startswith(strip_prefix):
        line = line[len(strip_prefix):]
    parts = _split(line.strip("\r\n"), delim)
    if delim != "ws":
        # tab/; inicial (PerkinElmer começa a linha com tab) e final (Cary/Excel) geram célula vazia
        while parts and parts[0] == "":
            parts = parts[1:]
        while parts and parts[-1] == "":
            parts = parts[:-1]
    vals = []
    for p in parts:
        v = _tofloat(p, dec)
        if v is None:
            return None
        vals.append(v)
    return vals if vals else None


def _sniff(lines, start=0):
    """Escolhe (delim, dec, ncol) que dá o maior bloco de linhas numéricas consistentes."""
    sample = lines[start:start + 3000]
    best = None
    for delim in (";", "\t", ",", "ws"):
        for dec in (".", ","):
            if delim == "," and dec == ",":
                continue
            rows = [_row(ln, delim, dec) for ln in sample]
            good = [r for r in rows if r is not None and len(r) >= 2]
            if len(good) < 3:
                continue
            ncol = max(set(len(r) for r in good), key=lambda n: sum(len(r) == n for r in good))
            score = sum(len(r) == ncol for r in good)
            key = (score, ncol, dec == "," and delim in (";", "\t", "ws"))
            if best is None or key > best[0]:
                best = (key, delim, dec, ncol)
    return best


def _data_block(lines, delim, dec, ncol=None, start=0, prefix=None):
    """Índices [i0, i1) do maior bloco contíguo de linhas numéricas (tolera linhas vazias)."""
    best, cur0, cnt = (0, 0, 0), None, 0
    blanks = 0
    for i in range(start, len(lines) + 1):
        ln = lines[i] if i < len(lines) else "x"
        r = _row(ln, delim, dec, prefix) if ln.strip() else "blank"
        if r == "blank":
            blanks += 1
            if blanks <= 1:
                continue
            r = None
        ok = r is not None and r != "blank" and len(r) >= 2 and (ncol is None or len(r) >= ncol - 1)
        if ok:
            blanks = 0
            if cur0 is None:
                cur0, cnt = i, 0
            cnt += 1
        else:
            if cur0 is not None and cnt > best[2]:
                best = (cur0, i, cnt)
            cur0 = None
    return best[0], best[1]


# ============================================================================ nomes e unidades
_DEG = re.compile(r"(?:°|º|›|�|\?|\xb0|deg\.?\s*|grad\s*)\s*([CFK])\b", re.I)


def _norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return s.lower().strip()


def norm_unit(u: str) -> str:
    u = (u or "").strip().strip("()[]{} ").replace("µ", "u").replace("μ", "u")
    u = _DEG.sub(lambda m: "°" + m.group(1).upper(), u)
    u = u.replace(" ", "")
    lu = u.lower()
    table = {"c": "°C", "°c": "°C", "degc": "°C", "k": "K", "s": "s", "sec": "s", "seconds": "s", "min": "min",
             "minutes": "min", "h": "h", "hr": "h", "mg": "mg", "ug": "ug", "g": "g", "%": "%", "wt%": "%",
             "mw": "mW", "w": "W", "uw": "uW", "w/g": "W/g", "mw/mg": "W/g", "mw/g": "mW/g", "j/(g*s)": "W/g",
             "uv": "uV", "uv/mg": "uV/mg", "mv": "mV", "%/min": "%/min", "%/°c": "%/°C", "%/k": "%/°C",
             "mg/min": "mg/min", "mg/°c": "mg/°C", "mg/k": "mg/°C", "%/s": "%/s", "mg/s": "mg/s",
             "ml/min": "mL/min", "bar": "bar", "°c/min": "K/min", "k/min": "K/min", "w/(g·°c)": "W/(g·K)"}
    return table.get(lu, u)


def split_name_unit(label: str):
    """'Temp./°C' -> ('Temp.', '°C'); 'DSC/(mW/mg)' -> ('DSC','mW/mg'); 'Weight (mg)'; 'Ts [°C]';
    'HeatFlow(mW)'; 'Deriv. Weight (%/°C)'; 'dm/dt' (sem unidade)."""
    s = (label or "").strip().strip('"').strip()
    m = re.match(r"^(.*?)\s*[\(\[]([^()\[\]]*(?:\([^()]*\)[^()\[\]]*)?)[\)\]]\s*$", s)
    if m and m.group(1):
        return m.group(1).strip(), norm_unit(m.group(2))
    m = re.match(r"^(.*?)/\s*\((.+)\)\s*$", s)          # DSC/(mW/mg)
    if m:
        return m.group(1).strip(), norm_unit(m.group(2))
    m = re.match(r"^(.*?[A-Za-z.\)])\s*/\s*([^/]+)$", s)   # Temp./°C, Time/min, Mass/%
    if m:
        u = norm_unit(m.group(2))
        if u in ("°C", "K", "s", "min", "h", "mg", "ug", "g", "%", "mW", "W", "uW", "W/g", "uV", "uV/mg", "mL/min",
                 "%/min", "mg/min", "bar", "mV", "°F") or re.search(r"\d", m.group(2)) is None and len(m.group(2)) <= 7:
            if not (m.group(1).lower() in ("dm", "dt") or m.group(2).lower() in ("dt", "dt ")):
                return m.group(1).strip(), u
    return s, ""


def _role(name: str, unit: str) -> Optional[str]:
    n = _norm(name)
    u = unit
    if re.search(r"std\s*dev|stddev|desvio|sigma|setpoint\s*pressure|balance\s*temp|baseline\s*weight", n):
        return None
    if re.search(r"^seg|segment|^step$|^etapa|^stage$", n) and "time" not in n:
        return "segment"
    if re.search(r"deriv|dtg|dm/dt|d\(?m|derivada", n) or u in ("%/min", "%/°C", "mg/min", "mg/°C", "%/s", "mg/s"):
        return "dtg"
    if re.search(r"^(t|time|zeit|tempo|time\s*\(.*)$", n) or (n.startswith("time") and "step" not in n) \
            or (u in ("s", "min", "h") and re.search(r"time|tempo|zeit|^t$", n)):
        return "time"
    if re.search(r"heat\s*flow|heatflow|^hf\b|^dsc|fluxo|warmestrom|wärmestrom|heat\s*flux|^q$", n) \
            or (u in ("mW", "W/g", "uW", "mW/g", "W") and not re.search(r"power|pot[eê]ncia|heater", n)):
        return "hf"
    if re.search(r"^dta|delta\s*t|^dt\b|Δt", n) or (u in ("uV", "uV/mg", "mV") and not re.search(r"temp", n)):
        return "hf"
    if re.search(r"weight|mass|^tg$|^tga|massa|peso|gewicht|masse|^m$", n) and not re.search(r"flow|purge|deriv", n):
        return "mass"
    if re.search(r"temp|^t[sr]?$|^t\s|^ts|^tr$|furnace\s*t|ofen|forno", n) or u in ("°C", "K", "°F"):
        if u not in ("°C", "K", "°F", "") and not re.search(r"temp", n):
            return None
        if re.search(r"program|ref|furnace|ofen|forno|setpoint|^tr$|block|oven|set\s*point", n):
            return "T_ref"
        return "T"
    return None


def _assign(names, units):
    """Mapeia colunas -> canais canônicos. Prioriza a melhor coluna por canal."""
    roles = {}
    score = {}
    for i, (nm, un) in enumerate(zip(names, units)):
        r = _role(nm, un)
        if r is None:
            continue
        nn = _norm(nm)
        s = 1.0
        if r == "T":
            s += 2 * bool(re.search(r"sample|amostra|probe|^ts$", nn)) + 1 * (nn in ("temp", "temp.", "temperature"))
        if r == "time":
            s += 2 * (nn in ("t", "time", "tempo", "zeit")) - 2 * ("step" in nn)
        if r == "mass":
            s += 1 * ("unsubtracted" in nn) + 0.5 * (un == "mg") - 1 * ("change" in nn)
        if r == "hf":
            s += 1 * ("normalized" in nn or un == "W/g") + 1 * bool(re.search(r"heat|dsc", nn))
        key = r
        if r == "mass" and un == "%":
            key = "mass_pct"
        if key not in roles or s > score[key]:
            roles[key], score[key] = i, s
    return roles


def _unit_convert(role, arr, unit):
    a = np.asarray(arr, float)
    u = unit
    if role in ("T", "T_ref"):
        if u == "K" or (u == "" and np.nanmedian(a) > 250 and np.nanmax(a) < 3000 and np.nanmin(a) > 150):
            return a - 273.15, "°C"
        if u == "°F":
            return (a - 32) / 1.8, "°C"
        return a, "°C"
    if role == "time":
        if u == "s":
            return a / 60.0, "min"
        if u == "h":
            return a * 60.0, "min"
        if u == "ms":
            return a / 60000.0, "min"
        return a, "min"
    if role == "mass":
        if u == "ug":
            return a / 1000.0, "mg"
        if u == "g":
            return a * 1000.0, "mg"
        if u == "%":
            return a, "%"
        return a, "mg"
    if role == "hf":
        if u == "W":
            return a * 1000.0, "mW"
        if u == "uW":
            return a / 1000.0, "mW"
        if u == "mW/g":
            return a / 1000.0, "W/g"
        if u == "mV":
            return a * 1000.0, "uV"
        return a, (u or "?")
    return a, u


def _build(names, units, cols, meta, title, source) -> Run:
    names = [n if n else f"col{i}" for i, n in enumerate(names)]
    units = [norm_unit(u) for u in units]
    roles = _assign(names, units)
    data, du, extra, eu = {}, {}, {}, {}
    used = set()
    for key, i in roles.items():
        role = "mass" if key == "mass_pct" else key
        arr, u = _unit_convert(role, cols[i], units[i])
        if key == "mass" and u == "%":
            key = "mass_pct"
        if key == "segment":
            arr = np.round(arr)
        data[key], du[key] = arr, u
        used.add(i)
    for i, (nm, un) in enumerate(zip(names, units)):
        if i in used:
            continue
        k = re.sub(r"[^a-z0-9]+", "_", _norm(nm)).strip("_") or f"col{i}"
        while k in extra:
            k += "_"
        extra[k], eu[k] = np.asarray(cols[i], float), un
    meta.setdefault("warnings", [])
    meta["columns"] = [f"{n} [{u}]" if u else n for n, u in zip(names, units)]
    meta["column_roles"] = {k: names[i] for k, i in roles.items()}
    return _finish(Run(data=data, units=du, extra=extra, extra_units=eu, title=title, source=source, meta=meta))


def _finish(r: Run) -> Run:
    """Completa canais derivados, tipo de ensaio e avisos."""
    m = r.meta
    w = m.setdefault("warnings", [])
    sm = r.sample_mass
    if r.has("mass") and r.units.get("mass") == "mg" and not m.get("is_blank"):
        mm = r.data["mass"]
        first = mm[np.isfinite(mm)][0] if np.isfinite(mm).any() else np.nan
        if m.get("mass_is_change"):             # NGB: variação de massa a partir de 0
            if sm:
                r.data["mass"] = mm + sm
            else:
                w.append("Massa é variação (Δm) e a massa inicial não foi informada: % relativa impossível.")
        mm = r.data["mass"]
        first = mm[np.isfinite(mm)][0]
        if not r.has("mass_pct"):
            ref = sm if (sm and abs(first - sm) < 0.05 * sm) else first
            if sm and abs(first - sm) >= 0.05 * sm:
                w.append(f"Massa inicial medida {first:.4g} mg difere da massa declarada {sm:.4g} mg; "
                         "% calculada sobre o primeiro ponto.")
            r.data["mass_pct"] = 100.0 * mm / ref
            r.units["mass_pct"] = "%"
        if not sm and first > 0:
            m["sample_mass_mg"] = float(first)
    elif r.has("mass_pct") and not r.has("mass") and sm:
        r.data["mass"] = r.data["mass_pct"] / 100.0 * sm
        r.units["mass"] = "mg"
    if r.has("hf"):
        u = r.units.get("hf", "")
        if u == "mW" and r.sample_mass and not m.get("hf_kept_absolute"):
            m["hf_mW_to_Wg_mass"] = r.sample_mass   # conversão feita sob demanda em ta_process
        if u in ("uV", "uV/mg", "K") and m.get("kind_hint") != "DSC":
            m.setdefault("kind_hint", "DTA" if u == "K" or "DTA" in json.dumps(m.get("column_roles", {})) else None)
    has_m, has_h = r.has("mass") or r.has("mass_pct"), r.has("hf")
    hint = m.get("kind_hint")
    if has_m and has_h:
        m["kind"] = "STA"
    elif has_m:
        m["kind"] = "TGA"
    elif has_h:
        m["kind"] = "DTA" if hint == "DTA" or r.units.get("hf") in ("K",) else "DSC"
    else:
        m["kind"] = "?"
        w.append("Nenhum canal de massa ou de fluxo de calor reconhecido; confira os nomes das colunas.")
    if not r.has("T"):
        if r.has("T_ref"):
            r.data["T"], r.units["T"] = r.data["T_ref"], "°C"
            w.append("Sem temperatura da amostra: usando a temperatura de programa/forno.")
        else:
            w.append("Sem coluna de temperatura reconhecida.")
    if not r.has("time"):
        w.append("Sem coluna de tempo: DTG em %/min e integração em tempo usarão a taxa nominal.")
    # direção endo declarada em rótulos
    if has_h and not m.get("endo"):
        txt = " ".join(m.get("columns", [])) + " " + " ".join(str(v) for v in m.get("header", {}).values())[:4000]
        t = _norm(txt)
        if re.search(r"endo\s*up|endo\s*↑|↑\s*endo|exo\s*down|exo\s*↓", t):
            m["endo"], m["endo_source"] = "up", "rótulo no arquivo"
        elif re.search(r"exo\s*up|exo\s*↑|↑\s*exo|\^\s*exo|endo\s*down|endo\s*↓", t):
            m["endo"], m["endo_source"] = "down", "rótulo no arquivo"
        elif m.get("vendor_endo"):
            m["endo"], m["endo_source"] = m["vendor_endo"], f"convenção típica {m.get('vendor', '')} (confirmar)"
    return r


# ============================================================================ leitores por fabricante
def _header_kv(lines, sep_re=r"\t|:\s+|:$|:(?=\S)"):
    kv = {}
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        p = re.split(sep_re, s, maxsplit=1)
        if len(p) == 2 and p[0].strip() and len(p[0]) < 60:
            k = p[0].strip().strip("#").strip()
            if k not in kv:
                kv[k] = p[1].strip()
    return kv


def _num_in(s, default=None):
    if s is None:
        return default
    m = re.search(r"[-+]?\d+(?:[.,]\d+)?(?:[Ee][-+]?\d+)?", str(s))
    return float(m.group(0).replace(",", ".")) if m else default


def read_netzsch_ascii(text: str, path: str) -> List[Run]:
    lines = text.splitlines()
    hdr = {}
    i_col = None
    for i, ln in enumerate(lines):
        if ln.startswith("##"):
            i_col = i
            break
        if ln.startswith("#"):
            k, _, v = ln[1:].partition(":")
            hdr[k.strip()] = v.strip().rstrip(";")
    if i_col is None:
        raise ValueError("cabeçalho de colunas '##' não encontrado")
    delim = {"SEMICOLON": ";", "COMMA": ",", "TAB": "\t"}.get(hdr.get("SEPARATOR", "SEMICOLON").upper(), ";")
    dec = "," if hdr.get("DECIMAL", "POINT").upper() == "COMMA" else "."
    labels = _split(lines[i_col][2:], delim)
    nu = [split_name_unit(x) for x in labels]
    rows = [_row(ln, delim, dec) for ln in lines[i_col + 1:]]
    rows = [r for r in rows if r is not None and len(r) == len(labels)]
    arr = np.array(rows, float)
    cols = [arr[:, j] for j in range(arr.shape[1])]
    names, units = [a for a, _ in nu], [b for _, b in nu]
    meta = {"format": "NETZSCH Proteus ASCII", "vendor": "NETZSCH", "header": hdr,
            "instrument": hdr.get("INSTRUMENT", ""), "operator": hdr.get("OPERATOR", ""),
            "date": hdr.get("DATE/TIME", ""), "crucible": hdr.get("TYPE OF CRUCIBLE", ""),
            "sample_mass_mg": _num_in(hdr.get("SAMPLE MASS /mg")),
            "program": [v for k, v in hdr.items() if k.startswith("SEG.")] or hdr.get("RANGE", ""),
            "atmosphere": "; ".join(v for k, v in hdr.items() if k.startswith(("PURGE", "PROTECTIVE", "ATMOSPHERE")) and v),
            "correction_file": hdr.get("CORR. FILE", ""), "netzsch_exo_flag": hdr.get("EXO"),
            "vendor_endo": "down", "kind_hint": {"DSC": "DSC", "TG": None, "DTA": "DTA", "STA": None}.get(
                hdr.get("MTYPE", "").upper()[:3].strip(), None)}
    rates = [float(x.replace(",", ".")) for x in re.findall(r"([\d.,]+)\s*\(K/min\)", json.dumps(meta["program"]))]
    if rates:
        meta["heating_rate"] = rates
    title = hdr.get("IDENTITY") or hdr.get("SAMPLE") or os.path.splitext(os.path.basename(path))[0]
    return [_build(names, units, cols, meta, title, path)]


def read_ta_ua(text: str, path: str) -> List[Run]:
    lines = text.splitlines()
    i0 = next(i for i, ln in enumerate(lines) if ln.strip() == "StartOfData")
    hdr = _header_kv(lines[:i0], r"\t")
    sigs = []
    for ln in lines[:i0]:
        m = re.match(r"^Sig(\d+)\t(.*)$", ln.strip())
        if m:
            sigs.append(m.group(2).strip())
    nu = [split_name_unit(s) for s in sigs]
    rows = [_row(ln, "\t", ".") for ln in lines[i0 + 1:]]
    rows = [r for r in rows if r is not None and len(r) == len(sigs)]
    arr = np.array(rows, float)
    # linhas-marcador do UA: tempo negativo (-3, -1) no início/fim
    good = arr[:, 0] >= 0 if nu and _role(nu[0][0], nu[0][1]) == "time" else np.ones(len(arr), bool)
    arr = arr[good]
    size = hdr.get("Size", "")
    sm = _num_in(size)
    if sm and "ug" in norm_unit(size.split()[-1] if size.split() else ""):
        sm /= 1000
    xc = [ln.strip().split("\t", 1)[1] for ln in lines[:i0] if ln.startswith("Xcomment\t") and "\t" in ln.strip()]
    prog = [ln.strip().split("\t", 1)[1] for ln in lines[:i0]
            if ln.startswith(("OrgMethod\t", "Method\t")) and "\t" in ln.strip()]
    meta = {"format": "TA Instruments Universal Analysis (export)", "vendor": "TA Instruments", "header": hdr,
            "instrument": hdr.get("Instrument", ""), "operator": hdr.get("Operator", ""),
            "date": f"{hdr.get('Date', '')} {hdr.get('Time', '')}".strip(), "sample_mass_mg": sm,
            "program": prog, "comments": xc, "atmosphere": "; ".join(x for x in xc if re.search(r"gas", x, re.I)),
            "crucible": "; ".join(x for x in xc if re.search(r"pan|crucible|cadinho", x, re.I)),
            "vendor_endo": "down", "kind_hint": "DSC" if "DSC" in hdr.get("Module", "") else None}
    rates = [float(x) for x in re.findall(r"Ramp\s+([\d.]+)\s*°?C/min", " ".join(prog))]
    if rates:
        meta["heating_rate"] = rates
    title = hdr.get("Sample") or os.path.splitext(os.path.basename(path))[0]
    return [_build([a for a, _ in nu], [b for _, b in nu], [arr[:, j] for j in range(arr.shape[1])], meta, title, path)]


def read_trios(text: str, path: str) -> List[Run]:
    lines = text.splitlines()
    first_step = next(i for i, ln in enumerate(lines) if ln.strip().lower() == "[step]")
    hdr = _header_kv(lines[:first_step], r"\t")
    blocks = [i for i, ln in enumerate(lines) if ln.strip().lower() == "[step]"] + [len(lines)]
    allcols, names, units, seg = None, None, None, []
    stepnames = []
    for b0, b1 in zip(blocks[:-1], blocks[1:]):
        blk = lines[b0 + 1:b1]
        # localizar linha de nomes: primeira linha com >= 2 células não numéricas seguida de unidades/dados
        j = 0
        stepname = blk[0].strip() if blk else ""
        while j < len(blk) and (len(blk[j].split("\t")) < 2 or _row(blk[j], "\t", ".", "Data point") is not None):
            j += 1
        if j >= len(blk):
            continue
        nm = [x.strip() for x in blk[j].split("\t")]
        un = [""] * len(nm)
        k = j + 1
        if k < len(blk) and _row(blk[k], "\t", ".", "Data point") is None:
            un = [x.strip() for x in blk[k].split("\t")] + [""] * len(nm)
            un = un[:len(nm)]
            k += 1
        rows = []
        for ln in blk[k:]:
            s = ln[len("Data point"):] if ln.startswith("Data point") else ln
            r = _row(s.lstrip("\t"), "\t", ".")
            if r is None:
                continue
            rows.append((r + [np.nan] * len(nm))[:len(nm)])
        if not rows:
            continue
        arr = np.array(rows, float)
        if names is None:
            names, units, allcols = nm, un, [arr]
        else:
            if nm != names:            # alinhar por nome
                idx = [nm.index(x) if x in nm else None for x in names]
                arr = np.column_stack([arr[:, i] if i is not None else np.full(len(arr), np.nan) for i in idx])
            allcols.append(arr)
        seg += [len(stepnames) + 1] * len(arr)
        stepnames.append(stepname)
    if names is None:
        raise ValueError("nenhum bloco [step] com dados")
    arr = np.vstack(allcols)
    names = list(names) + (["Segment"] if len(stepnames) > 1 else [])
    units = list(units) + ([""] if len(stepnames) > 1 else [])
    cols = [arr[:, j] for j in range(arr.shape[1])] + ([np.array(seg, float)] if len(stepnames) > 1 else [])
    sm = None
    for k, v in hdr.items():
        if re.search(r"^(size|sample\s*mass|mass|weight)\b", k, re.I):
            sm = _num_in(v)
    prog = hdr.get("proceduresegments", "")
    meta = {"format": "TA Instruments TRIOS (export)", "vendor": "TA Instruments", "header": hdr,
            "instrument": hdr.get("Instrument name", ""), "operator": hdr.get("Operator", hdr.get("operator", "")),
            "date": hdr.get("Run date", hdr.get("rundate", "")), "sample_mass_mg": sm, "program": prog,
            "steps": stepnames, "vendor_endo": "down"}
    rates = [float(x) for x in re.findall(r"Ramp\s+([\d.]+)\s*°?C/min", prog)]
    if rates:
        meta["heating_rate"] = rates
    title = hdr.get("Sample name") or hdr.get("Filename") or os.path.splitext(os.path.basename(path))[0]
    return [_build(names, units, cols, meta, title, path)]


def read_mettler(text: str, path: str) -> List[Run]:
    lines = text.splitlines()
    runs = []
    hdr_lines = []
    i = 0
    sample, sm, method = "", None, ""
    for j, ln in enumerate(lines[:200]):
        s = ln.strip()
        if s.startswith("SampleID:") and j + 1 < len(lines):
            v = lines[j + 1].strip()
            m = re.match(r"^(.*?),\s*([\d.,]+)\s*(mg|ug|µg|g)\b", v)
            if m:
                sample, sm = m.group(1).strip(), _num_in(m.group(2))
                if m.group(3) in ("ug", "µg"):
                    sm /= 1000
            else:
                sample = v
        if s.startswith("Method:") and j + 1 < len(lines):
            method = lines[j + 1].strip()
    while i < len(lines):
        if re.match(r"^\s*Index\s+", lines[i]):
            nm = lines[i].split()
            un = []
            k = i + 1
            if k < len(lines) and "[" in lines[k]:
                un = re.findall(r"\[([^\]]*)\]", lines[k])
                k += 1
            un = ([""] * (len(nm) - len(un)) + un) if len(un) < len(nm) else un[:len(nm)]
            rows = []
            while k < len(lines):
                r = _row(lines[k], "ws", ".")
                if r is None or len(r) != len(nm):
                    if lines[k].strip() == "":
                        k += 1
                        if k < len(lines) and _row(lines[k], "ws", ".") is not None:
                            continue
                    break
                rows.append(r)
                k += 1
            if rows:
                arr = np.array(rows, float)
                names = ["Index" if n == "Index" else {"t": "Time", "Ts": "Sample temperature",
                                                      "Tr": "Program temperature"}.get(n, n) for n in nm]
                vunit = norm_unit(un[nm.index("Value")]) if "Value" in nm else ""
                if "Value" in nm:
                    names[nm.index("Value")] = {"mg": "Weight", "%": "Weight", "mW": "Heat flow", "W/g": "Heat flow",
                                                "uV": "DTA"}.get(vunit, "Value")
                meta = {"format": "Mettler Toledo STARe (export)", "vendor": "Mettler Toledo",
                        "sample_mass_mg": sm, "program": method, "vendor_endo": "down",
                        "header": {"Method": method, "Sample": sample}}
                cols = [arr[:, j] for j in range(arr.shape[1])]
                keep = [j for j, n in enumerate(names) if n != "Index"]
                runs.append(_build([names[j] for j in keep], [un[j] for j in keep], [cols[j] for j in keep],
                                   meta, sample or os.path.splitext(os.path.basename(path))[0], path))
            i = k
        else:
            i += 1
    if not runs:
        raise ValueError("bloco 'Index t Ts ...' não encontrado")
    if sm is not None and sm <= 0 or re.search(r"\bblank\b|branco|leer", sample or "", re.I):
        for r in runs:
            r.meta["is_blank"] = True
            r.meta["sample_mass_mg"] = None
            r.data.pop("mass_pct", None)
            r.meta["warnings"].append("Curva de BRANCO (cadinho vazio): use-a para subtrair de uma amostra.")
            _finish(r)
    return merge_same_axis(runs)


def read_perkinelmer(text: str, path: str) -> List[Run]:
    lines = text.splitlines()
    hdr = _header_kv(lines[:60], r":\t|:\s+")
    sm = _num_in(hdr.get("Sample Weight"))
    blocks, i = [], 0
    while i < len(lines) - 2:
        a, b = lines[i], lines[i + 1]
        if re.search(r"\bTime\b", a) and "\t" in a and _row(b, "\t", ".") is None and _row(lines[i + 2], "\t", ".") is not None:
            n1 = [x.strip() for x in a.strip("\n").split("\t")]
            n2 = [x.strip() for x in b.strip("\n").split("\t")] + [""] * len(n1)
            nm = [(x + " " + y).strip() for x, y in zip(n1, n2[:len(n1)])]
            while nm and nm[0] == "":
                nm = nm[1:]
            step = lines[i - 1].strip() if i > 0 else ""
            k, rows = i + 2, []
            while k < len(lines):
                r = _row(lines[k], "\t", ".")
                if r is None:
                    if re.match(r"^\s*\d+\)\s+\S", lines[k]):        # nova etapa do método
                        blocks.append((step, nm, rows))
                        step, rows = lines[k].strip(), []
                        k += 1
                        continue
                    if not lines[k].strip():
                        k += 1
                        continue
                    break
                rows.append(r)
                k += 1
            blocks.append((step, nm, rows))
            i = k
        else:
            i += 1
    if not blocks:
        raise ValueError("tabela PerkinElmer não encontrada")
    nm = blocks[0][1]
    defaults = {"time": "min", "weight": "mg", "temperature": "°C", "temp": "°C", "heat flow": "mW", "flow": "mL/min"}
    units, names = [], []
    for x in nm:
        n, u = split_name_unit(x)
        if not u:
            lx = _norm(x)
            u = next((v for k, v in defaults.items() if k in lx), "")
            if "weight" in lx and "%" in x:
                u = "%"
        names.append(n)
        units.append(u)
    rows, seg = [], []
    for si, (step, _, rr) in enumerate(blocks):
        L = len(nm)
        rows += [(r + [np.nan] * L)[:L] for r in rr]
        seg += [si + 1] * len(rr)
    arr = np.array(rows, float)
    cols = [arr[:, j] for j in range(arr.shape[1])] + [np.array(seg, float)]
    prog = [ln.strip() for ln in lines if re.match(r"^\s*\d+\)\s+(Hold|Heat|Cool|Isothermal)", ln)]
    endo = "up" if any("endo up" in _norm(x) for x in nm) else "up"
    meta = {"format": "PerkinElmer Pyris (export)", "vendor": "PerkinElmer", "header": hdr, "sample_mass_mg": sm,
            "date": hdr.get("Data Collected", ""), "program": prog, "steps": [b[0] for b in blocks],
            "vendor_endo": endo}
    rates = [float(x) for x in re.findall(r"at\s+([\d.]+)\s*°C/min", " ".join(prog))]
    if rates:
        meta["heating_rate"] = rates
    title = hdr.get("Sample ID") or os.path.splitext(os.path.basename(path))[0]
    return [_build(names + ["Segment"], units + [""], cols, meta, title, path)]


def _split_header_ws(line: str, ncol: int):
    """Cabeçalho separado por espaços com nomes que contêm espaços ('Sample Temperature(°C)'):
    funde palavras sem unidade com a palavra seguinte até sobrar o número de colunas dos dados."""
    s = line.strip()
    p = [x for x in re.split(r"\s{2,}|\t", s) if x]
    if len(p) == ncol:
        return p
    tok = s.split()
    has_unit = lambda t: bool(re.search(r"[\(\[].*[\)\]]|/", t))
    idx_like = ("index", "#", "no", "no.", "nr", "point", "n")
    while len(tok) > ncol:
        cand = [i for i in range(len(tok) - 1) if not has_unit(tok[i]) and _norm(tok[i]) not in idx_like]
        pref = [i for i in cand if has_unit(tok[i + 1])]
        if not (pref or cand):
            break
        i = (pref or cand)[0]
        tok[i:i + 2] = [tok[i] + " " + tok[i + 1]]
    return tok


def parse_columns_spec(spec: str):
    """'time:s,T:K,mass:mg' ou 'T,mass_pct' -> [(papel, unidade)]. Papéis: time, T, T_ref, mass,
    mass_pct, hf, dtg, segment, skip."""
    out = []
    for item in spec.split(","):
        role, _, unit = item.strip().partition(":")
        out.append((role.strip(), norm_unit(unit) if unit else ""))
    return out


def read_generic(text: str, path: str, sheet_meta: Optional[dict] = None, columns=None) -> List[Run]:
    lines = text.splitlines()
    best = _sniff(lines)
    if best is None:
        raise ValueError(f"não encontrei tabela numérica em {path}")
    _, delim, dec, ncol = best
    i0, i1 = _data_block(lines, delim, dec, ncol)
    pre = [ln for ln in lines[:i0] if ln.strip()]
    names, units = [f"col{j}" for j in range(ncol)], [""] * ncol
    if pre:
        last = _split(pre[-1].strip(), delim)
        if delim != "ws":
            while last and last[0] == "" and len(last) > ncol:
                last = last[1:]
        prev = _split(pre[-2].strip(), delim) if len(pre) >= 2 else []
        unit_row = len(last) >= 2 and sum(bool(re.match(r"^[\[\(]?[%°ºa-zA-Zµ/·\-\^\d\(\) ]{0,12}[\]\)]?$", c))
                                          and norm_unit(c) in ("°C", "K", "s", "min", "h", "mg", "ug", "g", "%", "mW",
                                                               "W/g", "uV", "uV/mg", "W", "mL/min", "%/min", "%/°C",
                                                               "mg/min", "", "bar")
                                          for c in last) >= len(last) - 1 and any(c.strip() for c in last) \
            and len(prev) >= 2
        if unit_row:
            names = (prev + [""] * ncol)[:ncol]
            units = [norm_unit(c) for c in (last + [""] * ncol)[:ncol]]
        else:
            if delim == "ws" and len(last) != ncol:
                last = _split_header_ws(pre[-1], ncol)
            nu = [split_name_unit(c) for c in (last + [""] * ncol)[:ncol]]
            names, units = [a or f"col{j}" for j, (a, _) in enumerate(nu)], [b for _, b in nu]
    rows = [_row(ln, delim, dec) for ln in lines[i0:i1]]
    rows = [(r + [np.nan] * ncol)[:ncol] for r in rows if r is not None]
    arr = np.array(rows, float)
    hdr = _header_kv(pre[:-1] if pre else [], r"\t|:\s*|;|=")
    sm = None
    for k, v in hdr.items():
        if re.search(r"(sample\s*)?(mass|weight|massa|size|peso|initial mass|einwaage)", k, re.I) and _num_in(v):
            sm = _num_in(v)
            if re.search(r"\bug\b|µg", v):
                sm /= 1000
            if sm and sm > 0:
                break
    meta = {"format": "texto/CSV genérico", "delimiter": "espaços" if delim == "ws" else repr(delim), "decimal": dec,
            "header": hdr, "sample_mass_mg": sm}
    txt = " ".join(pre[:60]).lower()
    if "zone start time" in txt or "setsoft" in txt or "calisto" in txt:
        meta.update(format="Setaram Calisto/SETSOFT (export)", vendor="Setaram")
    elif "shimadzu" in txt or "ta-60" in txt:
        meta.update(format="Shimadzu (export)", vendor="Shimadzu")
    elif "linseis" in txt:
        meta.update(format="Linseis (export)", vendor="Linseis")
    elif "hitachi" in txt or "exstar" in txt or "nexta" in txt:
        meta.update(format="Hitachi/SII (export)", vendor="Hitachi")
    if sheet_meta:
        meta.update(sheet_meta)
    if columns:
        spec = parse_columns_spec(columns) if isinstance(columns, str) else columns
        defaults = {"time": "min", "T": "°C", "T_ref": "°C", "mass": "mg", "mass_pct": "%", "hf": "mW"}
        names = [{"time": "Time", "T": "Temperature", "T_ref": "Program temperature", "mass": "Mass",
                  "mass_pct": "Mass", "hf": "Heat flow", "dtg": "DTG", "segment": "Segment"}.get(r, "skip")
                 for r, _ in spec] + names[len(spec):]
        units = [(u or defaults.get(r, "")) if r != "mass_pct" else "%" for r, u in spec] + units[len(spec):]
        meta["columns_spec"] = columns
    keep = [j for j, n in enumerate(names) if _norm(n) not in ("index", "#", "no", "no.", "point", "nr", "skip")]
    title = os.path.splitext(os.path.basename(path))[0]
    return [_build([names[j] for j in keep], [units[j] for j in keep], [arr[:, j] for j in keep], meta, title, path)]


def read_excel(path: str) -> List[Run]:
    import pandas as pd
    try:
        sheets = pd.read_excel(path, sheet_name=None, header=None, dtype=str)
    except ImportError as e:
        raise ImportError(f"{e}. Para .xls instale xlrd (pip install xlrd); para .xlsx, openpyxl.")
    details = {}
    for name, df in sheets.items():
        if re.search(r"detail|info|param|header", str(name), re.I) and df.shape[1] >= 2:
            for _, row in df.iterrows():
                k, v = str(row.iloc[0]), str(row.iloc[1])
                if k and k != "nan":
                    details[k.strip()] = v.strip()
    runs = []
    for name, df in sheets.items():
        if name in details or re.search(r"detail|info|param", str(name), re.I):
            continue
        df = df.fillna("")
        lines = ["\t".join(str(v).strip() for v in row) for row in df.itertuples(index=False)]
        try:
            got = read_generic("\n".join(lines), path, {"sheet": str(name)})
        except ValueError:
            continue
        for r in got:
            r.meta["format"] = "Excel (" + r.meta["format"] + ")"
            r.meta["header"].update(details)
            for k, v in details.items():
                if re.search(r"^(size|sample\s*mass|sample\s*weight)", k, re.I) and _num_in(v) and not r.sample_mass:
                    r.meta["sample_mass_mg"] = _num_in(v)
                    _finish(r)
            if details.get("Sample name"):
                r.title = details["Sample name"] + (f" [{name}]" if len(sheets) > 2 else "")
            if any("TA" in v or "TRIOS" in v for v in details.values()) or "Instrument name" in details:
                r.meta.update(vendor="TA Instruments", vendor_endo="down")
        runs += got
    if not runs:
        raise ValueError("nenhuma planilha com dados numéricos")
    return concat_segments(runs)


def merge_same_axis(runs: List[Run]) -> List[Run]:
    """Funde curvas do mesmo ensaio (mesmo eixo de tempo) num Run só - p.ex. Mettler TGA/DSC
    exporta massa e fluxo de calor como duas 'Curve' -> um Run STA."""
    if len(runs) < 2:
        return runs
    base = runs[0]
    out = [base]
    for r in runs[1:]:
        if (base.has("time") and r.has("time") and len(base.time) == len(r.time)
                and np.allclose(base.time, r.time, equal_nan=True)):
            for k, v in r.data.items():
                if k not in base.data:
                    base.data[k], base.units[k] = v, r.units.get(k, "")
            for k in ("endo", "endo_source", "vendor_endo"):
                if r.meta.get(k) and not base.meta.get(k):
                    base.meta[k] = r.meta[k]
            base.meta["columns"] = base.meta.get("columns", []) + r.meta.get("columns", [])
            base.meta["column_roles"].update(r.meta.get("column_roles", {}))
            base.meta.pop("kind", None)
            _finish(base)
        else:
            out.append(r)
            base = r
    return out


def concat_segments(runs: List[Run]) -> List[Run]:
    """Concatena runs consecutivos com os mesmos canais e tempo crescente (planilhas/etapas
    exportadas separadamente) num único Run com canal 'segment'."""
    if len(runs) < 2:
        return runs
    groups, cur = [], [runs[0]]
    for r in runs[1:]:
        p = cur[-1]
        ok = set(p.data) == set(r.data) and p.has("time") and r.has("time") and \
            np.nanmin(r.time) >= np.nanmax(p.time) - 1e-6 and p.source == r.source
        if ok:
            cur.append(r)
        else:
            groups.append(cur)
            cur = [r]
    groups.append(cur)
    out = []
    for g in groups:
        if len(g) == 1:
            out.append(g[0])
            continue
        base = g[0].copy()
        for k in base.data:
            base.data[k] = np.concatenate([r.data[k] for r in g])
        base.extra = {k: np.concatenate([r.extra.get(k, np.full(r.n, np.nan)) for r in g]) for k in base.extra}
        base.data["segment"] = np.concatenate([np.full(r.n, i + 1.0) for i, r in enumerate(g)])
        base.units["segment"] = ""
        base.meta["steps"] = [r.meta.get("sheet", "") for r in g]
        if base.has("mass") and base.sample_mass:
            base.data["mass_pct"] = 100 * base.data["mass"] / base.sample_mass
        out.append(base)
    return out


# ============================================================================ NETZSCH binário (.ngb-*)
def _netzsch_sensitivity(T, c):
    z = (np.asarray(T, float) - c["p0"]) / c["p1"]
    return (c["p2"] + c["p3"] * z + c["p4"] * z ** 2 + c["p5"] * z ** 3) * np.exp(-z ** 2)


def read_ngb(path: str, baseline: Optional[str] = None, run: str = "auto", calibrate: bool = True) -> List[Run]:
    """NETZSCH Proteus .ngb-ss3 (amostra), .ngb-bs3 (correção), .ngb-ds3 (amostra+correção).
    Usa o pacote `pyngb` (MIT; https://github.com/GraysonBellamy/pyngb).
    - baseline: arquivo de correção (.ngb-bs3) para subtrair (mesmo programa de temperatura);
    - run='auto': em .ngb-ds3 devolve a curva corrigida (como o Proteus mostra);
    - calibrate: converte DSC µV -> mW pela sensibilidade gravada no arquivo."""
    try:
        import pyngb
    except ImportError:
        raise ImportError("Arquivos NETZSCH .ngb-* exigem o pacote pyngb: pip install pyngb "
                          "(ou exporte como ASCII no Proteus: Extras > Export data).")
    kw = {}
    if baseline:
        kw["baseline_file"] = baseline
    ext = os.path.splitext(path)[1].lower()
    if run == "auto":
        run = "corrected" if ext in (".ngb-ds3",) and not baseline else "sample"
    if run != "sample":
        kw["run"] = run
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        md, t = pyngb.read_ngb(path, return_metadata=True, **kw)
    cols = {c: t[c].to_numpy() for c in t.column_names}
    units = {c: (t.schema.field(c).metadata or {}).get(b"units", b"").decode() for c in t.column_names}
    w = []
    data, du, extra, eu = {}, {}, {}, {}
    if "time" in cols:
        data["time"], du["time"] = _unit_convert("time", cols["time"], norm_unit(units["time"]))
    if "sample_temperature" in cols:
        data["T"], du["T"] = cols["sample_temperature"].astype(float), "°C"
    if "furnace_temperature" in cols:
        data["T_ref"], du["T_ref"] = cols["furnace_temperature"].astype(float), "°C"
    sm = md.get("sample_mass")
    sm = float(sm) if sm is not None and sm > 0 else None
    if "mass" in cols and np.nanmax(np.abs(cols["mass"])) > 0:
        data["mass"], du["mass"] = cols["mass"].astype(float), "mg"
    endo, endo_src = "down", "convenção NETZSCH (calibração com entalpia de fusão negativa)"
    if "dsc_signal" in cols and np.nanmax(np.abs(cols["dsc_signal"])) > 0:
        sig = cols["dsc_signal"].astype(float)
        cal = md.get("calibration_constants") or {}
        sens = md.get("sensitivity_calibration") or {}
        fps = sens.get("fixpoints", []) if isinstance(sens, dict) else []
        if fps:
            areas = [fp.get("peak_area") for fp in fps if fp.get("peak_area") is not None]
            ents = [fp.get("enthalpy") for fp in fps if fp.get("enthalpy") is not None]
            if areas and ents and np.median(ents) < 0:
                endo = "down" if np.median(areas) < 0 else "up"
                endo_src = "sinal das áreas de fusão dos padrões de calibração do próprio arquivo"
        if calibrate and all(k in cal for k in ("p0", "p1", "p2", "p3", "p4", "p5")) and cal["p1"]:
            S = _netzsch_sensitivity(data["T"], cal)
            ok = S > 0.05 * np.nanmax(S)
            mw = np.where(ok, sig / np.where(ok, S, 1), np.nan)
            data["hf"], du["hf"] = mw, "mW"
            if (~ok).any():
                w.append(f"Sensibilidade DSC fora da faixa calibrada em {int((~ok).sum())} pontos (NaN).")
            extra["dsc_uV"], eu["dsc_uV"] = sig, "uV"
        else:
            data["hf"], du["hf"] = sig, "uV"
            w.append("DSC em µV (sem calibração de sensibilidade aplicada): áreas não são entalpias.")
    for c in cols:
        if c not in ("time", "sample_temperature", "furnace_temperature", "mass", "dsc_signal"):
            extra[c], eu[c] = cols[c].astype(float), norm_unit(units.get(c, ""))
    prog = md.get("temperature_program") or {}
    segs, rates = [], []
    for k in sorted(prog, key=lambda s: int(re.sub(r"\D", "", s) or 0)):
        st = prog[k]
        if st.get("stage_type") == 1:
            segs.append(f"{st.get('temperature'):.0f} °C, {st.get('heating_rate'):.1f} K/min, {st.get('time', 0) / 60:.1f} min")
            if st.get("heating_rate"):
                rates.append(float(st["heating_rate"]))
    gases = [f"{md.get(f'{g}_mfc_gas_formula', md.get(f'{g}_mfc_gas', ''))} {md.get(f'{g}_mfc_flow', 0):g} mL/min"
             for g in ("purge_1", "purge_2", "protective") if md.get(f"{g}_mfc_flow")]
    meta = {"format": f"NETZSCH binário {ext} (pyngb)", "vendor": "NETZSCH",
            "instrument": md.get("instrument_model") or md.get("instrument", ""), "operator": md.get("operator", ""),
            "date": md.get("date_performed", ""), "crucible": md.get("crucible_type", ""),
            "sample_mass_mg": sm, "program": segs, "heating_rate": rates, "atmosphere": "; ".join(gases),
            "measurement_type": md.get("measurement_type"), "run": run, "baseline_file": baseline or "",
            "correction_file": md.get("correction_file_path", ""), "mass_is_change": True,
            "endo": endo if "hf" in data else None, "endo_source": endo_src, "warnings": w,
            "header": {k: v for k, v in md.items() if isinstance(v, (str, int, float))}}
    if md.get("measurement_type") == "correction" or ext == ".ngb-bs3":
        w.append("Arquivo de CORREÇÃO (linha de base/branco), não de amostra.")
        meta["mass_is_change"] = True
        if not sm:
            meta["sample_mass_mg"] = None
    if baseline:
        w.append(f"Correção subtraída: {os.path.basename(baseline)}")
    elif run == "sample" and md.get("measurement_type") == "sample" and "mass" in data:
        cf = md.get("correction_file_path") or ""
        w.append("Curva SEM correção de empuxo/linha de base (arquivo de amostra puro). Para resultados "
                 "quantitativos passe o arquivo de correção (.ngb-bs3)" + (f": {cf.split(chr(92))[-1]}" if cf else "") + ".")
    title = md.get("sample_name") or md.get("sample_id") or os.path.splitext(os.path.basename(path))[0]
    r = Run(data=data, units=du, extra=extra, extra_units=eu, title=str(title), source=path, meta=meta)
    if "mass" in data and not sm:
        r.extra["mass_change"] = r.data.pop("mass")
        r.extra_units["mass_change"] = "mg"
        r.meta["is_blank"] = r.meta.get("measurement_type") == "correction"
        w.append("Sem massa inicial válida: só a variação de massa (mg) está disponível.")
    return [_finish(r)]


# ============================================================================ detecção e entrada
def detect_format(path: str) -> str:
    with open(path, "rb") as fh:
        head = fh.read(65536)
    ext = os.path.splitext(path)[1].lower()
    if head[:2] == b"PK":
        try:
            with zipfile.ZipFile(path) as z:
                if any(n.startswith("Streams/stream_") for n in z.namelist()):
                    return "ngb"
        except zipfile.BadZipFile:
            pass
        return "excel"
    if head[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" or ext in (".xls", ".xlsx", ".xlsm"):
        return "excel"
    if ext.startswith(".ngb"):
        return "ngb"
    if ext in (".tri", ".001", ".002", ".003") and b"StartOfData" not in head:
        raise ValueError(f"Binário TA ({ext}) não é legível fora do TRIOS/Universal Analysis: exporte como texto "
                         "(TRIOS: File > Export > 'Export to LIMS'/txt; UA: File > Export Data File > 'Signals Only').")
    if ext in (".stx", ".dsc", ".tga", ".tad", ".pyr", ".tg8d", ".dsc8d"):
        raise ValueError(f"Formato binário {ext} não suportado: exporte como TXT/CSV no software do equipamento.")
    t = _decode_bytes(head)
    if re.search(r"^#[A-Z .]+:", t, re.M) and re.search(r"^##", t, re.M):
        return "netzsch"
    if re.search(r"^StartOfData\s*$", t, re.M) or (re.search(r"^Nsig\t", t, re.M) and re.search(r"^Sig1\t", t, re.M)):
        return "ta_ua"
    if re.search(r"^\[step\]\s*$", t, re.M | re.I):
        return "trios"
    if re.search(r"STARe|^\s*Index\s+t\s+Ts", t, re.M):
        return "mettler"
    if re.search(r"PerkinElmer|Pyris|Unsubtracted", t):
        return "perkinelmer"
    return "generic"


def load_all(path: str, columns=None, **kw) -> List[Run]:
    """columns: força os papéis das colunas em arquivos sem cabeçalho, ex. 'time:min,mass:mg'."""
    fmt = detect_format(path)
    if columns:
        with open(path, "rb") as fh:
            return read_generic(_decode_bytes(fh.read()), path, columns=columns)
    if fmt == "ngb":
        return read_ngb(path, **{k: v for k, v in kw.items() if k in ("baseline", "run", "calibrate")})
    if fmt == "excel":
        return read_excel(path)
    with open(path, "rb") as fh:
        text = _decode_bytes(fh.read())
    reader = {"netzsch": read_netzsch_ascii, "ta_ua": read_ta_ua, "trios": read_trios, "mettler": read_mettler,
              "perkinelmer": read_perkinelmer, "generic": read_generic}[fmt]
    try:
        return reader(text, path)
    except (ValueError, StopIteration, IndexError) as e:
        if fmt != "generic":
            warnings.warn(f"leitor {fmt} falhou ({e}); tentando leitor genérico")
            return read_generic(text, path)
        raise


def load(path: str, index: int = 0, **kw) -> Run:
    runs = load_all(path, **kw)
    if len(runs) > 1:
        warnings.warn(f"{path}: {len(runs)} curvas no arquivo; usando a de índice {index}.")
    return runs[index]


def save_csv(r: Run, path: str, decimal: str = ".", extra: bool = False):
    sep = ";" if decimal == "," else ","
    keys = [k for k in ("time", "T", "T_ref", "mass", "mass_pct", "hf", "dtg", "segment") if r.has(k)]
    hdr = [f"{ {'time': 'Time', 'T': 'Temperature', 'T_ref': 'Program temperature', 'mass': 'Mass',
                  'mass_pct': 'Mass', 'hf': 'Heat flow' if r.kind != 'DTA' else 'DTA', 'dtg': 'DTG',
                  'segment': 'Segment'}[k]} ({r.units.get(k, '')})" for k in keys]
    cols = [r.data[k] for k in keys]
    if extra:
        for k, v in r.extra.items():
            hdr.append(f"{k} ({r.extra_units.get(k, '')})")
            cols.append(v)
    buf = io.StringIO()
    buf.write(sep.join(hdr) + "\n")
    for row in zip(*cols):
        cells = [f"{v:.6g}" for v in row]
        if decimal == ",":
            cells = [c.replace(".", ",") for c in cells]
        buf.write(sep.join(cells) + "\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(buf.getvalue())
