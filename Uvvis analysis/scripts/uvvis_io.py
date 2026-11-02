"""
uvvis_io.py - leitura de espectros UV-Vis(-NIR)

Formatos: JCAMP-DX (.jdx/.dx), CSV/TXT/PRN/ASC (inclusive vírgula decimal, exportações
Cary pareadas, Shimadzu UVProbe, Ocean Insight/OceanView), Excel (.xlsx/.xls),
Galactic GRAMS .spc (novo 0x4B e antigo 0x4D), PerkinElmer .sp (Lambda / Spectrum).

    from uvvis_io import load, load_all
    s = load("amostra.csv")        # um Spectrum
    specs = load_all("cary.csv")   # todos os espectros do arquivo

Todo Spectrum sai com x em COMPRIMENTO DE ONDA (nm), na ordem do arquivo, e `ykind`:
    'A'  absorbância        'T' / '%T'  transmitância      'R' / '%R'  reflectância
    'KM' Kubelka-Munk F(R)  'eps' absortividade molar       'logeps' log10(ε)
    'alpha' coeficiente de absorção (cm-1)   'counts' contagens brutas   'arb' desconhecido
`meta['ykind_guessed'] = True` quando o tipo foi inferido pela forma (não declarado):
confirme com o usuário - %T e %R têm a mesma faixa numérica.

Dependências: numpy (pandas+openpyxl só para Excel).
"""
from __future__ import annotations

import io
import os
import re
import struct
import warnings
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

__all__ = ["Spectrum", "load", "load_all", "read_jcamp", "read_csv", "read_spc", "read_sp",
           "read_excel", "save_csv", "HC_EV_NM"]

HC_EV_NM = 1239.84193  # h*c em eV·nm  (E[eV] = 1239.84 / λ[nm])


@dataclass
class Spectrum:
    x: np.ndarray                      # nm
    y: np.ndarray
    ykind: str = "arb"
    title: str = ""
    source: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def npoints(self) -> int:
        return int(self.x.size)

    @property
    def step(self) -> float:
        return float(np.median(np.abs(np.diff(self.x)))) if self.x.size > 1 else float("nan")

    @property
    def xrange(self):
        """(mín, máx) em nm."""
        return float(np.nanmin(self.x)), float(np.nanmax(self.x))

    @property
    def energy(self) -> np.ndarray:
        """Energia do fóton em eV para cada ponto."""
        return HC_EV_NM / self.x

    def copy(self, **kw) -> "Spectrum":
        d = dict(x=self.x.copy(), y=self.y.copy(), ykind=self.ykind, title=self.title,
                 source=self.source, meta=dict(self.meta))
        d.update(kw)
        return Spectrum(**d)

    def at(self, wl: float) -> float:
        """Valor interpolado em λ (nm)."""
        o = np.argsort(self.x)
        return float(np.interp(wl, self.x[o], self.y[o]))

    # ---- conversões de ordenada -------------------------------------------
    def to_absorbance(self, floor: float = 1e-6) -> "Spectrum":
        """T/%T -> A = -log10(T). R/%R -> pseudo-absorbância log10(1/R) (DRS).
        A, KM, eps, alpha, logeps, counts, arb passam inalterados."""
        k = self.ykind
        if k not in ("T", "%T", "R", "%R"):
            return self.copy()
        t = self.y / 100.0 if k in ("%T", "%R") else self.y.astype(float)
        nbad = int(np.sum(t <= 0))
        if nbad:
            warnings.warn(f"{nbad} pontos com T/R <= 0; recortados em {floor}.")
        out = self.copy(y=-np.log10(np.clip(t, floor, None)), ykind="A")
        out.meta["converted_from"] = k + (" (log 1/R)" if k in ("R", "%R") else "")
        return out

    def to_transmittance(self, percent: bool = True) -> "Spectrum":
        k = self.ykind
        if k in ("%T", "T"):
            t = self.y / 100.0 if k == "%T" else self.y
        elif k == "A":
            t = 10.0 ** (-self.y)
        else:
            raise ValueError(f"Não sei converter ykind={k!r} em transmitância.")
        out = self.copy(y=t * (100.0 if percent else 1.0), ykind="%T" if percent else "T")
        out.meta["converted_from"] = k
        return out

    def to_km(self, floor: float = 1e-4) -> "Spectrum":
        """Reflectância difusa -> Kubelka-Munk F(R) = (1-R)^2 / 2R  (R fracionária,
        medida contra padrão branco - BaSO4/Spectralon)."""
        k = self.ykind
        if k == "KM":
            return self.copy()
        if k not in ("R", "%R"):
            raise ValueError(f"Kubelka-Munk exige reflectância (R ou %R); ykind={k!r}. "
                             "Se o arquivo é %R mas foi lido como %T, ajuste s.ykind='%R'.")
        r = self.y / 100.0 if k == "%R" else self.y.astype(float)
        r = np.clip(r, floor, 1.0)
        out = self.copy(y=(1 - r) ** 2 / (2 * r), ykind="KM")
        out.meta["converted_from"] = k
        return out

    def sorted(self, ascending: bool = True) -> "Spectrum":
        idx = np.argsort(self.x)
        if not ascending:
            idx = idx[::-1]
        return self.copy(x=self.x[idx], y=self.y[idx])

    def crop(self, lo: float, hi: float) -> "Spectrum":
        lo, hi = min(lo, hi), max(lo, hi)
        m = (self.x >= lo) & (self.x <= hi)
        return self.copy(x=self.x[m], y=self.y[m])

    def summary(self) -> str:
        lo, hi = self.xrange
        lines = [
            f"Título      : {self.title or '(sem título)'}",
            f"Arquivo     : {self.source}",
            f"Formato     : {self.meta.get('format', '?')}",
            f"Pontos      : {self.npoints}   intervalo ~{self.step:.4g} nm",
            f"Faixa       : {lo:.1f} - {hi:.1f} nm  ({HC_EV_NM / hi:.2f} - {HC_EV_NM / lo:.2f} eV)",
            f"Ordenada    : {self.ykind}" + ("  (INFERIDO - confirmar)" if self.meta.get("ykind_guessed") else ""),
            f"y min/max   : {np.nanmin(self.y):.4g} / {np.nanmax(self.y):.4g}",
        ]
        for k in ("instrument", "date", "origin", "owner", "comment", "integration_time",
                  "xunits_original", "yunits_original", "x_label", "y_label", "column", "sheet"):
            if self.meta.get(k) not in (None, ""):
                lines.append(f"{k:<12}: {self.meta[k]}")
        for w in self.meta.get("warnings", []):
            lines.append(f"AVISO       : {w}")
        return "\n".join(lines)


# --------------------------------------------------------------------------
# Unidades
# --------------------------------------------------------------------------
def _x_to_nm(x: np.ndarray, unit: str):
    u = (unit or "").strip().upper().replace(" ", "")
    x = np.asarray(x, float)
    if u in ("", "NM", "NANOMETERS", "NANOMETER", "WAVELENGTH(NM)", "WAVELENGTHNM", "WAVELENGTH"):
        return x, "nm"
    if u in ("1/CM", "CM-1", "CM^-1", "WAVENUMBER", "WAVENUMBERS", "CM−1"):
        if np.nanmin(x) <= 0:   # interferograma / eixo não espectral: não converte
            return x, "cm-1 (não convertido: valores <= 0)"
        return 1e7 / x, "cm-1"
    if u in ("MICROMETERS", "MICROMETER", "UM", "µM", "MICRONS", "MICRON"):
        return x * 1000.0, "um"
    if u in ("EV", "ELECTRONVOLTS", "ELECTRONVOLT"):
        return HC_EV_NM / x, "eV"
    if u in ("A", "ANGSTROM", "ANGSTROMS", "Å"):
        return x / 10.0, "angstrom"
    if "NM" in u or "NANOM" in u:
        return x, "nm"
    if "CM" in u:
        return 1e7 / x, "cm-1"
    return x, u.lower() or "nm"


def _ykind_from_label(label: str) -> Optional[str]:
    if not label:
        return None
    u = label.strip().upper()
    u2 = u.replace(" ", "").replace(".", "")
    if "LOG" in u and ("EPS" in u or "ABSORPTIVITY" in u or "Ε" in u):
        return "logeps"
    if "KUBELKA" in u or u2 in ("KM", "K-M", "F(R)", "FR", "F(R∞)"):
        return "KM"
    if "%T" in u2 or "T%" in u2 or ("TRANSMIT" in u and "%" in u):
        return "%T"
    if "TRANSMIT" in u:
        return "T"
    if "%R" in u2 or "R%" in u2 or ("REFLECT" in u and "%" in u):
        return "%R"
    if "REFLECT" in u:
        return "R"
    if "LOG(1/R)" in u2 or "LOG1/R" in u2:
        return "A"
    if "ABSORB" in u or u2 in ("A", "ABS", "AU", "OD") or re.search(r"\bABS\b|\(A\)|\bA\.U\.|\bOD\b", u):
        return "A"
    if "EPSILON" in u or "ABSORPTIVITY" in u or "EXTINCTION" in u or "M-1CM-1" in u2 or "LMOL" in u2:
        return "eps"
    if "ALPHA" in u or "ABSORPTIONCOEF" in u2:
        return "alpha"
    if "COUNT" in u or "INTENSITY" in u or "SCOPE" in u or "SAMPLE" == u2 or "ENERGY" in u:
        return "counts"
    return None


def guess_ykind(y: np.ndarray) -> str:
    """Inferência pela faixa numérica. Contagens: valores > 150. Percentual (T ou R): 0-110
    com mediana alta. Fração 0-1.05 com mediana alta: T. Resto pequeno: absorbância."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return "arb"
    lo, hi, med = np.percentile(y, 1), np.percentile(y, 99), np.median(y)
    if hi > 150:
        return "counts"
    if 5 < hi <= 115 and lo > -5 and med > 3:
        return "%T"            # pode ser %R - o leitor marca como inferido
    if hi <= 1.05 and lo >= -0.02 and med > 0.3 and (hi - med) < (med - lo):
        return "T"
    if hi < 8:
        return "A"
    return "arb"


def _finish(spec: Spectrum, declared: Optional[str]) -> Spectrum:
    w = spec.meta.setdefault("warnings", [])
    if declared is None:
        spec.ykind = guess_ykind(spec.y)
        spec.meta["ykind_guessed"] = True
        if spec.ykind == "%T":
            w.append("Ordenada percentual sem rótulo: assumida %T. Se for reflectância difusa, use %R.")
    else:
        spec.ykind = declared
        mx = np.nanmax(spec.y)
        if declared == "T" and mx > 2.0:
            spec.ykind = "%T"
        elif declared == "%T" and mx <= 1.5:
            spec.ykind = "T"
        elif declared == "R" and mx > 2.0:
            spec.ykind = "%R"
        elif declared == "%R" and mx <= 1.5:
            spec.ykind = "R"
    if spec.ykind == "A":
        mx = float(np.nanmax(spec.y))
        if mx > 3.0:
            w.append(f"Absorbância máxima {mx:.2f}: acima de ~2.5-3 a maioria dos espectrofotômetros "
                     "perde linearidade (luz espúria) - dilua ou use caminho óptico menor.")
    if spec.ykind in ("T", "%T"):
        tmin = np.nanmin(spec.y) / (100.0 if spec.ykind == "%T" else 1.0)
        if tmin < 0.003:
            w.append(f"Transmitância mínima {tmin * 100:.2g}% - região saturada.")
    return spec


# --------------------------------------------------------------------------
# JCAMP-DX (mesmo decodificador validado na skill de FT-IR)
# --------------------------------------------------------------------------
_SQZ = {c: (str(i), 1) for i, c in enumerate("@ABCDEFGHI")}
_SQZ.update({c: ("-" + str(i), 1) for i, c in enumerate("abcdefghi", start=1)})
_DIF = {c: (str(i), 2) for i, c in enumerate("%JKLMNOPQR")}
_DIF.update({c: ("-" + str(i), 2) for i, c in enumerate("jklmnopqr", start=1)})
_DUP = {c: (str(i), 3) for i, c in enumerate("STUVWXYZs", start=1)}
_ASDF = {}
_ASDF.update(_SQZ)
_ASDF.update(_DIF)
_ASDF.update(_DUP)
_AFFN_RE = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[Ee][+-]?\d+)?")
_AFFN_ONLY = re.compile(r"^[\d\s.,+\-Ee\t]*$")


def _norm_label(lbl: str) -> str:
    return re.sub(r"[\s\-/_]", "", lbl).upper()


def _is_affn_exponents(line: str) -> bool:
    for m in re.finditer(r"[Ee]", line):
        nxt = line[m.end():m.end() + 2]
        before = line[:m.start()]
        if re.match(r"[+-]\d", nxt):
            continue
        if re.search(r"\d\.\d*$", before) and re.match(r"\d", nxt):
            continue
        return False
    return True


def _decode_asdf_line(line: str):
    """Decodifica uma linha (X++(Y..Y)) em AFFN/PAC/SQZ/DIF/DUP.
    Retorna (x_inicial, [y...], terminou_em_DIF)."""
    line = line.split("$$", 1)[0].strip()
    if not line:
        return None
    # AFFN/PAC puro (inclui expoentes). 'E'/'e' também são dígitos SQZ (+5/-5), então só
    # tratamos como expoente quando vem seguido de sinal ou de mantissa com ponto decimal.
    if _AFFN_ONLY.match(line) and _is_affn_exponents(line):
        nums = [float(t) for t in _AFFN_RE.findall(line.replace(",", " "))]
        if not nums:
            return None
        return nums[0], nums[1:], False
    # comprimido: tokeniza caractere a caractere
    tokens = []  # [modo, texto]  modo 0=AFFN/PAC abs, 1=SQZ abs, 2=DIF, 3=DUP
    for ch in line:
        if ch in _ASDF:
            s, mode = _ASDF[ch]
            tokens.append([mode, s])
        elif ch in "+-":
            tokens.append([0, ch])
        elif ch.isdigit() or ch == ".":
            if not tokens or tokens[-1][1] == "" or tokens[-1] is None:
                tokens.append([0, ch])
            else:
                tokens[-1][1] += ch
        elif ch in " \t,;":
            tokens.append([0, ""])  # separador
        else:
            continue
    tokens = [t for t in tokens if t[1] not in ("", "+", "-")]
    if not tokens:
        return None
    x0 = float(tokens[0][1])
    ys: List[float] = []
    last_mode = None
    last_dif = 0.0
    for mode, s in tokens[1:]:
        v = float(s)
        if mode in (0, 1):
            ys.append(v)
            last_mode = "abs"
        elif mode == 2:
            if not ys:
                raise ValueError("DIF sem valor absoluto anterior na linha: " + line[:60])
            ys.append(ys[-1] + v)
            last_dif = v
            last_mode = "dif"
        elif mode == 3:
            n = int(v) - 1
            for _ in range(n):
                if last_mode == "dif":
                    ys.append(ys[-1] + last_dif)
                else:
                    ys.append(ys[-1])
    ends_dif = last_mode == "dif"
    return x0, ys, ends_dif


def _parse_jcamp_block(lines: List[str]) -> dict:
    """Recebe as linhas de um bloco (##TITLE ... ##END) e devolve dict com header e dados."""
    hdr = {}
    data_kind = None
    data_lines: List[str] = []
    cur = None
    for raw in lines:
        line = raw.rstrip("\r\n")
        if line.startswith("##"):
            m = re.match(r"##(.*?)=(.*)", line)
            if not m:
                continue
            lbl, val = _norm_label(m.group(1)), m.group(2).split("$$", 1)[0].strip()
            if lbl in ("XYDATA", "XYPOINTS", "PEAKTABLE", "DATATABLE"):
                data_kind = (lbl, val.replace(" ", "").upper())
                cur = "data"
                continue
            cur = lbl
            if lbl == "END":
                break
            hdr[lbl] = val
        elif cur == "data":
            data_lines.append(line)
        elif cur and cur in hdr and not line.startswith("$$"):
            hdr[cur] += "\n" + line  # rótulo multilinha
    return {"hdr": hdr, "data_kind": data_kind, "data_lines": data_lines}


def _jcamp_block_to_spectrum(blk: dict, path: str) -> Optional[Spectrum]:
    hdr, kind, dl = blk["hdr"], blk["data_kind"], blk["data_lines"]
    if kind is None:
        return None
    lbl, form = kind
    f = lambda k, d=None: _tofloat(hdr.get(k), d)
    xfac, yfac = f("XFACTOR", 1.0), f("YFACTOR", 1.0)
    firstx, lastx, npts = f("FIRSTX"), f("LASTX"), hdr.get("NPOINTS")
    npts = int(float(npts)) if npts not in (None, "") else None
    warn = []

    if lbl in ("XYDATA", "DATATABLE") and "X++" in form:
        xs, ys = [], []
        prev_ended_dif = False
        deltax = None
        if firstx is not None and lastx is not None and npts and npts > 1:
            deltax = (lastx - firstx) / (npts - 1)
        elif f("DELTAX") is not None:
            deltax = f("DELTAX")
        for line in dl:
            dec = _decode_asdf_line(line)
            if dec is None:
                continue
            x0, yl, ends_dif = dec
            x0 *= xfac
            if prev_ended_dif and yl and ys:
                # primeiro Y é valor de checagem (repete o último da linha anterior)
                if abs(yl[0] - ys[-1] / yfac) > 1e-6 * max(1.0, abs(yl[0])):
                    warn.append(f"Checagem Y falhou perto de x={x0:g}")
                yl = yl[1:]
                start = 1
            else:
                start = 0
            for k, yv in enumerate(yl):
                xs.append(x0 + (k + start) * (deltax or 0.0))
                ys.append(yv * yfac)
            prev_ended_dif = ends_dif
        x, y = np.array(xs, float), np.array(ys, float)
        if npts and len(y) == npts and firstx is not None and lastx is not None:
            x = np.linspace(firstx, lastx, npts)
        elif npts and len(y) != npts:
            warn.append(f"NPOINTS={npts} mas {len(y)} valores decodificados.")
    elif lbl in ("XYPOINTS", "PEAKTABLE", "DATATABLE"):
        nums = []
        for line in dl:
            line = line.split("$$", 1)[0]
            nums += [float(t) for t in _AFFN_RE.findall(line.replace(";", " ").replace(",", " "))]
        if "XYW" in form:
            arr = np.array(nums[: len(nums) // 3 * 3]).reshape(-1, 3)
        else:
            arr = np.array(nums[: len(nums) // 2 * 2]).reshape(-1, 2)
        x, y = arr[:, 0] * xfac, arr[:, 1] * yfac
    else:
        warn.append(f"Forma de dados não suportada: {lbl}={form}")
        return None

    xunit = hdr.get("XUNITS", "NANOMETERS")
    x, xu = _x_to_nm(x, xunit)
    yunits = hdr.get("YUNITS", "")
    declared = _ykind_from_label(yunits)
    title = hdr.get("TITLE", "").strip()
    meta = {
        "format": "JCAMP-DX" + (" (peak table)" if lbl == "PEAKTABLE" else ""),
        "jcamp_version": hdr.get("JCAMPDX", ""),
        "data_type": hdr.get("DATATYPE", ""),
        "origin": hdr.get("ORIGIN", ""), "owner": hdr.get("OWNER", ""),
        "date": hdr.get("LONGDATE", hdr.get("DATE", "")),
        "instrument": hdr.get("SPECTROMETER/DATASYSTEM".replace("/", ""), hdr.get("SPECTROMETERDATASYSTEM", "")),
        "resolution": hdr.get("RESOLUTION", ""), "state": hdr.get("STATE", ""),
        "cas": hdr.get("CASREGISTRYNO", ""), "formula": hdr.get("MOLFORM", ""),
        "xunits_original": xunit, "yunits_original": yunits,
        "is_peak_table": lbl == "PEAKTABLE",
        "header": hdr, "warnings": warn,
    }
    sp = Spectrum(x=x, y=y, title=title, source=path, meta=meta)
    return _finish(sp, declared)


def _tofloat(v, default=None):
    if v in (None, ""):
        return default
    try:
        return float(str(v).split()[0].replace(",", "."))
    except ValueError:
        return default


def read_jcamp(path: str) -> List[Spectrum]:
    """Lê .jdx/.dx/.jcamp. Suporta XYDATA=(X++(Y..Y)) em AFFN, PAC, SQZ, DIF, DUP,
    XYPOINTS/PEAK TABLE, arquivos compostos (LINK) com vários blocos."""
    with open(path, "rb") as fh:
        raw = fh.read()
    text = _decode_bytes(raw)
    lines = text.splitlines()
    # quebra em blocos: cada ##TITLE abre um bloco; blocos LINK (sem dados) são ignorados
    blocks, cur = [], None
    for ln in lines:
        if re.match(r"##\s*TITLE\s*=", ln, re.I):
            if cur:
                blocks.append(cur)
            cur = [ln]
        elif cur is not None:
            cur.append(ln)
    if cur:
        blocks.append(cur)
    out = []
    for b in blocks:
        blk = _parse_jcamp_block(b)
        try:
            sp = _jcamp_block_to_spectrum(blk, path)
        except Exception as e:  # noqa
            warnings.warn(f"Bloco '{blk['hdr'].get('TITLE','?')}' ignorado: {e}")
            continue
        if sp is not None and sp.npoints > 0:
            out.append(sp)
    if not out:
        raise ValueError(f"Nenhum bloco de dados legível em {path} "
                         "(NTUPLES/espectros 2D não são suportados por este leitor).")
    return out


# --------------------------------------------------------------------------
# CSV / TXT / DPT / PRN
# --------------------------------------------------------------------------
def _decode_bytes(raw: bytes) -> str:
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    for enc in ("utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def _try_parse(lines, delim, dec):
    rows = []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        parts = re.split(r"\s+", s) if delim == "ws" else [p.strip() for p in s.split(delim)]
        parts = [p for p in parts if p != ""] if delim == "ws" else parts
        vals = []
        ok = True
        for p in parts:
            if p == "":
                vals.append(np.nan)
                continue
            q = p.strip('"').strip("'")
            if dec == ",":
                if "." in q and "," in q:
                    q = q.replace(".", "")  # milhar com ponto
                q = q.replace(",", ".")
            try:
                vals.append(float(q))
            except ValueError:
                ok = False
                break
        rows.append(vals if ok else None)
    return rows


def _has_decimal_comma(lines, delim) -> bool:
    """Há números do tipo 123,45 nos campos (separados por `delim`)?"""
    pat = re.compile(r"^\s*[+-]?\d+,\d+(?:[Ee][+-]?\d+)?\s*$")
    n = 0
    for ln in lines[:400]:
        parts = re.split(r"\s+", ln.strip()) if delim == "ws" else ln.split(delim)
        n += sum(bool(pat.match(p.strip('"'))) for p in parts)
        if n >= 3:
            return True
    return False




def _guess_xunit(x_raw, xlabel):
    xl = (xlabel or "").lower()
    if "ev" in xl.replace(" ", "") and "nm" not in xl:
        return "EV"
    if "cm" in xl and "nm" not in xl:
        return "1/CM"
    if "µm" in xl or "um)" in xl or "micro" in xl:
        return "MICROMETERS"
    if "nm" in xl or "wavelength" in xl or "comprimento" in xl or "lambda" in xl or "λ" in xl:
        return "NANOMETERS"
    mx, mn = np.nanmax(x_raw), np.nanmin(x_raw)
    if mx <= 7.0 and mn > 0.3:
        return "EV"
    if mx > 4000:
        return "1/CM"
    return "NANOMETERS"


def _read_wide(lines, path) -> Optional[List[Spectrum]]:
    """Layout 'largo': 1ª linha = eixo x, demais linhas = um espectro cada (1ª célula = nome)."""
    for delim in (";", "\t", ","):
        for dec in (".", ","):
            if delim == "," and dec == ",":
                continue
            rows = []
            for ln in lines:
                if not ln.strip():
                    continue
                parts = ln.split(delim)
                if len(parts) < 20:
                    break
                r = _try_parse([delim.join(parts[1:])], delim, dec)[0]
                if r is None:
                    break
                rows.append((parts[0].strip().strip('"'), np.array(r, float)))
            if len(rows) >= 2 and all(len(r[1]) == len(rows[0][1]) for r in rows):
                xr = rows[0][1]
                x, xu = _x_to_nm(xr, _guess_xunit(xr, rows[0][0]))
                out = []
                for name, y in rows[1:]:
                    meta = {"format": "CSV/texto (largo: espectros em linhas)", "delimiter": repr(delim),
                            "decimal": dec, "xunits_original": xu, "warnings": []}
                    out.append(_finish(Spectrum(x=x.copy(), y=y, title=name or os.path.basename(path),
                                                source=path, meta=meta), None))
                return out
    return None


_OCEAN_KEYS = {"Integration Time": "integration_time", "Spectrometer": "instrument", "Date": "date",
               "Scans to average": "scans", "Boxcar width": "boxcar", "Processing": "processing"}


def _parse_lines(lines, path, xcol=0, ycols=None, extra_meta=None) -> List[Spectrum]:
    best = None
    for delim in (";", "\t", ",", "ws"):
        for dec in (".", ","):
            if delim == "," and dec == ",":
                continue
            if dec == "," and not _has_decimal_comma(lines, delim):
                continue
            rows = _try_parse(lines, delim, dec)
            good = [r for r in rows if r is not None and sum(np.isfinite(r)) >= 2]
            if len(good) < 5:
                continue
            ncol = max(set(len(r) for r in good), key=lambda n: sum(len(r) == n for r in good))
            score = sum(len(r) == ncol for r in good)
            key = (score, ncol, dec == "," and delim in (";", "\t"))
            if best is None or key > best[0]:
                best = (key, delim, dec, ncol)
    wide = _read_wide(lines, path)
    if wide and (best is None or best[0][0] < 5 or best[3] > 50):
        return wide
    if best is None:
        raise ValueError(f"Não consegui achar colunas numéricas em {path}")
    _, delim, dec, ncol = best
    rows = _try_parse(lines, delim, dec)
    nonempty = [ln for ln in lines if ln.strip()]
    first_data = next((i for i, r in enumerate(rows) if r is not None and len(r) == ncol
                       and sum(np.isfinite(r)) >= 2), 0)
    pre = nonempty[:first_data]
    data = np.array([r for r in rows if r is not None and len(r) == ncol], float)
    # remove colunas totalmente vazias (vírgula final das exportações Cary)
    keep = [c for c in range(ncol) if np.isfinite(data[:, c]).any()]
    data = data[:, keep]
    ncol = data.shape[1]

    def split(line):
        if delim == "ws":
            return line.split()
        return [h.strip().strip('"') for h in line.split("\t" if delim == "\t" else delim)]

    names = names0 = None
    if pre:
        h = split(pre[-1])
        h = [h[c] if c < len(h) else "" for c in keep]
        names = h if len(h) == ncol else None
        if len(pre) >= 2:
            h0 = split(pre[-2])
            names0 = [h0[c] if c < len(h0) else "" for c in keep]
    meta_common = {"format": "CSV/texto", "delimiter": "whitespace" if delim == "ws" else repr(delim),
                   "decimal": dec, "header_lines": pre[:40], "warnings": []}
    for ln in pre:
        for k, mk in _OCEAN_KEYS.items():
            if ln.strip().startswith(k):
                meta_common[mk] = ln.split(":", 1)[-1].strip()
    if any(">>>>>Begin" in ln for ln in pre):
        meta_common["format"] = "Ocean Insight / OceanView"
    if extra_meta:
        meta_common.update(extra_meta)

    # ---- layout pareado (Cary: x,y,x,y,...) ?
    pairs = None
    if ncol >= 4 and ncol % 2 == 0 and xcol == 0 and ycols is None:
        xs = data[:, 0::2]
        lab_x = names and all(re.search(r"wave|nm|λ|compr", (names[c] or "").lower()) for c in range(0, ncol, 2))
        mono = all(np.all(np.diff(col[np.isfinite(col)]) > 0) or np.all(np.diff(col[np.isfinite(col)]) < 0)
                   for col in xs.T if np.isfinite(col).sum() > 3)
        if lab_x or (mono and np.nanmax(np.abs(xs - xs[:, :1])) < 1e-6 * np.nanmax(np.abs(xs)) + 1e-9):
            pairs = [(c, c + 1) for c in range(0, ncol, 2)]
    if pairs is None:
        cols = [c for c in range(ncol) if c != xcol] if ycols is None else list(ycols)
        pairs = [(xcol, c) for c in cols]
        meta_common["format"] = meta_common["format"] if "Ocean" in meta_common["format"] else "CSV/texto"
    else:
        meta_common["format"] = "CSV pareado (Cary/Varian)"
    out = []
    base = os.path.splitext(os.path.basename(path))[0]
    for cx, cy in pairs:
        xr, y = data[:, cx], data[:, cy]
        ok = np.isfinite(xr) & np.isfinite(y)
        xr, y = xr[ok], y[ok]
        if xr.size < 3:
            continue
        xlabel = names[cx] if names else ""
        ylabel = names[cy] if names else ""
        xu = _guess_xunit(xr, xlabel)
        x, xu_n = _x_to_nm(xr, xu)
        declared = _ykind_from_label(ylabel)
        if declared is None and pre:
            if "Ocean" in meta_common["format"]:
                # OceanView: o modo vem no nome do nó ("..._Absorbance__0__...") ou em 'Processing'
                hdr = " ".join(pre).upper()
                declared = ("A" if "ABSORB" in hdr else "%T" if "TRANSMI" in hdr else
                            "%R" if "REFLECT" in hdr else None)
            else:
                declared = _ykind_from_label(" ".join(pre[-3:]))
        title = base
        if len(pairs) > 1:
            cand = (names0[cx] if names0 and names0[cx] else "") or ylabel or f"col{cy}"
            title = cand
        meta = dict(meta_common, column=cy, x_label=xlabel, y_label=ylabel, xunits_original=xu_n,
                    warnings=list(meta_common["warnings"]))
        sp = Spectrum(x=x, y=y, title=title, source=path, meta=meta)
        out.append(_finish(sp, declared))
    if not out:
        raise ValueError(f"Nenhum espectro legível em {path}")
    return out


def read_csv(path: str, xcol: int = 0, ycols=None) -> List[Spectrum]:
    """Texto delimitado (CSV, TXT, PRN, ASC, DAT...). Detecta separador e vírgula decimal,
    cabeçalho/rodapé, colunas pareadas (Cary: x,y,x,y) com nomes na 1ª linha, cabeçalho do
    OceanView ('>>>>>Begin Spectral Data<<<<<') e layout largo. Cada coluna y vira um Spectrum."""
    with open(path, "rb") as fh:
        text = _decode_bytes(fh.read())
    return _parse_lines(text.splitlines(), path, xcol, ycols)


def read_excel(path: str, sheet=None) -> List[Spectrum]:
    """Excel (.xlsx/.xls): cada planilha é lida como uma tabela de texto (mesmas regras do CSV)."""
    import pandas as pd
    sheets = pd.read_excel(path, sheet_name=sheet if sheet is not None else None, header=None, dtype=str)
    if not isinstance(sheets, dict):
        sheets = {sheet: sheets}
    out = []
    for name, df in sheets.items():
        df = df.fillna("")
        lines = ["\t".join(str(v).strip() for v in row) for row in df.itertuples(index=False)]
        try:
            got = _parse_lines(lines, path, extra_meta={"sheet": name})
        except ValueError:
            continue
        hint = _ykind_from_label(str(name))
        for s in got:
            s.meta["format"] = "Excel (" + s.meta.get("format", "") + ")"
            if s.meta.get("ykind_guessed") and hint:
                s.ykind = hint
                s.meta["ykind_guessed"] = False
                s.meta["warnings"] = [w for w in s.meta["warnings"] if "percentual" not in w]
            if len(sheets) > 1 and len(got) == 1:
                s.title = str(name)
        out += got
    if not out:
        raise ValueError(f"Nenhuma tabela numérica em {path}")
    return out


# --------------------------------------------------------------------------
# Galactic GRAMS .spc
# --------------------------------------------------------------------------
_SPC_X = {0: "arb", 1: "1/CM", 2: "MICROMETERS", 3: "NANOMETERS", 13: "1/CM", 14: "EV"}
_SPC_Y = {2: "A", 3: "KM", 4: "counts", 10: "A", 11: None, 12: "counts", 13: "arb", 23: "eps",
          128: "T", 129: "R", 131: "counts"}
_SPC_NEW = struct.Struct("<BBBbiddiBBBBi9s9sH32s130s30siiBBHf48sfifB187s")
assert _SPC_NEW.size == 512


def _swap_words(raw: bytes) -> np.ndarray:
    a = np.frombuffer(raw, dtype="<u2").reshape(-1, 2)[:, ::-1].copy()
    return a.view("<i4").ravel()


def _roughness(y):
    y = np.asarray(y, float)
    if y.size < 5 or not np.isfinite(y).all() or np.ptp(y) == 0:
        return np.inf
    return float(np.median(np.abs(np.diff(y, 2))) / np.ptp(y))


def read_spc(path: str) -> List[Spectrum]:
    """Galactic/Thermo GRAMS .spc. Formato novo (0x4B; x uniforme ou explícito, multi-subarquivo,
    y float/int32/int16 com expoente) e antigo (0x4D, inteiros com palavras trocadas)."""
    with open(path, "rb") as fh:
        b = fh.read()
    if len(b) < 256:
        raise ValueError("arquivo .spc curto demais")
    ver = b[1]
    out = []
    base = os.path.splitext(os.path.basename(path))[0]
    if ver == 0x4B:
        (flg, _, _, fexp, npts, first, last, nsub, xt, yt, zt, post, fdate, fres, fsrc, peakpt, spare,
         cmnt, catxt, logoff, mods, procs, level, sampin, factor, method, zinc, wplanes, winc, wtype,
         resv) = _SPC_NEW.unpack_from(b, 0)
        TSPREC, TMULTI, TXYXYS, TXVALS = 0x01, 0x04, 0x40, 0x80
        pos = 512
        x_common = None
        if flg & TXVALS and not flg & TXYXYS:
            x_common = np.frombuffer(b, "<f4", npts, pos).astype(float)
            pos += 4 * npts
        elif not flg & TXYXYS:
            x_common = np.linspace(first, last, npts)
        nsub = nsub if (flg & TMULTI) else 1
        comment = cmnt.split(b"\x00")[0].decode("latin-1", "replace").strip()
        xunit = _SPC_X.get(xt, "NANOMETERS")
        ykind = _SPC_Y.get(yt, None)
        for k in range(max(1, nsub)):
            if pos + 32 > len(b):
                break
            sflg, sexp, sidx, stime, snext, snois, snpts, sscan, swl, _ = struct.unpack_from("<BbhfffiIf4s", b, pos)
            pos += 32
            if flg & TXYXYS:
                n = snpts
                xr = np.frombuffer(b, "<f4", n, pos).astype(float)
                pos += 4 * n
            else:
                n = npts
                xr = x_common
            exp = sexp if (flg & TMULTI) else fexp
            if exp == -128:
                y = np.frombuffer(b, "<f4", n, pos).astype(float)
                pos += 4 * n
            elif flg & TSPREC:
                y = np.frombuffer(b, "<i2", n, pos).astype(float) * 2.0 ** (exp - 16)
                pos += 2 * n
            else:
                y = np.frombuffer(b, "<i4", n, pos).astype(float) * 2.0 ** (exp - 32)
                pos += 4 * n
            x, xu = _x_to_nm(xr, xunit)
            meta = {"format": "Galactic SPC", "comment": comment, "spc_xtype": xt, "spc_ytype": yt,
                    "xunits_original": xu, "subfile": k, "subfile_z": stime, "warnings": []}
            if xu not in ("nm", "cm-1", "um", "eV"):
                meta["warnings"].append(f"Eixo x do .spc tem tipo {xt} (não é comprimento de onda/número de onda).")
            title = comment[:60] or base
            if nsub > 1:
                title = f"{base} #{k} (z={stime:g})"
            sp = Spectrum(x=x, y=y, title=title, source=path, meta=meta)
            out.append(_finish(sp, ykind))
        if logoff and logoff < len(b):
            try:
                sz, = struct.unpack_from("<i", b, logoff)
                txtoff, = struct.unpack_from("<i", b, logoff + 8)
                txt = b[logoff + txtoff:logoff + sz].split(b"\x00")[0].decode("latin-1", "replace")
                for s in out:
                    s.meta["spc_log"] = txt[:2000]
            except struct.error:
                pass
    elif ver == 0x4D:
        flg = b[0]
        oexp, = struct.unpack_from("<h", b, 2)
        onpts, ofirst, olast = struct.unpack_from("<fff", b, 4)
        oxt, oyt = b[16], b[17]
        cmnt = b[62:192].split(b"\x00")[0].decode("latin-1", "replace").strip()
        n = int(round(onpts))
        multi = bool(flg & 0x04)
        xunit = _SPC_X.get(oxt, "NANOMETERS")
        x, xu = _x_to_nm(np.linspace(ofirst, olast, n), xunit)
        pos = 256
        k = 0
        while pos + 4 * n <= len(b) and (k == 0 or multi):
            if k > 0 or multi:
                if k > 0:
                    if pos + 32 + 4 * n > len(b):
                        break
                    sexp = struct.unpack_from("<b", b, pos + 1)[0]
                    pos += 32
                else:
                    sexp = oexp
            else:
                sexp = oexp
            raw = b[pos:pos + 4 * n]
            pos += 4 * n
            cand = []
            if sexp == -128 or oexp == -128:
                cand.append(np.frombuffer(raw, "<f4").astype(float))
            cand.append(_swap_words(raw).astype(float) * 2.0 ** (sexp - 32))
            cand.append(np.frombuffer(raw, "<i4").astype(float) * 2.0 ** (sexp - 32))
            y = min(cand, key=_roughness)   # a ordem de palavras correta dá a curva suave
            meta = {"format": "Galactic SPC (antigo 0x4D)", "comment": cmnt, "xunits_original": xu,
                    "subfile": k, "warnings": ["Formato .spc antigo: ordem de palavras inferida pela suavidade da curva."]}
            title = cmnt[:60] or base
            if multi:
                title = f"{base} #{k}"
            out.append(_finish(Spectrum(x=x.copy(), y=y, title=title, source=path, meta=meta),
                               _SPC_Y.get(oyt)))
            k += 1
    else:
        raise ValueError(f"Versão de .spc desconhecida ({ver:#x}).")
    if not out:
        raise ValueError(f"Nenhum subarquivo legível em {path}")
    return out


# --------------------------------------------------------------------------
# PerkinElmer .sp (UV WinLab / Lambda e Spectrum - mesmo contêiner 'PEPE')
# --------------------------------------------------------------------------
_SP_ABSC_RANGE, _SP_ORD_RANGE, _SP_INTERVAL, _SP_NPTS = 35698, 35699, 35700, 35701
_SP_XLABEL, _SP_YLABEL, _SP_DATA, _SP_NAME, _SP_ALIAS = 35703, 35704, 35708, 35709, 35713
_T_RANGE, _T_DOUBLE, _T_INT, _T_ARRAY, _T_STRING = 29981, 29979, 29995, 29974, 29987


def _sp_strings(buf: bytes) -> List[str]:
    out, i = [], 0
    while True:
        i = buf.find(b"#u", i)
        if i < 0 or i + 4 > len(buf):
            break
        n = struct.unpack_from("<h", buf, i + 2)[0]
        if 0 < n < 512 and i + 4 + n <= len(buf):
            s = buf[i + 4:i + 4 + n].decode("latin-1", errors="replace").strip("\x00 ")
            if s:
                out.append(s)
            i += 4 + n
        else:
            i += 2
    return out


def read_sp(path: str) -> List[Spectrum]:
    """PerkinElmer .sp ('PEPE'). Limites do eixo INCLUSIVOS (n = |fim-início|/passo + 1);
    um valor de preenchimento ~0 antes dos dados é descartado se existir."""
    with open(path, "rb") as fh:
        b = fh.read()
    if b[:4] != b"PEPE":
        raise ValueError(f"{path} não começa com 'PEPE'.")
    descr = b[4:44].split(b"\x00")[0].decode("latin-1", errors="replace")
    start = None
    for i in range(44, len(b) - 24):
        if struct.unpack_from("<H", b, i)[0] == _SP_ABSC_RANGE and struct.unpack_from("<H", b, i + 6)[0] == _T_RANGE:
            start = i
            break
    if start is None:
        raise ValueError("Membro AbscissaRange não encontrado - variante de .sp desconhecida.")
    info, data, i = {}, None, start
    while i + 8 <= len(b):
        mid, size = struct.unpack_from("<Hi", b, i)
        if size < 0 or i + 6 + size > len(b):
            break
        p = b[i + 6:i + 6 + size]
        t = struct.unpack_from("<H", p, 0)[0] if len(p) >= 2 else None
        if mid == _SP_ABSC_RANGE and t == _T_RANGE:
            info["x_first"], info["x_last"] = struct.unpack_from("<dd", p, 2)
        elif mid == _SP_INTERVAL and t == _T_DOUBLE:
            info["interval"] = struct.unpack_from("<d", p, 2)[0]
        elif mid == _SP_NPTS and t == _T_INT:
            info["npoints"] = struct.unpack_from("<i", p, 2)[0]
        elif mid in (_SP_XLABEL, _SP_YLABEL, _SP_NAME, _SP_ALIAS) and t == _T_STRING:
            n = struct.unpack_from("<h", p, 2)[0]
            key = {_SP_XLABEL: "x_label", _SP_YLABEL: "y_label", _SP_NAME: "name", _SP_ALIAS: "alias"}[mid]
            info[key] = p[4:4 + n].decode("latin-1", errors="replace")
        elif mid == _SP_DATA and t == _T_ARRAY:
            n = struct.unpack_from("<I", p, 2)[0]
            data = np.frombuffer(p[6:6 + n], dtype="<f8").copy()
        i += 6 + size
    if data is None or "x_first" not in info:
        raise ValueError("Bloco de dados do .sp não encontrado.")
    warn = []
    npts = info.get("npoints", data.size)
    if data.size == npts + 1 and abs(data[0]) < 1e-12:
        data = data[1:]
        warn.append("Valor de preenchimento inicial descartado.")
    if data.size != npts:
        warn.append(f"npoints={npts} mas {data.size} valores; eixo reconstruído.")
        npts = data.size
    xr = np.linspace(info["x_first"], info["x_last"], npts)
    xlab = info.get("x_label", "nm")
    x, xu = _x_to_nm(xr, "1/CM" if "cm" in xlab.lower() else "NANOMETERS")
    strings = _sp_strings(b)
    meta = {"format": "PerkinElmer .sp", "description": descr, "x_label": xlab, "y_label": info.get("y_label", ""),
            "xunits_original": xu, "sp_header": info, "sp_strings": strings[:200], "warnings": warn}
    title = info.get("alias") or info.get("name") or os.path.splitext(os.path.basename(path))[0]
    title = os.path.splitext(os.path.basename(title.replace("\\", "/")))[0]
    return [_finish(Spectrum(x=x, y=data, title=title, source=path, meta=meta),
                    _ykind_from_label(info.get("y_label", "")))]


# --------------------------------------------------------------------------
# Entrada unificada
# --------------------------------------------------------------------------
_JCAMP_EXT = {".jdx", ".dx", ".jcamp", ".jcm"}


def detect_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    with open(path, "rb") as fh:
        head = fh.read(2048)
    if head[:4] == b"PEPE":
        return "sp"
    if len(head) > 2 and head[1] in (0x4B, 0x4D) and ext == ".spc":
        return "spc"
    if b"##TITLE" in head.upper() or b"##JCAMP" in head.upper():
        return "jcamp"
    if head[:4] == b"PK\x03\x04" or head[:8] == b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1" or ext in (".xlsx", ".xls", ".xlsm"):
        return "excel"
    if ext == ".sp":
        return "sp"
    if ext in _JCAMP_EXT:
        return "jcamp"
    if ext == ".spc":
        return "spc"
    if ext in (".vsp", ".dsw", ".bsw", ".ksd"):
        raise ValueError(f"Formato binário {ext} (Cary/Varian) não suportado: exporte como CSV.")
    if ext in (".vspd", ".uvd", ".pdx"):
        raise ValueError(f"Formato {ext} não suportado: exporte como CSV/TXT ou JCAMP no software do equipamento.")
    return "csv"


def load_all(path: str, **kw) -> List[Spectrum]:
    fmt = detect_format(path)
    if fmt == "sp":
        return read_sp(path)
    if fmt == "spc":
        return read_spc(path)
    if fmt == "jcamp":
        return read_jcamp(path)
    if fmt == "excel":
        return read_excel(path, **kw)
    return read_csv(path, **kw)


def load(path: str, index: int = 0, **kw) -> Spectrum:
    specs = load_all(path, **kw)
    if len(specs) > 1:
        warnings.warn(f"{path}: {len(specs)} espectros no arquivo; usando o índice {index}.")
    return specs[index]


def save_csv(spec: Spectrum, path: str, decimal: str = ".", sep: Optional[str] = None, energy: bool = False):
    sep = sep or (";" if decimal == "," else ",")
    ylab = {"A": "Absorbance", "%T": "Transmittance (%)", "T": "Transmittance", "%R": "Reflectance (%)",
            "R": "Reflectance", "KM": "F(R)", "eps": "epsilon (M-1 cm-1)", "alpha": "alpha (cm-1)",
            "counts": "Counts"}.get(spec.ykind, "Intensity")
    s = spec.sorted()
    buf = io.StringIO()
    buf.write(f"Wavelength (nm){sep}" + (f"Energy (eV){sep}" if energy else "") + f"{ylab}\n")
    for xv, yv in zip(s.x, s.y):
        cols = [f"{xv:.3f}"] + ([f"{HC_EV_NM / xv:.5f}"] if energy else []) + [f"{yv:.6g}"]
        if decimal == ",":
            cols = [c.replace(".", ",") for c in cols]
        buf.write(sep.join(cols) + "\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(buf.getvalue())
