"""
ftir_io.py - leitura de espectros FT-IR (JCAMP-DX, CSV/TXT/DPT, PerkinElmer .sp)

Uso rápido:
    from ftir_io import load, load_all
    s = load("amostra.sp")          # um Spectrum
    specs = load_all("dados.csv")   # lista (CSV multi-coluna, JCAMP com vários blocos)

Todo Spectrum sai com o eixo x em número de onda (cm-1), na ordem do arquivo,
e com `ykind` dizendo o que é o eixo y:
    'A'   absorbância            'T'   transmitância fracionária (0-1)
    '%T'  transmitância em %     'R' / '%R'  reflectância
    'KM'  Kubelka-Munk           'eps' absortividade (ex.: PNNL/NIST quantitativo)
    'arb' desconhecido / unidades arbitrárias
`meta['ykind_guessed'] = True` quando o tipo foi inferido pela forma dos dados
(e não declarado no arquivo) - confirme com o usuário antes de concluir algo.

Só depende de numpy.
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

__all__ = ["Spectrum", "load", "load_all", "read_jcamp", "read_csv", "read_sp",
           "save_csv", "save_jcamp"]


# --------------------------------------------------------------------------
# Modelo de dados
# --------------------------------------------------------------------------
@dataclass
class Spectrum:
    x: np.ndarray                      # cm-1
    y: np.ndarray
    ykind: str = "arb"
    title: str = ""
    source: str = ""
    meta: dict = field(default_factory=dict)

    # ---- propriedades úteis -------------------------------------------
    @property
    def npoints(self) -> int:
        return int(self.x.size)

    @property
    def step(self) -> float:
        """Intervalo médio entre pontos (cm-1, positivo)."""
        if self.x.size < 2:
            return float("nan")
        return float(np.median(np.abs(np.diff(self.x))))

    @property
    def xrange(self):
        return float(np.max(self.x)), float(np.min(self.x))

    def copy(self, **kw) -> "Spectrum":
        d = dict(x=self.x.copy(), y=self.y.copy(), ykind=self.ykind,
                 title=self.title, source=self.source, meta=dict(self.meta))
        d.update(kw)
        return Spectrum(**d)

    # ---- conversões de ordenada ----------------------------------------
    def to_absorbance(self, floor: float = 1e-5) -> "Spectrum":
        """Converte T/%T/R/%R para absorbância (-log10). A, KM, eps, arb passam direto."""
        k = self.ykind
        if k in ("A", "eps", "KM", "arb"):
            return self.copy()
        t = self.y / 100.0 if k in ("%T", "%R") else self.y.astype(float)
        nbad = int(np.sum(t <= 0))
        if nbad:
            warnings.warn(f"{nbad} pontos com T<=0 (saturação/ruído); recortados em {floor}.")
        a = -np.log10(np.clip(t, floor, None))
        newk = "A"
        out = self.copy(y=a, ykind=newk)
        out.meta["converted_from"] = k
        return out

    def to_transmittance(self, percent: bool = True) -> "Spectrum":
        k = self.ykind
        if k in ("%T", "T"):
            t = self.y / 100.0 if k == "%T" else self.y
        elif k == "A":
            t = 10.0 ** (-self.y)
        else:
            raise ValueError(f"Não sei converter ykind={k!r} para transmitância.")
        out = self.copy(y=t * (100.0 if percent else 1.0), ykind="%T" if percent else "T")
        out.meta["converted_from"] = k
        return out

    def sorted(self, descending: bool = True) -> "Spectrum":
        idx = np.argsort(self.x)
        if descending:
            idx = idx[::-1]
        return self.copy(x=self.x[idx], y=self.y[idx])

    def crop(self, hi: float, lo: float) -> "Spectrum":
        hi, lo = max(hi, lo), min(hi, lo)
        m = (self.x <= hi) & (self.x >= lo)
        return self.copy(x=self.x[m], y=self.y[m])

    def summary(self) -> str:
        hi, lo = self.xrange
        lines = [
            f"Título      : {self.title or '(sem título)'}",
            f"Arquivo     : {self.source}",
            f"Formato     : {self.meta.get('format', '?')}",
            f"Pontos      : {self.npoints}   intervalo ~{self.step:.4g} cm-1",
            f"Faixa       : {hi:.1f} -> {lo:.1f} cm-1",
            f"Ordenada    : {self.ykind}" + ("  (INFERIDO - confirmar)" if self.meta.get("ykind_guessed") else ""),
            f"y min/max   : {np.nanmin(self.y):.4g} / {np.nanmax(self.y):.4g}",
        ]
        for k in ("origin", "owner", "date", "instrument", "accessory", "resolution",
                  "scans", "state", "xunits_original", "yunits_original", "jcamp_version",
                  "x_label", "y_label", "column"):
            if self.meta.get(k) not in (None, ""):
                lines.append(f"{k:<12}: {self.meta[k]}")
        for w in self.meta.get("warnings", []):
            lines.append(f"AVISO       : {w}")
        return "\n".join(lines)


# --------------------------------------------------------------------------
# Helpers de unidade
# --------------------------------------------------------------------------
def _x_to_wavenumber(x: np.ndarray, unit: str):
    """Retorna (x_em_cm-1, nome_da_unidade_original_normalizada)."""
    u = (unit or "").strip().upper().replace(" ", "")
    if u in ("", "1/CM", "CM-1", "CM^-1", "CM**-1", "WAVENUMBER", "WAVENUMBERS", "1/CENTIMETERS", "CM−1"):
        return x.astype(float), "cm-1"
    if u in ("MICROMETERS", "MICROMETER", "UM", "µM", "MICRONS", "MICRON", "Μm".upper()):
        return 1e4 / x.astype(float), "um"
    if u in ("NANOMETERS", "NANOMETER", "NM"):
        return 1e7 / x.astype(float), "nm"
    return x.astype(float), u.lower() or "cm-1"


def _ykind_from_label(label: str) -> Optional[str]:
    if not label:
        return None
    u = label.strip().upper()
    u2 = u.replace(" ", "")
    if "KUBELKA" in u or u2 in ("KM", "K-M", "KUBELKA-MUNK"):
        return "KM"
    if "%T" in u2 or "PERCENTTRANS" in u2 or ("TRANSMITTANCE" in u and "%" in u):
        return "%T"
    if "TRANSMIT" in u:
        return "T"
    if "%R" in u2 or ("REFLECT" in u and "%" in u):
        return "%R"
    if "REFLECT" in u:
        return "R"
    if "ABSORB" in u or u2 in ("A", "ABS", "ABS.", "AU", "A.U.") or re.search(r"\bABS\b|\(A\)|\bA\.U\.", u):
        return "A"
    if "MOL" in u or "ABSORPTIVITY" in u or "EXTINCTION" in u or "(BASE10)" in u2:
        return "eps"
    if "ARBITRARY" in u or u2 in ("COUNTS", "INTENSITY", "SINGLEBEAM"):
        return "arb"
    return None


def guess_ykind(y: np.ndarray) -> str:
    """Infere A vs T vs %T pela forma: em absorbância a maior parte dos pontos fica
    perto do mínimo (linha de base) e as bandas apontam para cima; em transmitância
    a linha de base fica perto do máximo e as bandas apontam para baixo."""
    y = np.asarray(y, float)
    y = y[np.isfinite(y)]
    if y.size == 0:
        return "arb"
    lo, hi, med = np.percentile(y, 1), np.percentile(y, 99), np.median(y)
    span = hi - lo
    if span <= 0:
        return "arb"
    up = (hi - med) / span        # bandas para cima -> grande
    if up > 0.5:                  # mediana perto do mínimo
        return "A" if hi < 6 else "arb"
    # bandas para baixo (mediana perto do máximo) -> transmitância
    if 1.5 < hi <= 120 and med > 0.5 * hi:
        return "%T"
    if hi <= 1.2 and med > 0.5 * hi:
        return "T"
    return "A" if (hi < 6 and lo > -1) else "arb"


def _finish(spec: Spectrum, declared: Optional[str]) -> Spectrum:
    """Ajusta ykind: usa o declarado, corrige T vs %T pelo intervalo de valores."""
    w = spec.meta.setdefault("warnings", [])
    if declared is None:
        spec.ykind = guess_ykind(spec.y)
        spec.meta["ykind_guessed"] = True
    else:
        spec.ykind = declared
        mx = np.nanmax(spec.y)
        if declared == "T" and mx > 2.0:
            spec.ykind = "%T"
            w.append("Declarado TRANSMITTANCE mas valores até %.3g: tratado como %%T." % mx)
        elif declared == "%T" and mx <= 1.5:
            spec.ykind = "T"
            w.append("Declarado %%T mas valores <= 1.5: tratado como T fracionária.")
        elif declared == "R" and mx > 2.0:
            spec.ykind = "%R"
    if spec.ykind in ("T", "%T"):
        tmin = np.nanmin(spec.y) / (100.0 if spec.ykind == "%T" else 1.0)
        if tmin < 0.005:
            w.append(f"Transmitância mínima {tmin*100:.2g}% - bandas mais fortes podem estar saturadas.")
    if spec.ykind == "A" and np.nanmax(spec.y) > 2.5:
        w.append(f"Absorbância máxima {np.nanmax(spec.y):.2f} > 2.5 - possível saturação/não-linearidade.")
    return spec


# --------------------------------------------------------------------------
# JCAMP-DX
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

    xunit = hdr.get("XUNITS", "1/CM")
    x, xu = _x_to_wavenumber(x, xunit)
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


def _read_wide(lines, path) -> Optional[List[Spectrum]]:
    """Layout 'largo' (transposto): 1ª linha = eixo x, demais linhas = um espectro cada,
    1ª célula de cada linha = nome. Comum em exportações em lote / quimiometria."""
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
                x, xu = _x_to_wavenumber(xr, "MICROMETERS" if np.nanmax(xr) < 30 else "1/CM")
                out = []
                for name, y in rows[1:]:
                    meta = {"format": "CSV/texto (largo: espectros em linhas)", "delimiter": repr(delim),
                            "decimal": dec, "xunits_original": xu, "warnings": []}
                    sp = Spectrum(x=x.copy(), y=y, title=name or os.path.basename(path),
                                  source=path, meta=meta)
                    out.append(_finish(sp, None))
                return out
    return None


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


def read_csv(path: str, xcol: int = 0, ycols=None) -> List[Spectrum]:
    """Lê texto delimitado (CSV, TXT, DPT, PRN, XY, ASC). Detecta separador (; , tab espaço),
    vírgula decimal (padrão BR), linhas de cabeçalho e várias colunas de y.
    Cada coluna y vira um Spectrum."""
    with open(path, "rb") as fh:
        text = _decode_bytes(fh.read())
    lines = text.splitlines()
    best = None
    for delim in (";", "\t", ",", "ws"):
        for dec in (".", ","):
            if delim == "," and dec == ",":
                continue
            if dec == "," and not _has_decimal_comma(lines, delim):
                continue
            rows = _try_parse(lines, delim, dec)
            good = [r for r in rows if r is not None and len(r) >= 2]
            if len(good) < 5:
                continue
            ncol = max(set(len(r) for r in good), key=lambda n: sum(len(r) == n for r in good))
            score = sum(len(r) == ncol for r in good)
            # decimal ',' com delimitador ';' ou tab é muito comum no Brasil/Europa
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
    # cabeçalho = última linha não numérica antes do primeiro bloco de dados
    header, first_data = None, None
    nonempty = [ln for ln in lines if ln.strip()]
    for i, r in enumerate(rows):
        if r is not None and len(r) == ncol:
            first_data = i
            break
        if r is None:
            header = nonempty[i]
    pre = nonempty[:first_data] if first_data else []
    data = np.array([r for r in rows if r is not None and len(r) == ncol], float)
    names = None
    if header:
        sp = "\t" if delim == "\t" else (None if delim == "ws" else delim)
        names = [h.strip().strip('"') for h in (header.split(sp) if sp else header.split())]
        if len(names) != ncol:
            names = None
    x_raw = data[:, xcol]
    xlabel = names[xcol] if names else ""
    xunit = ""
    xl = xlabel.lower()
    if "µm" in xl or "um)" in xl or "micro" in xl:
        xunit = "um"
    elif "nm" in xl and "cm" not in xl:
        xunit = "nm"
    elif np.nanmax(x_raw) < 30 and np.nanmin(x_raw) > 1:
        xunit = "um"
    x, xu = _x_to_wavenumber(x_raw, {"um": "MICROMETERS", "nm": "NANOMETERS"}.get(xunit, "1/CM"))
    cols = [c for c in range(ncol) if c != xcol] if ycols is None else list(ycols)
    out = []
    for c in cols:
        y = data[:, c]
        ok = np.isfinite(x) & np.isfinite(y)
        label = names[c] if names else ""
        declared = _ykind_from_label(label)
        if declared is None:
            # procura pistas no texto antes dos dados (ex.: 'Y Units: %T')
            declared = _ykind_from_label(" ".join(pre)[-300:]) if pre else None
        base = os.path.splitext(os.path.basename(path))[0]
        title = label if (label and len(cols) > 1) else base
        meta = {"format": "CSV/texto", "delimiter": "whitespace" if delim == "ws" else repr(delim),
                "decimal": dec, "column": c, "x_label": xlabel, "y_label": label,
                "xunits_original": xu, "header_lines": pre[:20], "warnings": []}
        sp = Spectrum(x=x[ok], y=y[ok], title=title, source=path, meta=meta)
        out.append(_finish(sp, declared))
    return out


# --------------------------------------------------------------------------
# PerkinElmer Spectrum .sp (binário 'PEPE')
# --------------------------------------------------------------------------
# IDs de membros do bloco de dados (unsigned short) e os códigos de tipo esperados
_SP_ABSC_RANGE, _SP_ORD_RANGE, _SP_INTERVAL, _SP_NPTS = 35698, 35699, 35700, 35701
_SP_XLABEL, _SP_YLABEL, _SP_DATA, _SP_NAME, _SP_ALIAS = 35703, 35704, 35708, 35709, 35713
_T_RANGE, _T_DOUBLE, _T_INT, _T_ARRAY, _T_STRING = 29981, 29979, 29995, 29974, 29987


def _sp_strings(buf: bytes) -> List[str]:
    """Extrai strings '#u'+len(short)+texto de um trecho (metadados do instrumento)."""
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
    """Lê PerkinElmer Spectrum .sp ('PEPE'). Procura o membro AbscissaRange (35698, tipo
    29981) e percorre os membros seguintes (intervalo, nº de pontos, rótulos, dados)."""
    with open(path, "rb") as fh:
        b = fh.read()
    if b[:4] != b"PEPE":
        raise ValueError(f"{path} não começa com 'PEPE' - não é um .sp PerkinElmer.")
    descr = b[4:44].split(b"\x00")[0].decode("latin-1", errors="replace")
    info = {}
    start = None
    for i in range(44, len(b) - 24):
        if struct.unpack_from("<H", b, i)[0] == _SP_ABSC_RANGE and \
           struct.unpack_from("<H", b, i + 6)[0] == _T_RANGE:
            start = i
            break
    if start is None:
        raise ValueError("Membro AbscissaRange não encontrado - variante de .sp desconhecida.")
    i, data = start, None
    while i + 8 <= len(b):
        mid, size = struct.unpack_from("<Hi", b, i)
        if size < 0 or i + 6 + size > len(b):
            break
        p = b[i + 6:i + 6 + size]
        t = struct.unpack_from("<H", p, 0)[0] if len(p) >= 2 else None
        if mid == _SP_ABSC_RANGE and t == _T_RANGE:
            info["x_first"], info["x_last"] = struct.unpack_from("<dd", p, 2)
        elif mid == _SP_ORD_RANGE and t == _T_RANGE:
            info["y_min"], info["y_max"] = struct.unpack_from("<dd", p, 2)
        elif mid == _SP_INTERVAL and t == _T_DOUBLE:
            info["interval"] = struct.unpack_from("<d", p, 2)[0]
        elif mid == _SP_NPTS and t == _T_INT:
            info["npoints"] = struct.unpack_from("<i", p, 2)[0]
        elif mid in (_SP_XLABEL, _SP_YLABEL, _SP_NAME, _SP_ALIAS) and t == _T_STRING:
            n = struct.unpack_from("<h", p, 2)[0]
            key = {_SP_XLABEL: "x_label", _SP_YLABEL: "y_label",
                   _SP_NAME: "name", _SP_ALIAS: "alias"}[mid]
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
        warn.append("Valor de preenchimento (pad) inicial descartado.")
    if data.size != npts:
        warn.append(f"npoints={npts} mas {data.size} valores; eixo x reconstruído com {data.size}.")
        npts = data.size
    x = np.linspace(info["x_first"], info["x_last"], npts)  # limites INCLUSIVOS
    if "interval" in info and npts > 1:
        exp_step = (info["x_last"] - info["x_first"]) / (npts - 1)
        if abs(exp_step - info["interval"]) > 1e-6 * max(1, abs(info["interval"])):
            warn.append(f"Intervalo declarado {info['interval']} != calculado {exp_step:.6g}.")
    xlab = info.get("x_label", "cm-1")
    x, xu = _x_to_wavenumber(x, "NM" if xlab.strip().lower() == "nm" else "1/CM")
    # metadados textuais (instrumento, acessório, varreduras, resolução...)
    strings = _sp_strings(b)
    meta = {"format": "PerkinElmer .sp", "description": descr, "x_label": xlab,
            "y_label": info.get("y_label", ""), "xunits_original": xu,
            "sp_header": info, "sp_strings": strings[:200], "warnings": warn}
    kv = {}
    for k_, v_ in zip(strings, strings[1:]):
        kv.setdefault(k_, v_)
    for key, pats in (("accessory", ("Accessory", "IR Accessory", "Sampling")),
                      ("resolution", ("Resolution",)), ("scans", ("Scans", "Accumulations")),
                      ("instrument", ("Instrument", "Model"))):
        for s_ in strings:
            if any(s_.lower().startswith(p_.lower()) for p_ in pats):
                v_ = kv.get(s_)
                if v_:
                    meta[key] = v_
                    break
    title = info.get("alias") or info.get("name") or os.path.splitext(os.path.basename(path))[0]
    title = os.path.splitext(os.path.basename(title.replace("\\", "/")))[0]
    sp = Spectrum(x=x, y=data, title=title, source=path, meta=meta)
    return [_finish(sp, _ykind_from_label(info.get("y_label", "")))]


# --------------------------------------------------------------------------
# Entrada unificada
# --------------------------------------------------------------------------
_TEXT_EXT = {".csv", ".txt", ".dpt", ".prn", ".xy", ".asc", ".tsv", ".dat"}
_JCAMP_EXT = {".jdx", ".dx", ".jcamp", ".jcm"}


def detect_format(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    with open(path, "rb") as fh:
        head = fh.read(2048)
    if head[:4] == b"PEPE":
        return "sp"
    if b"##TITLE" in head.upper() or b"##JCAMP" in head.upper():
        return "jcamp"
    if ext == ".sp":
        return "sp"
    if ext in _JCAMP_EXT:
        return "jcamp"
    if ext in (".spa", ".spg"):
        raise ValueError("Thermo OMNIC .spa/.spg não suportado: exporte como CSV ou JCAMP.")
    if ext[:2] == ".0" or ext in (".opus",):
        raise ValueError("Bruker OPUS binário não suportado: exporte como .dpt/CSV ou JCAMP "
                         "(ou instale 'brukeropusreader').")
    return "csv"


def load_all(path: str, **kw) -> List[Spectrum]:
    fmt = detect_format(path)
    if fmt == "sp":
        return read_sp(path)
    if fmt == "jcamp":
        return read_jcamp(path)
    return read_csv(path, **kw)


def load(path: str, index: int = 0, **kw) -> Spectrum:
    specs = load_all(path, **kw)
    if len(specs) > 1:
        warnings.warn(f"{path}: {len(specs)} espectros no arquivo; usando o índice {index}.")
    return specs[index]


# --------------------------------------------------------------------------
# Saída
# --------------------------------------------------------------------------
def save_csv(spec: Spectrum, path: str, decimal: str = ".", sep: Optional[str] = None):
    sep = sep or (";" if decimal == "," else ",")
    ylab = {"A": "Absorbance", "%T": "Transmittance (%)", "T": "Transmittance",
            "KM": "Kubelka-Munk", "%R": "Reflectance (%)", "R": "Reflectance"}.get(spec.ykind, "Intensity")
    buf = io.StringIO()
    buf.write(f"Wavenumber (cm-1){sep}{ylab}\n")
    for xv, yv in zip(spec.x, spec.y):
        a, b_ = f"{xv:.4f}", f"{yv:.6g}"
        if decimal == ",":
            a, b_ = a.replace(".", ","), b_.replace(".", ",")
        buf.write(f"{a}{sep}{b_}\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(buf.getvalue())


def save_jcamp(spec: Spectrum, path: str):
    """Grava JCAMP-DX 4.24 em XYPOINTS AFFN (simples e lido por qualquer software)."""
    yun = {"A": "ABSORBANCE", "T": "TRANSMITTANCE", "%T": "TRANSMITTANCE",
           "KM": "KUBELKA-MUNK", "R": "REFLECTANCE", "%R": "REFLECTANCE"}.get(spec.ykind, "ARBITRARY UNITS")
    y = spec.y / 100.0 if spec.ykind in ("%T", "%R") else spec.y
    lines = [f"##TITLE={spec.title or 'spectrum'}", "##JCAMP-DX=4.24", "##DATA TYPE=INFRARED SPECTRUM",
             "##XUNITS=1/CM", f"##YUNITS={yun}", "##XFACTOR=1", "##YFACTOR=1",
             f"##FIRSTX={spec.x[0]:.6f}", f"##LASTX={spec.x[-1]:.6f}", f"##NPOINTS={spec.npoints}",
             f"##FIRSTY={y[0]:.8g}", "##XYPOINTS=(XY..XY)"]
    lines += [f"{xv:.6f}, {yv:.8g}" for xv, yv in zip(spec.x, y)]
    lines.append("##END=")
    with open(path, "w", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join(lines) + "\n")
