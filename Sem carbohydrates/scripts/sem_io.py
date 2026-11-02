"""
sem_io.py - leitura de imagens de MEV (TIF e exportações), metadados e escala

    from sem_io import load
    im = load("amido_01.tif")             # objeto SEMImage
    im.data          # float32 0..1, SEM a barra de informações (databar)
    im.px_um         # tamanho do pixel em µm (None se não calibrado)
    im.px_source     # de onde veio a escala (metadados, barra de escala + OCR, manual...)
    im.meta          # fabricante, kV, aumento, WD, detector, HFW, vácuo, data...
    im.summary()

Escala, em ordem de prioridade:
  1. manual: pixel_size="0.05um" ou scale="10um=235px" ou bar_value="10um" (mede a barra)
  2. metadados no TIF: FEI/Thermo (tag 34682, inclui Phenom/Quanta/Helios/Prisma), Zeiss
     (tag 34118, SmartSEM), TESCAN (tag 50431, [MAIN] PixelSizeX), ImageJ (unit + XResolution),
     Olympus SIS, OME
  3. arquivo .txt acompanhante com o mesmo nome: Hitachi ([SemImageFile] PixelSize em nm) e JEOL
     ($$SM_MICRON_BAR px + $$SM_MICRON_MARKER, em relação a $CM_FULL_SIZE)
  4. barra de escala desenhada na imagem: detecção da linha + OCR do valor (tesseract, se
     instalado). Sem OCR, informe bar_value e a barra é medida automaticamente.
Também aceita PNG/JPG/BMP (mesma lógica, sem metadados de fabricante).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np

UNIT_UM = {"nm": 1e-3, "um": 1.0, "µm": 1.0, "μm": 1.0, "micron": 1.0, "microns": 1.0, "mm": 1e3,
           "m": 1e6, "cm": 1e4, "pm": 1.0, "jim": 1.0, "ym": 1.0, "wm": 1.0, "pum": 1.0, "hm": 1.0}


@dataclass
class SEMImage:
    data: np.ndarray                       # float32 0..1, sem databar
    path: str = ""
    px_um: Optional[float] = None
    px_source: str = "não calibrado"
    meta: dict = field(default_factory=dict)
    raw: Optional[np.ndarray] = None       # imagem completa original (com databar), para conferência
    databar: Optional[tuple] = None        # (linha_inicial, altura) da barra removida
    scalebar: Optional[dict] = None        # barra de escala detectada (x0, x1, y, px, valor_um, texto)
    warnings: list = field(default_factory=list)

    @property
    def title(self):
        return os.path.splitext(os.path.basename(self.path))[0] if self.path else "imagem"

    @property
    def shape(self):
        return self.data.shape

    @property
    def field_um(self):
        if not self.px_um:
            return None
        return self.data.shape[1] * self.px_um, self.data.shape[0] * self.px_um

    def summary(self) -> str:
        m = self.meta
        L = [f"Arquivo     : {self.path}",
             f"Fabricante  : {m.get('vendor', '?')}  {m.get('instrument', '')}".rstrip(),
             f"Tamanho     : {self.data.shape[1]} × {self.data.shape[0]} px (original {m.get('raw_shape', '?')}, "
             f"{m.get('dtype', '?')})"]
        if self.databar:
            L.append(f"Databar     : {self.databar[1]} linhas removidas (a partir da linha {self.databar[0]})")
        if self.px_um:
            fw, fh = self.field_um
            L.append(f"Pixel       : {fmt_len(self.px_um)} /px  ({self.px_source})")
            L.append(f"Campo       : {fmt_len(fw)} × {fmt_len(fh)}")
        else:
            L.append("Pixel       : NÃO CALIBRADO - informe --bar-value 10um, --scale 10um=235px ou --pixel-size")
        if self.scalebar:
            sb = self.scalebar
            L.append(f"Barra       : {sb['px']:.0f} px" + (f" = {fmt_len(sb['value_um'])}" if sb.get('value_um') else "")
                     + (f"  (OCR: '{sb['text']}')" if sb.get("text") else ""))
        for k, lab in [("hv_kV", "Tensão"), ("mag", "Aumento"), ("wd_mm", "WD"), ("detector", "Detector"),
                       ("vacuum", "Vácuo"), ("pressure_Pa", "Pressão"), ("hfw_um", "HFW"), ("date", "Data"),
                       ("dwell_us", "Dwell"), ("tilt_deg", "Inclinação"), ("spot", "Spot"), ("current_pA", "Corrente")]:
            if m.get(k) not in (None, ""):
                v = m[k]
                unit = {"hv_kV": " kV", "wd_mm": " mm", "hfw_um": " µm", "pressure_Pa": " Pa", "dwell_us": " µs",
                        "tilt_deg": "°", "current_pA": " pA", "mag": "×"}.get(k, "")
                L.append(f"{lab:<12}: {v:g}{unit}" if isinstance(v, (int, float)) else f"{lab:<12}: {v}{unit}")
        for w in self.warnings:
            L.append(f"AVISO: {w}")
        return "\n".join(L)


# ----------------------------------------------------------------------------- utilidades

def fmt_len(um: float) -> str:
    """Formata comprimento em µm com a unidade mais legível."""
    if um is None or not np.isfinite(um):
        return "?"
    if um >= 1000:
        return f"{um / 1000:.3g} mm"
    if um >= 1:
        return f"{um:.4g} µm"
    return f"{um * 1000:.4g} nm"


def parse_length(s) -> float:
    """'10um', '10 µm', '500nm', '0.5 mm', 10 (µm) -> µm."""
    if isinstance(s, (int, float)):
        return float(s)
    t = str(s).strip().replace(",", ".").replace("μ", "µ")
    m = re.match(r"^\s*([0-9]*\.?[0-9]+(?:e-?\d+)?)\s*([a-zA-Zµ]*)\s*$", t)
    if not m:
        raise ValueError(f"comprimento inválido: {s!r} (use p.ex. 10um, 500nm, 0.2mm)")
    v, u = float(m.group(1)), (m.group(2) or "um").lower().replace("µ", "u")
    if u in ("u", "um", "micron", "microns", "mic"):
        return v
    if u == "nm":
        return v * 1e-3
    if u == "mm":
        return v * 1e3
    if u == "m":
        return v * 1e6
    raise ValueError(f"unidade desconhecida em {s!r}")


def parse_scale(s: str) -> float:
    """'10um=235px' -> µm/px."""
    a, _, b = str(s).partition("=")
    px = float(re.sub(r"[^0-9.]", "", b.replace(",", ".")))
    return parse_length(a) / px


def _to_float(a: np.ndarray):
    """Converte para cinza float32 0..1. RGB: luminância (ou o canal, se iguais)."""
    a = np.asarray(a)
    dtype = str(a.dtype)
    if a.ndim == 3 and a.shape[-1] in (3, 4):
        rgb = a[..., :3].astype(np.float64)
        if np.allclose(rgb[..., 0], rgb[..., 1]) and np.allclose(rgb[..., 1], rgb[..., 2]):
            g = rgb[..., 0]
        else:
            g = 0.299 * rgb[..., 0] + 0.587 * rgb[..., 1] + 0.114 * rgb[..., 2]
        a2 = g
    elif a.ndim == 3:            # pilha: primeira página
        a2 = a[0].astype(np.float64)
    else:
        a2 = a.astype(np.float64)
    if np.issubdtype(a.dtype, np.integer):
        mx = np.iinfo(a.dtype).max
        if a.dtype == np.uint16 and a2.max() <= 4095:         # 12 bits gravado em 16
            mx = 4095
        a2 = a2 / mx
    else:
        lo, hi = np.nanmin(a2), np.nanmax(a2)
        a2 = (a2 - lo) / (hi - lo) if hi > lo else a2 * 0
    return a2.astype(np.float32), dtype


# ----------------------------------------------------------------------------- metadados por fabricante

def _num(x):
    try:
        return float(str(x).replace(",", ".").split()[0])
    except Exception:
        return None


def _fei(md: dict, meta: dict):
    sc = md.get("Scan", {}) or md.get("EScan", {})
    px = sc.get("PixelWidth")
    eb = md.get("EBeam", {}) or {}
    bm = md.get("Beam", {}) or {}
    sysd = md.get("System", {}) or {}
    meta["vendor"] = "Thermo Fisher / FEI"
    meta["instrument"] = str(sysd.get("SystemType", "") or sysd.get("Type", "")).strip('"')
    hv = eb.get("HV", bm.get("HV"))
    if _num(hv):
        meta["hv_kV"] = round(_num(hv) / 1000, 3)
    wd = _num(eb.get("WD")) or _num((md.get("Stage", {}) or {}).get("WorkingDistance"))
    if wd and wd > 0:
        meta["wd_mm"] = round(wd * 1e3, 2)
    hfw = _num(eb.get("HFW")) or _num(sc.get("HorFieldsize"))
    if hfw:
        meta["hfw_um"] = round(hfw * 1e6, 4)
    det = md.get("Detectors", {}) or {}
    if det.get("Name"):
        meta["detector"] = f"{det.get('Name')} ({det.get('Mode', '')})".replace(" ()", "")
    vac = md.get("Vacuum", {}) or {}
    if vac.get("UserMode"):
        meta["vacuum"] = vac.get("UserMode")
    if _num(vac.get("ChPressure")) is not None:
        meta["pressure_Pa"] = _num(vac.get("ChPressure"))
    u = md.get("User", {}) or {}
    if u.get("Date"):
        meta["date"] = f"{u.get('Date')} {u.get('Time', '')}".strip()
    dw = _num(sc.get("Dwelltime")) or _num(sc.get("Dwell"))
    if dw:
        meta["dwell_us"] = round(dw * 1e6, 3)
    if _num(eb.get("BeamCurrent")):
        meta["current_pA"] = round(_num(eb["BeamCurrent"]) * 1e12, 2)
    img = md.get("Image", {}) or {}
    pf = md.get("PrivateFei", {}) or {}
    meta["fei_resy"] = img.get("ResolutionY")
    meta["fei_databar"] = pf.get("DatabarHeight")
    if hfw and _num(img.get("MagCanvasRealWidth")):
        meta["mag"] = round(_num(img["MagCanvasRealWidth"]) / hfw)
    px = _num(px)
    return px * 1e6 if px and px > 0 else None


def _zeiss(md: dict, meta: dict):
    def g(k):
        v = md.get(k)
        return v if isinstance(v, tuple) else None
    meta["vendor"] = "Zeiss"
    sv = md.get("sv_serial_number")
    if isinstance(sv, tuple):
        meta["instrument"] = str(sv[1])
    px = None
    for k in ("ap_image_pixel_size", "ap_pixel_size"):
        v = g(k)
        if v and len(v) >= 3 and _num(v[1]):
            px = _num(v[1]) * UNIT_UM.get(str(v[2]).replace("μ", "µ"), 1.0)
            if px:
                break
    for k, key, conv in [("ap_actualkv", "hv_kV", 1), ("ap_manualkv", "hv_kV", 1), ("ap_wd", "wd_mm", 1),
                         ("ap_mag", "mag", 1), ("ap_tilt_angle", "tilt_deg", 1)]:
        v = g(k)
        if v and key not in meta and _num(v[1]) is not None:
            val = _num(v[1])
            unit_s = (str(v[2]) if len(v) > 2 else str(v[1])).replace(" ", "")
            if key == "mag" and unit_s.upper().endswith("KX"):
                val *= 1000
            if k in ("ap_actualkv", "ap_manualkv") and len(v) > 2 and str(v[2]) == "V":
                val /= 1000
            meta[key] = val
    for k in ("dp_detector_channel", "dp_detector_type"):
        v = g(k)
        if v and "detector" not in meta:
            meta["detector"] = str(v[1])
    v = g("ap_date")
    if v:
        meta["date"] = str(v[1]) + (" " + str(g("ap_time")[1]) if g("ap_time") else "")
    v = g("dp_vp_mode") or g("ap_vp_mode")
    if v:
        meta["vacuum"] = str(v[1])
    v = g("ap_chamber_pressure")
    if v and _num(v[1]):
        meta["pressure_Pa"] = _num(v[1]) * (100.0 if "mbar" in str(v[-1]).lower() else 1.0)
    return px


def _ini_text(text: str) -> dict:
    out = {}
    for line in text.replace("\r", "\n").split("\n"):
        if "=" in line:
            k, _, v = line.partition("=")
            out[k.strip()] = v.strip()
    return out


def _tescan(text: str, meta: dict):
    d = _ini_text(text)
    meta["vendor"] = "TESCAN"
    meta["instrument"] = d.get("Device", d.get("DeviceModel", ""))
    if _num(d.get("HV")):
        meta["hv_kV"] = round(_num(d["HV"]) / 1000, 3)
    if _num(d.get("WD")):
        meta["wd_mm"] = round(_num(d["WD"]) * 1e3, 2)
    if d.get("Detector"):
        meta["detector"] = d.get("Detector")
    if d.get("Date"):
        meta["date"] = f"{d.get('Date')} {d.get('Time', '')}".strip()
    if d.get("VacMode"):
        meta["vacuum"] = d.get("VacMode")
    if _num(d.get("ChamberPressure")) is not None:
        meta["pressure_Pa"] = _num(d.get("ChamberPressure"))
    if _num(d.get("Magnification")):
        meta["mag"] = _num(d["Magnification"])
    px = _num(d.get("PixelSizeX"))
    return px * 1e6 if px else None


def _hitachi(text: str, meta: dict, img_w: int):
    d = _ini_text(text)
    meta["vendor"] = "Hitachi"
    meta["instrument"] = d.get("InstructName", "")
    v = d.get("AcceleratingVoltage", "")
    if _num(v) is not None:
        val = _num(v)
        meta["hv_kV"] = val / 1000 if ("volt" in v.lower() and "kv" not in v.lower()) or val > 100 else val
    v = d.get("WorkingDistance", "")
    if _num(v) is not None:
        val = _num(v)
        meta["wd_mm"] = val / 1000 if ("um" in v.lower() or val > 200) else val
    if _num(d.get("Magnification")):
        meta["mag"] = _num(d["Magnification"])
    for k in ("SignalName", "Signal"):
        if d.get(k):
            meta["detector"] = d[k]
            break
    if d.get("Date"):
        meta["date"] = f"{d.get('Date')} {d.get('Time', '')}".strip()
    if d.get("Vacuum"):
        meta["vacuum"] = d.get("Vacuum")
    px = _num(d.get("PixelSize"))
    if px:
        ds = d.get("DataSize", "")
        m = re.match(r"(\d+)\s*x\s*(\d+)", ds)
        if m and img_w and int(m.group(1)) != img_w:
            px *= int(m.group(1)) / img_w
        return px * 1e-3                  # nm -> µm
    return None


def _jeol(text: str, meta: dict, img_w: int):
    d = {}
    for line in text.replace("\r", "\n").split("\n"):
        line = line.strip()
        if line.startswith("$"):
            k, _, v = line.partition(" ")
            d[k.lstrip("$")] = v.strip()
    meta["vendor"] = "JEOL"
    meta["instrument"] = d.get("CM_INSTRUMENT", d.get("SM_INSTRUMENT", ""))
    if _num(d.get("CM_ACCEL_VOLT")) is not None:
        meta["hv_kV"] = _num(d["CM_ACCEL_VOLT"])
    if _num(d.get("CM_MAG")) is not None:
        meta["mag"] = _num(d["CM_MAG"])
    if _num(d.get("SM_WD")) is not None:
        meta["wd_mm"] = _num(d["SM_WD"])
    if d.get("CM_SIGNAL"):
        meta["detector"] = d["CM_SIGNAL"]
    if d.get("CM_DATE"):
        meta["date"] = f"{d.get('CM_DATE')} {d.get('CM_TIME', '')}".strip()
    full_w = None
    if d.get("CM_FULL_SIZE"):
        try:
            full_w = int(d["CM_FULL_SIZE"].split()[0])
        except Exception:
            pass
    bar = _num(d.get("SM_MICRON_BAR"))
    mk = d.get("SM_MICRON_MARKER")
    if bar and mk:
        try:
            px = parse_length(mk.replace(" ", "")) / bar
            if full_w and img_w and full_w != img_w:
                px *= full_w / img_w
            return px
        except ValueError:
            pass
    if _num(d.get("SM_PIXEL_SIZE")):                    # alguns modelos (nm)
        return _num(d["SM_PIXEL_SIZE"]) * 1e-3
    return None


def _sidecar_text(path: str):
    base = os.path.splitext(path)[0]
    for ext in (".txt", ".TXT", ".Txt", ".hdr"):
        p = base + ext
        if os.path.exists(p):
            for enc in ("utf-8", "latin-1", "utf-16"):
                try:
                    with open(p, encoding=enc) as f:
                        return f.read(), p
                except Exception:
                    continue
    return None, None


def _tif_metadata(path: str, meta: dict, warnings: list):
    import tifffile
    px, src = None, None
    with tifffile.TiffFile(path) as t:
        page = t.pages[0]
        try:
            arr = page.asarray()
        except ValueError as e:              # p.ex. TIF com compressão JPEG sem imagecodecs
            try:
                from PIL import Image
                arr = np.asarray(Image.open(path))
            except Exception:
                raise ValueError(f"{e} - instale com: pip install imagecodecs") from e
        h, w = page.shape[:2]
        # FEI/Thermo
        try:
            md = t.fei_metadata
        except Exception:
            md = None
        if md:
            px = _fei(md, meta)
            src = "metadados FEI/Thermo (Scan.PixelWidth)" if px else None
        if px is None:
            try:
                md = t.sem_metadata
            except Exception:
                md = None
            if md:
                px = _zeiss(md, meta)
                src = "metadados Zeiss (Image Pixel Size)" if px else None
        texts = []
        for tag in page.tags.values():
            v = tag.value
            if isinstance(v, bytes):
                try:
                    v = v.decode("latin-1")
                except Exception:
                    continue
            if isinstance(v, str) and len(v) > 20:
                texts.append((tag.code, v))
        if px is None:
            for code, s in texts:
                if code == 50431 or ("PixelSizeX=" in s and "[MAIN]" in s):
                    px = _tescan(s[s.find("["):] if "[" in s else s, meta)
                    src = "metadados TESCAN (PixelSizeX)" if px else None
                    break
        if px is None:
            for code, s in texts:                  # FEI em texto (Phenom/Quanta antigos não decodificados)
                if "PixelWidth=" in s:
                    d = _ini_text(s)
                    pw = _num(d.get("PixelWidth"))
                    if pw:
                        meta.setdefault("vendor", "Thermo Fisher / FEI / Phenom")
                        px, src = pw * 1e6, "metadados FEI/Phenom (PixelWidth, texto)"
                        if _num(d.get("HV")):
                            meta["hv_kV"] = _num(d["HV"]) / 1000
                        break
                if "[SemImageFile]" in s:
                    px = _hitachi(s, meta, w)
                    src = "metadados Hitachi embutidos (PixelSize)" if px else None
                    if px:
                        break
                if "$CM_" in s or "$$SM_" in s:
                    px = _jeol(s, meta, w)
                    src = "metadados JEOL embutidos (micron bar)" if px else None
                    if px:
                        break
        if px is None and t.is_imagej:
            ij = t.imagej_metadata or {}
            unit = str(ij.get("unit", "")).replace("\\u00B5", "µ")
            xr = page.tags.get(282)
            if unit and unit.lower() not in ("pixel", "pixels", "") and xr:
                num, den = xr.value
                if num:
                    try:
                        px = (den / num) * parse_length(f"1{unit.replace('micron', 'um')}")
                        src = f"calibração ImageJ ({unit})"
                        meta.setdefault("vendor", "ImageJ/Fiji")
                    except ValueError:
                        px = None
        if px is None and t.is_ome:
            m = re.search(r'PhysicalSizeX="([0-9.eE-]+)"(?:\s+PhysicalSizeXUnit="([^"]+)")?', t.ome_metadata or "")
            if m:
                px = float(m.group(1)) * UNIT_UM.get((m.group(2) or "µm").replace("μ", "µ"), 1.0)
                src = "metadados OME"
        if px is None and t.is_sis:
            sis = t.sis_metadata or {}
            if sis.get("pixelsizex"):
                px, src = float(sis["pixelsizex"]) * 1e6, "metadados Olympus SIS"
    return arr, px, src


# ----------------------------------------------------------------------------- databar e barra de escala

def detect_databar(a: np.ndarray, max_frac: float = 0.3, min_rows: int = 8):
    """Detecta a barra de informações (faixa uniforme com texto) na parte de baixo (ou de cima).
    Retorna (linha_inicial, altura, lado) ou None. a: float 0..1."""
    H, W = a.shape
    q = np.round(a * 255).astype(np.int32)

    def row_stats(r):
        vals, cnt = np.unique(q[r], return_counts=True)
        j = int(np.argmax(cnt))
        return vals[j], cnt[j] / W

    for side in ("bottom", "top"):
        order = list(range(H - 1, -1, -1)) if side == "bottom" else list(range(H))
        st = [row_stats(r) for r in order[:min(14, H)]]
        cand = [v for v, f in st[3:] if f >= 0.4]
        if not cand:
            continue
        vv, cc = np.unique(cand, return_counts=True)
        mode0 = int(vv[np.argmax(cc)])
        n = 0
        started = False
        for k, r in enumerate(order):
            v, f = row_stats(r)
            if f >= 0.4 and abs(int(v) - mode0) <= 3:
                n += 1
                started = True
                continue
            if f >= 0.9:                     # linha de borda/separadora
                n += 1
                if started and n >= min_rows:
                    break
                continue
            if not started and k < 3:        # borda de 1-3 linhas não uniforme
                n += 1
                continue
            break
        if mode0 is None or n < min_rows or n > max_frac * H:
            continue
        # a imagem logo acima/abaixo não pode ser também uniforme (evita fundo preto saturado)
        r_in = H - n - 3 if side == "bottom" else n + 2
        if 0 <= r_in < H:
            v, f = row_stats(r_in)
            if f > 0.6 and abs(int(v) - int(mode0)) <= 3:
                continue
        # tem que ter algum "texto" (pixels contrastantes) dentro da faixa
        band = q[H - n:] if side == "bottom" else q[:n]
        contrast = np.mean(np.abs(band - mode0) > 60)
        if contrast < 0.005:
            continue
        start = H - n if side == "bottom" else 0
        return start, n, side
    return None


def find_scalebar(band: np.ndarray, bg=None):
    """Procura a barra de escala (linha horizontal fina, possivelmente dividida pelo texto ou
    tracejada) numa região (float 0..1). Retorna dict(x0, x1, y, px, bbox) em coordenadas da região."""
    from scipy import ndimage as ndi
    H, W = band.shape
    q = np.round(band * 255).astype(np.int32)
    if bg is None:
        vals, cnt = np.unique(q, return_counts=True)
        bg = vals[np.argmax(cnt)]
    fg = np.abs(q - bg) > 60
    lab, n = ndi.label(fg, structure=np.ones((3, 3)))
    if n == 0:
        return None
    objs = ndi.find_objects(lab)
    cands = []
    for i, sl in enumerate(objs, start=1):
        ys, xs = sl
        h, w = ys.stop - ys.start, xs.stop - xs.start
        if w < 3 or w > 0.75 * W:
            continue
        m = lab[sl] == i
        runs = m.sum(axis=1)
        k = int(np.argmax(runs))
        run = runs[k]
        # linha principal: corrida longa e espessura pequena nas colunas da corrida
        thick = int((runs >= 0.7 * run).sum())
        if run < max(3, 0.8 * w) or thick > max(8, 0.25 * H):
            continue
        if xs.start <= 1 or xs.stop >= W - 1:
            continue
        # pixels fora da faixa da linha só podem ser ticks/pontas nas extremidades (não um glifo com "pé")
        band_rows = np.where(runs >= 0.7 * run)[0]
        off = m.copy()
        off[max(0, band_rows.min() - 1):band_rows.max() + 2] = False
        if off.any():
            cols = np.where(off.any(axis=0))[0]
            if np.any((cols > 4) & (cols < w - 5)):
                continue
        cands.append(dict(x0=xs.start, x1=xs.stop - 1, y=ys.start + k, w=w, thick=thick, h=h))
    if not cands:
        return None
    # agrupa segmentos colineares (barra dividida pelo texto, ou tracejada)
    cands.sort(key=lambda c: c["y"])
    groups = []
    for c in cands:
        for g in groups:
            if abs(g[0]["y"] - c["y"]) <= 2 and abs(g[0]["thick"] - c["thick"]) <= 2:
                g.append(c)
                break
        else:
            groups.append([c])
    best = None
    for g in groups:
        g.sort(key=lambda c: c["x0"])
        big = [c for c in g if c["w"] >= 10]
        small = [c for c in g if c["w"] < 10]
        # barra tracejada (Hitachi/JEOL antigos): >= 5 traços de largura parecida, espaçamento regular
        if len(small) >= 5 and not big:
            ws = np.array([c["w"] for c in small])
            gaps = np.diff([c["x0"] for c in small])
            if ws.std() <= 0.35 * ws.mean() + 1 and gaps.std() <= 0.35 * gaps.mean() + 1:
                cand = _mk(small)
                if best is None or cand["px"] > best["px"]:
                    best = cand
            continue
        if not big:
            continue
        segs = [big[0]]
        for c in big[1:] + [None]:
            if c is not None and c["x0"] - segs[-1]["x1"] <= 1.5 * max(s["w"] for s in segs + [c]):
                segs.append(c)
                continue
            cand = _mk(segs)
            if cand["px"] >= 15 and (best is None or cand["px"] > best["px"]):
                best = cand
            if c is not None:
                segs = [c]
    if best is None or best["px"] < 15:
        return None
    return best


def _is_bar_like(m: np.ndarray) -> bool:
    """Componente é um segmento de barra: linha horizontal longa e fina; pixels fora da linha só nas
    extremidades (ticks/pontas de seta)."""
    h, w = m.shape
    runs = m.sum(axis=1)
    run = runs.max()
    if run < max(3, 0.8 * w):
        return False
    band_rows = np.where(runs >= 0.7 * run)[0]
    if band_rows.max() - band_rows.min() + 1 > 8:
        return False
    off = m.copy()
    off[max(0, band_rows.min() - 1):band_rows.max() + 2] = False
    if off.any():
        cols = np.where(off.any(axis=0))[0]
        if np.any((cols > 4) & (cols < w - 5)):
            return False
    return True


def _mk(segs):
    x0, x1 = segs[0]["x0"], max(s["x1"] for s in segs)
    return dict(x0=x0, x1=x1, y=int(np.median([s["y"] for s in segs])), px=float(x1 - x0),
                nseg=len(segs), thick=segs[0]["thick"])


def _text_line_near_bar(band: np.ndarray, bar: dict):
    """Isola a linha de texto do rótulo da barra (componentes que não são a barra, perto dela)."""
    from scipy import ndimage as ndi
    H, W = band.shape
    q = np.round(band * 255).astype(np.int32)
    vals, cnt = np.unique(q, return_counts=True)
    bg = vals[np.argmax(cnt)]
    fg = np.abs(q - bg) > 60
    lab, n = ndi.label(fg, structure=np.ones((3, 3)))
    pad = 0.6 * bar["px"] + 10
    comps = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        ys, xs = sl
        h, w = ys.stop - ys.start, xs.stop - xs.start
        if h < 4 or h > 0.8 * H or w > 0.6 * bar["px"] + 20:
            continue
        xc = 0.5 * (xs.start + xs.stop)
        if xc < bar["x0"] - pad or xc > bar["x1"] + pad:
            continue
        xc = 0.5 * (xs.start + xs.stop)
        # descarta a própria barra (segmentos com corrida horizontal longa) e seus ticks
        mk = lab[sl] == i
        if w >= 10 and _is_bar_like(mk):
            continue
        if w <= 3 and (abs(xc - bar["x0"]) < 6 or abs(xc - bar["x1"]) < 6):
            continue
        comps.append((ys.start, ys.stop, xs.start, xs.stop, i))
    if not comps:
        return None
    # agrupa por linha (sobreposição vertical) e escolhe a mais próxima da barra
    comps.sort()
    lines = []
    for c in comps:
        for L in lines:
            if c[0] < L["y1"] - 1 and c[1] > L["y0"] + 1:
                L["c"].append(c)
                L["y0"], L["y1"] = min(L["y0"], c[0]), max(L["y1"], c[1])
                break
        else:
            lines.append(dict(y0=c[0], y1=c[1], c=[c]))
    lines = [L for L in lines if len(L["c"]) >= 2]
    if not lines:
        return None
    L = min(lines, key=lambda L: abs((L["y0"] + L["y1"]) / 2 - bar["y"]))
    # separa a linha em "palavras" (lacuna > ~1 altura de caractere) e fica com a mais próxima da barra
    cs = sorted(L["c"], key=lambda c: c[2])
    hchar = np.median([c[1] - c[0] for c in cs])
    words = [[cs[0]]]
    for c in cs[1:]:
        if c[2] - max(w[3] for w in words[-1]) > 1.1 * hchar:
            words.append([c])
        else:
            words[-1].append(c)

    def dist(w):
        x0, x1 = min(c[2] for c in w), max(c[3] for c in w)
        if x1 >= bar["x0"] and x0 <= bar["x1"]:
            return 0.0
        return float(min(abs(x0 - bar["x1"]), abs(bar["x0"] - x1)))
    # rótulo = palavra mais próxima + vizinha imediata (número e unidade às vezes separados)
    words.sort(key=lambda w: min(c[2] for c in w))
    k = int(np.argmin([dist(w) for w in words]))
    sel = list(words[k])
    for j in (k - 1, k + 1):
        if 0 <= j < len(words):
            gap = (min(c[2] for c in words[j]) - max(c[3] for c in words[k])) if j > k else \
                  (min(c[2] for c in words[k]) - max(c[3] for c in words[j]))
            if gap < 2.5 * hchar and dist(words[j]) <= dist(words[k]) + 4 * hchar:
                sel += words[j]
    x0 = min(c[2] for c in sel)
    x1 = max(c[3] for c in sel)
    mask = np.isin(lab[L["y0"]:L["y1"], x0:x1], [c[4] for c in sel])
    return mask


def _nice(v):
    m = v / 10 ** np.floor(np.log10(v))
    return any(abs(m - k) < 1e-6 for k in (1, 2, 2.5, 3, 4, 5))


def ocr_scalebar(band: np.ndarray, bar: dict):
    """Lê o valor da barra de escala por OCR (tesseract), com várias versões ampliadas do rótulo
    isolado e votação. Retorna (valor_um, texto) ou (None, motivo)."""
    try:
        import pytesseract
        from PIL import Image, ImageFilter, ImageOps
        pytesseract.get_tesseract_version()
    except Exception:
        return None, "OCR indisponível (instale tesseract + pytesseract) - use --bar-value"
    mask = _text_line_near_bar(band, bar)
    if mask is None:
        return None, "rótulo da barra não encontrado - use --bar-value"
    g = np.where(mask, 0, 255).astype(np.uint8)
    h = g.shape[0]
    votes, texts = {}, []
    for target in (40, 60, 90):
        sc = max(2, int(round(target / max(h, 1))))
        for blur in (0.3, 0.5):
            im = Image.fromarray(g).resize((g.shape[1] * sc, g.shape[0] * sc), Image.NEAREST)
            im = im.filter(ImageFilter.GaussianBlur(sc * blur)).point(lambda v: 0 if v < 140 else 255)
            im = ImageOps.expand(im, border=30, fill=255)
            try:
                txt = pytesseract.image_to_string(im, config="--psm 7").strip()
            except Exception as e:
                return None, f"OCR falhou: {e}"
            texts.append(txt)
            v = _parse_bar_text(txt)
            if v:
                key = round(v, 9)
                votes[key] = votes.get(key, 0) + (2 if _nice(v) else 1)
    if not votes:
        return None, "OCR não leu número + unidade (" + " | ".join(t for t in texts[:3] if t) + ") - use --bar-value"
    best = max(votes, key=votes.get)
    txt = next(t for t in texts if _parse_bar_text(t) and round(_parse_bar_text(t), 9) == best)
    return best, txt


def _parse_bar_text(txt: str):
    t = txt.replace("μ", "µ").replace(",", ".")
    t = re.sub(r"(?<=[0-9SOo])[Oo](?=[0-9Oo\s])|^[Oo]|(?<=\s)[Oo]", "0", t)
    t = re.sub(r"\bS(?=[0O0-9]{1,3}\b)|(?<=\b)S(?=\d)", "5", t)
    t = t.replace("O", "0").replace("l", "1").replace("I", "1")
    t = re.sub(r"(?<=\d)\s+(?=\d)", "", t)
    for m in re.finditer(r"(\d+(?:\.\d+)?)\s*(nm|µm|um|urn|pm|prn|pum|pur|jim|ym|wm|hm|mm|m\b)", t, re.I):
        v, u = float(m.group(1)), m.group(2).lower()
        if v <= 0:
            continue
        f = {"nm": 1e-3, "mm": 1e3}.get(u, 1.0)        # pm/prn/urn/ym/jim = µm mal lido; "m" sozinho -> µm
        return v * f
    m = re.search(r"(\d+(?:\.\d+)?)\s*(µ|u|y|p)\b", t)          # 'm' perdido: '200 y', '10 p'
    if m and float(m.group(1)) > 0:
        return float(m.group(1))
    return None


# ----------------------------------------------------------------------------- carga principal

def load(path: str, pixel_size=None, scale=None, bar_value=None, databar="auto", crop=None,
         ocr: bool = True) -> SEMImage:
    """Carrega uma imagem de MEV, separa a databar e determina a escala.
    pixel_size: '0.05um' | µm/px; scale: '10um=235px'; bar_value: '10um' (valor da barra desenhada);
    databar: 'auto' | 0 (não remover) | N (linhas da base a remover); crop: (x0, y0, x1, y1) px."""
    meta, warnings = {}, []
    ext = os.path.splitext(path)[1].lower()
    px, src = None, None
    if ext in (".tif", ".tiff"):
        arr, px, src = _tif_metadata(path, meta, warnings)
    else:
        from PIL import Image
        arr = np.asarray(Image.open(path))
        meta["vendor"] = "?"
    data, dtype = _to_float(arr)
    meta["dtype"] = dtype
    meta["raw_shape"] = f"{data.shape[1]}×{data.shape[0]}"
    H, W = data.shape
    if px is None:
        txt, p = _sidecar_text(path)
        if txt:
            if "[SemImageFile]" in txt or "PixelSize=" in txt:
                px = _hitachi(txt, meta, W)
                src = f"arquivo Hitachi {os.path.basename(p)} (PixelSize)" if px else None
            elif "$CM_" in txt or "$$SM_" in txt or "$SM_" in txt:
                px = _jeol(txt, meta, W)
                src = f"arquivo JEOL {os.path.basename(p)} (micron bar)" if px else None
    # databar
    db = None
    if databar == "auto":
        rh = meta.get("fei_databar")
        resy = meta.get("fei_resy")
        if resy and int(resy) < H:
            db = (int(resy), H - int(resy), "bottom")
        else:
            db = detect_databar(data)
    elif databar and int(databar) > 0:
        db = (H - int(databar), int(databar), "bottom")
    raw = data
    scalebar = None
    if db:
        start, n, side = db
        band = data[start:start + n] if side == "bottom" else data[:n]
        data = data[:start] if side == "bottom" else data[n:]
        bar = find_scalebar(band)
        if bar:
            scalebar = dict(bar)
            if bar_value is not None:
                scalebar["value_um"] = parse_length(bar_value)
            elif ocr:
                v, txt = ocr_scalebar(band, bar)
                scalebar["value_um"], scalebar["text"] = v, txt if v else None
                if v is None and px is None:
                    warnings.append(txt)
    elif bar_value is not None:
        # sem databar: procura a barra numa faixa inferior de 15% da imagem
        y0 = int(0.85 * H)
        bar = find_scalebar(data[y0:])
        if bar:
            scalebar = dict(bar, value_um=parse_length(bar_value))
            scalebar["y"] += y0
            warnings.append("barra de escala sobreposta à imagem (sem databar): recorte a região antes de medir")
    # escala final
    bar_px = (scalebar["value_um"] / scalebar["px"]) if scalebar and scalebar.get("value_um") else None
    if pixel_size is not None:
        px, src = parse_length(pixel_size), "manual (--pixel-size)"
    elif scale is not None:
        px, src = parse_scale(scale), f"manual (--scale {scale})"
    elif bar_value is not None and bar_px:
        px, src = bar_px, f"barra de escala medida ({scalebar['px']:.0f} px = {fmt_len(scalebar['value_um'])})"
    elif px is None and bar_px:
        px, src = bar_px, (f"barra de escala ({scalebar['px']:.0f} px) + OCR '{scalebar.get('text')}'"
                           " - CONFIRA o valor lido")
    if scalebar and scalebar.get("nseg", 1) >= 5:
        warnings.append("barra de escala tracejada: comprimento medido do 1º ao último traço - confira")
    if px and bar_px and src and "barra" not in src:
        dev = abs(bar_px - px) / px
        tol = 0.10 if scalebar.get("nseg", 1) >= 5 else 0.03
        if dev > tol:
            warnings.append(f"escala dos metadados ({fmt_len(px)}/px) difere {100 * dev:.0f}% da barra desenhada "
                            f"({fmt_len(bar_px)}/px) - confira (OCR ou imagem redimensionada?)")
    if px is None:
        warnings.append("imagem sem escala: medidas sairão em pixels")
    if crop:
        x0, y0, x1, y1 = [int(v) for v in crop]
        data = data[y0:y1, x0:x1]
    if px and meta.get("hfw_um") and "FEI" in str(meta.get("vendor")):
        exp = meta["hfw_um"] / W
        if abs(exp - px) / px > 0.02:
            warnings.append("HFW/largura não bate com PixelWidth - imagem redimensionada?")
    sat_hi = float(np.mean(data >= 0.995))
    sat_lo = float(np.mean(data <= 0.005))
    meta["sat_hi"], meta["sat_lo"] = sat_hi, sat_lo
    if sat_hi > 0.02:
        warnings.append(f"{100 * sat_hi:.1f}% dos pixels saturados no branco (carregamento? contraste alto?)")
    if sat_lo > 0.10:
        warnings.append(f"{100 * sat_lo:.1f}% dos pixels no preto absoluto (fundo saturado/contraste)")
    return SEMImage(data=np.ascontiguousarray(data), path=path, px_um=px, px_source=src or "não calibrado",
                    meta=meta, raw=raw, databar=(db[0], db[1]) if db else None, scalebar=scalebar,
                    warnings=warnings)


def save_tif(im: SEMImage, path: str):
    """Grava a imagem (sem databar) como TIF 8 bits calibrado no formato ImageJ (abre com escala no
    Fiji/ImageJ)."""
    import tifffile
    a = np.clip(im.data * 255 + 0.5, 0, 255).astype(np.uint8)
    if im.px_um:
        tifffile.imwrite(path, a, imagej=True, resolution=(1 / im.px_um, 1 / im.px_um),
                         metadata={"unit": "um"})
    else:
        tifffile.imwrite(path, a)
