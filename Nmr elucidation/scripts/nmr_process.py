"""
nmr_process.py - processamento 1D: FID -> espectro (filtro digital Bruker, apodização, zero-filling,
FFT, fase automática), linha de base, ruído/SNR, referência pelo solvente/TMS

    from nmr_process import prepare, check_reference, reference, baseline, noise
    s = prepare(s)                 # linha de base + checagem/ajuste de referência (o que os comandos usam)
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from nmr_data import SOLVENTS
from nmr_io import NMRSpectrum, _f, _s, _shift_reference


# ----------------------------------------------------------------------------- FID

def _phase(spec, p0, p1):
    n = spec.size
    return spec * np.exp(1j * np.deg2rad(p0 + p1 * np.arange(n) / n))


def _acme_score(sp, gamma=1000.0):
    r = sp.real
    d = np.abs(np.diff(r))
    p = d / (d.sum() + 1e-30)
    ent = -np.sum(p * np.log(p + 1e-30))
    neg = np.where(r < 0, r, 0)
    pen = gamma * np.sum(neg ** 2) / (np.sum(r ** 2) + 1e-30)
    return ent + pen


def autophase(spec, p0=0.0, p1=0.0, fit_p1=True):
    """Fase de ordem 0 e 1 por ACME (Chen et al., J. Magn. Reson. 158 (2002) 164) em versão decimada."""
    from scipy.optimize import minimize
    k = max(1, spec.size // 16384)
    sub = spec[::k]

    def f(x):
        return _acme_score(_phase(sub, x[0], x[1] if fit_p1 else 0.0))
    best = None
    starts = [(p0, p1), (0.0, 0.0)] + [(a, p1) for a in (90.0, 180.0, 270.0)]
    for st in starts:
        r = minimize(f, x0=np.array(st, float), method="Nelder-Mead",
                     options=dict(xatol=0.05, fatol=1e-7, maxiter=1500))
        if best is None or r.fun < best.fun:
            best = r
    a0, a1 = best.x[0], (best.x[1] if fit_p1 else 0.0)
    out = _phase(spec, a0, a1)
    if np.sum(out.real) < 0:                       # espectro invertido
        out = -out
        a0 += 180
    return out, float(a0 % 360), float(a1)


def fid_to_spectrum(fid, dic, nucleus="1H", title="", path="", solvent=None, experiment="1D", meta=None,
                    warnings=None, lb=None, zf=2, phase="auto", filter_removed=False) -> NMRSpectrum:
    import nmrglue as ng
    meta = dict(meta or {})
    warnings = list(warnings or [])
    fid = np.asarray(fid, complex)
    n = fid.size
    # filtro digital Bruker
    grp = _f(dic.get("$GRPDLY"))
    decim, dspfvs = _f(dic.get("$DECIM")), _f(dic.get("$DSPFVS"))
    if not filter_removed and (grp is not None or decim is not None):
        try:
            if grp is not None and grp > 0:
                fid = ng.bruker.rm_dig_filter(fid, int(decim or 1), int(dspfvs or 10), grp)
            elif decim and dspfvs:
                fid = ng.bruker.rm_dig_filter(fid, int(round(decim)), int(dspfvs))
        except Exception as e:
            warnings.append(f"não removi o filtro digital Bruker ({e}): a fase de 1ª ordem compensa")
    n = fid.size
    # largura espectral: eixo de tempo do JCAMP (mais confiável) > $SW_h/$SWH > $SW (ppm) × SF
    swh = None
    fl = _s(dic.get("FIRST")).split(",")
    ll = _s(dic.get("LAST")).split(",")
    nt = None
    for key in ("VARDIM", "VAR_DIM"):
        if dic.get(key):
            try:
                nt = int(float(_s(dic.get(key)).split(",")[0]))
            except ValueError:
                pass
    try:
        dw = (float(ll[0]) - float(fl[0])) / ((nt or n) - 1)
        if dw > 0:
            swh = 1.0 / dw
    except (ValueError, IndexError, ZeroDivisionError):
        pass
    if swh is None:
        swh = _f(dic.get("$SW_h")) or _f(dic.get("$SWH"))
    if swh is None and _f(dic.get("$SW")) and (_f(dic.get("$SF")) or _f(dic.get(".OBSERVEFREQUENCY"))):
        swh = _f(dic.get("$SW")) * (_f(dic.get("$SF")) or _f(dic.get(".OBSERVEFREQUENCY")))
    sf = _f(dic.get("$SF")) or _f(dic.get("$BF1")) or _f(dic.get(".OBSERVEFREQUENCY"))
    sfo1 = _f(dic.get("$SFO1")) or _f(dic.get(".OBSERVEFREQUENCY"))
    if lb is None:
        lb = 0.3 if nucleus == "1H" else 1.0
    t = np.arange(n) / swh
    fid = fid * np.exp(-np.pi * lb * t)
    fid[0] *= 0.5
    size = int(2 ** np.ceil(np.log2(n)) * zf)
    fid = np.concatenate([fid, np.zeros(size - n, complex)])
    spec = np.fft.fftshift(np.fft.fft(fid))[::-1]          # frequência decrescente (ppm alto à esquerda)
    if _s(dic.get("$REVERSE")).lower() in ("yes", "true", "1"):
        spec = spec[::-1]
    p0 = _f(dic.get("$PHC0"), 0.0)
    p1 = _f(dic.get("$PHC1"), 0.0)
    if phase == "auto":
        spec, a0, a1 = autophase(spec, 0.0, 0.0)
        meta["phase_desc"] = f"automática ACME (p0 = {a0:.1f}°, p1 = {a1:.1f}°)"
    elif isinstance(phase, (tuple, list)):
        spec = _phase(spec, *phase)
        meta["phase_desc"] = f"manual (p0 = {phase[0]}°, p1 = {phase[1]}°)"
    # eixo
    f = (np.arange(size) - size / 2) * swh / size
    f = f[::-1]                                              # Hz relativos à portadora, decrescente
    if "_varian" in dic:
        rfl, rfp, sw, sfrq = dic["_varian"]
        ppm = (f + sw / 2 - rfl + rfp) / sfrq if rfl else f / sfrq
        ppm = ppm - (ppm.max() - (sw - rfl + rfp) / sfrq) if rfl else ppm
        sf = sfrq
        axis = "Varian rfl/rfp"
    elif sf and sfo1 and abs(sfo1 - sf) > 1e-9 and ("$SF" in dic or "$BF1" in dic):
        center = (sfo1 - sf) / sf * 1e6
        ppm = center + f / sf
        axis = "portadora ($SFO1 − $SF)"
    elif _f(dic.get("$OFFSET")) is not None and sf:
        ppm = _f(dic.get("$OFFSET")) - (f[0] - f) / sf
        axis = "$OFFSET"
    else:
        ref = _shift_reference(dic)
        sf = sf or 400.0
        if ref is not None:
            pt, val = ref
            i = min(max(int(round(pt * size / n)) - 1, 0), size - 1)
            ppm = val + (f - f[i]) / sf
            axis = f".SHIFT REFERENCE ({val:g} ppm)"
        else:
            ppm = f / sf
            axis = "sem referência (centro = 0 ppm)"
            warnings.append("FID sem informação de referência: use --ref (ex.: --ref 7.31=7.26)")
    meta.update(from_fid=True, lb=lb, zf=zf, axis=axis, format=meta.get("format", "JCAMP-DX NMR FID"))
    s = NMRSpectrum(ppm=ppm, y=spec.real.copy(), yi=spec.imag.copy(), nucleus=nucleus, sf=float(sf), solvent=solvent,
                    title=title, path=path, experiment=experiment, meta=meta, warnings=warnings)
    return s


# ----------------------------------------------------------------------------- ruído e linha de base

def noise(y) -> float:
    """σ do ruído: MAD robusto do espectro (a maior parte dos pontos é linha de base). Usa também o MAD
    das diferenças (insensível a deriva) e fica com o maior - o das diferenças subestima σ quando há
    zero-filling (pontos vizinhos correlacionados)."""
    y = np.asarray(y, float)
    a = 1.4826 * np.median(np.abs(y - np.median(y)))
    d = np.diff(y)
    b = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    return float(max(a, b))


def _whittaker(y, w, lam):
    from scipy import sparse
    from scipy.sparse.linalg import spsolve
    n = y.size
    D = sparse.diags([1.0, -2.0, 1.0], [0, -1, -2], shape=(n, n - 2))
    W = sparse.diags(w)
    Z = W + lam * D.dot(D.T)
    return spsolve(Z.tocsc(), w * y)


def baseline(s: NMRSpectrum, lam=None, iters=4, k=3.0):
    """Linha de base por suavização de Whittaker ponderada (sinais excluídos iterativamente; funciona com
    sinais positivos e negativos, p.ex. DEPT). Retorna (y_corrigido, linha_de_base)."""
    y = np.asarray(s.y, float)
    n = y.size
    step = max(1, n // 8192)
    yd = ndi.uniform_filter1d(y, step)[::step] if step > 1 else y.copy()
    sig = noise(y) * (1 if step == 1 else 1 / np.sqrt(step))
    w = np.ones_like(yd)
    lam = lam or (1e6 if s.nucleus == "1H" else 1e7) / max(1, (yd.size / 8192) ** 4)
    base = np.full_like(yd, np.median(yd))
    pad = max(3, yd.size // 400)
    for _ in range(iters):
        r = yd - base
        sigr = 1.4826 * np.median(np.abs(r[w > 0] - np.median(r[w > 0]))) if (w > 0).any() else sig
        sigr = max(sigr, sig)
        sig_mask = np.abs(r) > k * sigr
        sig_mask = ndi.binary_dilation(sig_mask, iterations=pad)
        w = (~sig_mask).astype(float)
        if w.sum() < 10:
            break
        base = _whittaker(yd, w, lam)
    if step > 1:
        xs = np.arange(yd.size) * step + (step - 1) / 2
        base = np.interp(np.arange(n), xs, base)
    return y - base, base


# ----------------------------------------------------------------------------- referência

def find_peak_near(s: NMRSpectrum, target, tol):
    m = (s.ppm >= target - tol) & (s.ppm <= target + tol)
    if not m.any():
        return None
    idx = np.where(m)[0]
    j = idx[np.argmax(s.y[idx])]
    return float(s.ppm[j]), float(s.y[j])


def _local_max_near(s, x, tol_ppm):
    m = (s.ppm >= x - tol_ppm) & (s.ppm <= x + tol_ppm)
    if not m.any():
        return None
    idx = np.where(m)[0]
    j = idx[np.argmax(s.y[idx])]
    return float(s.ppm[j]), float(s.y[j])


def _match_pattern(s, center, J, rel, sigma, win):
    """Procura o padrão do solvente (multipleto por acoplamento com ²H) perto de 'center'. Retorna
    (posição_central, escore) do melhor candidato ou None."""
    m = (s.ppm >= center - win) & (s.ppm <= center + win)
    if not m.any():
        return None
    idx = np.where(m)[0]
    yy = s.y[idx]
    # candidatos: máximos locais acima de 5σ
    cand = [i for k, i in enumerate(idx[1:-1], start=1) if yy[k] >= yy[k - 1] and yy[k] >= yy[k + 1] and yy[k] > 5 * sigma]
    if not cand:
        return None
    n = len(rel)
    dp = J / s.sf
    rel = np.asarray(rel, float) / max(rel)
    res = []
    step = abs(s.ppm[1] - s.ppm[0])
    for i in cand:
        c = s.ppm[i]
        if n == 1 or dp < 2.5 * step:
            # singleto ou padrão não resolvido: precisa estar isolado (sem linhas de multipleto a 1,5-10 Hz)
            h = s.y[i]
            lo_, hi_ = 1.5 / s.sf, 10.0 / s.sf
            near = (np.abs(s.ppm - c) >= lo_) & (np.abs(s.ppm - c) <= hi_)
            if near.any():
                yy2 = s.y[near]
                # máximos locais nessa coroa
                loc = (yy2[1:-1] >= yy2[:-2]) & (yy2[1:-1] >= yy2[2:]) & (yy2[1:-1] > 0.3 * h)
                if loc.any():
                    continue
            res.append((c, 1.0 if n == 1 else 0.85, -h * np.exp(-((c - center) / 0.05) ** 2)))
            continue
        hs = []
        for k in range(n):
            x = c + (k - (n - 1) / 2) * dp
            p = _local_max_near(s, x, max(0.12 * dp, 1.5 * step))
            hs.append(p[1] if p else 0.0)
        hs = np.array(hs)
        if hs.max() <= 0 or np.argmax(hs) != (n - 1) // 2 and n > 1:
            continue
        if s.y[i] < 0.5 * yy.max():
            continue
        if (hs[rel > 0.3] < 3 * sigma).any():
            continue
        h = hs / hs.max()
        err = np.mean(np.abs(h - rel))
        if err < 0.25:
            res.append((c, 1 - err, abs(c - center)))
    if not res:
        return None
    res.sort(key=lambda r: (-round(r[1], 1), r[2]))
    return res[0][0], res[0][1]


def check_reference(s: NMRSpectrum, sigma=None) -> dict:
    """Confere a posição do sinal residual do solvente reconhecendo o padrão de acoplamento com ²H
    (CDCl3 ¹³C tripleto 1:1:1; DMSO-d6 septeto; CHD2 quintetos em ¹H...). Retorna dict com o desvio."""
    from nmr_data import SOLVENT_PATTERNS
    sigma = sigma or noise(s.y)
    res = dict(solvent=s.solvent, expected=None, found=None, delta=0.0, ok=True, matched=False, note="")
    pats = SOLVENT_PATTERNS.get(s.nucleus, {}).get(s.solvent or "")
    if s.experiment.startswith("DEPT"):
        res["note"] = "DEPT não mostra o solvente: use a referência do ¹³C correspondente"
        return res
    if not pats:
        res["note"] = ("solvente desconhecido" if not s.solvent else f"sem padrão de referência para {s.solvent} em "
                       f"{s.nucleus}") + ": referência do arquivo mantida"
        return res
    win = 0.12 if s.nucleus == "1H" else 1.2
    found = None
    for (c, J, rel) in pats:
        r = _match_pattern(s, c, J, rel, sigma, win)
        if r and (found is None or r[1] > found[2]):
            found = (c, r[0], r[1])
    if found is None:
        res["note"] = "sinal residual do solvente não reconhecido perto do esperado - referência do arquivo mantida"
        return res
    c, pos, score = found
    res.update(expected=c, found=pos, delta=pos - c, matched=score > 0.8)
    lim = 0.015 if s.nucleus == "1H" else 0.2
    if abs(res["delta"]) > lim:
        res["ok"] = False
        res["note"] = f"solvente em {pos:.3f} ppm (esperado {c:.2f}; desvio {res['delta']:+.3f} ppm)"
    else:
        res["note"] = f"solvente em {pos:.3f} ppm (esperado {c:.2f}) - referência OK"
    # TMS presente e no zero: o arquivo provavelmente está certo pelo TMS
    if s.nucleus in ("1H", "13C") and not res["ok"]:
        t = _local_max_near(s, 0.0, 0.03 if s.nucleus == "1H" else 0.3)
        if t and t[1] > 20 * sigma and abs(t[0]) < (0.01 if s.nucleus == "1H" else 0.1):
            res["tms"] = True
            res["note"] += "; TMS em 0,00 - mantida a referência pelo TMS"
            res["ok"] = True
            res["delta"] = 0.0
    return res


def reference(s: NMRSpectrum, delta: float) -> NMRSpectrum:
    out = s.copy(ppm=s.ppm - delta)
    out.meta["ref_shift"] = out.meta.get("ref_shift", 0.0) + delta
    return out


def parse_ref(spec: str):
    """'7.31=7.26' -> (7.31, 7.26); 'tms' -> ('tms', 0.0)."""
    t = str(spec).strip().lower()
    if t in ("tms", "dss", "tsp"):
        return t, 0.0
    a, _, b = t.partition("=")
    return float(a), float(b)


def phase_quality(y) -> float:
    """Fração da área dos pontos significativos (|y| > 6σ) que é negativa (≈0 em espectro bem fasado;
    DEPT/APT são exceção)."""
    y = np.asarray(y, float) - np.median(y)
    sig = noise(y)
    m = np.abs(y) > 6 * sig
    if not m.any():
        return 0.0
    yy = y[m]
    return float(-yy[yy < 0].sum() / (np.abs(yy).sum() + 1e-30))


def rephase(s: NMRSpectrum) -> NMRSpectrum:
    """Refasa um espectro já processado: usa a parte imaginária gravada ou reconstrói pela transformada
    de Hilbert e aplica fase automática (ACME)."""
    from scipy.signal import hilbert
    if s.yi is not None:
        z = s.y + 1j * s.yi
    else:
        z = hilbert(s.y)
    z2, a0, a1 = autophase(z)
    out = s.copy(y=z2.real.copy(), yi=z2.imag.copy())
    out.meta["phase_desc"] = f"refasado (ACME: p0 = {a0:.1f}°, p1 = {a1:.1f}°)"
    return out


def prepare(s: NMRSpectrum, do_baseline=True, ref="auto", rephase_auto=True) -> tuple[NMRSpectrum, dict]:
    """Linha de base e referência. ref: 'auto' (pelo solvente, se desviado), 'file' (manter), 'tms'
    (pico mais próximo de 0 vira 0,00) ou '7.31=7.26' (pico observado = valor)."""
    info = {}
    if s.kind != "spectrum":
        return s, info
    pq = phase_quality(s.y)
    info["phase_neg_frac"] = pq
    if rephase_auto and pq > 0.12 and not (s.experiment.startswith("DEPT") or s.experiment == "APT"):
        s2 = rephase(s)
        pq2 = phase_quality(s2.y)
        if pq2 < pq * 0.6:
            s = s2
            info["rephased"] = f"fase corrigida automaticamente (área negativa {100 * pq:.0f}% → {100 * pq2:.0f}%)"
        else:
            info["rephased"] = f"fase ruim (área negativa {100 * pq:.0f}%) e a correção automática não melhorou"
    if do_baseline:
        y, b = baseline(s)
        s = s.copy(y=y)
        info["baseline"] = "Whittaker ponderada"
    sig = noise(s.y)
    info["noise"] = sig
    if ref == "file":
        info["ref"] = check_reference(s, sig)
        info["ref"]["note"] += " (mantida: --ref file)"
        return s, info
    if ref in ("tms", "dss", "tsp"):
        p = find_peak_near(s, 0.0, 0.3)
        if p:
            s = reference(s, p[0])
            info["ref"] = dict(note=f"{ref.upper()} {p[0]:+.3f} → 0,00 ppm", delta=p[0], ok=True)
        else:
            info["ref"] = dict(note="TMS não encontrado perto de 0 ppm", delta=0.0, ok=False)
        return s, info
    if ref not in ("auto", None):
        obs, val = parse_ref(ref)
        p = find_peak_near(s, obs, 0.05 if s.nucleus == "1H" else 0.5)
        d = (p[0] if p else obs) - val
        s = reference(s, d)
        info["ref"] = dict(note=f"pico em {obs:g} → {val:g} ppm (deslocamento {d:+.3f})", delta=d, ok=True)
        return s, info
    chk = check_reference(s, sig)
    big = 0.1 if s.nucleus == "1H" else 1.0
    if not chk["ok"] and chk["expected"] is not None:
        if abs(chk["delta"]) <= big:
            s = reference(s, chk["delta"])
            chk["note"] += " → re-referenciado"
        else:
            chk["note"] += (" → NÃO corrigido automaticamente (desvio grande: confira o solvente ou use "
                            "--ref observado=esperado)")
    info["ref"] = chk
    return s, info
