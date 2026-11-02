"""
uvvis_bandgap.py - band gap óptico (Tauc), energia de Urbach e verificações cruzadas

    from uvvis_bandgap import tauc, urbach, derivative_edge, onset_wavelength
    r = tauc(s, transition="direct")      # s em A (solução/filme), %T, %R/R (-> Kubelka-Munk) ou α
    r["Eg"], r["Eg_baseline"], r["fit_range"], r["r2"]
    u = urbach(s, Eg=r["Eg"])            # energia de Urbach (meV)

Tauc: (F·hν)^p = B (hν − Eg), com F ∝ α (absorbância de filme/solução, F(R) de pó, ou α real).
    transition   p      (n de Tauc, F·hν ∝ (hν−Eg)^n)
    direct       2      n = 1/2  direta permitida
    indirect     1/2    n = 2    indireta permitida
    direct_forbidden   2/3  n = 3/2
    indirect_forbidden 1/3  n = 3
Eg não depende do fator de escala de F (espessura, concentração): só a forma importa.

Seleção automática da região linear: parte do ponto de MAIOR inclinação da borda e expande para
os dois lados enquanto o ajuste linear mantém R² ≥ r2_min. Também reporta Eg pelo método da
interseção com a linha de base sub-gap (Makuła, Pacia & Macyk, J. Phys. Chem. Lett. 2018, 9,
6814), mais robusto quando há absorção abaixo do gap (defeitos, dopantes, espalhamento).
"""
from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter

from uvvis_io import HC_EV_NM, Spectrum

POWERS = {"direct": 2.0, "indirect": 0.5, "direct_forbidden": 2.0 / 3.0, "indirect_forbidden": 1.0 / 3.0}
LABELS = {"direct": "(F·hν)²", "indirect": "(F·hν)^{1/2}", "direct_forbidden": "(F·hν)^{2/3}",
          "indirect_forbidden": "(F·hν)^{1/3}"}


def absorption_function(s: Spectrum, thickness_cm: float | None = None):
    """Escolhe F ∝ α conforme o tipo de dado. Retorna (E crescente, F, descrição)."""
    k = s.ykind
    if k in ("R", "%R", "KM"):
        f = s.to_km()
        desc = "F(R) Kubelka-Munk"
    elif k in ("T", "%T", "A"):
        f = s.to_absorbance()
        desc = "A"
        if thickness_cm:
            f = f.copy(y=2.302585 * f.y / thickness_cm)
            desc = f"α = 2.303·A/d (d = {thickness_cm:g} cm), cm⁻¹"
    elif k in ("alpha", "eps"):
        f, desc = s, k
    else:
        f, desc = s, f"{k} (assumido ∝ α)"
    E = HC_EV_NM / f.x
    o = np.argsort(E)
    return E[o], f.y[o].astype(float), desc


def _fit(E, Y):
    m, b = np.polyfit(E, Y, 1)
    pred = m * E + b
    ss = np.sum((Y - Y.mean()) ** 2)
    r2 = 1 - np.sum((Y - pred) ** 2) / ss if ss > 0 else 0.0
    return float(m), float(b), float(r2)


def _window_stats(E, Y, L):
    """R², inclinação e intercepto de TODAS as janelas contíguas de L pontos (somas acumuladas)."""
    c = lambda v: np.concatenate([[0.0], np.cumsum(v)])
    Sx, Sy, Sxx, Syy, Sxy = c(E), c(Y), c(E * E), c(Y * Y), c(E * Y)
    i = np.arange(0, len(E) - L + 1)
    j = i + L
    sx, sy = Sx[j] - Sx[i], Sy[j] - Sy[i]
    sxx, syy, sxy = Sxx[j] - Sxx[i], Syy[j] - Syy[i], Sxy[j] - Sxy[i]
    vx = sxx - sx * sx / L
    vy = syy - sy * sy / L
    cxy = sxy - sx * sy / L
    with np.errstate(invalid="ignore", divide="ignore"):
        m = cxy / vx
        r2 = cxy * cxy / (vx * vy)
        b = (sy - m * sx) / L
    return i, m, b, np.nan_to_num(r2)


def _best_linear(E, Y, r2_min, min_pts, slope_min, y_min, prefer="long"):
    """Maior segmento linear (R² ≥ r2_min) com inclinação ≥ slope_min e sinal médio ≥ y_min.
    Entre segmentos do mesmo tamanho, o de maior R². Se nenhum atinge r2_min, relaxa até 0.95."""
    n = len(E)
    for thr in (r2_min, 0.99, 0.98, 0.95):
        for L in range(n, min_pts - 1, -max(1, n // 200)):
            i, m, b, r2 = _window_stats(E, Y, L)
            meanY = (np.concatenate([[0.0], np.cumsum(Y)])[i + L] - np.concatenate([[0.0], np.cumsum(Y)])[i]) / L
            ok = (r2 >= thr) & (m >= slope_min) & (meanY >= y_min) & (m > 0)
            if ok.any():
                if prefer == "steep":
                    # todas as janelas válidas de tamanho >= min_pts: a mais íngreme
                    best = None
                    for L2 in range(min_pts, L + 1, max(1, (L - min_pts) // 40 or 1)):
                        i2, m2, b2, r22 = _window_stats(E, Y, L2)
                        ok2 = (r22 >= thr) & (m2 > 0)
                        if ok2.any():
                            k2 = np.where(ok2)[0][np.argmax(m2[ok2])]
                            if best is None or m2[k2] > best[0]:
                                best = (m2[k2], int(i2[k2]), int(i2[k2] + L2 - 1))
                    return best[1], best[2], thr
                k = np.where(ok)[0][np.argmax(r2[ok])]
                return int(i[k]), int(i[k] + L - 1), thr
    return None


def tauc(s: Spectrum, transition: str = "direct", thickness_cm: float | None = None,
         energy_window=None, fit_window=None, r2_min: float = 0.995, smooth_pts: int = 0) -> dict:
    """Band gap pelo gráfico de Tauc.
    energy_window: (E1, E2) eV para procurar a borda (exclui excitons/saturação/outra fase).
    fit_window: (E1, E2) eV para FORÇAR a região de ajuste (reprodução de um ajuste manual).
    Retorna Eg (interseção com y=0), Eg_baseline (interseção com a reta da região sub-gap),
    fit_range, slope, intercept, r2, Eg_spread (variação ao encurtar/alongar a janela ±25%),
    arrays E, Y para o gráfico, e avisos."""
    p = POWERS[transition]
    E, F, desc = absorption_function(s, thickness_cm)
    warn = []
    if np.mean(F < 0) > 0.05:
        warn.append("Mais de 5% dos pontos com F < 0 (linha de base negativa): corrija o offset antes.")
    Fp = np.clip(F, 0, None)
    Y = (Fp * E) ** p
    m_all = np.isfinite(Y)
    if energy_window is not None:
        m_all &= (E >= min(energy_window)) & (E <= max(energy_window))
    Ew, Yw = E[m_all], Y[m_all]
    if Ew.size < 15:
        raise ValueError("poucos pontos na janela de energia")
    dE = float(np.median(np.diff(Ew)))
    w = smooth_pts or max(5, int(round(0.04 / max(dE, 1e-6))) | 1)
    w = min(w, (Ew.size - 2) | 1)
    Ys = savgol_filter(Yw, w, 2)
    dY = np.gradient(Ys, Ew)
    min_pts = max(6, int(round(0.05 / dE)))
    used_thr = r2_min
    i_on = 0
    other_edges = []
    if fit_window is not None:
        mm = (Ew >= min(fit_window)) & (Ew <= max(fit_window))
        lo, hi = np.where(mm)[0][[0, -1]]
        auto = False
    else:
        # bordas candidatas = picos de dY/dE; cada uma vai do onset (dY < 10% do pico, à esquerda)
        # até onde dY volta a cair (à direita). A PRINCIPAL é a de maior subida em Y; as outras
        # (dopantes, defeitos, sensibilizadores, segunda fase) são reportadas.
        from scipy.signal import find_peaks as _fp
        w2 = min(max(7, int(round(0.10 / max(dE, 1e-6))) | 1), (Ew.size - 2) | 1)
        dYp = np.clip(np.gradient(savgol_filter(Yw, w2, 2), Ew), 0, None)   # derivada bem suavizada
        pk, _ = _fp(dYp, prominence=0.08 * float(np.nanmax(dYp)), distance=max(3, min_pts))
        if pk.size == 0:
            pk = np.array([int(np.argmax(dYp))])
        pk = np.sort(pk)
        # agrupa picos separados por vale raso (mesma borda com ruído/estrutura)
        groups = [[pk[0]]]
        for q in pk[1:]:
            pv = groups[-1][-1]
            valley = float(np.min(dYp[pv:q + 1]))
            if valley > 0.3 * min(dYp[pv], dYp[q]):
                groups[-1].append(q)
            else:
                groups.append([q])
        merged = []
        for gi, gq in enumerate(groups):
            qmax = max(gq, key=lambda z: dYp[z])
            if gi == 0:
                l = gq[0]
                while l > 0 and dYp[l] > 0.10 * dYp[qmax]:
                    l -= 1
            else:
                l = groups[gi - 1][-1] + int(np.argmin(dYp[groups[gi - 1][-1]:gq[0] + 1]))
            if gi == len(groups) - 1:
                r = gq[-1]
                while r < dYp.size - 1 and dYp[r] > 0.10 * dYp[qmax]:
                    r += 1
            else:
                r = gq[-1] + int(np.argmin(dYp[gq[-1]:groups[gi + 1][0] + 1]))
            merged.append((float(Ys[r] - Ys[l]), int(l), int(qmax), int(r)))
        main = max(merged, key=lambda e: e[0])
        # só interessam bordas ABAIXO da principal (absorção sub-gap); acima são transições normais/ruído
        other_edges = [float(Ew[e[2]]) for e in merged if e[2] < main[2] and e[0] > 0.03 * main[0]]
        i_on = main[1]
        base = float(Ys[i_on])
        loc = (Ew >= Ew[i_on]) & (Ew <= Ew[i_on] + 1.0)
        Ytop = float(np.nanmax(Ys[loc]))
        top = np.where((np.arange(Ys.size) > i_on) & (Ys >= base + 0.85 * (Ytop - base)))[0]
        i_top = int(top[0]) if top.size else Ys.size - 1
        a0, a1 = max(0, i_on - min_pts // 3), min(Ys.size - 1, i_top + min_pts // 3)
        Es, Yseg = Ew[a0:a1 + 1], Yw[a0:a1 + 1]
        found = None
        if Es.size >= min_pts:
            found = _best_linear(Es, Yseg, r2_min, max(min_pts, int(0.3 * Es.size)), 0.0, -np.inf,
                                 prefer="steep")
        if found is None:
            raise ValueError("borda de absorção não encontrada (sem trecho linear ascendente); "
                             "restrinja energy_window ou use fit_window")
        found = (found[0] + a0, found[1] + a0, found[2])
        lo, hi, used_thr = found
        auto = True
    m, b, r2 = _fit(Ew[lo:hi + 1], Yw[lo:hi + 1])
    if m <= 0:
        raise ValueError("inclinação não positiva na região de ajuste: borda não encontrada nessa janela")
    Eg = -b / m
    # sensibilidade: encurta/alonga 25% de cada lado
    n = hi - lo + 1
    alt = []
    for dl, dh in ((n // 4, 0), (0, -n // 4), (-n // 4, 0), (0, n // 4), (n // 4, -n // 4)):
        l2, h2 = max(0, lo + dl), min(Ew.size - 1, hi + dh)
        if h2 - l2 >= 4:
            mm_, bb_, _ = _fit(Ew[l2:h2 + 1], Yw[l2:h2 + 1])
            if mm_ > 0:
                alt.append(-bb_ / mm_)
    spread = float(np.ptp(alt)) if alt else np.nan
    # linha de base sub-gap: região logo abaixo da borda (entre Eg-0.5 e o início do ajuste)
    Eg_base, base_line = np.nan, None
    e_on = Ew[i_on] if auto else Ew[lo]
    e_hi = min(e_on, Ew[lo]) - 0.05
    sub = (Ew <= e_hi) & (Ew > e_hi - 0.5)
    if sub.sum() >= 6:
        mb, bb, r2b = _fit(Ew[sub], Yw[sub])
        if abs(m - mb) > 1e-12:
            Eg_base = (bb - b) / (m - mb)
            base_line = (mb, bb, (float(Ew[sub].min()), float(Ew[sub].max())))
            if Eg_base < Eg - 1e-9:
                warn.append("Reta sub-gap com inclinação incomum (interseção abaixo de Eg): confira a região.")
    # checagens
    if Ew[hi] >= Ew.max() - 2 * dE:
        warn.append("A região linear encosta no limite de alta energia dos dados: borda pode estar truncada "
                    "(solvente/cubeta/vidro cortando, ou saturação).")
    if s.ykind == "A":
        lam_fit = HC_EV_NM / Ew[lo:hi + 1]
        Amax = max(s.at(l) for l in (lam_fit.min(), lam_fit.max()))
        if Amax > 2.5:
            warn.append(f"Absorbância até {Amax:.2f} na região de ajuste: saturação distorce a borda.")
    if used_thr < r2_min and auto:
        warn.append(f"Nenhum trecho com R² ≥ {r2_min}; critério relaxado para {used_thr} (R² = {r2:.4f}).")
    if n < min_pts:
        warn.append("Região linear curta: Eg pouco confiável - inspecione o gráfico.")
    if other_edges:
        warn.append("Outra(s) borda(s)/banda(s) de absorção em " + ", ".join(f"{e:.2f}" for e in other_edges) +
                    " eV (inflexão): dopante/defeito/sensibilizador/segunda fase. O ajuste usou a borda principal; "
                    "para a outra, use energy_window.")
    if spread > 0.05:
        warn.append(f"Eg varia {spread:.2f} eV com a escolha da janela: reporte a incerteza ou fixe a janela.")
    return dict(Eg=float(Eg), Eg_baseline=float(Eg_base), transition=transition, power=p, source=desc,
                fit_range=(float(Ew[lo]), float(Ew[hi])), slope=m, intercept=b, r2=r2, n_fit=int(n),
                Eg_spread=spread, auto=auto, E=E, Y=Y, base_line=base_line, label=LABELS[transition],
                lambda_edge_nm=float(HC_EV_NM / Eg) if Eg > 0 else np.nan, other_edges_eV=other_edges,
                warnings=warn)


def urbach(s: Spectrum, Eg: float | None = None, thickness_cm: float | None = None,
           window=None, r2_min: float = 0.99) -> dict:
    """Energia de Urbach: na cauda exponencial abaixo do gap, α ∝ exp(hν/E_U) ⇒ ln F vs hν é
    linear com inclinação 1/E_U. Janela automática: entre Eg−0.5 eV e Eg, só pontos com F bem
    acima do ruído (F > 8σ), maior trecho linear com R² ≥ r2_min. Resultado em meV.
    Só tem significado se a cauda for medida acima da linha de base/espalhamento."""
    E, F, desc = absorption_function(s, thickness_cm)
    d2 = np.diff(F, 2)
    sig = float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6)) + 1e-12
    valid = F > 8 * sig
    E_all, L_all = E[F > 0], np.log(F[F > 0])
    if window is not None:
        sel = valid & (E >= min(window)) & (E <= max(window))
    elif Eg is not None:
        sel = valid & (E >= Eg - 0.5) & (E <= Eg)
    else:
        raise ValueError("informe Eg (de tauc) ou window=(E1, E2)")
    Ew, Lw = E[sel], np.log(F[sel])
    if Ew.size < 8:
        raise ValueError("poucos pontos acima do ruído na cauda de Urbach")
    if window is not None:
        lo, hi = 0, Ew.size - 1
    else:
        dE = float(np.median(np.diff(Ew)))
        w = min(max(5, int(round(0.03 / dE)) | 1), (Ew.size - 2) | 1)
        dd = np.gradient(savgol_filter(Lw, w, 2), Ew)
        found = _best_linear(Ew, Lw, r2_min, max(6, int(round(0.04 / dE))), 0.3 * float(np.nanmax(dd)), -np.inf)
        if found is None:
            raise ValueError("cauda exponencial não encontrada; informe window=(E1, E2)")
        lo, hi, _ = found
    sl, ic, r2 = _fit(Ew[lo:hi + 1], Lw[lo:hi + 1])
    return dict(Eu_meV=float(1000.0 / sl) if sl > 0 else np.nan, slope=sl, intercept=ic, r2=r2,
                fit_range=(float(Ew[lo]), float(Ew[hi])), E=E_all, lnF=L_all, source=desc)


def derivative_edge(s: Spectrum, thickness_cm: float | None = None, energy_window=None,
                    Eg: float | None = None) -> dict:
    """Verificação rápida: energia de máxima inclinação dF/dE da borda (≈ ponto de inflexão).
    Costuma ficar um pouco ACIMA do Eg de Tauc; serve para conferir ordem de grandeza.
    Com Eg (de tauc), procura só em Eg ± 0.4 eV (gaps indiretos sobem até o fim do espectro)."""
    E, F, _ = absorption_function(s, thickness_cm)
    if energy_window is None and Eg is not None:
        energy_window = (Eg - 0.4, Eg + 0.4)
    m = np.ones_like(E, bool) if energy_window is None else (E >= min(energy_window)) & (E <= max(energy_window))
    E, F = E[m], F[m]
    dE = float(np.median(np.diff(E)))
    w = min(max(5, int(round(0.05 / dE)) | 1), (E.size - 2) | 1)
    d = savgol_filter(F, w, 2, deriv=1, delta=dE)
    i = int(np.argmax(d))
    return dict(E_inflection=float(E[i]), lambda_nm=float(HC_EV_NM / E[i]))


def onset_wavelength(s: Spectrum, window_nm=None, Eg: float | None = None) -> dict:
    """Método do 'onset'/corte: tangente no ponto de inflexão da borda (A vs λ) intercepta a
    linha de base de comprimentos de onda longos (reta ajustada no platô após a borda).
    Eg ≈ 1240/λ_onset. Muito usado em química de materiais; tende a dar valores menores que Tauc."""
    a = s.to_absorbance() if s.ykind in ("T", "%T") else (s.to_km() if s.ykind in ("R", "%R") else s)
    if window_nm is None and Eg is not None:
        window_nm = (HC_EV_NM / (Eg + 0.6), HC_EV_NM / (Eg - 0.6))
    o = np.argsort(a.x)
    x, y = a.x[o], a.y[o]
    if window_nm is not None:
        m = (x >= min(window_nm)) & (x <= max(window_nm))
        x, y = x[m], y[m]
    st = float(np.median(np.diff(x)))
    w = min(max(5, int(round(8 / st)) | 1), (x.size - 2) | 1)
    d = savgol_filter(y, w, 2, deriv=1, delta=st)
    i = int(np.argmin(d))                         # descida mais íngreme (A cai com λ)
    k = max(3, int(round(5 / st)))
    sl, ic = np.polyfit(x[max(0, i - k):i + k + 1], y[max(0, i - k):i + k + 1], 1)
    # linha de base: trecho após a borda onde |d| é pequeno
    tail = np.where((np.arange(x.size) > i) & (np.abs(d) < 0.05 * abs(d[i])))[0]
    if tail.size >= 5:
        tail = tail[: max(5, tail.size // 2)]
        bs, bi = np.polyfit(x[tail], y[tail], 1)
    else:
        bs, bi = 0.0, float(np.min(y))
    lam = (bi - ic) / (sl - bs)
    return dict(lambda_onset_nm=float(lam), Eg_onset=float(HC_EV_NM / lam), inflection_nm=float(x[i]))
