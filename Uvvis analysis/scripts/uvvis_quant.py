"""
uvvis_quant.py - curva de calibração (Beer-Lambert), predição de amostras e cinética

    from uvvis_quant import calibrate, predict, kinetics
    cal = calibrate(conc, A)                      # dict com slope, intercept, R², LOD, LOQ, ...
    cal = calibrate(conc, A, path_cm=1, conc_unit="mol/L")   # também dá ε
    p = predict(cal, A_amostra, n_rep=3)           # concentração ± IC95%
    k = kinetics(t, A, model="pfo")               # pseudo-1ª ordem: k, t½, R², % degradação

Convenções (ICH Q2(R2)): LOD = 3.3·s(y/x)/b, LOQ = 10·s(y/x)/b, com s(y/x) o desvio-padrão
residual da regressão. Incerteza da predição: s_x0 = s(y/x)/b·√(1/n + 1/m + (ȳ0−ȳ)²/(b²Σ(xi−x̄)²)).
"""
from __future__ import annotations

import numpy as np
from scipy import stats
from scipy.optimize import curve_fit


def calibrate(conc, A, through_origin: bool = False, path_cm: float | None = None,
              conc_unit: str = "", A_max_linear: float = 1.5) -> dict:
    """Regressão linear A = b·c + a (ou A = b·c). Inclui:
      R², R² ajustado, s(y/x), erros-padrão, IC95% de b e a, LOD/LOQ, resíduos,
      teste de Mandel (quadrático vs linear) para não-linearidade,
      teste do intercepto (a é estatisticamente diferente de 0?),
      ε = b / caminho óptico, se conc em mol/L e path_cm dado.
    Pontos com A > A_max_linear geram aviso (desvios de Beer-Lambert/luz espúria)."""
    x = np.asarray(conc, float)
    y = np.asarray(A, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = x.size
    if n < 3:
        raise ValueError("são necessários pelo menos 3 padrões (recomendado ≥ 5-6)")
    warn = []
    if through_origin:
        b = float(np.sum(x * y) / np.sum(x * x))
        a = 0.0
        resid = y - b * x
        dof = n - 1
        syx = float(np.sqrt(np.sum(resid ** 2) / dof))
        se_b = syx / np.sqrt(np.sum(x * x))
        se_a = np.nan
        ss_tot = np.sum(y ** 2)          # R² não centrado (convenção para modelo sem intercepto)
    else:
        b, a = np.polyfit(x, y, 1)
        b, a = float(b), float(a)
        resid = y - (b * x + a)
        dof = n - 2
        syx = float(np.sqrt(np.sum(resid ** 2) / dof))
        sxx = np.sum((x - x.mean()) ** 2)
        se_b = syx / np.sqrt(sxx)
        se_a = syx * np.sqrt(np.sum(x ** 2) / (n * sxx))
        ss_tot = np.sum((y - y.mean()) ** 2)
    r2 = float(1 - np.sum(resid ** 2) / ss_tot)
    p = 1 if through_origin else 2
    r2_adj = float(1 - (1 - r2) * (n - 1) / (n - p)) if n > p else np.nan
    t = stats.t.ppf(0.975, dof)
    out = dict(slope=b, intercept=a, se_slope=float(se_b), se_intercept=float(se_a),
               ci_slope=(b - t * se_b, b + t * se_b),
               ci_intercept=(a - t * se_a, a + t * se_a) if not through_origin else (np.nan, np.nan),
               r2=r2, r2_adj=r2_adj, syx=syx, n=n, dof=dof, through_origin=through_origin,
               lod=3.3 * syx / b, loq=10 * syx / b, conc=x, A=y, residuals=resid,
               x_mean=float(x.mean()), y_mean=float(y.mean()), sxx=float(np.sum((x - x.mean()) ** 2)),
               conc_unit=conc_unit, range=(float(x.min()), float(x.max())))
    if not through_origin:
        tval = a / se_a if se_a > 0 else np.inf
        out["intercept_p"] = float(2 * (1 - stats.t.cdf(abs(tval), dof)))
        if out["intercept_p"] < 0.05:
            warn.append(f"Intercepto ({a:.4g}) significativamente ≠ 0 (p={out['intercept_p']:.3f}): "
                        "branco mal subtraído, interferente ou não-linearidade.")
    # teste de Mandel: o modelo quadrático reduz significativamente o resíduo?
    if n >= 5 and not through_origin:
        c2 = np.polyfit(x, y, 2)
        r_q = y - np.polyval(c2, x)
        ss1, ss2 = np.sum(resid ** 2), np.sum(r_q ** 2)
        F = (ss1 - ss2) / (ss2 / (n - 3)) if ss2 > 0 else np.inf
        pF = float(1 - stats.f.cdf(F, 1, n - 3))
        out["mandel_F"], out["mandel_p"] = float(F), pF
        if pF < 0.05:
            warn.append(f"Teste de Mandel indica curvatura (p={pF:.3f}): reduza a faixa ou use os pontos lineares.")
    if np.any(y > A_max_linear):
        warn.append(f"{int(np.sum(y > A_max_linear))} padrão(ões) com A > {A_max_linear}: confira a linearidade.")
    if r2 < 0.995:
        warn.append(f"R² = {r2:.4f} (< 0.995): verifique pontos discrepantes nos resíduos.")
    # resíduos com tendência (sinais em sequência) sugerem curvatura
    if n >= 6:
        o = np.argsort(x)
        signs = np.sign(resid[o])
        runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
        if runs <= 2:
            warn.append("Resíduos com tendência sistemática (poucos trocas de sinal): possível não-linearidade.")
    if path_cm:
        out["epsilon"] = b / path_cm
        out["epsilon_unit"] = f"(A)/({conc_unit or 'unid. conc.'}·cm)"
        if conc_unit.lower().replace(" ", "") in ("mol/l", "m", "mol·l-1", "moll-1"):
            out["epsilon_unit"] = "L mol⁻¹ cm⁻¹"
    out["warnings"] = warn
    return out


def predict(cal: dict, A0, n_rep: int | None = None, dilution: float = 1.0) -> dict:
    """Concentração de amostra(s) a partir da absorbância. A0 pode ser lista de réplicas
    (média usada, n_rep = len). Retorna conc, s_x0, IC95% e alerta se fora da faixa da curva.
    `dilution` multiplica o resultado (fator de diluição da amostra)."""
    A0 = np.atleast_1d(np.asarray(A0, float))
    m = n_rep or A0.size
    y0 = float(A0.mean())
    b, a = cal["slope"], cal["intercept"]
    c = (y0 - a) / b
    if cal["through_origin"]:
        s = cal["syx"] / abs(b) * np.sqrt(1 / m + y0 ** 2 / (b ** 2 * np.sum(cal["conc"] ** 2)))
    else:
        s = cal["syx"] / abs(b) * np.sqrt(1 / cal["n"] + 1 / m + (y0 - cal["y_mean"]) ** 2 / (b ** 2 * cal["sxx"]))
    t = stats.t.ppf(0.975, cal["dof"])
    lo, hi = cal["range"]
    flags = []
    if c < lo or c > hi:
        flags.append("FORA da faixa da curva - extrapolação; dilua/concentre a amostra")
    if c < cal["loq"]:
        flags.append("abaixo do LOQ" if c >= cal["lod"] else "abaixo do LOD")
    return dict(A=y0, n_rep=m, conc=c * dilution, s_conc=float(s) * dilution,
                ci95=((c - t * s) * dilution, (c + t * s) * dilution), dilution=dilution, flags=flags)


# ------------------------------------------------------------------ cinética
def kinetics(t, A, model: str = "pfo", A0: float | None = None, A_inf: float = 0.0, fit: str = "linear") -> dict:
    """Cinética a partir de A(t) (ou C(t): C/C0 = (A−A∞)/(A0−A∞) se Beer-Lambert vale).
    model: 'pfo' pseudo-1ª ordem   ln(C0/C) = k·t          (Langmuir-Hinshelwood a baixa conc.)
           'pso' pseudo-2ª ordem   1/C − 1/C0 = k'·t       (k' em 1/(unid.A·tempo) se usar A)
           'zero' ordem zero       C0 − C = k0·t
    fit: 'linear' (forma linearizada, como se publica) ou 'nonlinear' (ajuste direto de C(t),
    menos sensível a pontos finais ruidosos). A0 padrão = valor em t=0 (início da reação, após
    o equilíbrio de adsorção no escuro quando houver). Pontos com t<0 (fase escura) são ignorados
    no ajuste mas reportados."""
    t = np.asarray(t, float)
    A = np.asarray(A, float)
    o = np.argsort(t)
    t, A = t[o], A[o]
    dark = t < 0
    if A0 is None:
        if not np.any(t == 0):
            raise ValueError("forneça o ponto t=0 ou A0")
        A0 = float(A[t == 0][0])
    tt, AA = t[~dark], A[~dark]
    C = (AA - A_inf) / (A0 - A_inf)          # C/C0
    C = np.clip(C, 1e-6, None)
    res = dict(model=model, t=tt, C_rel=C, A0=A0, A_inf=A_inf, dark_points=int(dark.sum()))
    if model == "pfo":
        Y = np.log(1 / C)
        if fit == "nonlinear":
            (k,), cov = curve_fit(lambda x, k: np.exp(-k * x), tt, C, p0=[1 / max(tt.max(), 1e-9)])
            pred = np.exp(-k * tt)
            se = float(np.sqrt(cov[0, 0]))
            r2 = 1 - np.sum((C - pred) ** 2) / np.sum((C - C.mean()) ** 2)
            b0 = 0.0
        else:
            k, b0 = np.polyfit(tt, Y, 1)
            pred_Y = k * tt + b0
            r2 = 1 - np.sum((Y - pred_Y) ** 2) / np.sum((Y - Y.mean()) ** 2)
            n = tt.size
            se = float(np.sqrt(np.sum((Y - pred_Y) ** 2) / max(n - 2, 1)) / np.sqrt(np.sum((tt - tt.mean()) ** 2)))
        res.update(k=float(k), se_k=se, intercept=float(b0), r2=float(r2), half_life=float(np.log(2) / k),
                   linear_y=Y, linear_label="ln(C₀/C)")
    elif model == "pso":
        Y = 1 / (C * (A0 - A_inf)) - 1 / (A0 - A_inf)
        k, b0 = np.polyfit(tt, Y, 1)
        pred_Y = k * tt + b0
        r2 = 1 - np.sum((Y - pred_Y) ** 2) / np.sum((Y - Y.mean()) ** 2)
        res.update(k=float(k), intercept=float(b0), r2=float(r2), half_life=float(1 / (k * (A0 - A_inf))),
                   linear_y=Y, linear_label="1/C − 1/C₀")
    elif model == "zero":
        Y = (1 - C) * (A0 - A_inf)
        k, b0 = np.polyfit(tt, Y, 1)
        pred_Y = k * tt + b0
        r2 = 1 - np.sum((Y - pred_Y) ** 2) / np.sum((Y - Y.mean()) ** 2)
        res.update(k=float(k), intercept=float(b0), r2=float(r2), half_life=float((A0 - A_inf) / (2 * k)),
                   linear_y=Y, linear_label="C₀ − C")
    else:
        raise ValueError(model)
    res["degradation_pct"] = float(100 * (1 - C[-1]))
    if dark.any():
        # A[0] = início da fase escura (menor t); A0 = fim do escuro (t = 0)
        res["dark_adsorption_pct"] = float(100 * (1 - (A0 - A_inf) / (A[0] - A_inf))) if A[0] != A_inf else np.nan
    return res


def compare_kinetic_models(t, A, **kw) -> list:
    """Ajusta ordem zero, pseudo-1ª e pseudo-2ª; devolve lista ordenada por R²."""
    out = []
    for m in ("zero", "pfo", "pso"):
        try:
            r = kinetics(t, A, model=m, **kw)
            out.append((m, r["r2"], r["k"]))
        except Exception as e:  # noqa
            out.append((m, np.nan, str(e)))
    return sorted(out, key=lambda z: -np.nan_to_num(z[1], nan=-1))
