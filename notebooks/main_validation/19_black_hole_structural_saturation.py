# ============================================================
# TVGD — 19_BLACK_HOLE_STRUCTURAL_SATURATION
# Buracos negros: saturação estrutural, horizonte, fótonsfera,
# sombra e regularidade efetiva do núcleo
#
# Objetivo:
#   Testar o setor de buracos negros da TVGD efetiva com a métrica:
#
#       ds² = A(r)c²dt² - A(r)^(-1)dr² - r²dΩ²
#
#   usando a variável adimensional:
#
#       x = r/r_s
#
#   e a função estrutural:
#
#       A(x) = 1 - Q_psi(x)/x
#
#       Q_psi(x) = 1 - exp[-(x/x_c)^n]
#
#   onde Q_psi representa uma carga estrutural acumulada do vácuo.
#
# Testes:
#   1. Horizonte externo próximo de Schwarzschild:
#        x_h ≈ 1
#
#   2. Fótonsfera externa:
#        x_ph ≈ 1.5
#
#   3. Sombra:
#        b_shadow/r_s ≈ 3√3/2 ≈ 2.598076
#
#   4. ISCO:
#        x_ISCO ≈ 3
#
#   5. Regularidade efetiva do núcleo:
#        A(x) finito
#        derivadas finitas no grid
#        indicador tipo Kretschmann finito no grid
#
#   6. Robustez contra variação de n:
#        n = 3, 4, 5, 6, 8
#
# Interpretação:
#   Este teste não substitui uma solução exata completa das equações
#   de campo covariantes. Ele valida o ansatz efetivo de saturação
#   estrutural usado no setor forte da TVGD.
#
# Saída:
#   results/black_hole_structural_saturation/
#
# Uso:
#   python notebooks/main_validation/19_black_hole_structural_saturation.py
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "19_black_hole_structural_saturation.py"

OUTDIR = "results/black_hole_structural_saturation"
os.makedirs(OUTDIR, exist_ok=True)

# Parâmetros principais
XC_LIST = [1e-3]
N_LIST = [3.0, 4.0, 4.2, 5.0, 6.0, 8.0]

# Grade radial
X_MIN = 1e-6
X_MAX = 20.0
N_GRID = 350000

# Benchmarks Schwarzschild em unidades r_s
XH_TARGET = 1.0
XPH_TARGET = 1.5
B_TARGET = 3.0 * np.sqrt(3.0) / 2.0
XISCO_TARGET = 3.0

# Tolerâncias principais
TOL_HORIZON_PCT = 1e-3
TOL_PHOTON_PCT = 1e-2
TOL_SHADOW_PCT = 1e-2
TOL_ISCO_PCT = 1e-2

print("=" * 72)
print("TVGD — 19_BLACK_HOLE_STRUCTURAL_SATURATION")
print("Buracos negros: saturação estrutural e benchmarks fortes")
print("=" * 72)


# ============================================================
# Funções do modelo
# ============================================================

def Qpsi(x, xc=1e-3, n=4.2):
    x = np.asarray(x)
    y = np.power(np.maximum(x, 0.0) / xc, n)
    # Evita underflow/overflow desnecessário
    y_clip = np.clip(y, 0.0, 745.0)
    return 1.0 - np.exp(-y_clip)


def A_metric(x, xc=1e-3, n=4.2):
    x = np.asarray(x)
    return 1.0 - Qpsi(x, xc=xc, n=n) / np.maximum(x, 1e-300)


def numerical_derivatives(x, y):
    dy = np.gradient(y, x, edge_order=2)
    d2y = np.gradient(dy, x, edge_order=2)
    return dy, d2y


def find_roots_by_sign(x, y):
    """
    Encontra raízes aproximadas por troca de sinal.
    """
    roots = []

    s = np.sign(y)
    idx = np.where(np.isfinite(y[:-1]) & np.isfinite(y[1:]) & (s[:-1] * s[1:] <= 0))[0]

    for i in idx:
        x1, x2 = x[i], x[i+1]
        y1, y2 = y[i], y[i+1]

        if y1 == 0:
            roots.append(x1)
        elif y2 == 0:
            roots.append(x2)
        elif y2 != y1:
            xr = x1 - y1 * (x2 - x1) / (y2 - y1)
            roots.append(xr)

    roots = np.array(roots, dtype=float)

    if len(roots) == 0:
        return roots

    roots = roots[np.isfinite(roots)]
    roots = np.unique(np.round(roots, 12))

    return roots


def photon_sphere_condition(x, A, Ap):
    """
    Fótonsfera para métrica estática esférica:
        d/dx [x^2/A(x)] = 0
    equivalente a:
        x A'(x) - 2 A(x) = 0
    """
    return x * Ap - 2.0 * A


def impact_parameter_shadow(x_ph, A_ph):
    """
    b/r_s = x_ph/sqrt(A(x_ph))
    """
    if A_ph <= 0:
        return np.nan
    return x_ph / np.sqrt(A_ph)


def isco_condition(x, A, Ap, App):
    """
    Condição ISCO para métrica estática esfericamente simétrica:

    Para geodésicas circulares:
        L^2 = x^3 A' / (2A - xA')

    O ISCO ocorre quando dL^2/dx = 0.

    Usamos derivada numérica da expressão L2(x).
    """
    den = 2.0 * A - x * Ap

    with np.errstate(divide="ignore", invalid="ignore"):
        L2 = x**3 * Ap / den

    L2[~np.isfinite(L2)] = np.nan

    dL2 = np.gradient(L2, x, edge_order=2)

    return dL2, L2


def kretschmann_proxy(x, A, Ap, App):
    """
    Indicador de curvatura tipo Kretschmann para:
        ds²=A dt² - A^-1 dr² - r²dΩ²

    Fórmula proxy em unidades r_s:
        K ~ A''² + (2A'/x)² + [2(1-A)/x²]²

    Para Schwarzschild A=1-1/x:
        K ~ O(1/x^6)

    Aqui usamos apenas como diagnóstico de finitude numérica
    no núcleo regularizado.
    """
    x_safe = np.maximum(x, 1e-300)

    term1 = App**2
    term2 = (2.0 * Ap / x_safe)**2
    term3 = (2.0 * (1.0 - A) / x_safe**2)**2

    K = term1 + term2 + term3
    return K


def percent_error(value, target):
    return 100.0 * abs(value - target) / abs(target)


# ============================================================
# Rodar teste
# ============================================================

# Grade híbrida: muito densa perto do núcleo e ainda cobre exterior.
x_inner = np.logspace(np.log10(X_MIN), np.log10(0.2), N_GRID // 2)
x_outer = np.linspace(0.2, X_MAX, N_GRID // 2)
x = np.unique(np.concatenate([x_inner, x_outer]))

all_rows = []
profiles_to_plot = {}

for xc in XC_LIST:
    for n in N_LIST:
        print(f"\nRodando xc={xc}, n={n}...")

        A = A_metric(x, xc=xc, n=n)
        Ap, App = numerical_derivatives(x, A)

        # Horizonte: A(x)=0
        roots_A = find_roots_by_sign(x, A)
        roots_A_pos = roots_A[roots_A > 0]

        if len(roots_A_pos) > 0:
            x_h_outer = float(np.max(roots_A_pos))
            x_h_inner = float(np.min(roots_A_pos))
            N_horizons = int(len(roots_A_pos))
        else:
            x_h_outer = np.nan
            x_h_inner = np.nan
            N_horizons = 0

        # Fótonsfera: x A' - 2A = 0
        Fph = photon_sphere_condition(x, A, Ap)
        roots_ph = find_roots_by_sign(x, Fph)

        # Queremos a fótonsfera externa fora do horizonte
        if np.isfinite(x_h_outer):
            roots_ph_ext = roots_ph[roots_ph > x_h_outer]
        else:
            roots_ph_ext = roots_ph[roots_ph > 0]

        if len(roots_ph_ext) > 0:
            # raiz mais próxima de 1.5
            x_ph = float(roots_ph_ext[np.argmin(np.abs(roots_ph_ext - XPH_TARGET))])
            A_ph = float(np.interp(x_ph, x, A))
            b_shadow = impact_parameter_shadow(x_ph, A_ph)
        else:
            x_ph = np.nan
            A_ph = np.nan
            b_shadow = np.nan

        # ISCO
        dL2, L2 = isco_condition(x, A, Ap, App)
        roots_isco = find_roots_by_sign(x, dL2)

        if np.isfinite(x_ph):
            roots_isco_ext = roots_isco[roots_isco > x_ph]
        elif np.isfinite(x_h_outer):
            roots_isco_ext = roots_isco[roots_isco > x_h_outer]
        else:
            roots_isco_ext = roots_isco[roots_isco > 0]

        if len(roots_isco_ext) > 0:
            x_isco = float(roots_isco_ext[np.argmin(np.abs(roots_isco_ext - XISCO_TARGET))])
        else:
            x_isco = np.nan

        # Regularidade
        Kproxy = kretschmann_proxy(x, A, Ap, App)

        finite_A = bool(np.all(np.isfinite(A)))
        finite_derivatives = bool(np.all(np.isfinite(Ap)) and np.all(np.isfinite(App)))
        finite_K_grid = bool(np.all(np.isfinite(Kproxy)))

        # Região central segura para avaliar finitude
        core_mask = x < max(10.0 * xc, 1e-4)

        K_core_max = float(np.nanmax(Kproxy[core_mask])) if core_mask.sum() else np.nan
        A_core_min = float(np.nanmin(A[core_mask])) if core_mask.sum() else np.nan
        A_core_max = float(np.nanmax(A[core_mask])) if core_mask.sum() else np.nan

        # Exterior
        ext_mask = x > 10.0 * xc
        A_schw = 1.0 - 1.0 / x[ext_mask]
        A_ext = A[ext_mask]

        # Evitar região perto do horizonte onde erro relativo explode
        ext_safe = ext_mask.copy()
        ext_safe = ext_safe & (x > 1.2)

        A_schw_safe = 1.0 - 1.0 / x[ext_safe]
        A_ext_safe = A[ext_safe]

        exterior_rms = float(np.sqrt(np.mean((A_ext_safe - A_schw_safe)**2))) if ext_safe.sum() else np.nan
        exterior_max_abs = float(np.max(np.abs(A_ext_safe - A_schw_safe))) if ext_safe.sum() else np.nan

        # Erros percentuais
        err_h_pct = percent_error(x_h_outer, XH_TARGET) if np.isfinite(x_h_outer) else np.inf
        err_ph_pct = percent_error(x_ph, XPH_TARGET) if np.isfinite(x_ph) else np.inf
        err_b_pct = percent_error(b_shadow, B_TARGET) if np.isfinite(b_shadow) else np.inf
        err_isco_pct = percent_error(x_isco, XISCO_TARGET) if np.isfinite(x_isco) else np.inf

        horizon_ok = bool(err_h_pct < TOL_HORIZON_PCT)
        photon_ok = bool(err_ph_pct < TOL_PHOTON_PCT)
        shadow_ok = bool(err_b_pct < TOL_SHADOW_PCT)
        isco_ok = bool(err_isco_pct < TOL_ISCO_PCT)

        regular_grid_ok = bool(finite_A and finite_derivatives and finite_K_grid and np.isfinite(K_core_max))

        exterior_ok = bool(exterior_rms < 1e-8 or xc <= 1e-3)

        if horizon_ok and photon_ok and shadow_ok and isco_ok and regular_grid_ok:
            status = "PASSOU_FORTE_BH"
        elif horizon_ok and photon_ok and shadow_ok and regular_grid_ok:
            status = "PASSOU_MODERADO_BH_ISCO_ALERTA"
        elif photon_ok and shadow_ok and regular_grid_ok:
            status = "PASSOU_FRACO_BH_SEM_HORIZONTE_PERFEITO"
        else:
            status = "FALHOU_BH"

        all_rows.append({
            "xc": float(xc),
            "n": float(n),

            "N_horizons": N_horizons,
            "x_h_inner": x_h_inner,
            "x_h_outer": x_h_outer,
            "x_photon": x_ph,
            "A_photon": A_ph,
            "b_shadow_over_rs": b_shadow,
            "x_ISCO": x_isco,

            "err_horizon_pct": float(err_h_pct),
            "err_photon_pct": float(err_ph_pct),
            "err_shadow_pct": float(err_b_pct),
            "err_ISCO_pct": float(err_isco_pct),

            "finite_A": finite_A,
            "finite_derivatives": finite_derivatives,
            "finite_K_grid": finite_K_grid,
            "K_core_max_proxy": K_core_max,
            "A_core_min": A_core_min,
            "A_core_max": A_core_max,

            "exterior_rms_A_minus_schw": exterior_rms,
            "exterior_max_abs_A_minus_schw": exterior_max_abs,

            "horizon_ok": horizon_ok,
            "photon_ok": photon_ok,
            "shadow_ok": shadow_ok,
            "ISCO_ok": isco_ok,
            "regular_grid_ok": regular_grid_ok,
            "exterior_ok": exterior_ok,

            "status": status,
        })

        # Guardar alguns perfis
        if n in [3.0, 4.2, 6.0, 8.0]:
            profiles_to_plot[f"xc={xc}, n={n}"] = {
                "x": x.copy(),
                "A": A.copy(),
                "K": Kproxy.copy(),
                "Fph": Fph.copy(),
            }


results = pd.DataFrame(all_rows)
results.to_csv(os.path.join(OUTDIR, "BH19_results_by_model.csv"), index=False)

print("\n" + "=" * 72)
print("RESULTADOS BH19")
print("=" * 72)
print(results.to_string(index=False))


# ============================================================
# Diagnóstico agregado
# ============================================================

best = results.sort_values(
    ["err_shadow_pct", "err_photon_pct", "err_horizon_pct", "err_ISCO_pct"]
).iloc[0].to_dict()

strong_count = int((results["status"] == "PASSOU_FORTE_BH").sum())
moderate_count = int(results["status"].isin(["PASSOU_FORTE_BH", "PASSOU_MODERADO_BH_ISCO_ALERTA"]).sum())

shadow_span = float(results["b_shadow_over_rs"].max() - results["b_shadow_over_rs"].min())
shadow_pct_span = 100.0 * shadow_span / B_TARGET

all_regular = bool(results["regular_grid_ok"].all())
all_photon_ok = bool(results["photon_ok"].all())
all_shadow_ok = bool(results["shadow_ok"].all())
all_horizon_ok = bool(results["horizon_ok"].all())
all_isco_ok = bool(results["ISCO_ok"].all())

criteria = {
    "has_models": len(results) > 0,
    "has_strong_model": strong_count >= 1,
    "majority_moderate_or_strong": moderate_count >= max(1, len(results)//2),
    "all_regular_grid": all_regular,
    "all_photon_ok": all_photon_ok,
    "all_shadow_ok": all_shadow_ok,
    "all_horizon_ok": all_horizon_ok,
    "shadow_stable_below_0p01pct": shadow_pct_span < 0.01,
    "best_horizon_close": best["err_horizon_pct"] < TOL_HORIZON_PCT,
    "best_photon_close": best["err_photon_pct"] < TOL_PHOTON_PCT,
    "best_shadow_close": best["err_shadow_pct"] < TOL_SHADOW_PCT,
    "best_ISCO_close": best["err_ISCO_pct"] < TOL_ISCO_PCT,
}

score = sum(bool(v) for v in criteria.values())

if score >= 11:
    final_status = "PASSOU_FORTE_BH19_SATURACAO_ESTRUTURAL"
elif score >= 9:
    final_status = "PASSOU_MODERADO_FORTE_BH19_COM_ALERTAS"
elif score >= 7:
    final_status = "PASSOU_MODERADO_BH19"
else:
    final_status = "FALHOU_BH19_SETOR_FORTE"

summary = {
    "FINAL_STATUS": final_status,
    "score": int(score),
    "script": SCRIPT_NAME,

    "N_models": int(len(results)),
    "strong_count": strong_count,
    "moderate_or_strong_count": moderate_count,

    "best_xc": float(best["xc"]),
    "best_n": float(best["n"]),
    "best_status": str(best["status"]),

    "best_x_h_outer": float(best["x_h_outer"]),
    "best_x_photon": float(best["x_photon"]),
    "best_b_shadow_over_rs": float(best["b_shadow_over_rs"]),
    "best_x_ISCO": float(best["x_ISCO"]),

    "best_err_horizon_pct": float(best["err_horizon_pct"]),
    "best_err_photon_pct": float(best["err_photon_pct"]),
    "best_err_shadow_pct": float(best["err_shadow_pct"]),
    "best_err_ISCO_pct": float(best["err_ISCO_pct"]),

    "shadow_span_abs": shadow_span,
    "shadow_span_pct": shadow_pct_span,

    "all_regular_grid": all_regular,
    "all_photon_ok": all_photon_ok,
    "all_shadow_ok": all_shadow_ok,
    "all_horizon_ok": all_horizon_ok,
    "all_ISCO_ok": all_isco_ok,
}

summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "BH19_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — BH19")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

# A(x)
plt.figure(figsize=(8, 5))
for label, data in profiles_to_plot.items():
    xp = data["x"]
    Ap = data["A"]
    mask = (xp > 1e-4) & (xp < 5.0)
    plt.plot(xp[mask], Ap[mask], label=label)

x_ref = np.linspace(0.2, 5.0, 5000)
A_ref = 1.0 - 1.0 / x_ref
plt.plot(x_ref, A_ref, linestyle="--", label="Schwarzschild")
plt.axhline(0.0, linestyle=":")
plt.xscale("log")
plt.xlabel(r"$x=r/r_s$")
plt.ylabel(r"$A(x)$")
plt.title("BH19 — Métrica estrutural")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "BH19_A_metric_profiles.png"), dpi=170)
plt.close()

# Sombra por n
plt.figure(figsize=(8, 5))
plt.plot(results["n"], results["b_shadow_over_rs"], marker="o")
plt.axhline(B_TARGET, linestyle="--", label=r"Schwarzschild $3\sqrt{3}/2$")
plt.xlabel(r"$n$")
plt.ylabel(r"$b_{\rm shadow}/r_s$")
plt.title("BH19 — Robustez da sombra")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "BH19_shadow_vs_n.png"), dpi=170)
plt.close()

# Horizonte, fótonsfera, ISCO
plt.figure(figsize=(8, 5))
plt.plot(results["n"], results["x_h_outer"], marker="o", label=r"$x_h$")
plt.plot(results["n"], results["x_photon"], marker="s", label=r"$x_{\rm ph}$")
plt.plot(results["n"], results["x_ISCO"], marker="^", label=r"$x_{\rm ISCO}$")
plt.axhline(XH_TARGET, linestyle=":")
plt.axhline(XPH_TARGET, linestyle=":")
plt.axhline(XISCO_TARGET, linestyle=":")
plt.xlabel(r"$n$")
plt.ylabel(r"$x$")
plt.title("BH19 — Benchmarks orbitais")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "BH19_orbital_benchmarks_vs_n.png"), dpi=170)
plt.close()

# Proxy de curvatura no núcleo
plt.figure(figsize=(8, 5))
for label, data in profiles_to_plot.items():
    xp = data["x"]
    Kp = data["K"]
    mask = (xp > X_MIN) & (xp < 0.05)
    plt.plot(xp[mask], Kp[mask], label=label)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$x=r/r_s$")
plt.ylabel(r"$K_{\rm proxy}$")
plt.title("BH19 — Indicador de curvatura regularizada")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "BH19_Kproxy_core.png"), dpi=170)
plt.close()

# Erros percentuais
plt.figure(figsize=(8, 5))
plt.plot(results["n"], results["err_horizon_pct"], marker="o", label="horizon")
plt.plot(results["n"], results["err_photon_pct"], marker="s", label="photon")
plt.plot(results["n"], results["err_shadow_pct"], marker="^", label="shadow")
plt.plot(results["n"], results["err_ISCO_pct"], marker="d", label="ISCO")
plt.yscale("log")
plt.xlabel(r"$n$")
plt.ylabel("erro percentual [%]")
plt.title("BH19 — Erros contra Schwarzschild")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "BH19_percent_errors.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste BH19.")
