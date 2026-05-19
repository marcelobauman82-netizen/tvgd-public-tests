# ============================================================
# TVGD PUBLIC TESTS
# TEST 08 — SINGULARITY AVOIDANCE / CURVATURE REGULARIZATION
# ============================================================
#
# Goal:
# Compare the Schwarzschild curvature singularity with a TVGD
# saturation-regularized black-hole core.
#
# The test checks:
#
# 1. Schwarzschild Kretschmann scalar diverges as K ~ 12/x^6.
# 2. TVGD regularized metric keeps K finite near x -> 0.
# 3. The exterior horizon, photon sphere and shadow are preserved.
# 4. The regularization affects the core but not the external regime.
#
# Metric form:
#
#     ds^2 = A(x) dt^2 - A(x)^(-1) dx^2 - x^2 dOmega^2
#
# Schwarzschild:
#
#     A_Schw(x) = 1 - 1/x
#
# TVGD regularized:
#
#     A_TVGD(x) = 1 - Q(x)/x
#
# with
#
#     Q(x) = 1 - exp[-(x/xc)^n].
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import brentq

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test08_singularity_avoidance_curvature"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Numerical setup
# ============================================================

X_MIN = 1e-6
X_MAX = 20.0
N_GRID = 20000

x_grid = np.logspace(np.log10(X_MIN), np.log10(X_MAX), N_GRID)

# TVGD core parameters to audit
CASES = [
    {"name": "TVGD_n6_xc0p5", "n": 6, "xc": 0.5},
    {"name": "TVGD_n8_xc0p5", "n": 8, "xc": 0.5},
    {"name": "TVGD_n6_xc0p3", "n": 6, "xc": 0.3},
    {"name": "TVGD_n8_xc0p3", "n": 8, "xc": 0.3},
]

# Schwarzschild benchmarks in units x = r/rs
SCHW_HORIZON = 1.0
SCHW_PHOTON = 1.5
SCHW_SHADOW_B = 3.0 * np.sqrt(3.0) / 2.0  # 2.598076...

# Criteria
MAX_HORIZON_REL_ERR = 1e-3
MAX_PHOTON_REL_ERR = 1e-3
MAX_SHADOW_REL_ERR = 1e-3
MIN_SCHW_VS_TVGD_CORE_RATIO = 1e12

print("=" * 72)
print("TEST 08 — SINGULARITY AVOIDANCE / CURVATURE REGULARIZATION")
print("=" * 72)

# ============================================================
# Schwarzschild functions
# ============================================================

def A_schwarzschild(x):
    x = np.asarray(x, dtype=float)
    return 1.0 - 1.0 / x


def K_schwarzschild(x):
    """
    Dimensionless Kretschmann scalar for Schwarzschild
    in units rs = 1:

        K = 12 / x^6
    """
    x = np.asarray(x, dtype=float)
    return 12.0 / x**6


# ============================================================
# TVGD regularized core functions
# ============================================================

def Q_tvgd(x, n, xc):
    """
    Saturated structural charge:

        Q(x) = 1 - exp[-(x/xc)^n]
    """
    x = np.asarray(x, dtype=float)
    y = (x / xc)**n
    return 1.0 - np.exp(-y)


def Q_derivatives_tvgd(x, n, xc):
    """
    Analytic Q, Q', Q'' for:

        Q(x) = 1 - exp[-q(x)]
        q(x) = (x/xc)^n

    Q'  = exp(-q) q'
    Q'' = exp(-q) [q'' - (q')^2]
    """
    x = np.asarray(x, dtype=float)

    q = (x / xc)**n
    expq = np.exp(-q)

    qp = n * q / x
    qpp = n * (n - 1.0) * q / x**2

    Q = 1.0 - expq
    Qp = expq * qp
    Qpp = expq * (qpp - qp**2)

    return Q, Qp, Qpp


def A_tvgd(x, n, xc):
    x = np.asarray(x, dtype=float)
    Q = Q_tvgd(x, n=n, xc=xc)
    return 1.0 - Q / x


def A_derivatives_tvgd(x, n, xc):
    """
    A(x) = 1 - Q(x)/x

    A'  = -Q'/x + Q/x^2
    A'' = -Q''/x + 2Q'/x^2 - 2Q/x^3
    """
    x = np.asarray(x, dtype=float)
    Q, Qp, Qpp = Q_derivatives_tvgd(x, n=n, xc=xc)

    A = 1.0 - Q / x
    Ap = -Qp / x + Q / x**2
    App = -Qpp / x + 2.0 * Qp / x**2 - 2.0 * Q / x**3

    return A, Ap, App


def K_from_A(x, A, Ap, App):
    """
    Kretschmann scalar for metric:

        ds^2 = A dt^2 - A^-1 dx^2 - x^2 dOmega^2

    Dimensionless formula:

        K = A''^2 + 4(A'/x)^2 + 4((A-1)/x^2)^2

    For A = 1 - 1/x, this gives K = 12/x^6.
    """
    x = np.asarray(x, dtype=float)
    return App**2 + 4.0 * (Ap / x)**2 + 4.0 * ((A - 1.0) / x**2)**2


def K_tvgd(x, n, xc):
    A, Ap, App = A_derivatives_tvgd(x, n=n, xc=xc)
    return K_from_A(x, A, Ap, App)


# ============================================================
# Horizon / photon sphere / shadow diagnostics
# ============================================================

def find_roots_by_scan(func, x_min=1e-6, x_max=20.0, n_scan=50000):
    xs = np.logspace(np.log10(x_min), np.log10(x_max), n_scan)
    ys = func(xs)

    roots = []

    for i in range(len(xs) - 1):
        y1 = ys[i]
        y2 = ys[i + 1]

        if not np.isfinite(y1) or not np.isfinite(y2):
            continue

        if y1 == 0:
            roots.append(xs[i])

        if y1 * y2 < 0:
            try:
                root = brentq(func, xs[i], xs[i + 1], maxiter=200)
                roots.append(root)
            except Exception:
                pass

    # unique roots
    roots_unique = []
    for r in roots:
        if all(abs(r - rr) > 1e-5 for rr in roots_unique):
            roots_unique.append(r)

    return sorted(roots_unique)


def photon_condition_tvgd(x, n, xc):
    """
    Photon sphere condition for static spherical metric:

        x A'(x) - 2 A(x) = 0
    """
    A, Ap, _ = A_derivatives_tvgd(x, n=n, xc=xc)
    return x * Ap - 2.0 * A


def find_outer_photon_sphere(n, xc):
    roots = find_roots_by_scan(
        lambda xx: photon_condition_tvgd(xx, n=n, xc=xc),
        x_min=1e-4,
        x_max=20.0,
        n_scan=50000
    )

    if len(roots) == 0:
        return np.nan

    # external photon sphere is the largest positive root
    return float(max(roots))


def shadow_impact_parameter(x_ph, n, xc):
    A = float(A_tvgd(x_ph, n=n, xc=xc))

    if A <= 0 or not np.isfinite(A):
        return np.nan

    return float(x_ph / np.sqrt(A))


# ============================================================
# Run audit
# ============================================================

K_schw = K_schwarzschild(x_grid)

rows = []

for case in CASES:
    name = case["name"]
    n = case["n"]
    xc = case["xc"]

    print("\n" + "-" * 72)
    print("Case:", name, "n =", n, "xc =", xc)

    A_vals, Ap_vals, App_vals = A_derivatives_tvgd(x_grid, n=n, xc=xc)
    K_vals = K_from_A(x_grid, A_vals, Ap_vals, App_vals)

    finite_A = bool(np.all(np.isfinite(A_vals)))
    finite_K = bool(np.all(np.isfinite(K_vals)))

    K_core = float(K_vals[0])
    K_core_window_max = float(np.max(K_vals[x_grid < 1e-3]))
    K_global_max = float(np.max(K_vals))

    K_schw_core = float(K_schw[0])
    core_regularization_ratio = float(K_schw_core / max(K_core_window_max, 1e-300))

    # Horizons
    horizon_roots = find_roots_by_scan(
        lambda xx: A_tvgd(xx, n=n, xc=xc),
        x_min=1e-6,
        x_max=20.0,
        n_scan=50000
    )

    if len(horizon_roots) > 0:
        outer_horizon = float(max(horizon_roots))
    else:
        outer_horizon = np.nan

    horizon_rel_err = (
        abs(outer_horizon - SCHW_HORIZON) / SCHW_HORIZON
        if np.isfinite(outer_horizon)
        else np.nan
    )

    # Photon sphere and shadow
    x_ph = find_outer_photon_sphere(n=n, xc=xc)
    b_shadow = shadow_impact_parameter(x_ph, n=n, xc=xc) if np.isfinite(x_ph) else np.nan

    photon_rel_err = (
        abs(x_ph - SCHW_PHOTON) / SCHW_PHOTON
        if np.isfinite(x_ph)
        else np.nan
    )

    shadow_rel_err = (
        abs(b_shadow - SCHW_SHADOW_B) / SCHW_SHADOW_B
        if np.isfinite(b_shadow)
        else np.nan
    )

    # Exterior agreement with Schwarzschild
    x_ext = np.array([2.0, 3.0, 5.0, 10.0])
    A_ext_tvgd = A_tvgd(x_ext, n=n, xc=xc)
    A_ext_schw = A_schwarzschild(x_ext)

    exterior_A_abs_err_max = float(np.max(np.abs(A_ext_tvgd - A_ext_schw)))
    exterior_A_rel_err_max = float(
        np.max(np.abs(A_ext_tvgd - A_ext_schw) / np.maximum(np.abs(A_ext_schw), 1e-300))
    )

    # Criteria
    regular_core_ok = (
        finite_A
        and finite_K
        and np.isfinite(K_core_window_max)
        and core_regularization_ratio > MIN_SCHW_VS_TVGD_CORE_RATIO
    )

    exterior_horizon_ok = (
        np.isfinite(outer_horizon)
        and horizon_rel_err < MAX_HORIZON_REL_ERR
    )

    photon_ok = (
        np.isfinite(x_ph)
        and photon_rel_err < MAX_PHOTON_REL_ERR
    )

    shadow_ok = (
        np.isfinite(b_shadow)
        and shadow_rel_err < MAX_SHADOW_REL_ERR
    )

    exterior_ok = exterior_A_abs_err_max < 1e-6

    passed = (
        regular_core_ok
        and exterior_horizon_ok
        and photon_ok
        and shadow_ok
        and exterior_ok
    )

    if passed:
        status = "PASSOU_SINGULARIDADE_REGULAR_EXTERIOR_PRESERVADO"
    elif regular_core_ok and exterior_horizon_ok:
        status = "PASSOU_REGULARIDADE_MAS_EXTERIOR_TENSIONADO"
    elif regular_core_ok:
        status = "PASSOU_REGULARIDADE_NUCLEAR_APENAS"
    else:
        status = "FALHOU_REGULARIDADE"

    print("K_schw_core:", K_schw_core)
    print("K_tvgd_core:", K_core)
    print("K_tvgd_core_window_max:", K_core_window_max)
    print("core_regularization_ratio:", core_regularization_ratio)
    print("outer_horizon:", outer_horizon)
    print("x_ph:", x_ph)
    print("b_shadow:", b_shadow)
    print("status:", status)

    rows.append({
        "case": name,
        "n": n,
        "xc": xc,

        "finite_A": finite_A,
        "finite_K": finite_K,

        "K_schw_at_xmin": K_schw_core,
        "K_tvgd_at_xmin": K_core,
        "K_tvgd_core_window_max_x_lt_1e_minus_3": K_core_window_max,
        "K_tvgd_global_max": K_global_max,
        "core_regularization_ratio_Kschw_over_Ktvgd": core_regularization_ratio,

        "n_horizons_found": len(horizon_roots),
        "horizon_roots": str(horizon_roots),
        "outer_horizon": outer_horizon,
        "horizon_rel_err_vs_schwarzschild": horizon_rel_err,

        "photon_sphere": x_ph,
        "photon_rel_err_vs_schwarzschild": photon_rel_err,

        "shadow_b": b_shadow,
        "shadow_rel_err_vs_schwarzschild": shadow_rel_err,

        "exterior_A_abs_err_max_x_2_3_5_10": exterior_A_abs_err_max,
        "exterior_A_rel_err_max_x_2_3_5_10": exterior_A_rel_err_max,

        "regular_core_ok": bool(regular_core_ok),
        "exterior_horizon_ok": bool(exterior_horizon_ok),
        "photon_ok": bool(photon_ok),
        "shadow_ok": bool(shadow_ok),
        "exterior_ok": bool(exterior_ok),
        "passed": bool(passed),
        "status": status,
    })

results_df = pd.DataFrame(rows)

N_total = len(results_df)
N_passed = int(results_df["passed"].sum())
frac_passed = N_passed / N_total

if N_passed == N_total:
    final_status = "PASSOU_FORTE_SINGULARIDADE_TVGD_REGULARIZADA"
elif N_passed >= 0.5 * N_total:
    final_status = "PASSOU_MODERADO_SINGULARIDADE_TVGD_REGULARIZADA"
elif results_df["regular_core_ok"].all():
    final_status = "PASSOU_REGULARIDADE_NUCLEAR_MAS_EXTERIOR_TENSIONADO"
else:
    final_status = "TENSIONADO_SINGULARIDADE_TVGD"

summary = {
    "N_total": N_total,
    "N_passed": N_passed,
    "frac_passed": frac_passed,

    "Schwarzschild_K_at_xmin": float(K_schw[0]),
    "Schwarzschild_K_at_x_1e_minus_5": float(K_schwarzschild(1e-5)),
    "Schwarzschild_K_at_x_1e_minus_4": float(K_schwarzschild(1e-4)),
    "Schwarzschild_K_at_x_1e_minus_3": float(K_schwarzschild(1e-3)),

    "benchmark_horizon": SCHW_HORIZON,
    "benchmark_photon_sphere": SCHW_PHOTON,
    "benchmark_shadow_b": SCHW_SHADOW_B,

    "min_core_regularization_ratio": float(
        results_df["core_regularization_ratio_Kschw_over_Ktvgd"].min()
    ),
    "max_horizon_rel_err": float(
        results_df["horizon_rel_err_vs_schwarzschild"].max()
    ),
    "max_photon_rel_err": float(
        results_df["photon_rel_err_vs_schwarzschild"].max()
    ),
    "max_shadow_rel_err": float(
        results_df["shadow_rel_err_vs_schwarzschild"].max()
    ),
    "max_exterior_A_abs_err": float(
        results_df["exterior_A_abs_err_max_x_2_3_5_10"].max()
    ),

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\n" + "=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
display(results_df)
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

results_path = os.path.join(OUTDIR, "test08_singularity_case_results.csv")
summary_path = os.path.join(OUTDIR, "test08_singularity_summary.csv")

results_df.to_csv(results_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

# Plot 1: Kretschmann comparison
plt.figure(figsize=(8, 5))
plt.loglog(x_grid, K_schw, label="Schwarzschild K = 12/x^6")

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    K_vals = K_tvgd(x_grid, n=n, xc=xc)
    plt.loglog(x_grid, K_vals, label=name)

plt.axvline(1.0, linestyle="--", color="black", label="x=1 horizon scale")
plt.xlabel("x = r/rs")
plt.ylabel("dimensionless Kretschmann K")
plt.title("Curvature singularity vs TVGD regularized core")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_kretschmann_comparison.png"), dpi=220)
plt.show()

# Plot 2: Core zoom
plt.figure(figsize=(8, 5))
core_mask = x_grid < 1e-2

plt.loglog(x_grid[core_mask], K_schw[core_mask], label="Schwarzschild")

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    K_vals = K_tvgd(x_grid, n=n, xc=xc)
    plt.loglog(x_grid[core_mask], K_vals[core_mask], label=name)

plt.xlabel("x = r/rs")
plt.ylabel("dimensionless Kretschmann K")
plt.title("Core curvature zoom")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_core_curvature_zoom.png"), dpi=220)
plt.show()

# Plot 3: Metric function
plt.figure(figsize=(8, 5))
x_metric = np.linspace(0.001, 5.0, 3000)
plt.plot(x_metric, A_schwarzschild(x_metric), label="Schwarzschild A(x)")

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    plt.plot(x_metric, A_tvgd(x_metric, n=n, xc=xc), label=name)

plt.axhline(0.0, color="black", linewidth=1)
plt.axvline(1.0, color="gray", linestyle="--", label="Schwarzschild horizon")
plt.xlabel("x = r/rs")
plt.ylabel("A(x)")
plt.title("Metric function A(x)")
plt.ylim(-2, 2)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_metric_function.png"), dpi=220)
plt.show()

# Plot 4: Horizon/photon/shadow relative errors
plt.figure(figsize=(8, 5))
labels = results_df["case"].values
xpos = np.arange(len(labels))

plt.scatter(xpos, results_df["horizon_rel_err_vs_schwarzschild"], label="horizon rel. err")
plt.scatter(xpos, results_df["photon_rel_err_vs_schwarzschild"], label="photon rel. err")
plt.scatter(xpos, results_df["shadow_rel_err_vs_schwarzschild"], label="shadow rel. err")

plt.axhline(MAX_HORIZON_REL_ERR, linestyle="--", color="black", label="1e-3 criterion")
plt.xticks(xpos, labels, rotation=30, ha="right")
plt.yscale("log")
plt.ylabel("relative error")
plt.title("External Schwarzschild benchmark preservation")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_external_benchmark_errors.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(results_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_SINGULARIDADE_TVGD_REGULARIZADA":
    print("TVGD regularization removes the curvature singularity while preserving the external Schwarzschild benchmarks.")
elif "PASSOU_MODERADO" in final_status:
    print("TVGD regularization removes the singularity in most tested cases, with some external benchmark tension.")
else:
    print("The tested TVGD regularization requires further refinement.")

print("\nDONE.")
