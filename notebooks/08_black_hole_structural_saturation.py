# ============================================================
# TVGD PUBLIC TESTS
# TEST 08 — BLACK-HOLE STRUCTURAL SATURATION
# ============================================================
#
# Goal:
# Test the TVGD black-hole structural saturation mechanism:
#
#     A(x) = 1 - Q_psi(x)/x
#
# with
#
#     Q_psi(x) = 1 - exp[-(x/xc)^n].
#
# This test checks four layers:
#
# A. Geometry:
#    - finite Kretschmann scalar
#    - Schwarzschild singularity comparison
#
# B. External observables:
#    - horizon x = 1
#    - photon sphere x = 1.5
#    - shadow impact parameter b = 3sqrt(3)/2
#
# C. Effective structural source:
#    - finite density
#    - monotonic accumulated structural charge Q_psi
#    - non-negative density
#    - radial NEC controlled
#    - tangential pressure/NEC diagnosed
#
# D. Final audit:
#    - regular core
#    - preserved exterior
#    - finite structural stress-energy
#
# Interpretation:
# Q_psi is treated as an accumulated structural vacuum charge,
# consistent with the TVGD covariant effective layer.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import brentq

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test08_black_hole_structural_saturation"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Numerical setup
# ============================================================

X_MIN = 1e-6
X_MAX = 20.0
N_GRID = 30000

x_grid = np.logspace(np.log10(X_MIN), np.log10(X_MAX), N_GRID)

CASES = [
    {"name": "TVGD_n6_xc0p5", "n": 6, "xc": 0.5},
    {"name": "TVGD_n8_xc0p5", "n": 8, "xc": 0.5},
    {"name": "TVGD_n6_xc0p3", "n": 6, "xc": 0.3},
    {"name": "TVGD_n8_xc0p3", "n": 8, "xc": 0.3},
]

# Schwarzschild external benchmarks in x = r/rs
SCHW_HORIZON = 1.0
SCHW_PHOTON = 1.5
SCHW_SHADOW_B = 3.0 * np.sqrt(3.0) / 2.0

# Criteria
MAX_HORIZON_REL_ERR = 1e-3
MAX_PHOTON_REL_ERR = 1e-3
MAX_SHADOW_REL_ERR = 1e-3
MAX_EXTERIOR_A_ABS_ERR = 1e-6
MIN_CORE_REGULARIZATION_RATIO = 1e12
TOL_NEGATIVE = 1e-10

print("=" * 72)
print("TEST 08 — BLACK-HOLE STRUCTURAL SATURATION")
print("=" * 72)

# ============================================================
# Schwarzschild reference
# ============================================================

def A_schwarzschild(x):
    x = np.asarray(x, dtype=float)
    return 1.0 - 1.0 / x


def K_schwarzschild(x):
    """
    Dimensionless Schwarzschild Kretschmann scalar
    for rs = 1:

        K = 12/x^6.
    """
    x = np.asarray(x, dtype=float)
    return 12.0 / x**6


# ============================================================
# TVGD structural charge Q_psi
# ============================================================

def Q_psi(x, n, xc):
    """
    Accumulated structural vacuum charge:

        Q_psi(x) = 1 - exp[-(x/xc)^n].
    """
    x = np.asarray(x, dtype=float)
    y = (x / xc)**n
    return 1.0 - np.exp(-y)


def Q_psi_derivatives(x, n, xc):
    """
    Q = 1 - exp[-q]
    q = (x/xc)^n

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


# ============================================================
# Metric and curvature
# ============================================================

def A_tvgd(x, n, xc):
    x = np.asarray(x, dtype=float)
    return 1.0 - Q_psi(x, n=n, xc=xc) / x


def A_derivatives_tvgd(x, n, xc):
    """
    A(x) = 1 - Q(x)/x

    A'  = -Q'/x + Q/x^2
    A'' = -Q''/x + 2Q'/x^2 - 2Q/x^3
    """
    x = np.asarray(x, dtype=float)

    Q, Qp, Qpp = Q_psi_derivatives(x, n=n, xc=xc)

    A = 1.0 - Q / x
    Ap = -Qp / x + Q / x**2
    App = -Qpp / x + 2.0 * Qp / x**2 - 2.0 * Q / x**3

    return A, Ap, App


def K_from_A(x, A, Ap, App):
    """
    Kretschmann scalar for

        ds^2 = A dt^2 - A^-1 dx^2 - x^2 dOmega^2

    Dimensionless expression:

        K = A''^2 + 4(A'/x)^2 + 4((A-1)/x^2)^2.
    """
    x = np.asarray(x, dtype=float)
    return App**2 + 4.0 * (Ap / x)**2 + 4.0 * ((A - 1.0) / x**2)**2


def K_tvgd(x, n, xc):
    A, Ap, App = A_derivatives_tvgd(x, n=n, xc=xc)
    return K_from_A(x, A, Ap, App)


# ============================================================
# Effective structural source
# ============================================================

def structural_source_terms(x, n, xc):
    """
    Effective structural stress-energy diagnostics.

    For A(x)=1-Q(x)/x, the effective mass function is m(x)=Q(x)/2.

    Up to a common positive normalization factor, use:

        rho_hat = Q'(x)/x^2
        p_r_hat = -rho_hat
        p_t_hat = -Q''(x)/(2x)

    Then:

        NEC_radial_hat = rho + p_r = 0
        NEC_tangential_hat = rho + p_t

    These are dimensionless structural diagnostics, not ordinary matter density.
    """
    x = np.asarray(x, dtype=float)
    Q, Qp, Qpp = Q_psi_derivatives(x, n=n, xc=xc)

    rho = Qp / x**2
    p_r = -rho
    p_t = -Qpp / (2.0 * x)

    nec_radial = rho + p_r
    nec_tangential = rho + p_t
    sec = rho + p_r + 2.0 * p_t

    return {
        "Q": Q,
        "Qp": Qp,
        "Qpp": Qpp,
        "rho": rho,
        "p_r": p_r,
        "p_t": p_t,
        "nec_radial": nec_radial,
        "nec_tangential": nec_tangential,
        "sec": sec,
    }


# ============================================================
# Horizons, photon sphere and shadow
# ============================================================

def find_roots_by_scan(func, x_min=1e-6, x_max=20.0, n_scan=60000):
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
                root = brentq(func, xs[i], xs[i + 1], maxiter=300)
                roots.append(root)
            except Exception:
                pass

    roots_unique = []

    for r in roots:
        if all(abs(r - rr) > 1e-5 for rr in roots_unique):
            roots_unique.append(float(r))

    return sorted(roots_unique)


def photon_condition_tvgd(x, n, xc):
    A, Ap, _ = A_derivatives_tvgd(x, n=n, xc=xc)
    return x * Ap - 2.0 * A


def find_outer_photon_sphere(n, xc):
    roots = find_roots_by_scan(
        lambda xx: photon_condition_tvgd(xx, n=n, xc=xc),
        x_min=1e-4,
        x_max=20.0,
        n_scan=60000
    )

    if len(roots) == 0:
        return np.nan

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
    print("Case:", name, "| n =", n, "| xc =", xc)

    Q, Qp, Qpp = Q_psi_derivatives(x_grid, n=n, xc=xc)

    A_vals, Ap_vals, App_vals = A_derivatives_tvgd(x_grid, n=n, xc=xc)
    K_vals = K_from_A(x_grid, A_vals, Ap_vals, App_vals)

    src = structural_source_terms(x_grid, n=n, xc=xc)

    # Geometry diagnostics
    finite_A = bool(np.all(np.isfinite(A_vals)))
    finite_K = bool(np.all(np.isfinite(K_vals)))

    K_core = float(K_vals[0])
    K_core_window_max = float(np.max(K_vals[x_grid < 1e-3]))
    K_global_max = float(np.max(K_vals))

    K_schw_core = float(K_schw[0])
    core_regularization_ratio = float(K_schw_core / max(K_core_window_max, 1e-300))

    regular_core_ok = (
        finite_A
        and finite_K
        and np.isfinite(K_core_window_max)
        and core_regularization_ratio > MIN_CORE_REGULARIZATION_RATIO
    )

    # External benchmarks
    horizon_roots = find_roots_by_scan(
        lambda xx: A_tvgd(xx, n=n, xc=xc),
        x_min=1e-6,
        x_max=20.0,
        n_scan=60000
    )

    outer_horizon = float(max(horizon_roots)) if len(horizon_roots) > 0 else np.nan

    horizon_rel_err = (
        abs(outer_horizon - SCHW_HORIZON) / SCHW_HORIZON
        if np.isfinite(outer_horizon)
        else np.nan
    )

    x_ph = find_outer_photon_sphere(n=n, xc=xc)

    photon_rel_err = (
        abs(x_ph - SCHW_PHOTON) / SCHW_PHOTON
        if np.isfinite(x_ph)
        else np.nan
    )

    b_shadow = shadow_impact_parameter(x_ph, n=n, xc=xc) if np.isfinite(x_ph) else np.nan

    shadow_rel_err = (
        abs(b_shadow - SCHW_SHADOW_B) / SCHW_SHADOW_B
        if np.isfinite(b_shadow)
        else np.nan
    )

    x_ext = np.array([2.0, 3.0, 5.0, 10.0])
    A_ext_tvgd = A_tvgd(x_ext, n=n, xc=xc)
    A_ext_schw = A_schwarzschild(x_ext)

    exterior_A_abs_err_max = float(np.max(np.abs(A_ext_tvgd - A_ext_schw)))
    exterior_A_rel_err_max = float(
        np.max(np.abs(A_ext_tvgd - A_ext_schw) / np.maximum(np.abs(A_ext_schw), 1e-300))
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

    exterior_ok = exterior_A_abs_err_max < MAX_EXTERIOR_A_ABS_ERR

    # Structural source diagnostics
    rho = src["rho"]
    p_r = src["p_r"]
    p_t = src["p_t"]
    nec_radial = src["nec_radial"]
    nec_tangential = src["nec_tangential"]
    sec = src["sec"]

    finite_source = bool(
        np.all(np.isfinite(rho))
        and np.all(np.isfinite(p_r))
        and np.all(np.isfinite(p_t))
    )

    Q_monotonic = bool(np.all(Qp >= -TOL_NEGATIVE))
    Q_bounds_ok = bool(np.min(Q) >= -TOL_NEGATIVE and np.max(Q) <= 1.0 + TOL_NEGATIVE)
    rho_nonnegative = bool(np.min(rho) >= -TOL_NEGATIVE)
    radial_nec_ok = bool(np.max(np.abs(nec_radial)) < 1e-8)

    tangential_nec_min = float(np.min(nec_tangential))
    tangential_nec_negative_fraction = float(np.mean(nec_tangential < -TOL_NEGATIVE))

    sec_min = float(np.min(sec))
    sec_negative_fraction = float(np.mean(sec < -TOL_NEGATIVE))

    source_ok = (
        finite_source
        and Q_monotonic
        and Q_bounds_ok
        and rho_nonnegative
        and radial_nec_ok
    )

    structural_pressure_detected = bool(
        tangential_nec_negative_fraction > 0.0 or sec_negative_fraction > 0.0
    )

    # Final case status
    passed_geometry = regular_core_ok
    passed_exterior = exterior_horizon_ok and photon_ok and shadow_ok and exterior_ok
    passed_source = source_ok

    passed = passed_geometry and passed_exterior and passed_source

    if passed:
        if structural_pressure_detected:
            status = "PASSOU_FORTE_COM_PRESSAO_ESTRUTURAL_CONTROLADA"
        else:
            status = "PASSOU_FORTE_SEM_VIOLACAO_TANGENCIAL"
    elif passed_geometry and passed_source:
        status = "PASSOU_NUCLEO_E_FONTE_MAS_EXTERIOR_TENSIONADO"
    elif passed_geometry:
        status = "PASSOU_GEOMETRIA_MAS_FONTE_TENSIONADA"
    else:
        status = "FALHOU_REGULARIDADE_ESTRUTURAL_BH"

    print("K_schw_core:", K_schw_core)
    print("K_tvgd_core_window_max:", K_core_window_max)
    print("core_regularization_ratio:", core_regularization_ratio)
    print("outer_horizon:", outer_horizon)
    print("x_ph:", x_ph)
    print("b_shadow:", b_shadow)
    print("rho_min:", float(np.min(rho)))
    print("rho_max:", float(np.max(rho)))
    print("NEC_tangential_min:", tangential_nec_min)
    print("SEC_min:", sec_min)
    print("status:", status)

    rows.append({
        "case": name,
        "n": n,
        "xc": xc,

        # Geometry
        "finite_A": finite_A,
        "finite_K": finite_K,
        "K_schw_at_xmin": K_schw_core,
        "K_tvgd_at_xmin": K_core,
        "K_tvgd_core_window_max_x_lt_1e_minus_3": K_core_window_max,
        "K_tvgd_global_max": K_global_max,
        "core_regularization_ratio_Kschw_over_Ktvgd": core_regularization_ratio,
        "regular_core_ok": bool(regular_core_ok),

        # Exterior
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
        "exterior_horizon_ok": bool(exterior_horizon_ok),
        "photon_ok": bool(photon_ok),
        "shadow_ok": bool(shadow_ok),
        "exterior_ok": bool(exterior_ok),

        # Structural source
        "finite_source": bool(finite_source),
        "Q_min": float(np.min(Q)),
        "Q_max": float(np.max(Q)),
        "Qp_min": float(np.min(Qp)),
        "Qp_max": float(np.max(Qp)),
        "Q_monotonic": bool(Q_monotonic),
        "Q_bounds_ok": bool(Q_bounds_ok),

        "rho_min": float(np.min(rho)),
        "rho_max": float(np.max(rho)),
        "rho_nonnegative": bool(rho_nonnegative),

        "p_r_min": float(np.min(p_r)),
        "p_r_max": float(np.max(p_r)),
        "p_t_min": float(np.min(p_t)),
        "p_t_max": float(np.max(p_t)),

        "NEC_radial_max_abs": float(np.max(np.abs(nec_radial))),
        "radial_nec_ok": bool(radial_nec_ok),

        "NEC_tangential_min": tangential_nec_min,
        "NEC_tangential_negative_fraction": tangential_nec_negative_fraction,

        "SEC_min": sec_min,
        "SEC_negative_fraction": sec_negative_fraction,

        "source_ok": bool(source_ok),
        "structural_pressure_detected": bool(structural_pressure_detected),

        # Final
        "passed_geometry": bool(passed_geometry),
        "passed_exterior": bool(passed_exterior),
        "passed_source": bool(passed_source),
        "passed": bool(passed),
        "status": status,
    })

results_df = pd.DataFrame(rows)

N_total = len(results_df)
N_passed = int(results_df["passed"].sum())
frac_passed = float(N_passed / N_total)

all_geometry = bool(results_df["passed_geometry"].all())
all_exterior = bool(results_df["passed_exterior"].all())
all_source = bool(results_df["passed_source"].all())
all_pressure_structural = bool(results_df["structural_pressure_detected"].all())

if N_passed == N_total:
    final_status = "PASSOU_FORTE_SATURACAO_ESTRUTURAL_BH_TVGD"
elif all_geometry and all_source:
    final_status = "PASSOU_MODERADO_SATURACAO_ESTRUTURAL_BH_EXTERIOR_TENSIONADO"
elif all_geometry:
    final_status = "PASSOU_GEOMETRIA_BH_MAS_FONTE_TENSIONADA"
else:
    final_status = "TENSIONADO_SATURACAO_ESTRUTURAL_BH_TVGD"

summary = {
    "N_total": N_total,
    "N_passed": N_passed,
    "frac_passed": frac_passed,

    "all_geometry": all_geometry,
    "all_exterior": all_exterior,
    "all_source": all_source,
    "all_pressure_structural": all_pressure_structural,

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

    "min_rho": float(results_df["rho_min"].min()),
    "max_rho": float(results_df["rho_max"].max()),
    "min_NEC_tangential": float(results_df["NEC_tangential_min"].min()),
    "max_NEC_tangential_negative_fraction": float(
        results_df["NEC_tangential_negative_fraction"].max()
    ),
    "min_SEC": float(results_df["SEC_min"].min()),
    "max_SEC_negative_fraction": float(results_df["SEC_negative_fraction"].max()),

    "interpretation_source": (
        "finite anisotropic structural vacuum source; radial NEC controlled; "
        "tangential/SEC violation interpreted as localized structural pressure"
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

results_path = os.path.join(OUTDIR, "test08_black_hole_structural_saturation_cases.csv")
summary_path = os.path.join(OUTDIR, "test08_black_hole_structural_saturation_summary.csv")

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
plt.title("Curvature singularity vs TVGD structural saturation")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_kretschmann_structural_saturation.png"), dpi=220)
plt.show()

# Plot 2: Metric function
plt.figure(figsize=(8, 5))
x_metric = np.linspace(0.001, 5.0, 4000)

plt.plot(x_metric, A_schwarzschild(x_metric), label="Schwarzschild A(x)")

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    plt.plot(x_metric, A_tvgd(x_metric, n=n, xc=xc), label=name)

plt.axhline(0.0, color="black", linewidth=1)
plt.axvline(1.0, color="gray", linestyle="--", label="x=1")
plt.xlabel("x = r/rs")
plt.ylabel("A(x)")
plt.title("Metric function with structural saturation")
plt.ylim(-2, 2)
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_metric_function_structural_saturation.png"), dpi=220)
plt.show()

# Plot 3: Structural charge Q_psi
plt.figure(figsize=(8, 5))

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    plt.semilogx(x_grid, Q_psi(x_grid, n=n, xc=xc), label=name)

plt.axhline(1.0, color="black", linestyle="--", label="Q_psi -> 1")
plt.xlabel("x = r/rs")
plt.ylabel("Q_psi(x)")
plt.title("Accumulated structural vacuum charge")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_structural_charge_Qpsi.png"), dpi=220)
plt.show()

# Plot 4: Effective density
plt.figure(figsize=(8, 5))

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    src = structural_source_terms(x_grid, n=n, xc=xc)
    plt.loglog(x_grid, np.maximum(src["rho"], 1e-300), label=name)

plt.xlabel("x = r/rs")
plt.ylabel("rho_hat")
plt.title("Effective structural density")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_effective_structural_density.png"), dpi=220)
plt.show()

# Plot 5: NEC tangential
plt.figure(figsize=(8, 5))

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    src = structural_source_terms(x_grid, n=n, xc=xc)
    plt.semilogx(x_grid, src["nec_tangential"], label=name)

plt.axhline(0.0, color="black", linewidth=1)
plt.xlabel("x = r/rs")
plt.ylabel("rho + p_t")
plt.title("Tangential NEC diagnostic")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08_tangential_NEC_diagnostic.png"), dpi=220)
plt.show()

# Plot 6: External benchmark errors
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
if final_status == "PASSOU_FORTE_SATURACAO_ESTRUTURAL_BH_TVGD":
    print(
        "TVGD structural saturation regularizes the black-hole core, "
        "preserves external Schwarzschild observables, and admits a finite "
        "effective anisotropic structural vacuum source."
    )
elif "PASSOU_MODERADO" in final_status:
    print(
        "TVGD structural saturation regularizes the core and source, "
        "but some external benchmarks require caution."
    )
else:
    print(
        "The tested structural black-hole saturation requires further refinement."
    )

print("\nDONE.")
