
# ============================================================
# TVGD PUBLIC TESTS
# TEST 08B v2 — ROTATING BLACK-HOLE STRUCTURAL SATURATION
# ============================================================
#
# Goal:
# Test a preliminary rotating extension of the TVGD black-hole
# structural saturation mechanism.
#
# Static TVGD:
#
#     A(x) = 1 - Q_psi(x)/x
#
# Rotating Kerr-like proxy:
#
#     Delta_TVGD(x) = x^2 - Q_psi(x) x + alpha^2
#
# where:
#
#     x = r / r_s
#     alpha = a / r_s
#     0 <= alpha <= 0.5 for Kerr black holes.
#
# For Q_psi -> 1:
#
#     Delta_TVGD -> x^2 - x + alpha^2,
#
# which recovers the Kerr external Delta in units r_s = 2M.
#
# This corrected version:
#
# 1. Treats infinite regularization ratio as a strong pass.
# 2. Separates moderate rotation alpha <= 0.4 from near-extremal alpha = 0.49.
# 3. Does not require the near-extremal case to pass for the moderate-rotation status.
# 4. Reports horizon loss near extremality as a physical tension/constraint on xc.
#
# Important:
# This is not a full exact Kerr-TVGD solution.
# It is a controlled rotating consistency test.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test08B_rotating_black_hole_structural_saturation"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Numerical setup
# ============================================================

X_MIN = 1e-6
X_MAX = 20.0
N_GRID = 50000

x_grid = np.logspace(np.log10(X_MIN), np.log10(X_MAX), N_GRID)

CASES = [
    {"name": "TVGD_n6_xc0p5", "n": 6, "xc": 0.5},
    {"name": "TVGD_n8_xc0p5", "n": 8, "xc": 0.5},
    {"name": "TVGD_n6_xc0p3", "n": 6, "xc": 0.3},
    {"name": "TVGD_n8_xc0p3", "n": 8, "xc": 0.3},
]

# alpha = a/r_s. In Kerr with r_s=2M, alpha <= 0.5.
SPINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.49]

ALPHA_MODERATE_MAX = 0.4
ALPHA_NEAR_EXTREMAL = 0.49

# Criteria
MAX_HORIZON_REL_ERR = 1e-3
MAX_ERGO_EQ_REL_ERR = 1e-3
MAX_EXTERNAL_DELTA_ABS_ERR = 1e-6
MAX_EXTERNAL_GTT_ABS_ERR = 1e-6
MIN_RING_PROXY_REGULARIZATION_RATIO = 1e12

print("=" * 72)
print("TEST 08B v2 — ROTATING BLACK-HOLE STRUCTURAL SATURATION")
print("=" * 72)

# ============================================================
# Structural charge
# ============================================================

def Q_psi(x, n, xc):
    x = np.asarray(x, dtype=float)
    return 1.0 - np.exp(-(x / xc)**n)


def Q_psi_derivatives(x, n, xc):
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
# Static and rotating metric proxies
# ============================================================

def A_static_tvgd(x, n, xc):
    x = np.asarray(x, dtype=float)
    return 1.0 - Q_psi(x, n=n, xc=xc) / x


def Delta_kerr(x, alpha):
    x = np.asarray(x, dtype=float)
    return x**2 - x + alpha**2


def Delta_tvgd_rotating(x, alpha, n, xc):
    x = np.asarray(x, dtype=float)
    return x**2 - Q_psi(x, n=n, xc=xc) * x + alpha**2


def kerr_outer_horizon(alpha):
    if alpha > 0.5:
        return np.nan
    return 0.5 * (1.0 + np.sqrt(max(0.0, 1.0 - 4.0 * alpha**2)))


def kerr_inner_horizon(alpha):
    if alpha > 0.5:
        return np.nan
    return 0.5 * (1.0 - np.sqrt(max(0.0, 1.0 - 4.0 * alpha**2)))

# ============================================================
# Ergosurface proxy
# ============================================================

def gtt_factor_tvgd(x, theta, alpha, n, xc):
    """
    Kerr-like g_tt proxy:

        g_tt factor = 1 - Q_psi(x) x / Sigma

    where:

        Sigma = x^2 + alpha^2 cos^2(theta).
    """
    x = np.asarray(x, dtype=float)
    Sigma = x**2 + alpha**2 * np.cos(theta)**2
    return 1.0 - Q_psi(x, n=n, xc=xc) * x / Sigma


def gtt_factor_kerr(x, theta, alpha):
    x = np.asarray(x, dtype=float)
    Sigma = x**2 + alpha**2 * np.cos(theta)**2
    return 1.0 - x / Sigma


def kerr_ergo_radius(theta, alpha):
    val = 1.0 - 4.0 * alpha**2 * np.cos(theta)**2
    if val < 0:
        return np.nan
    return 0.5 * (1.0 + np.sqrt(val))

# ============================================================
# Ring curvature proxy
# ============================================================

def kerr_ring_curvature_proxy(x, theta, alpha):
    """
    Positive Kerr-like ring-divergence proxy:

        K_proxy ~ 12 / Sigma^3.

    At the equator theta = pi/2, Sigma = x^2,
    so K_proxy ~ 12/x^6.
    """
    x = np.asarray(x, dtype=float)
    Sigma = x**2 + alpha**2 * np.cos(theta)**2
    return 12.0 / np.maximum(Sigma, 1e-300)**3


def tvgd_ring_curvature_proxy(x, theta, alpha, n, xc):
    """
    Regularized proxy:

        K_proxy_TVGD ~ 12 Q_psi(x)^2 / Sigma^3.

    Since Q_psi ~ x^n near x=0, at the equator:

        K_proxy_TVGD ~ x^(2n - 6).

    For n=6 or n=8, the ring proxy is strongly suppressed.
    """
    x = np.asarray(x, dtype=float)
    Sigma = x**2 + alpha**2 * np.cos(theta)**2
    Q = Q_psi(x, n=n, xc=xc)
    return 12.0 * Q**2 / np.maximum(Sigma, 1e-300)**3

# ============================================================
# Root finder
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

# ============================================================
# Audit
# ============================================================

rows = []

theta_equator = np.pi / 2.0
theta_pole = 0.0

x_external = np.array([2.0, 3.0, 5.0, 10.0])
x_core = x_grid[x_grid < 1e-3]

for case in CASES:
    name = case["name"]
    n = case["n"]
    xc = case["xc"]

    for alpha in SPINS:
        print("\n" + "-" * 72)
        print("Case:", name, "| n =", n, "| xc =", xc, "| alpha =", alpha)

        # Kerr benchmark
        xh_kerr_outer = kerr_outer_horizon(alpha)
        xh_kerr_inner = kerr_inner_horizon(alpha)

        # TVGD horizon roots from Delta_TVGD = 0
        roots = find_roots_by_scan(
            lambda xx: Delta_tvgd_rotating(xx, alpha=alpha, n=n, xc=xc),
            x_min=1e-6,
            x_max=20.0,
            n_scan=60000
        )

        outer_horizon = float(max(roots)) if len(roots) > 0 else np.nan
        inner_horizon = float(min(roots)) if len(roots) > 1 else np.nan

        horizon_rel_err = (
            abs(outer_horizon - xh_kerr_outer) / max(abs(xh_kerr_outer), 1e-300)
            if np.isfinite(outer_horizon) and np.isfinite(xh_kerr_outer)
            else np.nan
        )

        # Equatorial ergosurface
        ergo_roots_eq = find_roots_by_scan(
            lambda xx: gtt_factor_tvgd(
                xx, theta=theta_equator, alpha=alpha, n=n, xc=xc
            ),
            x_min=1e-6,
            x_max=20.0,
            n_scan=60000
        )

        outer_ergo_eq = float(max(ergo_roots_eq)) if len(ergo_roots_eq) > 0 else np.nan
        kerr_ergo_eq = kerr_ergo_radius(theta_equator, alpha)

        ergo_eq_rel_err = (
            abs(outer_ergo_eq - kerr_ergo_eq) / max(abs(kerr_ergo_eq), 1e-300)
            if np.isfinite(outer_ergo_eq) and np.isfinite(kerr_ergo_eq)
            else np.nan
        )

        # Polar ergosurface
        ergo_roots_pole = find_roots_by_scan(
            lambda xx: gtt_factor_tvgd(
                xx, theta=theta_pole, alpha=alpha, n=n, xc=xc
            ),
            x_min=1e-6,
            x_max=20.0,
            n_scan=60000
        )

        outer_ergo_pole = float(max(ergo_roots_pole)) if len(ergo_roots_pole) > 0 else np.nan
        kerr_ergo_pole = kerr_ergo_radius(theta_pole, alpha)

        # External Delta and g_tt agreement
        delta_t_ext = Delta_tvgd_rotating(x_external, alpha=alpha, n=n, xc=xc)
        delta_k_ext = Delta_kerr(x_external, alpha=alpha)
        external_delta_abs_err = float(np.max(np.abs(delta_t_ext - delta_k_ext)))

        gtt_t_eq = gtt_factor_tvgd(x_external, theta_equator, alpha, n, xc)
        gtt_k_eq = gtt_factor_kerr(x_external, theta_equator, alpha)

        gtt_t_pole = gtt_factor_tvgd(x_external, theta_pole, alpha, n, xc)
        gtt_k_pole = gtt_factor_kerr(x_external, theta_pole, alpha)

        external_gtt_abs_err = float(
            max(
                np.max(np.abs(gtt_t_eq - gtt_k_eq)),
                np.max(np.abs(gtt_t_pole - gtt_k_pole)),
            )
        )

        # Ring curvature proxy at equator
        K_kerr_core = kerr_ring_curvature_proxy(
            x_core, theta_equator, alpha=max(alpha, 1e-12)
        )
        K_tvgd_core = tvgd_ring_curvature_proxy(
            x_core, theta_equator, alpha=max(alpha, 1e-12), n=n, xc=xc
        )

        K_kerr_core_max = float(np.max(K_kerr_core))
        K_tvgd_core_max = float(np.max(K_tvgd_core))

        if K_tvgd_core_max == 0.0:
            ring_proxy_regularization_ratio = np.inf
        else:
            ring_proxy_regularization_ratio = float(K_kerr_core_max / K_tvgd_core_max)

        ring_proxy_ratio_is_infinite = bool(np.isinf(ring_proxy_regularization_ratio))

        # alpha = 0 static recovery check
        if alpha == 0.0:
            static_roots = find_roots_by_scan(
                lambda xx: A_static_tvgd(xx, n=n, xc=xc),
                x_min=1e-6,
                x_max=20.0,
                n_scan=60000
            )

            static_outer = float(max(static_roots)) if len(static_roots) > 0 else np.nan
            static_recovery_err = (
                abs(static_outer - outer_horizon)
                if np.isfinite(static_outer) and np.isfinite(outer_horizon)
                else np.nan
            )
        else:
            static_outer = np.nan
            static_recovery_err = np.nan

        # Criteria
        horizon_exists = bool(np.isfinite(outer_horizon))
        horizon_ok = bool(horizon_exists and horizon_rel_err < MAX_HORIZON_REL_ERR)

        ergo_eq_ok = bool(
            np.isfinite(outer_ergo_eq)
            and ergo_eq_rel_err < MAX_ERGO_EQ_REL_ERR
        )

        external_ok = bool(
            external_delta_abs_err < MAX_EXTERNAL_DELTA_ABS_ERR
            and external_gtt_abs_err < MAX_EXTERNAL_GTT_ABS_ERR
        )

        # Corrected: infinite ratio means extremely strong regularization
        ring_proxy_ok = bool(
            np.isinf(ring_proxy_regularization_ratio)
            or ring_proxy_regularization_ratio > MIN_RING_PROXY_REGULARIZATION_RATIO
        )

        if alpha == 0.0:
            static_recovery_ok = bool(
                np.isfinite(static_recovery_err)
                and static_recovery_err < 1e-8
            )
        else:
            static_recovery_ok = True

        is_moderate = bool(alpha <= ALPHA_MODERATE_MAX)
        is_near_extremal = bool(alpha >= ALPHA_NEAR_EXTREMAL)

        passed = (
            horizon_ok
            and ergo_eq_ok
            and external_ok
            and ring_proxy_ok
            and static_recovery_ok
        )

        # Physical classification
        if passed and is_moderate:
            status = "PASSOU_ROTACAO_MODERADA_SATURACAO_ESTRUTURAL"
        elif passed and is_near_extremal:
            status = "PASSOU_QUASE_EXTREMAL_SATURACAO_ESTRUTURAL"
        elif (not horizon_ok) and is_near_extremal and external_ok and ring_proxy_ok:
            status = "TENSAO_HORIZONTE_QUASE_EXTREMAL_RESTRICAO_XC"
        elif external_ok and ring_proxy_ok:
            status = "PASSOU_EXTERIOR_E_REGULARIZACAO_PROXY_MAS_HORIZONTE_TENSIONADO"
        else:
            status = "TENSIONADO_ROTATING_BH_TVGD"

        print("Kerr outer horizon:", xh_kerr_outer)
        print("TVGD outer horizon:", outer_horizon)
        print("horizon_rel_err:", horizon_rel_err)
        print("outer_ergo_eq:", outer_ergo_eq)
        print("ergo_eq_rel_err:", ergo_eq_rel_err)
        print("external_delta_abs_err:", external_delta_abs_err)
        print("external_gtt_abs_err:", external_gtt_abs_err)
        print("ring_proxy_regularization_ratio:", ring_proxy_regularization_ratio)
        print("ring_proxy_ok:", ring_proxy_ok)
        print("status:", status)

        rows.append({
            "case": name,
            "n": n,
            "xc": xc,
            "alpha_a_over_rs": alpha,
            "is_moderate_alpha_le_0p4": is_moderate,
            "is_near_extremal_alpha_0p49": is_near_extremal,

            "kerr_outer_horizon": xh_kerr_outer,
            "kerr_inner_horizon": xh_kerr_inner,
            "tvgd_outer_horizon": outer_horizon,
            "tvgd_inner_horizon": inner_horizon,
            "n_horizon_roots": len(roots),
            "horizon_roots": str(roots),
            "horizon_rel_err_vs_kerr": horizon_rel_err,
            "horizon_exists": horizon_exists,
            "horizon_ok": horizon_ok,

            "kerr_ergo_equator": kerr_ergo_eq,
            "tvgd_ergo_equator": outer_ergo_eq,
            "ergo_equator_rel_err_vs_kerr": ergo_eq_rel_err,
            "ergo_equator_ok": ergo_eq_ok,

            "kerr_ergo_pole": kerr_ergo_pole,
            "tvgd_ergo_pole": outer_ergo_pole,

            "external_delta_abs_err_x_2_3_5_10": external_delta_abs_err,
            "external_gtt_abs_err_x_2_3_5_10": external_gtt_abs_err,
            "external_ok": external_ok,

            "K_kerr_ring_proxy_core_max": K_kerr_core_max,
            "K_tvgd_ring_proxy_core_max": K_tvgd_core_max,
            "ring_proxy_regularization_ratio": ring_proxy_regularization_ratio,
            "ring_proxy_ratio_is_infinite": ring_proxy_ratio_is_infinite,
            "ring_proxy_ok": ring_proxy_ok,

            "static_outer_horizon_alpha0": static_outer,
            "static_recovery_err_alpha0": static_recovery_err,
            "static_recovery_ok": static_recovery_ok,

            "passed": passed,
            "status": status,
        })

results_df = pd.DataFrame(rows)

# ============================================================
# Final audit
# ============================================================

N_total = len(results_df)
N_passed = int(results_df["passed"].sum())
frac_passed = float(N_passed / N_total)

moderate_df = results_df[results_df["is_moderate_alpha_le_0p4"]].copy()
near_extreme_df = results_df[results_df["is_near_extremal_alpha_0p49"]].copy()

N_moderate = len(moderate_df)
N_moderate_passed = int(moderate_df["passed"].sum())
frac_moderate_passed = float(N_moderate_passed / N_moderate)

N_near_extreme = len(near_extreme_df)
N_near_extreme_passed = int(near_extreme_df["passed"].sum())
frac_near_extreme_passed = float(N_near_extreme_passed / N_near_extreme)

N_near_extreme_horizon_tension = int(
    ((near_extreme_df["horizon_ok"] == False)
     & (near_extreme_df["external_ok"] == True)
     & (near_extreme_df["ring_proxy_ok"] == True)).sum()
)

all_external = bool(results_df["external_ok"].all())
all_ring_proxy = bool(results_df["ring_proxy_ok"].all())
all_static_recovery = bool(results_df["static_recovery_ok"].all())

all_moderate_horizon = bool(moderate_df["horizon_ok"].all())
all_moderate_ergo = bool(moderate_df["ergo_equator_ok"].all())
all_moderate_external = bool(moderate_df["external_ok"].all())
all_moderate_ring_proxy = bool(moderate_df["ring_proxy_ok"].all())

# Robust min ratio ignoring infinities for reporting
finite_ratios = results_df[
    np.isfinite(results_df["ring_proxy_regularization_ratio"])
]["ring_proxy_regularization_ratio"].values

if len(finite_ratios) > 0:
    min_finite_ring_ratio = float(np.min(finite_ratios))
else:
    min_finite_ring_ratio = np.inf

if (
    frac_moderate_passed == 1.0
    and all_moderate_horizon
    and all_moderate_ergo
    and all_moderate_external
    and all_moderate_ring_proxy
):
    if frac_near_extreme_passed == 1.0:
        final_status = "PASSOU_FORTE_ROTATING_BH_STRUCTURAL_SATURATION_TVGD"
    else:
        final_status = "PASSOU_MODERADO_ROTACAO_MODERADA_COM_TENSAO_QUASE_EXTREMAL"
elif all_external and all_ring_proxy:
    final_status = "PASSOU_FRACO_EXTERIOR_E_PROXY_REGULARIZADO_ROTACAO"
else:
    final_status = "TENSIONADO_ROTATING_BH_STRUCTURAL_SATURATION_TVGD"

summary = {
    "N_total": N_total,
    "N_passed": N_passed,
    "frac_passed": frac_passed,

    "N_moderate_alpha_le_0p4": N_moderate,
    "N_moderate_passed": N_moderate_passed,
    "frac_moderate_passed": frac_moderate_passed,

    "N_near_extreme_alpha_0p49": N_near_extreme,
    "N_near_extreme_passed": N_near_extreme_passed,
    "frac_near_extreme_passed": frac_near_extreme_passed,
    "N_near_extreme_horizon_tension": N_near_extreme_horizon_tension,

    "all_external_all_spins": all_external,
    "all_ring_proxy_regularized_all_spins": all_ring_proxy,
    "all_static_recovery_alpha0": all_static_recovery,

    "all_moderate_horizon": all_moderate_horizon,
    "all_moderate_ergo_equator": all_moderate_ergo,
    "all_moderate_external": all_moderate_external,
    "all_moderate_ring_proxy": all_moderate_ring_proxy,

    "max_horizon_rel_err_all_finite": float(
        results_df["horizon_rel_err_vs_kerr"].replace([np.inf, -np.inf], np.nan).max()
    ),
    "max_ergo_equator_rel_err": float(
        results_df["ergo_equator_rel_err_vs_kerr"].replace([np.inf, -np.inf], np.nan).max()
    ),
    "max_external_delta_abs_err": float(results_df["external_delta_abs_err_x_2_3_5_10"].max()),
    "max_external_gtt_abs_err": float(results_df["external_gtt_abs_err_x_2_3_5_10"].max()),

    "min_finite_ring_proxy_regularization_ratio": min_finite_ring_ratio,
    "N_infinite_ring_regularization_ratio": int(results_df["ring_proxy_ratio_is_infinite"].sum()),

    "max_horizon_rel_err_moderate": float(
        moderate_df["horizon_rel_err_vs_kerr"].replace([np.inf, -np.inf], np.nan).max()
    ),
    "max_ergo_equator_rel_err_moderate": float(
        moderate_df["ergo_equator_rel_err_vs_kerr"].replace([np.inf, -np.inf], np.nan).max()
    ),
    "min_finite_ring_proxy_ratio_moderate": float(
        moderate_df[np.isfinite(moderate_df["ring_proxy_regularization_ratio"])]
        ["ring_proxy_regularization_ratio"].min()
    ),

    "interpretation": (
        "preliminary Kerr-like rotating consistency test; moderate rotation "
        "is evaluated separately from near-extremal alpha=0.49; infinite "
        "ring regularization ratio is treated as strong regularization"
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

results_path = os.path.join(
    OUTDIR,
    "test08B_rotating_black_hole_structural_saturation_cases.csv"
)

summary_path = os.path.join(
    OUTDIR,
    "test08B_rotating_black_hole_structural_saturation_summary.csv"
)

results_df.to_csv(results_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

# Plot 1: horizon vs spin
plt.figure(figsize=(8, 5))

for case in CASES:
    name = case["name"]
    sub = results_df[results_df["case"] == name]
    plt.plot(
        sub["alpha_a_over_rs"],
        sub["tvgd_outer_horizon"],
        marker="o",
        label=name
    )

alpha_plot = np.linspace(0.0, 0.5, 300)
plt.plot(
    alpha_plot,
    [kerr_outer_horizon(a) for a in alpha_plot],
    color="black",
    linestyle="--",
    label="Kerr outer horizon"
)

plt.xlabel("alpha = a/r_s")
plt.ylabel("outer horizon x_+")
plt.title("Rotating TVGD horizon vs Kerr")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08B_horizon_vs_spin.png"), dpi=220)
plt.show()

# Plot 2: equatorial ergosurface
plt.figure(figsize=(8, 5))

for case in CASES:
    name = case["name"]
    sub = results_df[results_df["case"] == name]
    plt.plot(
        sub["alpha_a_over_rs"],
        sub["tvgd_ergo_equator"],
        marker="o",
        label=name
    )

plt.axhline(1.0, color="black", linestyle="--", label="Kerr equatorial ergosurface")
plt.xlabel("alpha = a/r_s")
plt.ylabel("equatorial ergosurface x")
plt.title("Equatorial ergosurface")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08B_equatorial_ergosphere.png"), dpi=220)
plt.show()

# Plot 3: Delta comparison at alpha=0.4
plt.figure(figsize=(8, 5))

alpha_show = 0.4
x_plot = np.linspace(0.001, 3.0, 3000)

plt.plot(
    x_plot,
    Delta_kerr(x_plot, alpha_show),
    color="black",
    linestyle="--",
    label="Kerr Delta"
)

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    plt.plot(
        x_plot,
        Delta_tvgd_rotating(x_plot, alpha_show, n=n, xc=xc),
        label=name
    )

plt.axhline(0.0, color="black", linewidth=1)
plt.xlabel("x = r/r_s")
plt.ylabel("Delta(x)")
plt.title("Delta function at alpha=0.4")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08B_delta_alpha0p4.png"), dpi=220)
plt.show()

# Plot 4: ring curvature proxy
plt.figure(figsize=(8, 5))

alpha_show = 0.4
theta_show = np.pi / 2.0
x_core_plot = np.logspace(-6, -1, 1000)

plt.loglog(
    x_core_plot,
    kerr_ring_curvature_proxy(x_core_plot, theta_show, alpha_show),
    color="black",
    linestyle="--",
    label="Kerr ring proxy"
)

for case in CASES:
    n = case["n"]
    xc = case["xc"]
    name = case["name"]
    plt.loglog(
        x_core_plot,
        tvgd_ring_curvature_proxy(x_core_plot, theta_show, alpha_show, n=n, xc=xc),
        label=name
    )

plt.xlabel("x = r/r_s")
plt.ylabel("curvature proxy")
plt.title("Equatorial ring curvature proxy regularization")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08B_ring_curvature_proxy.png"), dpi=220)
plt.show()

# Plot 5: horizon relative error
plt.figure(figsize=(8, 5))

for case in CASES:
    name = case["name"]
    sub = results_df[results_df["case"] == name]
    plt.plot(
        sub["alpha_a_over_rs"],
        sub["horizon_rel_err_vs_kerr"],
        marker="o",
        label=f"{name}"
    )

plt.axhline(MAX_HORIZON_REL_ERR, color="black", linestyle="--", label="criterion")
plt.axvline(ALPHA_MODERATE_MAX, color="gray", linestyle=":", label="moderate limit")
plt.yscale("log")
plt.xlabel("alpha = a/r_s")
plt.ylabel("relative error")
plt.title("Outer horizon relative error vs Kerr")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test08B_horizon_relative_error.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(results_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_ROTATING_BH_STRUCTURAL_SATURATION_TVGD":
    print(
        "The preliminary rotating TVGD structural saturation preserves the Kerr external limit, "
        "keeps the horizon/ergoregion structure, and regularizes the equatorial ring proxy "
        "for all tested spins."
    )
elif final_status == "PASSOU_MODERADO_ROTACAO_MODERADA_COM_TENSAO_QUASE_EXTREMAL":
    print(
        "The preliminary rotating TVGD structural saturation passes for moderate rotation "
        "and regularizes the ring proxy, while near-extremal rotation shows a horizon constraint "
        "for some core scales."
    )
elif "PASSOU_FRACO" in final_status:
    print(
        "The rotating TVGD proxy preserves the external Kerr limit and regularizes the ring proxy, "
        "but horizon/ergoregion diagnostics require caution."
    )
else:
    print(
        "The rotating TVGD structural saturation proxy is tensioned and requires refinement."
    )

print("\nDONE.")
