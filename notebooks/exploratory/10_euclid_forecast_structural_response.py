
# ============================================================
# TVGD PUBLIC TESTS
# TEST 10 — EUCLID FORECAST / STRUCTURAL RESPONSE DETECTABILITY
# ============================================================
#
# Goal:
# Forecast whether a Euclid-like survey can distinguish a TVGD-like
# structural response evolution from a constant/no-evolution model.
#
# This is not a claim of Euclid detection.
# It is a falsifiability forecast:
#
#     What redshift baseline, sample size, and observational scatter
#     are needed to test the TVGD structural response?
#
# Structural model:
#
#     a0(z) = a0_0 * (1+z)^beta
#
# TVGD reference:
#
#     beta = 1
#
# Null model:
#
#     beta = 0
#
# Forecast observable:
#
#     y(z) = log10[a0(z)/a0_0] = beta * log10(1+z)
#
# Detection significance:
#
#     sigma_beta = sigma_y / sqrt(sum((x_i - mean(x))^2))
#
# where:
#
#     x_i = log10(1+z_i)
#
# A broad redshift baseline is more important than sample size alone.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test10_euclid_forecast_structural_response"
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 72)
print("TEST 10 — EUCLID FORECAST / STRUCTURAL RESPONSE DETECTABILITY")
print("=" * 72)

# ============================================================
# Forecast setup
# ============================================================

BETA_TVGD = 1.0
BETA_NULL = 0.0

# Effective observational scatter in log10(a0) or equivalent structural proxy.
# Conservative values: 0.75 dex resembles noisy galaxy samples.
# Better future-quality samples may reach 0.30, 0.20, or 0.10 dex.
SIGMA_DEX_VALUES = [0.75, 0.50, 0.30, 0.20, 0.10]

# Euclid-like and extended redshift baselines.
REDSHIFT_SCENARIOS = [
    {"scenario": "low_z_0p1_1p0", "z_min": 0.1, "z_max": 1.0},
    {"scenario": "euclid_core_0p5_2p0", "z_min": 0.5, "z_max": 2.0},
    {"scenario": "wide_0p1_2p0", "z_min": 0.1, "z_max": 2.0},
    {"scenario": "wide_0p1_3p0", "z_min": 0.1, "z_max": 3.0},
    {"scenario": "deep_0p5_3p0", "z_min": 0.5, "z_max": 3.0},
    {"scenario": "very_wide_0p1_5p0", "z_min": 0.1, "z_max": 5.0},
]

N_VALUES = [50, 100, 200, 400, 800, 1200, 2000, 5000, 10000]

# Detection thresholds
THRESHOLD_3SIGMA = 3.0
THRESHOLD_5SIGMA = 5.0

# Monte Carlo setup for one validation block
MC_SEED = 12345
MC_N_REALIZATIONS = 1000

rng = np.random.default_rng(MC_SEED)

# ============================================================
# Helper functions
# ============================================================

def x_of_z(z):
    return np.log10(1.0 + np.asarray(z, dtype=float))


def y_tvgd(z, beta=BETA_TVGD):
    return beta * x_of_z(z)


def forecast_sigma_beta(z_values, sigma_dex):
    """
    Linear model:

        y = beta x + intercept

    with constant Gaussian scatter sigma_dex.

    The uncertainty on beta is:

        sigma_beta = sigma / sqrt(sum((x - mean(x))^2))
    """
    x = x_of_z(z_values)
    sxx = np.sum((x - np.mean(x))**2)

    if sxx <= 0:
        return np.inf

    return sigma_dex / np.sqrt(sxx)


def detection_significance(beta_true, beta_null, sigma_beta):
    if not np.isfinite(sigma_beta) or sigma_beta <= 0:
        return 0.0
    return abs(beta_true - beta_null) / sigma_beta


def generate_redshifts(z_min, z_max, n, distribution="uniform"):
    """
    Simple forecast redshift distributions.

    uniform:
        z uniformly distributed.

    volume_like:
        approximate increasing counts with z by sampling z from
        a distribution weighted toward high z.
    """
    if distribution == "uniform":
        return rng.uniform(z_min, z_max, n)

    if distribution == "volume_like":
        u = rng.uniform(0.0, 1.0, n)
        return z_min + (z_max - z_min) * np.sqrt(u)

    raise ValueError("Unknown distribution")


def fit_beta_ols(z_values, y_values):
    """
    Fit y = intercept + beta x.
    """
    x = x_of_z(z_values)
    X = np.vstack([np.ones_like(x), x]).T
    coeff, *_ = np.linalg.lstsq(X, y_values, rcond=None)
    intercept, beta_hat = coeff
    return float(intercept), float(beta_hat)


# ============================================================
# Analytic forecast scan
# ============================================================

rows = []

for scenario in REDSHIFT_SCENARIOS:
    z_min = scenario["z_min"]
    z_max = scenario["z_max"]
    scenario_name = scenario["scenario"]

    for n in N_VALUES:
        # deterministic evenly spaced redshifts for analytic forecast
        z_values = np.linspace(z_min, z_max, n)

        x_span = float(x_of_z(z_max) - x_of_z(z_min))
        x_std = float(np.std(x_of_z(z_values)))

        for sigma_dex in SIGMA_DEX_VALUES:
            sigma_beta = forecast_sigma_beta(z_values, sigma_dex)
            sig = detection_significance(BETA_TVGD, BETA_NULL, sigma_beta)

            detectable_3sigma = sig >= THRESHOLD_3SIGMA
            detectable_5sigma = sig >= THRESHOLD_5SIGMA

            rows.append({
                "scenario": scenario_name,
                "z_min": z_min,
                "z_max": z_max,
                "N": n,
                "sigma_dex": sigma_dex,

                "beta_true_TVGD": BETA_TVGD,
                "beta_null": BETA_NULL,

                "x_span_log10_1pz": x_span,
                "x_std_log10_1pz": x_std,

                "sigma_beta_forecast": float(sigma_beta),
                "significance_beta1_vs_beta0": float(sig),

                "detectable_3sigma": bool(detectable_3sigma),
                "detectable_5sigma": bool(detectable_5sigma),
            })

forecast_df = pd.DataFrame(rows)

# ============================================================
# Minimum N table
# ============================================================

min_rows = []

for scenario in REDSHIFT_SCENARIOS:
    scenario_name = scenario["scenario"]
    for sigma_dex in SIGMA_DEX_VALUES:
        sub = forecast_df[
            (forecast_df["scenario"] == scenario_name)
            & (forecast_df["sigma_dex"] == sigma_dex)
        ].copy()

        sub3 = sub[sub["detectable_3sigma"]]
        sub5 = sub[sub["detectable_5sigma"]]

        N_3sigma = int(sub3["N"].min()) if len(sub3) > 0 else np.nan
        N_5sigma = int(sub5["N"].min()) if len(sub5) > 0 else np.nan

        max_sig_at_maxN = float(sub["significance_beta1_vs_beta0"].max())

        min_rows.append({
            "scenario": scenario_name,
            "sigma_dex": sigma_dex,
            "N_required_3sigma": N_3sigma,
            "N_required_5sigma": N_5sigma,
            "max_significance_in_grid": max_sig_at_maxN,
        })

minN_df = pd.DataFrame(min_rows)

# ============================================================
# Monte Carlo validation for representative cases
# ============================================================

mc_cases = [
    {"scenario": "euclid_core_0p5_2p0", "z_min": 0.5, "z_max": 2.0, "N": 800, "sigma_dex": 0.30},
    {"scenario": "wide_0p1_3p0", "z_min": 0.1, "z_max": 3.0, "N": 400, "sigma_dex": 0.30},
    {"scenario": "very_wide_0p1_5p0", "z_min": 0.1, "z_max": 5.0, "N": 200, "sigma_dex": 0.30},
    {"scenario": "euclid_core_0p5_2p0", "z_min": 0.5, "z_max": 2.0, "N": 2000, "sigma_dex": 0.20},
]

mc_rows = []

for case in mc_cases:
    scenario_name = case["scenario"]
    z_min = case["z_min"]
    z_max = case["z_max"]
    N = case["N"]
    sigma_dex = case["sigma_dex"]

    beta_hats = []

    for _ in range(MC_N_REALIZATIONS):
        z_values = generate_redshifts(z_min, z_max, N, distribution="uniform")
        y_true = y_tvgd(z_values, beta=BETA_TVGD)
        y_obs = y_true + rng.normal(0.0, sigma_dex, size=N)

        _, beta_hat = fit_beta_ols(z_values, y_obs)
        beta_hats.append(beta_hat)

    beta_hats = np.array(beta_hats)

    beta_mean = float(np.mean(beta_hats))
    beta_std = float(np.std(beta_hats, ddof=1))
    p_beta_gt_0 = float(np.mean(beta_hats > 0.0))
    p_beta_gt_0p5 = float(np.mean(beta_hats > 0.5))
    p_beta_between_0p5_1p5 = float(np.mean((beta_hats > 0.5) & (beta_hats < 1.5)))

    sig_mc = abs(beta_mean - BETA_NULL) / beta_std if beta_std > 0 else np.inf

    mc_rows.append({
        "scenario": scenario_name,
        "z_min": z_min,
        "z_max": z_max,
        "N": N,
        "sigma_dex": sigma_dex,
        "MC_realizations": MC_N_REALIZATIONS,
        "beta_hat_mean": beta_mean,
        "beta_hat_std": beta_std,
        "MC_significance_vs_beta0": float(sig_mc),
        "P_beta_gt_0": p_beta_gt_0,
        "P_beta_gt_0p5": p_beta_gt_0p5,
        "P_beta_between_0p5_1p5": p_beta_between_0p5_1p5,
    })

mc_df = pd.DataFrame(mc_rows)

# ============================================================
# Classification
# ============================================================

# Strong condition:
# at least one Euclid-like or wide scenario reaches 5 sigma with N <= 2000
# at sigma <= 0.30 dex.
strong_candidates = minN_df[
    (minN_df["sigma_dex"] <= 0.30)
    & (minN_df["N_required_5sigma"].notna())
    & (minN_df["N_required_5sigma"] <= 2000)
].copy()

# Moderate condition:
# at least one scenario reaches 3 sigma with N <= 2000 at sigma <= 0.50 dex.
moderate_candidates = minN_df[
    (minN_df["sigma_dex"] <= 0.50)
    & (minN_df["N_required_3sigma"].notna())
    & (minN_df["N_required_3sigma"] <= 2000)
].copy()

# Euclid-core specific condition
euclid_core = minN_df[minN_df["scenario"] == "euclid_core_0p5_2p0"].copy()

euclid_core_3sigma_sigma03 = euclid_core[
    (euclid_core["sigma_dex"] == 0.30)
    & (euclid_core["N_required_3sigma"].notna())
]

euclid_core_5sigma_sigma03 = euclid_core[
    (euclid_core["sigma_dex"] == 0.30)
    & (euclid_core["N_required_5sigma"].notna())
]

if len(strong_candidates) > 0:
    final_status = "PASSOU_FORTE_FORECAST_EUCLID_TVGD_DETECTAVEL"
elif len(moderate_candidates) > 0:
    final_status = "PASSOU_MODERADO_FORECAST_EUCLID_TVGD_TESTAVEL"
else:
    final_status = "TENSIONADO_FORECAST_EUCLID_TVGD_DIFICIL"

best_row = forecast_df.sort_values(
    "significance_beta1_vs_beta0",
    ascending=False
).iloc[0]

best_minN_3 = minN_df[minN_df["N_required_3sigma"].notna()].sort_values(
    "N_required_3sigma"
).iloc[0]

best_minN_5 = minN_df[minN_df["N_required_5sigma"].notna()].sort_values(
    "N_required_5sigma"
).iloc[0]

summary = {
    "beta_TVGD_reference": BETA_TVGD,
    "beta_null_reference": BETA_NULL,

    "N_forecast_rows": len(forecast_df),
    "N_minN_rows": len(minN_df),
    "N_MC_cases": len(mc_df),

    "N_strong_candidates": len(strong_candidates),
    "N_moderate_candidates": len(moderate_candidates),

    "best_scenario": best_row["scenario"],
    "best_N": int(best_row["N"]),
    "best_sigma_dex": float(best_row["sigma_dex"]),
    "best_significance": float(best_row["significance_beta1_vs_beta0"]),

    "best_3sigma_scenario": best_minN_3["scenario"],
    "best_3sigma_sigma_dex": float(best_minN_3["sigma_dex"]),
    "best_3sigma_N_required": int(best_minN_3["N_required_3sigma"]),

    "best_5sigma_scenario": best_minN_5["scenario"],
    "best_5sigma_sigma_dex": float(best_minN_5["sigma_dex"]),
    "best_5sigma_N_required": int(best_minN_5["N_required_5sigma"]),

    "euclid_core_sigma0p30_N_required_3sigma": (
        int(euclid_core_3sigma_sigma03["N_required_3sigma"].iloc[0])
        if len(euclid_core_3sigma_sigma03) > 0
        else np.nan
    ),
    "euclid_core_sigma0p30_N_required_5sigma": (
        int(euclid_core_5sigma_sigma03["N_required_5sigma"].iloc[0])
        if len(euclid_core_5sigma_sigma03) > 0
        else np.nan
    ),

    "interpretation": (
        "Euclid-like redshift baselines can test TVGD structural evolution if "
        "the effective scatter is sufficiently controlled. Broad redshift coverage "
        "is more important than sample size alone."
    ),

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

# ============================================================
# Display
# ============================================================

print("\nForecast results:")
display(forecast_df)

print("\nMinimum N requirements:")
display(minN_df)

print("\nMonte Carlo validation:")
display(mc_df)

print("\nSummary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

forecast_path = os.path.join(OUTDIR, "test10_euclid_forecast_grid.csv")
minN_path = os.path.join(OUTDIR, "test10_euclid_minimum_sample_requirements.csv")
mc_path = os.path.join(OUTDIR, "test10_euclid_monte_carlo_validation.csv")
summary_path = os.path.join(OUTDIR, "test10_euclid_forecast_summary.csv")

forecast_df.to_csv(forecast_path, index=False)
minN_df.to_csv(minN_path, index=False)
mc_df.to_csv(mc_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

# Plot 1: significance vs N for sigma=0.30
plt.figure(figsize=(8, 5))

for scenario in REDSHIFT_SCENARIOS:
    scenario_name = scenario["scenario"]
    sub = forecast_df[
        (forecast_df["scenario"] == scenario_name)
        & (forecast_df["sigma_dex"] == 0.30)
    ].copy()

    plt.plot(
        sub["N"],
        sub["significance_beta1_vs_beta0"],
        marker="o",
        label=scenario_name
    )

plt.axhline(3.0, color="black", linestyle="--", label="3σ")
plt.axhline(5.0, color="gray", linestyle="--", label="5σ")
plt.xscale("log")
plt.xlabel("Sample size N")
plt.ylabel("Detection significance for beta=1 vs beta=0")
plt.title("TVGD structural evolution forecast, sigma=0.30 dex")
plt.grid(True, which="both", alpha=0.3)
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test10_significance_vs_N_sigma0p30.png"), dpi=220)
plt.show()

# Plot 2: N required for 3 sigma
plt.figure(figsize=(9, 5))

plot_df = minN_df[minN_df["sigma_dex"] == 0.30].copy()
plt.bar(plot_df["scenario"], plot_df["N_required_3sigma"])
plt.xticks(rotation=45, ha="right")
plt.ylabel("N required for 3σ")
plt.title("Minimum sample size for 3σ detection, sigma=0.30 dex")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test10_N_required_3sigma_sigma0p30.png"), dpi=220)
plt.show()

# Plot 3: N required for 5 sigma
plt.figure(figsize=(9, 5))

plot_df = minN_df[minN_df["sigma_dex"] == 0.30].copy()
plt.bar(plot_df["scenario"], plot_df["N_required_5sigma"])
plt.xticks(rotation=45, ha="right")
plt.ylabel("N required for 5σ")
plt.title("Minimum sample size for 5σ detection, sigma=0.30 dex")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test10_N_required_5sigma_sigma0p30.png"), dpi=220)
plt.show()

# Plot 4: structural evolution curves
plt.figure(figsize=(8, 5))

z_plot = np.linspace(0.0, 5.0, 500)

for beta in [0.0, 0.5, 1.0, 1.5]:
    plt.plot(
        z_plot,
        y_tvgd(z_plot, beta=beta),
        label=f"beta={beta}"
    )

plt.xlabel("z")
plt.ylabel("log10[a0(z)/a0(0)]")
plt.title("Structural evolution forecast signatures")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test10_structural_evolution_curves.png"), dpi=220)
plt.show()

# Plot 5: Monte Carlo beta recovery
plt.figure(figsize=(8, 5))

plt.errorbar(
    mc_df["scenario"] + "\nN=" + mc_df["N"].astype(str),
    mc_df["beta_hat_mean"],
    yerr=mc_df["beta_hat_std"],
    fmt="o",
    capsize=4
)

plt.axhline(1.0, color="black", linestyle="--", label="TVGD beta=1")
plt.axhline(0.0, color="gray", linestyle="--", label="null beta=0")
plt.xticks(rotation=30, ha="right")
plt.ylabel("Recovered beta")
plt.title("Monte Carlo recovery of TVGD structural evolution")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test10_monte_carlo_beta_recovery.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(forecast_path)
print(minN_path)
print(mc_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_FORECAST_EUCLID_TVGD_DETECTAVEL":
    print(
        "The TVGD structural evolution is forecast to be detectable by Euclid-like "
        "or broad-redshift surveys, provided scatter is controlled. The prediction "
        "is falsifiable: beta=1 can be tested against beta=0."
    )
elif final_status == "PASSOU_MODERADO_FORECAST_EUCLID_TVGD_TESTAVEL":
    print(
        "The TVGD structural evolution is testable, but requires favorable scatter "
        "and/or broad redshift coverage."
    )
else:
    print(
        "The forecast indicates that the current Euclid-like configuration may have "
        "difficulty testing the TVGD structural evolution without improved precision "
        "or broader redshift coverage."
    )

print("\nDONE.")
