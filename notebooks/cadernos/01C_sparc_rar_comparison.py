# ============================================================
# TVGD PUBLIC TESTS
# TEST 01C — SPARC/RAR: TVGD vs EMPIRICAL RAR / MOND CONTROLS
# Public reproducibility script
# ============================================================

"""
Test 01C — SPARC/RAR: TVGD vs empirical RAR/MOND controls

This script compares the effective TVGD acceleration law against:

1. baryonic-only acceleration;
2. empirical RAR with fixed g_dagger;
3. empirical RAR with optimized g_dagger;
4. MOND simple interpolating function;
5. MOND standard interpolating function.

The main goal is not to claim that TVGD "defeats" RAR/MOND in this test.
The goal is to show that TVGD reproduces the empirical RAR behavior as
an effective realization of the same acceleration relation, while being
embedded in the broader TVGD theoretical program.

Author: Marcelo Bauman
"""

import os
import zipfile
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuration
# ============================================================

OUTDIR = "Resultados/test01C_sparc_rar_comparison"
DATADIR = os.path.join(OUTDIR, "data")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

SPARC_URLS = [
    "http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip",
    "https://zenodo.org/records/16284118/files/Rotmod_LTG.zip?download=1",
]

ZIP_PATH = os.path.join(DATADIR, "Rotmod_LTG.zip")

# Constants
A0_TVGD = 1.2e-10
GDAG_RAR = 1.2e-10

KPC_TO_M = 3.0856775814913673e19
KM_TO_M = 1000.0

# Fiducial stellar mass-to-light ratios
UPS_DISK = 0.5
UPS_BULGE = 0.7

MIN_POINTS_GALAXY = 3

# Grid for empirical RAR optimization
GDAG_GRID = np.linspace(0.6e-10, 2.0e-10, 80)


# ============================================================
# Model functions
# ============================================================

def P_tvgd(g_bar, a0=A0_TVGD):
    """
    TVGD structural response function.

    P(u) = 1 - exp(-u)

    with

    u = sqrt(g_bar / a0).
    """
    g_bar = np.asarray(g_bar, dtype=float)
    u = np.sqrt(np.maximum(g_bar, 0.0) / a0)
    return 1.0 - np.exp(-u)


def g_tvgd(g_bar, a0=A0_TVGD):
    """
    Effective TVGD acceleration.

    g_TVGD = g_bar / [1 - exp(-sqrt(g_bar/a0))]
    """
    g_bar = np.asarray(g_bar, dtype=float)
    P = P_tvgd(g_bar, a0=a0)
    return g_bar / np.maximum(P, 1e-300)


def g_rar_empirical(g_bar, gdag=GDAG_RAR):
    """
    Empirical RAR form:

    g_obs = g_bar / [1 - exp(-sqrt(g_bar/gdag))]

    This is mathematically identical to the tested TVGD response when
    gdag = a0.
    """
    g_bar = np.asarray(g_bar, dtype=float)
    x = np.sqrt(np.maximum(g_bar, 0.0) / gdag)
    denom = 1.0 - np.exp(-x)
    return g_bar / np.maximum(denom, 1e-300)


def g_mond_simple(g_bar, a0=A0_TVGD):
    """
    MOND control with simple interpolating function.

    Approximately solves:

    g * mu(g/a0) = g_bar

    with

    mu(x) = x/(1+x)

    giving

    g = 0.5 * (g_bar + sqrt(g_bar^2 + 4 g_bar a0)).
    """
    g_bar = np.asarray(g_bar, dtype=float)
    return 0.5 * (g_bar + np.sqrt(g_bar**2 + 4.0 * g_bar * a0))


def g_mond_standard(g_bar, a0=A0_TVGD):
    """
    MOND control with standard interpolating function:

    mu(x) = x / sqrt(1+x^2)

    Algebraic approximation:

    g = g_bar * sqrt(0.5 + 0.5 * sqrt(1 + 4(a0/g_bar)^2)).
    """
    g_bar = np.asarray(g_bar, dtype=float)
    x = np.maximum(g_bar, 1e-300)
    return x * np.sqrt(0.5 + 0.5 * np.sqrt(1.0 + 4.0 * (a0 / x)**2))


def log_residual(g_obs, g_model):
    """
    Logarithmic residual:

    log10(g_model) - log10(g_obs)
    """
    g_obs = np.asarray(g_obs, dtype=float)
    g_model = np.asarray(g_model, dtype=float)

    out = np.full_like(g_obs, np.nan, dtype=float)

    mask = (
        np.isfinite(g_obs)
        & np.isfinite(g_model)
        & (g_obs > 0)
        & (g_model > 0)
    )

    out[mask] = np.log10(g_model[mask]) - np.log10(g_obs[mask])
    return out


def metrics(g_obs, g_model):
    """
    Basic log-space metrics.
    """
    r = log_residual(g_obs, g_model)
    r = r[np.isfinite(r)]

    if len(r) == 0:
        return {
            "N": 0,
            "mean_abs_log": np.nan,
            "median_abs_log": np.nan,
            "rmse_log": np.nan,
            "bias_log": np.nan,
        }

    return {
        "N": int(len(r)),
        "mean_abs_log": float(np.mean(np.abs(r))),
        "median_abs_log": float(np.median(np.abs(r))),
        "rmse_log": float(np.sqrt(np.mean(r**2))),
        "bias_log": float(np.mean(r)),
    }


# ============================================================
# SPARC parser
# ============================================================

def read_rotmod_file(path):
    """
    Reads a SPARC *_rotmod.dat file.

    Typical format:

        Rad  Vobs  errV  Vgas  Vdisk  Vbul  SBdisk  SBbul
    """
    rows = []

    with open(path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#") or line.startswith(";"):
                continue

            parts = line.replace(",", " ").split()

            vals = []
            for p in parts:
                try:
                    vals.append(float(p))
                except Exception:
                    pass

            if len(vals) >= 6:
                vals = vals[:8] if len(vals) >= 8 else vals + [np.nan] * (8 - len(vals))
                rows.append(vals)

    if len(rows) == 0:
        return None

    df = pd.DataFrame(
        rows,
        columns=[
            "R_kpc",
            "Vobs_km_s",
            "eVobs_km_s",
            "Vgas_km_s",
            "Vdisk_km_s",
            "Vbul_km_s",
            "SBdisk",
            "SBbul",
        ],
    )

    name = os.path.basename(path).replace("_rotmod.dat", "")
    df["galaxy"] = name

    return df


# ============================================================
# Download and load SPARC
# ============================================================

print("=" * 72)
print("TEST 01C — SPARC/RAR: TVGD vs RAR/MOND CONTROLS")
print("=" * 72)

if not os.path.exists(ZIP_PATH):
    ok = False
    last_error = None

    for url in SPARC_URLS:
        try:
            print(f"Downloading SPARC from: {url}")
            urllib.request.urlretrieve(url, ZIP_PATH)
            ok = True
            print("Download completed.")
            break
        except Exception as e:
            last_error = e
            print("Download failed:", repr(e))

    if not ok:
        raise RuntimeError(f"Could not download Rotmod_LTG.zip. Error: {last_error}")
else:
    print("ZIP file already exists:", ZIP_PATH)

with zipfile.ZipFile(ZIP_PATH, "r") as zf:
    zf.extractall(DATADIR)

rotmod_files = []

for root, dirs, files in os.walk(DATADIR):
    for f in files:
        if f.lower().endswith("_rotmod.dat"):
            rotmod_files.append(os.path.join(root, f))

rotmod_files = sorted(rotmod_files)

print("Rotmod files found:", len(rotmod_files))

if len(rotmod_files) == 0:
    raise RuntimeError("No *_rotmod.dat files found.")


# ============================================================
# Build RAR table
# ============================================================

all_dfs = []
bad_files = []

for path in rotmod_files:
    tmp = read_rotmod_file(path)

    if tmp is None:
        bad_files.append(path)
    else:
        all_dfs.append(tmp)

raw = pd.concat(all_dfs, ignore_index=True)

print("Raw points:", len(raw))
print("Raw galaxies:", raw["galaxy"].nunique())
print("Bad files:", len(bad_files))

df = raw.copy()
df = df.replace([np.inf, -np.inf], np.nan)

df = df.dropna(
    subset=[
        "R_kpc",
        "Vobs_km_s",
        "Vgas_km_s",
        "Vdisk_km_s",
        "Vbul_km_s",
    ]
)

df = df[
    (df["R_kpc"] > 0)
    & (df["Vobs_km_s"] > 0)
].copy()

R_m = df["R_kpc"].values * KPC_TO_M

Vobs_m_s = df["Vobs_km_s"].values * KM_TO_M
Vgas_m_s = df["Vgas_km_s"].values * KM_TO_M
Vdisk_m_s = df["Vdisk_km_s"].values * KM_TO_M
Vbul_m_s = df["Vbul_km_s"].values * KM_TO_M

df["g_obs"] = Vobs_m_s**2 / R_m

g_gas = Vgas_m_s * np.abs(Vgas_m_s) / R_m
g_disk = UPS_DISK * Vdisk_m_s * np.abs(Vdisk_m_s) / R_m
g_bul = UPS_BULGE * Vbul_m_s * np.abs(Vbul_m_s) / R_m

df["g_bar"] = g_gas + g_disk + g_bul

df = df[
    np.isfinite(df["g_obs"])
    & np.isfinite(df["g_bar"])
    & (df["g_obs"] > 0)
    & (df["g_bar"] > 0)
].copy()

counts = df.groupby("galaxy").size()
good_galaxies = counts[counts >= MIN_POINTS_GALAXY].index
df = df[df["galaxy"].isin(good_galaxies)].copy()

print("Valid points:", len(df))
print("Valid galaxies:", df["galaxy"].nunique())


# ============================================================
# Models
# ============================================================

g_obs = df["g_obs"].values
g_bar = df["g_bar"].values

df["g_baryonic"] = g_bar
df["g_tvgd"] = g_tvgd(g_bar, a0=A0_TVGD)
df["g_rar_fixed"] = g_rar_empirical(g_bar, gdag=GDAG_RAR)
df["g_mond_simple"] = g_mond_simple(g_bar, a0=A0_TVGD)
df["g_mond_standard"] = g_mond_standard(g_bar, a0=A0_TVGD)

# Optimize empirical RAR gdag
grid_rows = []

for gdag in GDAG_GRID:
    g_model = g_rar_empirical(g_bar, gdag=gdag)
    m = metrics(g_obs, g_model)

    grid_rows.append({
        "gdag": gdag,
        "mean_abs_log": m["mean_abs_log"],
        "rmse_log": m["rmse_log"],
        "bias_log": m["bias_log"],
    })

gdag_grid_df = pd.DataFrame(grid_rows)

best_mae_row = gdag_grid_df.sort_values("mean_abs_log").iloc[0]
best_rmse_row = gdag_grid_df.sort_values("rmse_log").iloc[0]

BEST_GDAG_MAE = float(best_mae_row["gdag"])
BEST_GDAG_RMSE = float(best_rmse_row["gdag"])

df["g_rar_best_mae"] = g_rar_empirical(g_bar, gdag=BEST_GDAG_MAE)
df["g_rar_best_rmse"] = g_rar_empirical(g_bar, gdag=BEST_GDAG_RMSE)


# ============================================================
# Metrics
# ============================================================

models = {
    "baryonic_only": "g_baryonic",
    "TVGD_fixed_a0": "g_tvgd",
    "RAR_fixed_gdag": "g_rar_fixed",
    "RAR_best_MAE_gdag": "g_rar_best_mae",
    "RAR_best_RMSE_gdag": "g_rar_best_rmse",
    "MOND_simple": "g_mond_simple",
    "MOND_standard": "g_mond_standard",
}

metric_rows = []

for name, col in models.items():
    m = metrics(g_obs, df[col].values)

    metric_rows.append({
        "model": name,
        "N": m["N"],
        "mean_abs_log": m["mean_abs_log"],
        "median_abs_log": m["median_abs_log"],
        "rmse_log": m["rmse_log"],
        "bias_log": m["bias_log"],
    })

metrics_df = pd.DataFrame(metric_rows)

bar_mae = metrics_df.loc[
    metrics_df["model"] == "baryonic_only",
    "mean_abs_log",
].iloc[0]

bar_rmse = metrics_df.loc[
    metrics_df["model"] == "baryonic_only",
    "rmse_log",
].iloc[0]

metrics_df["improvement_mae_vs_baryonic_percent"] = 100.0 * (
    bar_mae - metrics_df["mean_abs_log"]
) / bar_mae

metrics_df["improvement_rmse_vs_baryonic_percent"] = 100.0 * (
    bar_rmse - metrics_df["rmse_log"]
) / bar_rmse

max_abs_tvgd_rar_fixed = float(
    np.max(np.abs(df["g_tvgd"].values - df["g_rar_fixed"].values))
)

rel_tvgd_rar_fixed = float(
    np.max(
        np.abs(df["g_tvgd"].values - df["g_rar_fixed"].values)
        / np.maximum(df["g_rar_fixed"].values, 1e-300)
    )
)

tvgd_mae = metrics_df.loc[
    metrics_df["model"] == "TVGD_fixed_a0",
    "mean_abs_log",
].iloc[0]

tvgd_rmse = metrics_df.loc[
    metrics_df["model"] == "TVGD_fixed_a0",
    "rmse_log",
].iloc[0]

rar_best_mae = metrics_df.loc[
    metrics_df["model"] == "RAR_best_MAE_gdag",
    "mean_abs_log",
].iloc[0]

rar_best_rmse = metrics_df.loc[
    metrics_df["model"] == "RAR_best_RMSE_gdag",
    "rmse_log",
].iloc[0]

tvgd_vs_best_mae_percent = 100.0 * (tvgd_mae - rar_best_mae) / rar_best_mae
tvgd_vs_best_rmse_percent = 100.0 * (tvgd_rmse - rar_best_rmse) / rar_best_rmse


# ============================================================
# Acceleration regimes
# ============================================================

df["u"] = np.sqrt(df["g_bar"].values / A0_TVGD)

regimes = {
    "low_u_u<1": df[df["u"] < 1],
    "transition_1<=u<3": df[(df["u"] >= 1) & (df["u"] < 3)],
    "high_u_u>=3": df[df["u"] >= 3],
}

regime_rows = []

for regime_name, sub in regimes.items():
    if len(sub) < 10:
        continue

    for model_name, col in models.items():
        m = metrics(sub["g_obs"].values, sub[col].values)

        regime_rows.append({
            "regime": regime_name,
            "model": model_name,
            "N": len(sub),
            "galaxies": sub["galaxy"].nunique(),
            "mean_abs_log": m["mean_abs_log"],
            "rmse_log": m["rmse_log"],
            "bias_log": m["bias_log"],
        })

regime_df = pd.DataFrame(regime_rows)


# ============================================================
# Automatic diagnostic
# ============================================================

tvgd_improvement_mae = metrics_df.loc[
    metrics_df["model"] == "TVGD_fixed_a0",
    "improvement_mae_vs_baryonic_percent",
].iloc[0]

tvgd_improvement_rmse = metrics_df.loc[
    metrics_df["model"] == "TVGD_fixed_a0",
    "improvement_rmse_vs_baryonic_percent",
].iloc[0]

equivalence_ok = rel_tvgd_rar_fixed < 1e-12

near_best_rar_ok = (
    abs(tvgd_vs_best_mae_percent) < 10.0
    and abs(tvgd_vs_best_rmse_percent) < 10.0
)

strong_baryonic_improvement_ok = (
    tvgd_improvement_mae > 50.0
    and tvgd_improvement_rmse > 50.0
)

if strong_baryonic_improvement_ok and equivalence_ok and near_best_rar_ok:
    status = "PASSOU_FORTE_COMO_REALIZACAO_TVGD_DA_RAR"
elif strong_baryonic_improvement_ok and near_best_rar_ok:
    status = "PASSOU_MODERADO_COMO_CONTROLE_RAR"
elif strong_baryonic_improvement_ok:
    status = "PASSOU_FRACO_CONTRA_BARIO"
else:
    status = "TENSIONADO_CONTROLE_RAR"

summary = {
    "N_points": int(len(df)),
    "N_galaxies": int(df["galaxy"].nunique()),
    "UPS_DISK": UPS_DISK,
    "UPS_BULGE": UPS_BULGE,
    "A0_TVGD": A0_TVGD,
    "GDAG_RAR_FIXED": GDAG_RAR,
    "BEST_GDAG_MAE": BEST_GDAG_MAE,
    "BEST_GDAG_RMSE": BEST_GDAG_RMSE,
    "TVGD_MAE_log": float(tvgd_mae),
    "TVGD_RMSE_log": float(tvgd_rmse),
    "RAR_best_MAE_log": float(rar_best_mae),
    "RAR_best_RMSE_log": float(rar_best_rmse),
    "TVGD_vs_best_RAR_MAE_percent": float(tvgd_vs_best_mae_percent),
    "TVGD_vs_best_RAR_RMSE_percent": float(tvgd_vs_best_rmse_percent),
    "TVGD_improvement_MAE_vs_baryonic_percent": float(tvgd_improvement_mae),
    "TVGD_improvement_RMSE_vs_baryonic_percent": float(tvgd_improvement_rmse),
    "max_abs_TVGD_minus_RAR_fixed": max_abs_tvgd_rar_fixed,
    "max_relative_TVGD_minus_RAR_fixed": rel_tvgd_rar_fixed,
    "equivalence_ok": bool(equivalence_ok),
    "near_best_rar_ok": bool(near_best_rar_ok),
    "strong_baryonic_improvement_ok": bool(strong_baryonic_improvement_ok),
    "status": status,
}


# ============================================================
# Print results
# ============================================================

print()
print("Global model metrics:")
print(metrics_df.to_string(index=False))

print()
print("Empirical RAR optimization:")
print("Best gdag by MAE:", BEST_GDAG_MAE)
print("Best gdag by RMSE:", BEST_GDAG_RMSE)

print()
print("TVGD vs fixed RAR equivalence:")
print("max_abs_TVGD_minus_RAR_fixed:", max_abs_tvgd_rar_fixed)
print("max_relative_TVGD_minus_RAR_fixed:", rel_tvgd_rar_fixed)

print()
print("TVGD vs optimized RAR:")
print("TVGD_vs_best_RAR_MAE_percent:", tvgd_vs_best_mae_percent)
print("TVGD_vs_best_RAR_RMSE_percent:", tvgd_vs_best_rmse_percent)

print()
print("Automatic diagnostic:")
for k, v in summary.items():
    print(f"{k}: {v}")

print()
print("Metrics by acceleration regime:")
print(regime_df.to_string(index=False))


# ============================================================
# Save results
# ============================================================

df.to_csv(os.path.join(OUTDIR, "test01C_sparc_rar_points.csv"), index=False)
metrics_df.to_csv(os.path.join(OUTDIR, "test01C_model_metrics.csv"), index=False)
regime_df.to_csv(os.path.join(OUTDIR, "test01C_metrics_by_regime.csv"), index=False)
gdag_grid_df.to_csv(os.path.join(OUTDIR, "test01C_gdag_grid.csv"), index=False)
pd.DataFrame([summary]).to_csv(os.path.join(OUTDIR, "test01C_summary.csv"), index=False)


# ============================================================
# Figures
# ============================================================

plt.figure(figsize=(7, 6))
plt.scatter(df["g_bar"], df["g_obs"], s=5, alpha=0.25, label="SPARC points")

x = np.logspace(
    np.log10(df["g_bar"].min()),
    np.log10(df["g_bar"].max()),
    500,
)

plt.plot(x, x, "--", label="Baryonic/Newtonian")
plt.plot(x, g_tvgd(x, A0_TVGD), label="TVGD fixed a0")
plt.plot(x, g_mond_simple(x, A0_TVGD), label="MOND simple", alpha=0.8)
plt.plot(x, g_mond_standard(x, A0_TVGD), label="MOND standard", alpha=0.8)

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$g_{\rm bar}$ [m/s$^2$]")
plt.ylabel(r"$g_{\rm obs}$ or model [m/s$^2$]")
plt.title("Test 01C — SPARC/RAR: TVGD vs RAR/MOND controls")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01C_rar_models.png"), dpi=220)
plt.show()

plot_df = metrics_df.sort_values("mean_abs_log")

plt.figure(figsize=(10, 5))
plt.bar(plot_df["model"], plot_df["mean_abs_log"])
plt.xticks(rotation=30, ha="right")
plt.ylabel("Mean absolute log residual")
plt.title("Model comparison — MAE log")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01C_model_mae_comparison.png"), dpi=220)
plt.show()

plot_df2 = metrics_df.sort_values("rmse_log")

plt.figure(figsize=(10, 5))
plt.bar(plot_df2["model"], plot_df2["rmse_log"])
plt.xticks(rotation=30, ha="right")
plt.ylabel("RMSE log")
plt.title("Model comparison — RMSE log")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01C_model_rmse_comparison.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(gdag_grid_df["gdag"], gdag_grid_df["mean_abs_log"], label="MAE log")
plt.plot(gdag_grid_df["gdag"], gdag_grid_df["rmse_log"], label="RMSE log")
plt.axvline(A0_TVGD, linestyle="--", label="TVGD a0")
plt.axvline(BEST_GDAG_MAE, linestyle=":", label="Best gdag MAE")
plt.axvline(BEST_GDAG_RMSE, linestyle="-.", label="Best gdag RMSE")
plt.xlabel(r"$g_\dagger$ [m/s$^2$]")
plt.ylabel("Residual metric")
plt.title("Empirical RAR gdag scan")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01C_gdag_scan.png"), dpi=220)
plt.show()

print()
print("Files saved in:")
print(OUTDIR)

print()
print("Compact summary:")
print(pd.DataFrame([summary]).to_string(index=False))

print()
print("DONE.")
