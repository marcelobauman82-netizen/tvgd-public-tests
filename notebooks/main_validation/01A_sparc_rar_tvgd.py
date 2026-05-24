# ============================================================
# TVGD PUBLIC TESTS
# TEST 01A — SPARC/RAR GALACTIC TEST
# Public reproducibility script
# ============================================================

"""
Test 01A — SPARC/RAR Galactic Test

This script tests the effective TVGD acceleration law against the SPARC
galaxy rotation-curve data.

The tested TVGD response is:

    P(u) = 1 - exp(-u)

with

    u = sqrt(g_bar / a0)

and

    g_TVGD = g_bar / P(u)

The script compares:

1. baryonic-only acceleration;
2. TVGD effective acceleration.

It uses fixed stellar mass-to-light ratios:

    Upsilon_disk  = 0.5
    Upsilon_bulge = 0.7

This first public version is intended as a reproducibility test, not as a
complete Bayesian model-selection analysis.

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

OUTDIR = "Resultados/test01_sparc_rar_tvgd"
DATADIR = os.path.join(OUTDIR, "data")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

SPARC_URLS = [
    "http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip",
    "https://zenodo.org/records/16284118/files/Rotmod_LTG.zip?download=1",
]

ZIP_PATH = os.path.join(DATADIR, "Rotmod_LTG.zip")

# Constants
A0 = 1.2e-10  # m/s^2
KPC_TO_M = 3.0856775814913673e19
KM_TO_M = 1000.0

# Fiducial stellar mass-to-light ratios
UPS_DISK = 0.5
UPS_BULGE = 0.7

# Minimum cuts
MIN_POINTS_GALAXY = 3
MIN_VALID_POINTS_TOTAL = 100


# ============================================================
# TVGD functions
# ============================================================

def P_tvgd(g_bar, a0=A0):
    """
    TVGD structural response function.

    P(u) = 1 - exp(-u)

    where

    u = sqrt(g_bar / a0).
    """
    g_bar = np.asarray(g_bar, dtype=float)
    u = np.sqrt(np.maximum(g_bar, 0.0) / a0)
    return 1.0 - np.exp(-u)


def g_tvgd(g_bar, a0=A0):
    """
    Effective TVGD acceleration.

    g_TVGD = g_bar / P(u)
    """
    g_bar = np.asarray(g_bar, dtype=float)
    P = P_tvgd(g_bar, a0=a0)
    return g_bar / np.maximum(P, 1e-300)


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

    Some files may contain comment lines or non-numeric entries.
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
# Download and extract SPARC
# ============================================================

print("=" * 72)
print("TEST 01A — SPARC/RAR GALACTIC TEST")
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
    raise RuntimeError("No *_rotmod.dat files found after extraction.")


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

# Unit conversion
R_m = df["R_kpc"].values * KPC_TO_M

Vobs_m_s = df["Vobs_km_s"].values * KM_TO_M
Vgas_m_s = df["Vgas_km_s"].values * KM_TO_M
Vdisk_m_s = df["Vdisk_km_s"].values * KM_TO_M
Vbul_m_s = df["Vbul_km_s"].values * KM_TO_M

# Observed acceleration
df["g_obs"] = Vobs_m_s**2 / R_m

# Baryonic acceleration
# Vgas may be negative in some SPARC entries, so V*abs(V) preserves sign.
g_gas = Vgas_m_s * np.abs(Vgas_m_s) / R_m
g_disk = UPS_DISK * Vdisk_m_s * np.abs(Vdisk_m_s) / R_m
g_bul = UPS_BULGE * Vbul_m_s * np.abs(Vbul_m_s) / R_m

df["g_bar"] = g_gas + g_disk + g_bul

# Physical cuts
df = df[
    np.isfinite(df["g_obs"])
    & np.isfinite(df["g_bar"])
    & (df["g_obs"] > 0)
    & (df["g_bar"] > 0)
].copy()

# Minimum number of points per galaxy
counts = df.groupby("galaxy").size()
good_galaxies = counts[counts >= MIN_POINTS_GALAXY].index

df = df[df["galaxy"].isin(good_galaxies)].copy()

print("Valid points:", len(df))
print("Valid galaxies:", df["galaxy"].nunique())

if len(df) < MIN_VALID_POINTS_TOTAL:
    raise RuntimeError("Too few valid points. Check parser and cuts.")


# ============================================================
# Models and residuals
# ============================================================

df["g_baryonic_model"] = df["g_bar"]
df["g_tvgd"] = g_tvgd(df["g_bar"].values, a0=A0)

df["res_baryonic_log"] = log_residual(
    df["g_obs"].values,
    df["g_baryonic_model"].values,
)

df["res_tvgd_log"] = log_residual(
    df["g_obs"].values,
    df["g_tvgd"].values,
)

df["u"] = np.sqrt(df["g_bar"].values / A0)

df["log_g_bar"] = np.log10(df["g_bar"].values)
df["log_g_obs"] = np.log10(df["g_obs"].values)
df["log_g_tvgd"] = np.log10(df["g_tvgd"].values)


# ============================================================
# Global metrics
# ============================================================

m_bar = metrics(df["g_obs"].values, df["g_baryonic_model"].values)
m_tvgd = metrics(df["g_obs"].values, df["g_tvgd"].values)

improvement_mae = 100.0 * (
    m_bar["mean_abs_log"] - m_tvgd["mean_abs_log"]
) / m_bar["mean_abs_log"]

improvement_rmse = 100.0 * (
    m_bar["rmse_log"] - m_tvgd["rmse_log"]
) / m_bar["rmse_log"]

print()
print("Global metrics:")
print("Baryonic only:", m_bar)
print("TVGD:", m_tvgd)
print(f"MAE log improvement: {improvement_mae:.6f}%")
print(f"RMSE log improvement: {improvement_rmse:.6f}%")


# ============================================================
# Metrics by acceleration regime
# ============================================================

regimes = {
    "low_u_u<1": df[df["u"] < 1],
    "transition_1<=u<3": df[(df["u"] >= 1) & (df["u"] < 3)],
    "high_u_u>=3": df[df["u"] >= 3],
}

rows_regime = []

for name, sub in regimes.items():
    if len(sub) < 10:
        continue

    mb = metrics(sub["g_obs"].values, sub["g_baryonic_model"].values)
    mt = metrics(sub["g_obs"].values, sub["g_tvgd"].values)

    imp_mae = 100.0 * (
        mb["mean_abs_log"] - mt["mean_abs_log"]
    ) / mb["mean_abs_log"]

    imp_rmse = 100.0 * (
        mb["rmse_log"] - mt["rmse_log"]
    ) / mb["rmse_log"]

    rows_regime.append({
        "regime": name,
        "N": len(sub),
        "galaxies": sub["galaxy"].nunique(),
        "baryonic_mean_abs_log": mb["mean_abs_log"],
        "tvgd_mean_abs_log": mt["mean_abs_log"],
        "baryonic_rmse_log": mb["rmse_log"],
        "tvgd_rmse_log": mt["rmse_log"],
        "improvement_mae_percent": imp_mae,
        "improvement_rmse_percent": imp_rmse,
        "median_u": float(np.median(sub["u"])),
    })

regime_df = pd.DataFrame(rows_regime)

print()
print("Metrics by acceleration regime:")
print(regime_df.to_string(index=False))


# ============================================================
# Automatic preliminary diagnostic
# ============================================================

global_ok = improvement_mae > 10.0 and improvement_rmse > 10.0

low_regime_ok = False
if "low_u_u<1" in regime_df["regime"].values:
    low_imp = regime_df.loc[
        regime_df["regime"] == "low_u_u<1",
        "improvement_mae_percent",
    ].iloc[0]
    low_regime_ok = low_imp > 20.0

transition_ok = False
if "transition_1<=u<3" in regime_df["regime"].values:
    trans_imp = regime_df.loc[
        regime_df["regime"] == "transition_1<=u<3",
        "improvement_mae_percent",
    ].iloc[0]
    transition_ok = trans_imp > 5.0

if global_ok and (low_regime_ok or transition_ok):
    status = "PASSOU_FORTE_SPARC_RAR_PRELIMINAR"
elif improvement_mae > 0 and improvement_rmse > 0:
    status = "PASSOU_MODERADO_SPARC_RAR_PRELIMINAR"
else:
    status = "TENSIONADO_SPARC_RAR_PRELIMINAR"

print()
print("Automatic diagnostic:")
print("global_ok:", global_ok)
print("low_regime_ok:", low_regime_ok)
print("transition_ok:", transition_ok)
print("status:", status)


# ============================================================
# Save results
# ============================================================

points_path = os.path.join(OUTDIR, "test01A_sparc_rar_points.csv")
regime_path = os.path.join(OUTDIR, "test01A_sparc_rar_by_regime.csv")
summary_path = os.path.join(OUTDIR, "test01A_sparc_rar_summary.csv")

df.to_csv(points_path, index=False)
regime_df.to_csv(regime_path, index=False)

summary = pd.DataFrame([{
    "N_points": len(df),
    "N_galaxies": df["galaxy"].nunique(),
    "UPS_DISK": UPS_DISK,
    "UPS_BULGE": UPS_BULGE,
    "A0": A0,
    "baryonic_mean_abs_log": m_bar["mean_abs_log"],
    "tvgd_mean_abs_log": m_tvgd["mean_abs_log"],
    "baryonic_rmse_log": m_bar["rmse_log"],
    "tvgd_rmse_log": m_tvgd["rmse_log"],
    "improvement_mae_percent": improvement_mae,
    "improvement_rmse_percent": improvement_rmse,
    "status": status,
}])

summary.to_csv(summary_path, index=False)


# ============================================================
# Figures
# ============================================================

# 1. RAR relation
plt.figure(figsize=(7, 6))
plt.scatter(
    df["g_bar"],
    df["g_obs"],
    s=5,
    alpha=0.25,
    label="SPARC points",
)

x = np.logspace(
    np.log10(df["g_bar"].min()),
    np.log10(df["g_bar"].max()),
    500,
)

plt.plot(x, x, "--", label="Baryonic/Newtonian")
plt.plot(x, g_tvgd(x), label="TVGD")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$g_{\rm bar}$ [m/s$^2$]")
plt.ylabel(r"$g_{\rm obs}$ or model [m/s$^2$]")
plt.title("Test 01A — SPARC/RAR: TVGD vs baryonic")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01A_sparc_rar_relation.png"), dpi=220)
plt.show()

# 2. Residual distributions
plt.figure(figsize=(8, 5))
plt.hist(
    df["res_baryonic_log"].dropna(),
    bins=80,
    alpha=0.5,
    label="Baryonic",
)
plt.hist(
    df["res_tvgd_log"].dropna(),
    bins=80,
    alpha=0.5,
    label="TVGD",
)
plt.axvline(0, color="black", linewidth=1)
plt.xlabel(r"$\log_{10}(g_{\rm model})-\log_{10}(g_{\rm obs})$")
plt.ylabel("N")
plt.title("Residual distribution")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01A_sparc_rar_residuals.png"), dpi=220)
plt.show()

# 3. Improvement by acceleration regime
if len(regime_df) > 0:
    plt.figure(figsize=(8, 5))
    plt.bar(
        regime_df["regime"],
        regime_df["improvement_mae_percent"],
    )
    plt.axhline(0, color="black", linewidth=1)
    plt.ylabel("MAE log improvement [%]")
    plt.title("TVGD improvement over baryonic model by acceleration regime")
    plt.xticks(rotation=20, ha="right")
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    plt.savefig(
        os.path.join(OUTDIR, "test01A_sparc_rar_improvement_by_regime.png"),
        dpi=220,
    )
    plt.show()

print()
print("Files saved in:")
print(OUTDIR)
print(points_path)
print(regime_path)
print(summary_path)

print()
print("Compact summary:")
print(summary.to_string(index=False))

print()
print("DONE.")
