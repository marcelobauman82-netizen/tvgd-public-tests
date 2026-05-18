# ============================================================
# TVGD PUBLIC TESTS
# TEST 01B — SPARC/RAR MASS-TO-LIGHT ROBUSTNESS
# Public reproducibility script
# ============================================================

"""
Test 01B — SPARC/RAR Mass-to-Light Robustness

This script tests whether the SPARC/RAR performance of the effective
TVGD acceleration law depends on a single stellar mass-to-light choice.

The tested TVGD response is:

    P(u) = 1 - exp(-u)

with

    u = sqrt(g_bar / a0)

and

    g_TVGD = g_bar / P(u).

The script scans a grid of stellar mass-to-light ratios:

    Upsilon_disk  = 0.30 ... 0.80
    Upsilon_bulge = 0.50 ... 1.00

and compares TVGD against the baryonic-only model for each grid point.

This is a robustness test, not a full Bayesian model-selection analysis.

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

OUTDIR = "Resultados/test01B_sparc_ml_robustness"
DATADIR = os.path.join(OUTDIR, "data")

os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

SPARC_URLS = [
    "http://astroweb.cwru.edu/SPARC/Rotmod_LTG.zip",
    "https://zenodo.org/records/16284118/files/Rotmod_LTG.zip?download=1",
]

ZIP_PATH = os.path.join(DATADIR, "Rotmod_LTG.zip")

# Constants
A0 = 1.2e-10
KPC_TO_M = 3.0856775814913673e19
KM_TO_M = 1000.0

# Mass-to-light grid
UPS_DISK_GRID = np.round(np.arange(0.3, 0.81, 0.05), 2)
UPS_BULGE_GRID = np.round(np.arange(0.5, 1.01, 0.10), 2)

MIN_POINTS_GALAXY = 3

# Preliminary robustness criteria
STRONG_IMPROVEMENT_MAE = 30.0
STRONG_IMPROVEMENT_RMSE = 25.0


# ============================================================
# TVGD functions
# ============================================================

def P_tvgd(g_bar, a0=A0):
    """
    TVGD structural response function.

    P(u) = 1 - exp(-u)
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
print("TEST 01B — SPARC/RAR MASS-TO-LIGHT ROBUSTNESS")
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
# Build base table independent of M/L
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

base = raw.copy()
base = base.replace([np.inf, -np.inf], np.nan)

base = base.dropna(
    subset=[
        "R_kpc",
        "Vobs_km_s",
        "Vgas_km_s",
        "Vdisk_km_s",
        "Vbul_km_s",
    ]
)

base = base[
    (base["R_kpc"] > 0)
    & (base["Vobs_km_s"] > 0)
].copy()

R_m = base["R_kpc"].values * KPC_TO_M

base["g_obs"] = (base["Vobs_km_s"].values * KM_TO_M) ** 2 / R_m

Vgas_m_s = base["Vgas_km_s"].values * KM_TO_M
Vdisk_m_s = base["Vdisk_km_s"].values * KM_TO_M
Vbul_m_s = base["Vbul_km_s"].values * KM_TO_M

base["g_gas_unit"] = Vgas_m_s * np.abs(Vgas_m_s) / R_m
base["g_disk_unit"] = Vdisk_m_s * np.abs(Vdisk_m_s) / R_m
base["g_bul_unit"] = Vbul_m_s * np.abs(Vbul_m_s) / R_m

base = base[
    np.isfinite(base["g_obs"])
    & (base["g_obs"] > 0)
].copy()

counts = base.groupby("galaxy").size()
good_galaxies = counts[counts >= MIN_POINTS_GALAXY].index
base = base[base["galaxy"].isin(good_galaxies)].copy()

print("Clean base points:", len(base))
print("Clean base galaxies:", base["galaxy"].nunique())


# ============================================================
# Run M/L grid
# ============================================================

results = []

for ups_disk in UPS_DISK_GRID:
    for ups_bulge in UPS_BULGE_GRID:

        g_bar = (
            base["g_gas_unit"].values
            + ups_disk * base["g_disk_unit"].values
            + ups_bulge * base["g_bul_unit"].values
        )

        mask = (
            np.isfinite(g_bar)
            & (g_bar > 0)
            & np.isfinite(base["g_obs"].values)
            & (base["g_obs"].values > 0)
        )

        sub_g_obs = base["g_obs"].values[mask]
        sub_g_bar = g_bar[mask]

        if len(sub_g_obs) < 100:
            continue

        sub_g_tvgd = g_tvgd(sub_g_bar, a0=A0)

        m_bar = metrics(sub_g_obs, sub_g_bar)
        m_tvgd = metrics(sub_g_obs, sub_g_tvgd)

        improvement_mae = 100.0 * (
            m_bar["mean_abs_log"] - m_tvgd["mean_abs_log"]
        ) / m_bar["mean_abs_log"]

        improvement_rmse = 100.0 * (
            m_bar["rmse_log"] - m_tvgd["rmse_log"]
        ) / m_bar["rmse_log"]

        u = np.sqrt(sub_g_bar / A0)

        low = u < 1
        trans = (u >= 1) & (u < 3)
        high = u >= 3

        def regime_improvement(mask_regime):
            if np.sum(mask_regime) < 10:
                return np.nan, np.nan, int(np.sum(mask_regime))

            mb = metrics(sub_g_obs[mask_regime], sub_g_bar[mask_regime])
            mt = metrics(sub_g_obs[mask_regime], sub_g_tvgd[mask_regime])

            imp_mae = 100.0 * (
                mb["mean_abs_log"] - mt["mean_abs_log"]
            ) / mb["mean_abs_log"]

            imp_rmse = 100.0 * (
                mb["rmse_log"] - mt["rmse_log"]
            ) / mb["rmse_log"]

            return float(imp_mae), float(imp_rmse), int(np.sum(mask_regime))

        low_mae, low_rmse, low_N = regime_improvement(low)
        trans_mae, trans_rmse, trans_N = regime_improvement(trans)
        high_mae, high_rmse, high_N = regime_improvement(high)

        results.append({
            "UPS_DISK": ups_disk,
            "UPS_BULGE": ups_bulge,
            "N_points": int(len(sub_g_obs)),
            "N_galaxies": int(base.loc[mask, "galaxy"].nunique()),
            "baryonic_mean_abs_log": m_bar["mean_abs_log"],
            "tvgd_mean_abs_log": m_tvgd["mean_abs_log"],
            "baryonic_rmse_log": m_bar["rmse_log"],
            "tvgd_rmse_log": m_tvgd["rmse_log"],
            "improvement_mae_percent": float(improvement_mae),
            "improvement_rmse_percent": float(improvement_rmse),
            "bias_tvgd_log": m_tvgd["bias_log"],
            "low_N": low_N,
            "low_improvement_mae_percent": low_mae,
            "low_improvement_rmse_percent": low_rmse,
            "transition_N": trans_N,
            "transition_improvement_mae_percent": trans_mae,
            "transition_improvement_rmse_percent": trans_rmse,
            "high_N": high_N,
            "high_improvement_mae_percent": high_mae,
            "high_improvement_rmse_percent": high_rmse,
        })

grid = pd.DataFrame(results)

if len(grid) == 0:
    raise RuntimeError("The M/L grid produced no valid results.")


# ============================================================
# Robustness diagnostic
# ============================================================

grid["pass_global_strong"] = (
    (grid["improvement_mae_percent"] > STRONG_IMPROVEMENT_MAE)
    & (grid["improvement_rmse_percent"] > STRONG_IMPROVEMENT_RMSE)
)

grid["pass_low_strong"] = grid["low_improvement_mae_percent"] > 40.0
grid["pass_transition_positive"] = grid["transition_improvement_mae_percent"] > 0.0

grid["pass_overall"] = (
    grid["pass_global_strong"]
    & grid["pass_low_strong"]
    & grid["pass_transition_positive"]
)

best_mae = grid.sort_values("tvgd_mean_abs_log").iloc[0]
best_rmse = grid.sort_values("tvgd_rmse_log").iloc[0]

frac_pass_global = grid["pass_global_strong"].mean()
frac_pass_overall = grid["pass_overall"].mean()

summary = {
    "N_grid": int(len(grid)),
    "UPS_DISK_min": float(grid["UPS_DISK"].min()),
    "UPS_DISK_max": float(grid["UPS_DISK"].max()),
    "UPS_BULGE_min": float(grid["UPS_BULGE"].min()),
    "UPS_BULGE_max": float(grid["UPS_BULGE"].max()),
    "frac_pass_global_strong": float(frac_pass_global),
    "frac_pass_overall": float(frac_pass_overall),
    "median_improvement_mae_percent": float(grid["improvement_mae_percent"].median()),
    "median_improvement_rmse_percent": float(grid["improvement_rmse_percent"].median()),
    "min_improvement_mae_percent": float(grid["improvement_mae_percent"].min()),
    "min_improvement_rmse_percent": float(grid["improvement_rmse_percent"].min()),
    "max_improvement_mae_percent": float(grid["improvement_mae_percent"].max()),
    "max_improvement_rmse_percent": float(grid["improvement_rmse_percent"].max()),
    "best_mae_UPS_DISK": float(best_mae["UPS_DISK"]),
    "best_mae_UPS_BULGE": float(best_mae["UPS_BULGE"]),
    "best_mae_tvgd_mean_abs_log": float(best_mae["tvgd_mean_abs_log"]),
    "best_rmse_UPS_DISK": float(best_rmse["UPS_DISK"]),
    "best_rmse_UPS_BULGE": float(best_rmse["UPS_BULGE"]),
    "best_rmse_tvgd_rmse_log": float(best_rmse["tvgd_rmse_log"]),
}

if frac_pass_overall > 0.75:
    status = "PASSOU_FORTE_ROBUSTEZ_ML"
elif frac_pass_global > 0.75:
    status = "PASSOU_MODERADO_ROBUSTEZ_ML"
elif grid["improvement_mae_percent"].median() > 0:
    status = "PASSOU_FRACO_ROBUSTEZ_ML"
else:
    status = "TENSIONADO_ROBUSTEZ_ML"

summary["status"] = status


# ============================================================
# Print results
# ============================================================

print()
print("Robustness summary:")
for k, v in summary.items():
    print(f"{k}: {v}")

print()
print("Best cases:")
print("Best MAE:")
print(best_mae[[
    "UPS_DISK",
    "UPS_BULGE",
    "tvgd_mean_abs_log",
    "improvement_mae_percent",
    "improvement_rmse_percent",
    "low_improvement_mae_percent",
    "transition_improvement_mae_percent",
]].to_string())

print()
print("Best RMSE:")
print(best_rmse[[
    "UPS_DISK",
    "UPS_BULGE",
    "tvgd_rmse_log",
    "improvement_mae_percent",
    "improvement_rmse_percent",
    "low_improvement_mae_percent",
    "transition_improvement_mae_percent",
]].to_string())


# ============================================================
# Save results
# ============================================================

grid_path = os.path.join(OUTDIR, "test01B_sparc_ml_robustness_grid.csv")
summary_path = os.path.join(OUTDIR, "test01B_sparc_ml_robustness_summary.csv")

grid.to_csv(grid_path, index=False)
pd.DataFrame([summary]).to_csv(summary_path, index=False)


# ============================================================
# Figures
# ============================================================

pivot_mae = grid.pivot(
    index="UPS_BULGE",
    columns="UPS_DISK",
    values="improvement_mae_percent",
)

plt.figure(figsize=(9, 6))
plt.imshow(
    pivot_mae.values,
    origin="lower",
    aspect="auto",
    extent=[
        grid["UPS_DISK"].min(),
        grid["UPS_DISK"].max(),
        grid["UPS_BULGE"].min(),
        grid["UPS_BULGE"].max(),
    ],
)
plt.colorbar(label="MAE log improvement [%]")
plt.xlabel(r"$\Upsilon_{\rm disk}$")
plt.ylabel(r"$\Upsilon_{\rm bulge}$")
plt.title("Test 01B — SPARC/RAR TVGD robustness in M/L")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01B_heatmap_mae_improvement.png"), dpi=220)
plt.show()

pivot_rmse = grid.pivot(
    index="UPS_BULGE",
    columns="UPS_DISK",
    values="improvement_rmse_percent",
)

plt.figure(figsize=(9, 6))
plt.imshow(
    pivot_rmse.values,
    origin="lower",
    aspect="auto",
    extent=[
        grid["UPS_DISK"].min(),
        grid["UPS_DISK"].max(),
        grid["UPS_BULGE"].min(),
        grid["UPS_BULGE"].max(),
    ],
)
plt.colorbar(label="RMSE log improvement [%]")
plt.xlabel(r"$\Upsilon_{\rm disk}$")
plt.ylabel(r"$\Upsilon_{\rm bulge}$")
plt.title("RMSE improvement across M/L grid")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01B_heatmap_rmse_improvement.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(grid["improvement_mae_percent"], bins=20, alpha=0.7, label="MAE")
plt.hist(grid["improvement_rmse_percent"], bins=20, alpha=0.7, label="RMSE")
plt.axvline(0, color="black", linewidth=1)
plt.xlabel("Improvement over baryonic model [%]")
plt.ylabel("N grid points")
plt.title("Distribution of TVGD improvement across M/L grid")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test01B_improvement_distribution.png"), dpi=220)
plt.show()

print()
print("Files saved in:")
print(OUTDIR)
print(grid_path)
print(summary_path)

print()
print("Compact summary:")
print(pd.DataFrame([summary]).to_string(index=False))

print()
print("DONE.")
