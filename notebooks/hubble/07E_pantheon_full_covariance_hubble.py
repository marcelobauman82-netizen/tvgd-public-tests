
# ============================================================
# TVGD PUBLIC TESTS
# TEST 07E — PANTHEON+SH0ES FULL COVARIANCE HUBBLE TEST
# ============================================================
#
# Goal:
# Test the TVGD Hubble structural transition using Pantheon+SH0ES
# with the full statistical+systematic covariance matrix.
#
# This is stronger than diagonal-error tests because it uses:
#
#     chi2 = r^T C^{-1} r
#
# Data source:
# PantheonPlusSH0ES/DataRelease
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib.request

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test07E_pantheon_full_covariance"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Public data URLs
# ============================================================

DATA_URLS = [
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat",
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat",
]

COV_URLS = [
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov",
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES_STAT+SYS.cov",
]

DATA_PATH = os.path.join(OUTDIR, "PantheonSH0ES.dat")
COV_PATH = os.path.join(OUTDIR, "PantheonSH0ES_STAT_SYS.cov")

# ============================================================
# Reproducibility
# ============================================================

SEED = 20260519
rng = np.random.default_rng(SEED)

# ============================================================
# Reference cosmology and TVGD parameters
# ============================================================

C_KM_S = 299792.458

H0_CMB = 67.4
H0_LOCAL = 73.04
OMEGA_M = 0.315

# TVGD Hubble structural transition
A_H = H0_LOCAL / H0_CMB - 1.0
S_NATURAL = 1.5

# Fixed value found in 07B
THETA0_07B = 28.311753

# Theta0 scan
THETA0_GRID = np.linspace(1.0, 80.0, 500)

# Main redshift cut
Z_MIN = 0.01
Z_MAX = 0.20

# BAO/CMB safety
Z_BAO_MIN = 0.106
Z_CMB = 1089.92
MAX_BAO_SHIFT_PERCENT = 0.5

# Placebo null
N_NULL = 300

# Automatic status thresholds
MIN_DELTA_CHI2_STRONG = 10.0
MIN_DELTA_CHI2_MODERATE = 2.0
MAX_P_NULL_STRONG = 0.05

# ============================================================
# Download helpers
# ============================================================

print("=" * 72)
print("TEST 07E — PANTHEON+SH0ES FULL COVARIANCE HUBBLE TEST")
print("=" * 72)

def download_first_available(urls, path):
    if os.path.exists(path):
        print("File already exists:", path)
        return

    last_err = None

    for url in urls:
        try:
            print("\nTrying download:")
            print(url)
            urllib.request.urlretrieve(url, path)
            print("Download completed:", path)
            return
        except Exception as e:
            last_err = e
            print("Failed:", e)

    raise RuntimeError(f"Could not download {path}. Last error: {last_err}")

download_first_available(DATA_URLS, DATA_PATH)
download_first_available(COV_URLS, COV_PATH)

# ============================================================
# Read data
# ============================================================

df_raw = pd.read_csv(DATA_PATH, sep=r"\s+", comment="#")
df_raw = df_raw.reset_index(drop=True)
df_raw["original_index"] = np.arange(len(df_raw))

required_cols = ["zHD", "MU_SH0ES", "MU_SH0ES_ERR_DIAG"]

for col in required_cols:
    if col not in df_raw.columns:
        raise ValueError(f"Required column missing: {col}")

df_all = df_raw[
    np.isfinite(df_raw["zHD"])
    & np.isfinite(df_raw["MU_SH0ES"])
    & np.isfinite(df_raw["MU_SH0ES_ERR_DIAG"])
    & (df_raw["MU_SH0ES_ERR_DIAG"] > 0)
].copy()

df = df_all[(df_all["zHD"] >= Z_MIN) & (df_all["zHD"] <= Z_MAX)].copy()
df = df.sort_values("zHD").reset_index(drop=True)

print("\nData:")
print("N total dat =", len(df_raw))
print("N after cuts =", len(df))
print("z_min =", df["zHD"].min())
print("z_max =", df["zHD"].max())

if len(df) < 20:
    raise RuntimeError("Sample too small after cuts.")

# ============================================================
# Read full covariance
# ============================================================

def read_pantheon_covariance(path):
    """
    Pantheon .cov files usually have:
        first entry = N
        then N*N flattened covariance entries.

    This function also accepts a raw flattened matrix without the first N.
    """
    vals = np.loadtxt(path).astype(float)

    n0 = int(round(vals[0]))
    rest = vals[1:]

    if len(rest) == n0 * n0:
        return rest.reshape((n0, n0))

    n = int(round(np.sqrt(len(vals))))
    if n * n == len(vals):
        return vals.reshape((n, n))

    raise ValueError("Unrecognized covariance format.")

C_full = read_pantheon_covariance(COV_PATH)

print("\nCovariance:")
print("shape =", C_full.shape)

if C_full.shape[0] != len(df_raw):
    raise RuntimeError(
        f"Covariance size ({C_full.shape[0]}) differs from dat size ({len(df_raw)})."
    )

idx = df["original_index"].values.astype(int)
C = C_full[np.ix_(idx, idx)]

# Numerical symmetrization
C = 0.5 * (C + C.T)

eigvals = np.linalg.eigvalsh(C)
eig_min = float(np.min(eigvals))
eig_max = float(np.max(eigvals))
cond_number = float(eig_max / max(eig_min, 1e-300))

print("eig_min =", eig_min)
print("eig_max =", eig_max)
print("condition_number =", cond_number)

jitter_added = 0.0

if eig_min <= 0:
    jitter_added = abs(eig_min) + 1e-12 * eig_max
    C = C + np.eye(C.shape[0]) * jitter_added
    print("Jitter added:", jitter_added)

Cinv = np.linalg.pinv(C, rcond=1e-12)

# ============================================================
# Cosmological distance functions
# ============================================================

def E_lcdm(z, omega_m=OMEGA_M):
    z = np.asarray(z, dtype=float)
    return np.sqrt(omega_m * (1.0 + z)**3 + (1.0 - omega_m))


def comoving_distance_mpc_single(z, H0=H0_CMB, omega_m=OMEGA_M, n_grid=2500):
    z = float(z)

    if z <= 0:
        return 0.0

    zz = np.linspace(0.0, z, n_grid)
    yy = 1.0 / E_lcdm(zz, omega_m=omega_m)
    integral = np.trapezoid(yy, zz)

    return (C_KM_S / H0) * integral


def mu_lcdm(z_array, H0=H0_CMB, omega_m=OMEGA_M):
    z_array = np.asarray(z_array, dtype=float)
    out = []

    for zi in z_array:
        DM = comoving_distance_mpc_single(zi, H0=H0, omega_m=omega_m)
        DL = (1.0 + zi) * DM
        out.append(5.0 * np.log10(max(DL, 1e-300)) + 25.0)

    return np.array(out)


# ============================================================
# TVGD structural Hubble transition
# ============================================================

def delta_theta_z(z, theta0, s=S_NATURAL):
    z = np.asarray(z, dtype=float)
    return theta0 * ((1.0 + np.maximum(z, 0.0))**s - 1.0)


def Xi_TVGD(z, theta0, s=S_NATURAL):
    dth = delta_theta_z(z, theta0=theta0, s=s)
    return 1.0 + A_H * np.exp(-np.maximum(dth, 0.0))


def delta_mu_tvgd(z, theta0):
    """
    If the effective inferred H0 increases by Xi(z),
    luminosity distance decreases by approximately 1/Xi(z).

    Therefore:

        delta_mu = -5 log10[Xi(z)]
    """
    xi = Xi_TVGD(z, theta0=theta0)
    return -5.0 * np.log10(xi)


def mu_tvgd_transition(z, theta0):
    return mu_lcdm(z, H0=H0_CMB) + delta_mu_tvgd(z, theta0)


def shift_percent(z, theta0):
    return 100.0 * (Xi_TVGD(z, theta0=theta0) - 1.0)


def bao_shift_max_percent(theta0):
    z_grid = np.linspace(Z_BAO_MIN, 2.0, 800)
    return float(np.max(np.abs(shift_percent(z_grid, theta0))))


def cmb_shift_abs(theta0):
    return float(abs(Xi_TVGD(Z_CMB, theta0=theta0) - 1.0))


# ============================================================
# Full covariance statistics
# ============================================================

def chi2_full(mu_obs, mu_model):
    r = mu_obs - mu_model
    return float(r @ Cinv @ r)


def best_offset_full(mu_obs, mu_model):
    """
    Best global offset with full covariance.
    """
    one = np.ones_like(mu_obs)
    r = mu_obs - mu_model

    num = one @ Cinv @ r
    den = one @ Cinv @ one

    return float(num / den)


def chi2_full_with_offset(mu_obs, mu_model):
    off = best_offset_full(mu_obs, mu_model)
    return chi2_full(mu_obs, mu_model + off), off


def aic(chi2, k):
    return chi2 + 2.0 * k


def bic(chi2, k, n):
    return chi2 + k * np.log(n)


# ============================================================
# Model evaluation
# ============================================================

z = df["zHD"].values
mu_obs = df["MU_SH0ES"].values

mu_cmb = mu_lcdm(z, H0=H0_CMB)
mu_local = mu_lcdm(z, H0=H0_LOCAL)
mu_07b = mu_tvgd_transition(z, THETA0_07B)

chi2_cmb = chi2_full(mu_obs, mu_cmb)
chi2_local = chi2_full(mu_obs, mu_local)
chi2_07b = chi2_full(mu_obs, mu_07b)

chi2_cmb_off, off_cmb = chi2_full_with_offset(mu_obs, mu_cmb)
chi2_local_off, off_local = chi2_full_with_offset(mu_obs, mu_local)
chi2_07b_off, off_07b = chi2_full_with_offset(mu_obs, mu_07b)

# ============================================================
# Theta0 scan
# ============================================================

scan_rows = []

for theta0 in THETA0_GRID:
    mu_t = mu_tvgd_transition(z, theta0)

    c2 = chi2_full(mu_obs, mu_t)
    c2_off, off = chi2_full_with_offset(mu_obs, mu_t)

    bao = bao_shift_max_percent(theta0)
    cmb = cmb_shift_abs(theta0)

    scan_rows.append({
        "theta0": theta0,
        "chi2_full": c2,
        "chi2_full_offset": c2_off,
        "offset_best": off,
        "bao_shift_max_percent": bao,
        "cmb_shift_abs": cmb,
        "shift_z0p01_percent": float(shift_percent(0.01, theta0)),
        "shift_z0p03_percent": float(shift_percent(0.03, theta0)),
        "shift_z0p05_percent": float(shift_percent(0.05, theta0)),
        "shift_z0p10_percent": float(shift_percent(0.10, theta0)),
    })

scan_df = pd.DataFrame(scan_rows)
scan_ok = scan_df[scan_df["bao_shift_max_percent"] < MAX_BAO_SHIFT_PERCENT].copy()

if len(scan_ok) > 0:
    best = scan_ok.sort_values("chi2_full").iloc[0]
    best_off = scan_ok.sort_values("chi2_full_offset").iloc[0]
else:
    best = scan_df.sort_values("chi2_full").iloc[0]
    best_off = scan_df.sort_values("chi2_full_offset").iloc[0]

theta_best = float(best["theta0"])
theta_best_off = float(best_off["theta0"])

mu_best = mu_tvgd_transition(z, theta_best)
mu_best_off_model = mu_tvgd_transition(z, theta_best_off)

chi2_best = chi2_full(mu_obs, mu_best)
chi2_best_off, off_best = chi2_full_with_offset(mu_obs, mu_best_off_model)

delta_chi2_best_vs_cmb = chi2_cmb - chi2_best
delta_chi2_07b_vs_cmb = chi2_cmb - chi2_07b
delta_chi2_best_vs_local = chi2_local - chi2_best

n = len(df)

aic_cmb = aic(chi2_cmb, k=0)
aic_tvgd = aic(chi2_best, k=1)

bic_cmb = bic(chi2_cmb, k=0, n=n)
bic_tvgd = bic(chi2_best, k=1, n=n)

delta_aic_cmb_minus_tvgd = aic_cmb - aic_tvgd
delta_bic_cmb_minus_tvgd = bic_cmb - bic_tvgd

# ============================================================
# Placebo null with full covariance
# ============================================================

delta_mu_grid = np.array([delta_mu_tvgd(z, th) for th in THETA0_GRID])
bao_grid = np.array([bao_shift_max_percent(th) for th in THETA0_GRID])
ok_theta = bao_grid < MAX_BAO_SHIFT_PERCENT

chi2_cmb_real = chi2_cmb
delta_chi2_real = delta_chi2_best_vs_cmb

null_delta = []
null_theta = []

for j in range(N_NULL):
    perm = rng.permutation(n)

    best_null_chi2 = np.inf
    best_null_theta = np.nan

    for k, theta0 in enumerate(THETA0_GRID):
        if not ok_theta[k]:
            continue

        mu_perm = mu_cmb + delta_mu_grid[k, perm]
        c2 = chi2_full(mu_obs, mu_perm)

        if c2 < best_null_chi2:
            best_null_chi2 = c2
            best_null_theta = theta0

    null_delta.append(chi2_cmb_real - best_null_chi2)
    null_theta.append(best_null_theta)

null_delta = np.array(null_delta)
null_theta = np.array(null_theta)

p_null_positive = float((np.sum(null_delta >= delta_chi2_real) + 1) / (N_NULL + 1))
p_null_two_sided_like = float((np.sum(np.abs(null_delta) >= abs(delta_chi2_real)) + 1) / (N_NULL + 1))

null_df = pd.DataFrame({
    "null_id": np.arange(N_NULL),
    "delta_chi2_null": null_delta,
    "theta_best_null": null_theta,
})

# ============================================================
# Automatic diagnosis
# ============================================================

bao_ok = float(best["bao_shift_max_percent"]) < MAX_BAO_SHIFT_PERCENT
aic_ok = delta_aic_cmb_minus_tvgd > 0
bic_ok = delta_bic_cmb_minus_tvgd > 0
null_ok = p_null_positive <= MAX_P_NULL_STRONG

if (
    delta_chi2_best_vs_cmb > MIN_DELTA_CHI2_STRONG
    and bao_ok
    and aic_ok
    and bic_ok
    and null_ok
):
    final_status = "PASSOU_FORTE_PANTHEON_FULL_COV_TVGD_HUBBLE"
elif (
    delta_chi2_best_vs_cmb > MIN_DELTA_CHI2_MODERATE
    and bao_ok
    and aic_ok
):
    final_status = "PASSOU_MODERADO_PANTHEON_FULL_COV_TVGD_HUBBLE"
elif delta_chi2_best_vs_cmb > 0 and bao_ok:
    final_status = "PASSOU_FRACO_PANTHEON_FULL_COV_TVGD_HUBBLE"
else:
    final_status = "TENSIONADO_PANTHEON_FULL_COV_TVGD_HUBBLE"

summary = {
    "N_sne": n,
    "z_min": float(df["zHD"].min()),
    "z_max": float(df["zHD"].max()),

    "H0_CMB": H0_CMB,
    "H0_LOCAL": H0_LOCAL,
    "A_H": A_H,
    "s_natural": S_NATURAL,
    "theta0_07B": THETA0_07B,

    "cov_eig_min_original": eig_min,
    "cov_eig_max_original": eig_max,
    "cov_condition_original": cond_number,
    "cov_jitter_added": jitter_added,

    "chi2_cmb_full": chi2_cmb,
    "chi2_local_full": chi2_local,
    "chi2_tvgd_07b_full": chi2_07b,
    "chi2_tvgd_best_full": chi2_best,

    "delta_chi2_07b_vs_cmb": delta_chi2_07b_vs_cmb,
    "delta_chi2_best_vs_cmb": delta_chi2_best_vs_cmb,
    "delta_chi2_best_vs_local": delta_chi2_best_vs_local,

    "theta0_best_full": theta_best,
    "theta0_best_offset_full": theta_best_off,

    "best_bao_shift_max_percent": float(best["bao_shift_max_percent"]),
    "best_cmb_shift_abs": float(best["cmb_shift_abs"]),
    "best_shift_z0p01_percent": float(best["shift_z0p01_percent"]),
    "best_shift_z0p03_percent": float(best["shift_z0p03_percent"]),
    "best_shift_z0p05_percent": float(best["shift_z0p05_percent"]),
    "best_shift_z0p10_percent": float(best["shift_z0p10_percent"]),

    "chi2_cmb_full_offset": chi2_cmb_off,
    "chi2_local_full_offset": chi2_local_off,
    "chi2_tvgd_07b_full_offset": chi2_07b_off,
    "chi2_tvgd_best_full_offset": chi2_best_off,

    "offset_cmb": off_cmb,
    "offset_local": off_local,
    "offset_tvgd_07b": off_07b,
    "offset_tvgd_best": off_best,

    "delta_chi2_best_offset_vs_cmb_offset": chi2_cmb_off - chi2_best_off,

    "aic_cmb": aic_cmb,
    "aic_tvgd": aic_tvgd,
    "bic_cmb": bic_cmb,
    "bic_tvgd": bic_tvgd,
    "delta_aic_cmb_minus_tvgd": delta_aic_cmb_minus_tvgd,
    "delta_bic_cmb_minus_tvgd": delta_bic_cmb_minus_tvgd,

    "N_null": N_NULL,
    "null_delta_mean": float(np.mean(null_delta)),
    "null_delta_median": float(np.median(null_delta)),
    "null_delta_std": float(np.std(null_delta, ddof=1)),
    "null_delta_p95": float(np.quantile(null_delta, 0.95)),
    "p_null_positive": p_null_positive,
    "p_null_two_sided_like": p_null_two_sided_like,

    "bao_ok": bool(bao_ok),
    "aic_ok": bool(aic_ok),
    "bic_ok": bool(bic_ok),
    "null_ok": bool(null_ok),

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nFinal summary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

scan_path = os.path.join(OUTDIR, "test07E_theta0_full_cov_scan.csv")
null_path = os.path.join(OUTDIR, "test07E_full_cov_null_placebo.csv")
summary_path = os.path.join(OUTDIR, "test07E_full_cov_summary.csv")

scan_df.to_csv(scan_path, index=False)
null_df.to_csv(null_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(scan_df["theta0"], scan_df["chi2_full"], label="TVGD full-cov chi2")
plt.axhline(chi2_cmb, color="black", linestyle="--", label="CMB-like")
plt.axhline(chi2_local, color="gray", linestyle=":", label="Local H0")
plt.axvline(THETA0_07B, color="red", linestyle="--", label="theta0 07B")
plt.axvline(theta_best, color="blue", linestyle="--", label="best theta0")
plt.xlabel("theta0")
plt.ylabel("chi2 full covariance")
plt.title("Pantheon+SH0ES full covariance TVGD scan")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07E_full_cov_chi2_scan.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(scan_df["theta0"], scan_df["bao_shift_max_percent"])
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO criterion")
plt.axvline(theta_best, color="blue", linestyle="--", label="best theta0")
plt.xlabel("theta0")
plt.ylabel("max BAO shift [%]")
plt.title("BAO safety")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07E_bao_safety.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(null_delta, bins=35, alpha=0.8, label="placebo null")
plt.axvline(delta_chi2_real, color="red", linestyle="--", label="real")
plt.axvline(0, color="black", linewidth=1)
plt.xlabel("Delta chi2 = chi2_CMB - chi2_TVGD")
plt.ylabel("N")
plt.title("Full covariance placebo null")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07E_null_delta_chi2.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.scatter(z, mu_obs - mu_cmb, s=8, alpha=0.35, label="Residual vs CMB")
plt.scatter(z, mu_obs - mu_best, s=8, alpha=0.35, label="Residual vs TVGD best")
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("zHD")
plt.ylabel("mu_obs - mu_model")
plt.title("Residuals with full covariance best fit")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07E_residuals.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(scan_path)
print(null_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_PANTHEON_FULL_COV_TVGD_HUBBLE":
    print("The TVGD structural transition improves Pantheon+SH0ES with full covariance, preserves BAO, and passes the placebo control.")
elif "PASSOU_MODERADO" in final_status:
    print("The TVGD structural transition improves Pantheon+SH0ES with full covariance and preserves BAO, but one or more controls still require caution.")
elif "PASSOU_FRACO" in final_status:
    print("The TVGD structural transition improves weakly with full covariance, but the Hubble sector is not closed.")
else:
    print("The TVGD structural transition is tensioned in the full-covariance Pantheon+SH0ES test.")

print("\nDONE.")
