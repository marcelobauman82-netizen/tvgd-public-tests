
# ============================================================
# TVGD PUBLIC TESTS
# TEST 07G — PANTHEON+ NON-CALIBRATOR SHAPE-ONLY HUBBLE TEST
# ============================================================
#
# Goal:
# Remove SH0ES calibrators and test whether the TVGD redshift-shape
# signal remains beyond a constant offset, using full covariance.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import urllib.request

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test07G_pantheon_noncalibrator_shape"
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

A_H = H0_LOCAL / H0_CMB - 1.0
S_NATURAL = 1.5

THETA0_GRID = np.linspace(1.0, 80.0, 500)

Z_MIN = 0.01
Z_MAX = 0.20

Z_BAO_MIN = 0.106
Z_CMB = 1089.92
MAX_BAO_SHIFT_PERCENT = 0.5

N_NULL = 300

MIN_SHAPE_GAIN = 2.0
MAX_P_SHAPE_NULL = 0.05

# ============================================================
# Download helpers
# ============================================================

print("=" * 72)
print("TEST 07G — PANTHEON+ NON-CALIBRATOR SHAPE-ONLY HUBBLE TEST")
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

# ============================================================
# Remove SH0ES calibrators if the column exists
# ============================================================

if "IS_CALIBRATOR" in df_all.columns:
    df_all["IS_CALIBRATOR"] = pd.to_numeric(
        df_all["IS_CALIBRATOR"], errors="coerce"
    ).fillna(0)
    df_all = df_all[df_all["IS_CALIBRATOR"] == 0].copy()
    calibrator_cut_used = True
else:
    calibrator_cut_used = False

df = df_all[(df_all["zHD"] >= Z_MIN) & (df_all["zHD"] <= Z_MAX)].copy()
df = df.sort_values("zHD").reset_index(drop=True)

print("\nData:")
print("Calibrator cut used:", calibrator_cut_used)
print("N total dat =", len(df_raw))
print("N after cuts =", len(df))
print("z_min =", df["zHD"].min())
print("z_max =", df["zHD"].max())

if len(df) < 50:
    raise RuntimeError("Sample too small after calibrator/redshift cuts.")

# ============================================================
# Read covariance
# ============================================================

def read_pantheon_covariance(path):
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

if C_full.shape[0] != len(df_raw):
    raise RuntimeError("Covariance and .dat file have different sizes.")

idx = df["original_index"].values.astype(int)

C = C_full[np.ix_(idx, idx)]
C = 0.5 * (C + C.T)

eigvals = np.linalg.eigvalsh(C)
eig_min = float(np.min(eigvals))
eig_max = float(np.max(eigvals))
cond_number = float(eig_max / max(eig_min, 1e-300))

jitter_added = 0.0

if eig_min <= 0:
    jitter_added = abs(eig_min) + 1e-12 * eig_max
    C = C + np.eye(C.shape[0]) * jitter_added

Cinv = np.linalg.pinv(C, rcond=1e-12)

print("\nCovariance:")
print("shape =", C.shape)
print("eig_min =", eig_min)
print("eig_max =", eig_max)
print("condition_number =", cond_number)
print("jitter_added =", jitter_added)

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
    xi = Xi_TVGD(z, theta0=theta0)
    return -5.0 * np.log10(xi)


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
    one = np.ones_like(mu_obs)
    r = mu_obs - mu_model

    num = one @ Cinv @ r
    den = one @ Cinv @ one

    return float(num / den)


def chi2_full_with_offset(mu_obs, mu_model):
    off = best_offset_full(mu_obs, mu_model)
    return chi2_full(mu_obs, mu_model + off), off


def project_out_constant(vec):
    one = np.ones_like(vec)

    alpha = float((one @ Cinv @ vec) / (one @ Cinv @ one))

    return vec - alpha * one, alpha


def best_amplitude_for_template(mu_obs, mu_base, template):
    r = mu_obs - mu_base
    den = template @ Cinv @ template

    if abs(den) < 1e-300:
        return np.nan, np.inf

    B = float((template @ Cinv @ r) / den)
    mu_model = mu_base + B * template

    return B, chi2_full(mu_obs, mu_model)


def aic(chi2, k):
    return chi2 + 2.0 * k


def bic(chi2, k, n):
    return chi2 + k * np.log(n)

# ============================================================
# Main data vectors
# ============================================================

z = df["zHD"].values
mu_obs = df["MU_SH0ES"].values
n = len(df)

mu_cmb = mu_lcdm(z, H0=H0_CMB)
mu_local = mu_lcdm(z, H0=H0_LOCAL)

chi2_cmb_abs = chi2_full(mu_obs, mu_cmb)
chi2_local_abs = chi2_full(mu_obs, mu_local)

chi2_cmb_offset, off_cmb = chi2_full_with_offset(mu_obs, mu_cmb)
chi2_local_offset, off_local = chi2_full_with_offset(mu_obs, mu_local)

mu_base_offset = mu_cmb + off_cmb

# ============================================================
# Shape-only scan
# ============================================================

scan_rows = []

for theta0 in THETA0_GRID:
    dmu_raw = delta_mu_tvgd(z, theta0)
    dmu_shape, removed_const = project_out_constant(dmu_raw)

    B_shape, chi2_shape = best_amplitude_for_template(
        mu_obs, mu_base_offset, dmu_shape
    )

    mu_abs = mu_cmb + dmu_raw
    chi2_abs = chi2_full(mu_obs, mu_abs)

    bao = bao_shift_max_percent(theta0)
    cmb = cmb_shift_abs(theta0)

    scan_rows.append({
        "theta0": theta0,
        "chi2_abs": chi2_abs,
        "chi2_shape_only": chi2_shape,
        "B_shape_best": B_shape,
        "removed_constant_from_dmu": removed_const,
        "bao_shift_max_percent": bao,
        "cmb_shift_abs": cmb,
        "shift_z0p01_percent": float(shift_percent(0.01, theta0)),
        "shift_z0p03_percent": float(shift_percent(0.03, theta0)),
        "shift_z0p05_percent": float(shift_percent(0.05, theta0)),
        "shift_z0p10_percent": float(shift_percent(0.10, theta0)),
    })

scan_df = pd.DataFrame(scan_rows)

scan_ok = scan_df[scan_df["bao_shift_max_percent"] < MAX_BAO_SHIFT_PERCENT].copy()

if len(scan_ok) == 0:
    scan_ok = scan_df.copy()

best_shape = scan_ok.sort_values("chi2_shape_only").iloc[0]
best_abs = scan_ok.sort_values("chi2_abs").iloc[0]

theta_shape_best = float(best_shape["theta0"])
theta_abs_best = float(best_abs["theta0"])

chi2_shape_best = float(best_shape["chi2_shape_only"])
chi2_abs_best = float(best_abs["chi2_abs"])

shape_gain_vs_cmb_offset = chi2_cmb_offset - chi2_shape_best
absolute_gain_vs_cmb = chi2_cmb_abs - chi2_abs_best

delta_aic_shape = aic(chi2_cmb_offset, k=1) - aic(chi2_shape_best, k=2)
delta_bic_shape = bic(chi2_cmb_offset, k=1, n=n) - bic(chi2_shape_best, k=2, n=n)

delta_aic_abs = aic(chi2_cmb_abs, k=0) - aic(chi2_abs_best, k=1)
delta_bic_abs = bic(chi2_cmb_abs, k=0, n=n) - bic(chi2_abs_best, k=1, n=n)

# ============================================================
# Shape-only placebo null
# ============================================================

shape_templates = []
theta_templates = []

for _, row in scan_ok.iterrows():
    theta0 = float(row["theta0"])
    dmu_raw = delta_mu_tvgd(z, theta0)
    dmu_shape, _ = project_out_constant(dmu_raw)

    shape_templates.append(dmu_shape)
    theta_templates.append(theta0)

shape_templates = np.array(shape_templates)
theta_templates = np.array(theta_templates)

null_shape_gain = []
null_theta = []
null_B = []

for j in range(N_NULL):
    perm = rng.permutation(n)

    best_c2 = np.inf
    best_th = np.nan
    best_b = np.nan

    for k in range(len(shape_templates)):
        template_perm = shape_templates[k, perm]
        B, c2 = best_amplitude_for_template(mu_obs, mu_base_offset, template_perm)

        if c2 < best_c2:
            best_c2 = c2
            best_th = theta_templates[k]
            best_b = B

    null_shape_gain.append(chi2_cmb_offset - best_c2)
    null_theta.append(best_th)
    null_B.append(best_b)

null_shape_gain = np.array(null_shape_gain)

p_shape_null_positive = float(
    (np.sum(null_shape_gain >= shape_gain_vs_cmb_offset) + 1) / (N_NULL + 1)
)

null_df = pd.DataFrame({
    "null_id": np.arange(N_NULL),
    "shape_gain_null": null_shape_gain,
    "theta_best_null": null_theta,
    "B_best_null": null_B,
})

# ============================================================
# Automatic diagnosis
# ============================================================

shape_gain_ok = shape_gain_vs_cmb_offset > MIN_SHAPE_GAIN
shape_null_ok = p_shape_null_positive <= MAX_P_SHAPE_NULL
aic_shape_ok = delta_aic_shape > 0
bic_shape_ok = delta_bic_shape > 0

if shape_gain_ok and shape_null_ok and aic_shape_ok and bic_shape_ok:
    final_status = "PASSOU_FORTE_SHAPE_ONLY_NONCALIBRATOR_TVGD_HUBBLE"
elif shape_gain_ok and aic_shape_ok:
    final_status = "PASSOU_MODERADO_SHAPE_ONLY_NONCALIBRATOR_TVGD_HUBBLE"
elif shape_gain_ok:
    final_status = "PASSOU_FRACO_SHAPE_ONLY_NONCALIBRATOR_TVGD_HUBBLE"
else:
    final_status = "TENSIONADO_SHAPE_ONLY_NONCALIBRATOR_TVGD_HUBBLE"

summary = {
    "N_sne": n,
    "z_min": float(np.min(z)),
    "z_max": float(np.max(z)),
    "calibrator_cut_used": calibrator_cut_used,

    "H0_CMB": H0_CMB,
    "H0_LOCAL": H0_LOCAL,
    "A_H": A_H,
    "s_natural": S_NATURAL,

    "cov_eig_min_original": eig_min,
    "cov_eig_max_original": eig_max,
    "cov_condition_original": cond_number,
    "cov_jitter_added": jitter_added,

    "chi2_cmb_abs": chi2_cmb_abs,
    "chi2_local_abs": chi2_local_abs,
    "chi2_cmb_offset": chi2_cmb_offset,
    "chi2_local_offset": chi2_local_offset,

    "off_cmb": off_cmb,
    "off_local": off_local,

    "theta_abs_best": theta_abs_best,
    "theta_shape_best": theta_shape_best,

    "chi2_abs_best": chi2_abs_best,
    "chi2_shape_best": chi2_shape_best,

    "absolute_gain_vs_cmb": absolute_gain_vs_cmb,
    "shape_gain_vs_cmb_offset": shape_gain_vs_cmb_offset,

    "best_shape_B": float(best_shape["B_shape_best"]),
    "best_shape_bao_shift_percent": float(best_shape["bao_shift_max_percent"]),
    "best_shape_cmb_shift_abs": float(best_shape["cmb_shift_abs"]),

    "best_shape_shift_z0p01_percent": float(best_shape["shift_z0p01_percent"]),
    "best_shape_shift_z0p03_percent": float(best_shape["shift_z0p03_percent"]),
    "best_shape_shift_z0p05_percent": float(best_shape["shift_z0p05_percent"]),
    "best_shape_shift_z0p10_percent": float(best_shape["shift_z0p10_percent"]),

    "delta_aic_abs_cmb_minus_tvgd": delta_aic_abs,
    "delta_bic_abs_cmb_minus_tvgd": delta_bic_abs,
    "delta_aic_shape_cmboffset_minus_tvgdshape": delta_aic_shape,
    "delta_bic_shape_cmboffset_minus_tvgdshape": delta_bic_shape,

    "N_null": N_NULL,
    "null_shape_gain_mean": float(np.mean(null_shape_gain)),
    "null_shape_gain_median": float(np.median(null_shape_gain)),
    "null_shape_gain_std": float(np.std(null_shape_gain, ddof=1)),
    "null_shape_gain_p95": float(np.quantile(null_shape_gain, 0.95)),
    "p_shape_null_positive": p_shape_null_positive,

    "shape_gain_ok": bool(shape_gain_ok),
    "shape_null_ok": bool(shape_null_ok),
    "aic_shape_ok": bool(aic_shape_ok),
    "bic_shape_ok": bool(bic_shape_ok),

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nFinal summary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

scan_path = os.path.join(OUTDIR, "test07G_noncalibrator_shape_scan.csv")
null_path = os.path.join(OUTDIR, "test07G_noncalibrator_shape_null.csv")
summary_path = os.path.join(OUTDIR, "test07G_noncalibrator_shape_summary.csv")

scan_df.to_csv(scan_path, index=False)
null_df.to_csv(null_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(scan_df["theta0"], scan_df["chi2_shape_only"], label="TVGD shape-only")
plt.axhline(chi2_cmb_offset, color="black", linestyle="--", label="CMB + offset")
plt.axvline(theta_shape_best, color="red", linestyle="--", label="best theta0")
plt.xlabel("theta0")
plt.ylabel("chi2 full covariance")
plt.title("TVGD shape-only non-calibrator scan")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07G_shape_chi2_scan.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(null_shape_gain, bins=35, alpha=0.8, label="shape null")
plt.axvline(shape_gain_vs_cmb_offset, color="red", linestyle="--", label="real shape gain")
plt.axvline(0, color="black", linewidth=1)
plt.xlabel("Shape gain over CMB+offset")
plt.ylabel("N")
plt.title("Non-calibrator shape-only placebo null")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07G_shape_null.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(scan_df["theta0"], scan_df["bao_shift_max_percent"])
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO criterion")
plt.axvline(theta_shape_best, color="blue", linestyle="--", label="best theta0")
plt.xlabel("theta0")
plt.ylabel("max BAO shift [%]")
plt.title("BAO safety")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07G_bao_safety.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(scan_path)
print(null_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_SHAPE_ONLY_NONCALIBRATOR_TVGD_HUBBLE":
    print("The TVGD redshift-shape survives without SH0ES calibrators and passes offset, AIC/BIC, and placebo controls.")
elif "PASSOU_MODERADO" in final_status:
    print("The TVGD redshift-shape shows shape-only gain without SH0ES calibrators, but still requires statistical caution.")
elif "PASSOU_FRACO" in final_status:
    print("The TVGD redshift-shape has shape-only gain, but does not pass enough controls.")
else:
    print("The TVGD redshift-shape was not isolated in the non-calibrator sample.")

print("\nDONE.")
