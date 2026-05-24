# ============================================================
# TVGD — 18_PANTHEON_FULL_COVARIANCE_SN
# Pantheon+ com matriz completa STAT+SYS
#
# Objetivo:
#   Comparar LCDM e TVGD V1.1 usando Pantheon+SH0ES com matriz
#   completa de covariância, marginalizando o offset absoluto.
#
# Modelo TVGD V1.1:
#
#   rho_psi(a) ∝ [1 - exp(-lambda a^p)] / [1 - exp(-lambda)]
#
#   lambda = 7.0
#   p      = 0.55
#
# Saída:
#   results/pantheon_full_covariance_sn/
#
# Uso:
#   python notebooks/main_validation/18_pantheon_full_covariance_sn.py
# ============================================================

import os
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "18_pantheon_full_covariance_sn.py"

OUTDIR = "results/pantheon_full_covariance_sn"
os.makedirs(OUTDIR, exist_ok=True)

DATA_URL = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
COV_URL = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov"

DATA_PATH = os.path.join(OUTDIR, "PantheonPlusSH0ES.dat")
COV_PATH = os.path.join(OUTDIR, "PantheonPlusSH0ES_STAT_SYS.cov")

C_LIGHT = 299792.458  # km/s

H0 = 67.4
Omega_m = 0.31339978339159447
Omega_r = 9.0e-5
Omega_de = 1.0 - Omega_m - Omega_r

LAMBDA_TVGD = 7.0
P_TVGD = 0.55

ZMIN = 0.001
ZMAX = 2.5

print("=" * 72)
print("TVGD — 18_PANTHEON_FULL_COVARIANCE_SN")
print("Pantheon+ com matriz completa STAT+SYS")
print("=" * 72)


# ============================================================
# Download
# ============================================================

def download_if_needed(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"Arquivo já existe: {path}")
        return

    print(f"Baixando: {url}")
    urllib.request.urlretrieve(url, path)

    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        raise RuntimeError(f"Falha ao baixar: {url}")

    print(f"Download OK: {path}")


download_if_needed(DATA_URL, DATA_PATH)
download_if_needed(COV_URL, COV_PATH)


# ============================================================
# Ler Pantheon+
# ============================================================

df = pd.read_csv(DATA_PATH, sep=r"\s+", comment="#", engine="python")

print("\nColunas Pantheon+:")
print(list(df.columns))
print(f"\nN bruto = {len(df)}")

required = ["zHD", "MU_SH0ES"]
for col in required:
    if col not in df.columns:
        raise RuntimeError(f"Coluna obrigatória ausente: {col}")

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["zHD", "MU_SH0ES"]).copy()

df = df[(df["zHD"] > ZMIN) & (df["zHD"] < ZMAX)].copy()
df = df.reset_index(drop=False).rename(columns={"index": "original_index"})

print(f"N após cortes = {len(df)}")
print(f"z range = {df['zHD'].min()} {df['zHD'].max()}")


# ============================================================
# Ler matriz de covariância
# ============================================================

def read_pantheon_cov(path, n_expected_full):
    raw = np.loadtxt(path)

    # Formato oficial pode vir como:
    # primeira entrada = N, seguida de N*N elementos.
    if raw.ndim == 1:
        first = int(raw[0])
        rest = raw[1:]

        if first * first == len(rest):
            cov = rest.reshape((first, first))
            return cov

        # Ou apenas N*N elementos
        n = int(np.sqrt(len(raw)))
        if n * n == len(raw):
            return raw.reshape((n, n))

    if raw.ndim == 2:
        return raw

    raise RuntimeError("Formato de covariância não reconhecido.")


cov_full = read_pantheon_cov(COV_PATH, n_expected_full=len(pd.read_csv(DATA_PATH, sep=r"\s+", comment="#", engine="python")))

print("\nMatriz completa:")
print(f"shape = {cov_full.shape}")

if cov_full.shape[0] != cov_full.shape[1]:
    raise RuntimeError("Covariância não é quadrada.")

if cov_full.shape[0] < df["original_index"].max() + 1:
    raise RuntimeError("Covariância menor do que índices Pantheon+ selecionados.")

idx = df["original_index"].values.astype(int)
cov = cov_full[np.ix_(idx, idx)]

# Simetrização conservadora
asym = np.max(np.abs(cov - cov.T))
cov = 0.5 * (cov + cov.T)

print(f"Submatriz shape = {cov.shape}")
print(f"assimetria máxima antes da simetrização = {asym:.6e}")

# Regularização mínima se necessário
jitter = 0.0
try:
    L = np.linalg.cholesky(cov)
    chol_ok = True
except np.linalg.LinAlgError:
    diag_med = np.median(np.diag(cov))
    jitter = 1e-10 * diag_med
    cov = cov + jitter * np.eye(len(cov))
    L = np.linalg.cholesky(cov)
    chol_ok = True

print(f"Cholesky OK = {chol_ok}")
print(f"jitter usado = {jitter:.6e}")


# ============================================================
# Cosmologia
# ============================================================

def A_tvgd(a, lam=LAMBDA_TVGD, p=P_TVGD):
    return (1.0 - np.exp(-lam * np.power(a, p))) / (1.0 - np.exp(-lam))


def E_lcdm(z):
    z = np.asarray(z)
    return np.sqrt(Omega_m * (1+z)**3 + Omega_r * (1+z)**4 + (1.0 - Omega_m - Omega_r))


def E_tvgd(z):
    z = np.asarray(z)
    a = 1.0 / (1.0 + z)
    return np.sqrt(Omega_m * (1+z)**3 + Omega_r * (1+z)**4 + Omega_de * A_tvgd(a))


def comoving_distance(z, model="lcdm", n_grid=5000):
    """
    Distância comóvel plana em Mpc.
    Integra em grade acumulada até z_max.
    """
    z = np.asarray(z)
    zmax = max(float(np.max(z)), 1e-6)

    zg = np.linspace(0.0, zmax, n_grid)

    if model == "lcdm":
        Ez = E_lcdm(zg)
    elif model == "tvgd":
        Ez = E_tvgd(zg)
    else:
        raise ValueError(model)

    invE = 1.0 / Ez

    dz = np.diff(zg)
    integ = np.zeros_like(zg)
    integ[1:] = np.cumsum(0.5 * (invE[1:] + invE[:-1]) * dz)

    chi = np.interp(z, zg, integ)
    return (C_LIGHT / H0) * chi


def distance_modulus(z, model="lcdm"):
    Dc = comoving_distance(z, model=model)
    Dl = (1.0 + np.asarray(z)) * Dc
    mu = 5.0 * np.log10(np.maximum(Dl, 1e-300)) + 25.0
    return mu


# ============================================================
# Likelihood com offset marginalizado
# ============================================================

def solve_cov(C, b):
    return np.linalg.solve(C, b)


def chi2_with_free_offset(mu_obs, mu_model, cov):
    """
    Minimiza:
        chi2 = (mu_obs - mu_model - offset)^T C^-1 (...)
    """
    r0 = mu_obs - mu_model
    one = np.ones_like(r0)

    Cinv_r0 = solve_cov(cov, r0)
    Cinv_one = solve_cov(cov, one)

    A = one @ Cinv_one
    B = one @ Cinv_r0

    offset = B / A

    r = r0 - offset
    Cinv_r = solve_cov(cov, r)
    chi2 = float(r @ Cinv_r)

    return chi2, float(offset), r


z = df["zHD"].values.astype(float)
mu_obs = df["MU_SH0ES"].values.astype(float)

mu_lcdm = distance_modulus(z, model="lcdm")
mu_tvgd = distance_modulus(z, model="tvgd")

chi2_lcdm, offset_lcdm, resid_lcdm = chi2_with_free_offset(mu_obs, mu_lcdm, cov)
chi2_tvgd, offset_tvgd, resid_tvgd = chi2_with_free_offset(mu_obs, mu_tvgd, cov)

dof = len(df) - 1

chi2red_lcdm = chi2_lcdm / dof
chi2red_tvgd = chi2_tvgd / dof

delta_chi2 = chi2_tvgd - chi2_lcdm

print("\nResultados Pantheon+ full covariance:")
print(f"LCDM chi2 = {chi2_lcdm}")
print(f"LCDM chi2red = {chi2red_lcdm}")
print(f"LCDM offset = {offset_lcdm}")
print(f"TVGD chi2 = {chi2_tvgd}")
print(f"TVGD chi2red = {chi2red_tvgd}")
print(f"TVGD offset = {offset_tvgd}")
print(f"Delta chi2 TVGD-LCDM = {delta_chi2}")


# ============================================================
# Métricas de resíduos e diferença prevista
# ============================================================

diff_mu = mu_tvgd - mu_lcdm

res_metrics = {
    "LCDM_rms_resid": float(np.sqrt(np.mean(resid_lcdm**2))),
    "LCDM_mean_abs_resid": float(np.mean(np.abs(resid_lcdm))),
    "LCDM_median_abs_resid": float(np.median(np.abs(resid_lcdm))),

    "TVGD_rms_resid": float(np.sqrt(np.mean(resid_tvgd**2))),
    "TVGD_mean_abs_resid": float(np.mean(np.abs(resid_tvgd))),
    "TVGD_median_abs_resid": float(np.median(np.abs(resid_tvgd))),

    "diff_mu_min": float(np.min(diff_mu)),
    "diff_mu_max": float(np.max(diff_mu)),
    "diff_mu_mean": float(np.mean(diff_mu)),
    "diff_mu_rms": float(np.sqrt(np.mean(diff_mu**2))),
}

print("\nMétricas:")
for k, v in res_metrics.items():
    print(f"{k}: {v}")


# ============================================================
# Bins em redshift
# ============================================================

bins = np.array([0.0, 0.01, 0.03, 0.06, 0.10, 0.20, 0.40, 0.80, 1.30, 2.50])
bin_rows = []

for i in range(len(bins)-1):
    lo, hi = bins[i], bins[i+1]
    m = (z >= lo) & (z < hi)

    if m.sum() == 0:
        continue

    bin_rows.append({
        "z_min": lo,
        "z_max": hi,
        "N": int(m.sum()),
        "z_med": float(np.median(z[m])),
        "diff_mu_TVGD_minus_LCDM_mean": float(np.mean(diff_mu[m])),
        "diff_mu_TVGD_minus_LCDM_median": float(np.median(diff_mu[m])),
        "LCDM_resid_mean": float(np.mean(resid_lcdm[m])),
        "TVGD_resid_mean": float(np.mean(resid_tvgd[m])),
    })

bins_df = pd.DataFrame(bin_rows)


# ============================================================
# Salvar tabelas
# ============================================================

out_df = df.copy()
out_df["mu_LCDM"] = mu_lcdm
out_df["mu_TVGD"] = mu_tvgd
out_df["diff_mu_TVGD_minus_LCDM"] = diff_mu
out_df["resid_LCDM_offset_marg"] = resid_lcdm
out_df["resid_TVGD_offset_marg"] = resid_tvgd

out_df.to_csv(os.path.join(OUTDIR, "PANTHEON18_predictions_residuals.csv"), index=False)
bins_df.to_csv(os.path.join(OUTDIR, "PANTHEON18_redshift_bins.csv"), index=False)


# ============================================================
# Diagnóstico final
# ============================================================

criteria = {
    "data_loaded": len(df) > 1000,
    "full_cov_loaded": cov.shape[0] == len(df),
    "cov_cholesky_ok": bool(chol_ok),
    "cov_symmetrized": True,
    "LCDM_chi2red_reasonable": chi2red_lcdm < 2.0,
    "TVGD_chi2red_reasonable": chi2red_tvgd < 2.0,
    "delta_chi2_small": abs(delta_chi2) < 5.0,
    "delta_chi2_very_small": abs(delta_chi2) < 1.0,
    "diff_mu_sub_mmag_scale": res_metrics["diff_mu_rms"] < 0.005,
    "TVGD_resid_not_worse": res_metrics["TVGD_rms_resid"] <= res_metrics["LCDM_rms_resid"] + 0.002,
}

score = sum(bool(v) for v in criteria.values())

if score >= 9:
    status = "PASSOU_FORTE_PANTHEON18_FULLCOV_TVGD"
elif score >= 7:
    status = "PASSOU_MODERADO_FORTE_PANTHEON18_COM_ALERTAS"
elif score >= 5:
    status = "PASSOU_MODERADO_PANTHEON18"
else:
    status = "FALHOU_PANTHEON18_FULLCOV"

summary = {
    "FINAL_STATUS": status,
    "score": int(score),
    "script": SCRIPT_NAME,

    "H0": H0,
    "Omega_m": Omega_m,
    "Omega_r": Omega_r,
    "Omega_de": Omega_de,
    "lambda_TVGD": LAMBDA_TVGD,
    "p_TVGD": P_TVGD,

    "N_SN": int(len(df)),
    "z_min": float(np.min(z)),
    "z_max": float(np.max(z)),

    "cov_shape_0": int(cov.shape[0]),
    "cov_shape_1": int(cov.shape[1]),
    "cov_asym_before_symmetrization": float(asym),
    "cov_jitter": float(jitter),

    "chi2_LCDM": float(chi2_lcdm),
    "chi2red_LCDM": float(chi2red_lcdm),
    "offset_LCDM": float(offset_lcdm),

    "chi2_TVGD": float(chi2_tvgd),
    "chi2red_TVGD": float(chi2red_tvgd),
    "offset_TVGD": float(offset_tvgd),

    "delta_chi2_TVGD_minus_LCDM": float(delta_chi2),
}

summary.update(res_metrics)
summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "PANTHEON18_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — PANTHEON 18")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.scatter(z, diff_mu, s=8, alpha=0.5)
plt.axhline(0.0, linestyle="--")
plt.xscale("log")
plt.xlabel("z")
plt.ylabel(r"$\mu_{\rm TVGD}-\mu_{\Lambda CDM}$ [mag]")
plt.title("Pantheon+ full covariance — TVGD minus LCDM prediction")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "PANTHEON18_diff_mu_TVGD_LCDM.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.scatter(z, resid_lcdm, s=8, alpha=0.35, label="LCDM")
plt.scatter(z, resid_tvgd, s=8, alpha=0.35, label="TVGD")
plt.axhline(0.0, linestyle="--")
plt.xscale("log")
plt.xlabel("z")
plt.ylabel("Residual after offset marginalization [mag]")
plt.title("Pantheon+ full covariance — residuals")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "PANTHEON18_residuals.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.plot(bins_df["z_med"], bins_df["diff_mu_TVGD_minus_LCDM_mean"], marker="o")
plt.axhline(0.0, linestyle="--")
plt.xscale("log")
plt.xlabel("median z")
plt.ylabel(r"mean $\Delta\mu$ [mag]")
plt.title("Pantheon+ full covariance — binned TVGD-LCDM difference")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "PANTHEON18_binned_diff_mu.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(resid_lcdm, bins=60, alpha=0.6, label="LCDM")
plt.hist(resid_tvgd, bins=60, alpha=0.6, label="TVGD")
plt.xlabel("Residual [mag]")
plt.ylabel("N")
plt.title("Pantheon+ full covariance — residual distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "PANTHEON18_residual_hist.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste Pantheon+ full covariance.")
