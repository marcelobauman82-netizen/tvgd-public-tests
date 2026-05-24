
# ============================================================
# TVGD - TESTE 29C
# Pantheon+ full-covariance bulk anisotropy test
#
# Objetivo:
# Testar uma modulacao direcional tipo bulk na taxa local inferida
# a partir de supernovas Pantheon+SH0ES, usando a matriz completa
# STAT+SYS de covariancia.
#
# Direcao principal:
#   bulk: (l,b) = (280 deg, 10 deg)
#
# Modelo linear:
#   mu_obs = mu_LCDM(H0_DYN) + offset + A_bulk * X_bulk
#
# onde:
#   X_bulk = -(5/ln 10) * (n.d_bulk)/H0_DYN
#
# Assim, A_bulk sai em km/s/Mpc.
#
# Saida principal:
#   FINAL_STATUS_29C
# ============================================================

import os
import sys
import subprocess
import warnings
warnings.filterwarnings("ignore")

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import requests
    from scipy.stats import chi2
except Exception:
    subprocess.check_call([sys.executable, "-m", "pip", "-q", "install",
                           "numpy", "pandas", "matplotlib", "requests", "scipy"])
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import requests
    from scipy.stats import chi2


# ============================================================
# CONFIGURACAO GERAL
# ============================================================

OUTDIR = "/content/tvgd_test_29C_pantheon_fullcov_bulk"
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 100)
print("TVGD - TESTE 29C: PANTHEON+ FULL-COVARIANCE BULK ANISOTROPY")
print("=" * 100)
print("Saida:", OUTDIR)

# Cosmologia de referencia usada no teste
C_LIGHT = 299792.458  # km/s
H0_DYN = 67.4
OMEGA_M = 0.315
OMEGA_L = 1.0 - OMEGA_M

# Direcao bulk principal
BULK_L = 280.0
BULK_B = 10.0

# Controle X-ray aproximado
XRAY_L = 280.0
XRAY_B = -20.0

# Controle oposto
CONTROL_L = 100.0
CONTROL_B = -10.0

# Bootstrap/null
N_BOOT = 800
N_NULL = 800
SEED = 290300
rng = np.random.default_rng(SEED)

# URLs Pantheon+SH0ES
URL_DAT = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES.dat"
URL_COV = "https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/Pantheon%2BSH0ES_STAT%2BSYS.cov"

DAT_PATH = os.path.join(OUTDIR, "Pantheon+SH0ES.dat")
COV_PATH = os.path.join(OUTDIR, "Pantheon+SH0ES_STAT+SYS.cov")

print("H0_DYN =", H0_DYN)
print("Omega_m =", OMEGA_M)
print("Bulk direction l,b =", BULK_L, BULK_B)


# ============================================================
# DOWNLOAD DOS DADOS
# ============================================================

def download_file(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print("Arquivo ja existe:", path)
        return

    print("Baixando:", url)
    r = requests.get(url, timeout=120)
    r.raise_for_status()

    with open(path, "wb") as f:
        f.write(r.content)

    print("Salvo:", path, "tamanho:", os.path.getsize(path))


download_file(URL_DAT, DAT_PATH)
download_file(URL_COV, COV_PATH)


# ============================================================
# LEITURA PANTHEON+ E COVARIANCIA
# ============================================================

sn = pd.read_csv(DAT_PATH, sep=r"\s+", comment="#")

print("\nPantheon+ carregado.")
print("N bruto:", len(sn))
print("Colunas:")
print(list(sn.columns))

required_cols = ["zHD", "MU_SH0ES", "MU_SH0ES_ERR_DIAG", "RA", "DEC"]
for c in required_cols:
    if c not in sn.columns:
        raise RuntimeError(f"Coluna obrigatoria ausente: {c}")

def load_pantheon_covariance(path, n_expected):
    """
    Formato Pantheon+:
    primeira entrada geralmente contem N,
    seguida por N*N valores da matriz achatada.
    """
    vals = np.loadtxt(path)

    if vals.size == n_expected * n_expected + 1:
        n = int(vals[0])
        arr = vals[1:]
    elif vals.size == n_expected * n_expected:
        n = n_expected
        arr = vals
    else:
        raise RuntimeError(
            f"Tamanho inesperado da covariancia: {vals.size}. "
            f"Esperado {n_expected*n_expected} ou {n_expected*n_expected+1}."
        )

    if n != n_expected:
        raise RuntimeError(f"N da covariancia {n} diferente de N dos dados {n_expected}")

    return arr.reshape((n_expected, n_expected))


COV_FULL = load_pantheon_covariance(COV_PATH, len(sn))
print("\nCovariancia carregada:", COV_FULL.shape)


# ============================================================
# FUNCOES COSMOLOGICAS
# ============================================================

def E_z(z):
    z = np.asarray(z, dtype=float)
    return np.sqrt(OMEGA_M * (1.0 + z)**3 + OMEGA_L)

def luminosity_distance_mpc(z):
    """
    Distancia de luminosidade plana usando integracao numerica simples.
    Suficiente para z baixo do teste.
    """
    z = np.asarray(z, dtype=float)
    out = np.zeros_like(z, dtype=float)

    for i, zi in enumerate(z):
        if zi <= 0:
            out[i] = np.nan
            continue

        grid = np.linspace(0.0, zi, 800)
        integ = np.trapz(1.0 / E_z(grid), grid)
        dc = (C_LIGHT / H0_DYN) * integ
        dl = (1.0 + zi) * dc
        out[i] = dl

    return out

def mu_lcdm(z):
    dl = luminosity_distance_mpc(z)
    return 5.0 * np.log10(dl) + 25.0


# ============================================================
# FUNCOES DE CEU
# ============================================================

def sky_unit_vectors_equatorial(ra_deg, dec_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))

    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    return np.vstack([x, y, z]).T

def galactic_to_equatorial(l_deg, b_deg):
    l = np.deg2rad(l_deg)
    b = np.deg2rad(b_deg)

    xg = np.cos(b) * np.cos(l)
    yg = np.cos(b) * np.sin(l)
    zg = np.sin(b)
    vg = np.array([xg, yg, zg])

    R_eq_to_gal = np.array([
        [-0.0548755604, -0.8734370902, -0.4838350155],
        [ 0.4941094279, -0.4448296300,  0.7469822445],
        [-0.8676661490, -0.1980763734,  0.4559837762]
    ])

    veq = R_eq_to_gal.T @ vg
    x, y, z = veq

    ra = np.arctan2(y, x)
    ra = np.mod(ra, 2.0 * np.pi)
    dec = np.arcsin(z / np.linalg.norm(veq))

    return np.rad2deg(ra), np.rad2deg(dec)

def unit_from_galactic(l_deg, b_deg):
    ra, dec = galactic_to_equatorial(l_deg, b_deg)
    return sky_unit_vectors_equatorial([ra], [dec])[0]

def add_projection(df, l_deg, b_deg, name):
    d = df.copy()
    v = unit_from_galactic(l_deg, b_deg)
    u = d[["ux", "uy", "uz"]].values
    d[name] = u @ v
    return d


# ============================================================
# PREPARAR DADOS
# ============================================================

sn = sn.copy()
sn["row_index_original"] = np.arange(len(sn))

sn["z"] = pd.to_numeric(sn["zHD"], errors="coerce")
sn["mu_obs"] = pd.to_numeric(sn["MU_SH0ES"], errors="coerce")
sn["mu_err_diag"] = pd.to_numeric(sn["MU_SH0ES_ERR_DIAG"], errors="coerce")
sn["RA"] = pd.to_numeric(sn["RA"], errors="coerce")
sn["DEC"] = pd.to_numeric(sn["DEC"], errors="coerce")

sn = sn.replace([np.inf, -np.inf], np.nan)
sn = sn.dropna(subset=["z", "mu_obs", "mu_err_diag", "RA", "DEC"]).copy()

uv = sky_unit_vectors_equatorial(sn["RA"].values, sn["DEC"].values)
sn["ux"] = uv[:, 0]
sn["uy"] = uv[:, 1]
sn["uz"] = uv[:, 2]

sn["mu_model"] = mu_lcdm(sn["z"].values)
sn["res_mu"] = sn["mu_obs"] - sn["mu_model"]

print("\nDados preparados:")
print("N apos limpeza:", len(sn))
print("z min/med/max:", sn["z"].min(), sn["z"].median(), sn["z"].max())


# ============================================================
# GLS COM COVARIANCIA COMPLETA
# ============================================================

def gls_fit(y, X, C):
    """
    Ajuste GLS:
    beta = (X^T C^-1 X)^-1 X^T C^-1 y
    chi2 = (y-Xb)^T C^-1 (y-Xb)
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    C = np.asarray(C, dtype=float)

    # Pequena regularizacao numerica
    diag_med = np.nanmedian(np.diag(C))
    C_reg = C + np.eye(C.shape[0]) * diag_med * 1e-10

    try:
        L = np.linalg.cholesky(C_reg)
        CiX = np.linalg.solve(L.T, np.linalg.solve(L, X))
        Ciy = np.linalg.solve(L.T, np.linalg.solve(L, y))
    except Exception:
        Cinv = np.linalg.pinv(C_reg)
        CiX = Cinv @ X
        Ciy = Cinv @ y

    Fisher = X.T @ CiX
    cov_beta = np.linalg.pinv(Fisher)
    beta = cov_beta @ (X.T @ Ciy)

    resid = y - X @ beta

    try:
        chi2_val = float(resid.T @ np.linalg.solve(L.T, np.linalg.solve(L, resid)))
    except Exception:
        chi2_val = float(resid.T @ np.linalg.pinv(C_reg) @ resid)

    dof = len(y) - X.shape[1]

    return beta, cov_beta, chi2_val, dof

def build_design(df, l_deg=None, b_deg=None, include_bulk=True):
    d = df.copy()

    X_parts = []
    names = []

    # Offset global
    X_parts.append(np.ones(len(d)))
    names.append("offset")

    if include_bulk:
        d = add_projection(d, l_deg, b_deg, "proj_bulk")
        proj = d["proj_bulk"].values

        # dmu/dA_H0 = -(5/ln10)*(proj/H0)
        x_bulk = -(5.0 / np.log(10.0)) * proj / H0_DYN

        X_parts.append(x_bulk)
        names.append("A_bulk_km_s_Mpc")

    X = np.vstack(X_parts).T
    return X, names

def fit_fixed_direction_fullcov(df_window, C_window, l_deg, b_deg, label):
    y = df_window["res_mu"].values

    X0, names0 = build_design(df_window, include_bulk=False)
    X1, names1 = build_design(df_window, l_deg=l_deg, b_deg=b_deg, include_bulk=True)

    beta0, cov0, chi20, dof0 = gls_fit(y, X0, C_window)
    beta1, cov1, chi21, dof1 = gls_fit(y, X1, C_window)

    dchi = chi20 - chi21
    pval = chi2.sf(dchi, 1)

    A = beta1[1]
    sigma_A = np.sqrt(max(cov1[1, 1], 0.0))
    snr = A / sigma_A if sigma_A > 0 else np.nan

    return {
        "label": label,
        "l_fixed": float(l_deg),
        "b_fixed": float(b_deg),
        "N": int(len(df_window)),
        "z_min": float(df_window["z"].min()),
        "z_max": float(df_window["z"].max()),
        "A_bulk": float(A),
        "sigma_A": float(sigma_A),
        "A_over_sigma": float(snr),
        "offset": float(beta1[0]),
        "chi2_null": float(chi20),
        "chi2_bulk": float(chi21),
        "delta_chi2": float(dchi),
        "p_chi2": float(pval),
        "dof_null": int(dof0),
        "dof_bulk": int(dof1)
    }


# ============================================================
# RODAR JANELAS PRINCIPAIS
# ============================================================

WINDOWS = {
    "z001_015_main": (0.01, 0.15),
    "z001_020": (0.01, 0.20),
    "z001_002_low": (0.01, 0.02),
    "z002_015_noverylow": (0.02, 0.15),
    "z0023_015": (0.023, 0.15),
    "z003_014": (0.03, 0.14),
}

rows = []

for win, (zmin, zmax) in WINDOWS.items():
    mask = (sn["z"] >= zmin) & (sn["z"] < zmax)
    sub = sn[mask].copy().reset_index(drop=True)

    if len(sub) < 20:
        continue

    idx = sn[mask]["row_index_original"].values.astype(int)
    Csub = COV_FULL[np.ix_(idx, idx)]

    print("\n" + "=" * 100)
    print("Janela:", win, "N =", len(sub), "z =", zmin, zmax)
    print("=" * 100)

    for label, lval, bval in [
        ("bulk_29C", BULK_L, BULK_B),
        ("xray_control", XRAY_L, XRAY_B),
        ("opposite_control", CONTROL_L, CONTROL_B)
    ]:
        res = fit_fixed_direction_fullcov(sub, Csub, lval, bval, label)
        res["window"] = win
        res["z_window_min"] = zmin
        res["z_window_max"] = zmax
        rows.append(res)

df_results = pd.DataFrame(rows)
df_results.to_csv(os.path.join(OUTDIR, "teste29C_resultados_janelas.csv"), index=False)

print("\n" + "=" * 100)
print("29C - RESULTADOS POR JANELA")
print("=" * 100)
print(df_results)


# ============================================================
# JANELA PRINCIPAL, BOOTSTRAP DIAGONAL E NULL FULLCOV
# ============================================================

MAIN_WINDOW = "z001_015_main"
main_row = df_results[(df_results["window"] == MAIN_WINDOW) &
                      (df_results["label"] == "bulk_29C")].iloc[0].to_dict()

zmin = main_row["z_window_min"]
zmax = main_row["z_window_max"]

mask_main = (sn["z"] >= zmin) & (sn["z"] < zmax)
main_df = sn[mask_main].copy().reset_index(drop=True)
main_idx = sn[mask_main]["row_index_original"].values.astype(int)
Cmain = COV_FULL[np.ix_(main_idx, main_idx)]

print("\nJanela principal:", MAIN_WINDOW)
print("N_main =", len(main_df))

# Bootstrap diagonal para estabilidade de amplitude
def fit_fixed_direction_diag(df_window, l_deg, b_deg, label):
    d = df_window.copy()
    d = add_projection(d, l_deg, b_deg, "proj_bulk")

    y = d["res_mu"].values
    sigma = d["mu_err_diag"].values

    x_bulk = -(5.0 / np.log(10.0)) * d["proj_bulk"].values / H0_DYN

    X0 = np.vstack([np.ones(len(d))]).T
    X1 = np.vstack([np.ones(len(d)), x_bulk]).T

    Cdiag = np.diag(sigma**2)

    beta0, cov0, chi20, dof0 = gls_fit(y, X0, Cdiag)
    beta1, cov1, chi21, dof1 = gls_fit(y, X1, Cdiag)

    dchi = chi20 - chi21
    pval = chi2.sf(dchi, 1)

    A = beta1[1]
    sigma_A = np.sqrt(max(cov1[1, 1], 0.0))

    return {
        "label": label,
        "N": int(len(d)),
        "A_bulk": float(A),
        "sigma_A": float(sigma_A),
        "A_over_sigma": float(A / sigma_A if sigma_A > 0 else np.nan),
        "delta_chi2": float(dchi),
        "p_chi2": float(pval)
    }

boot_rows = []

for i in range(N_BOOT):
    idx = rng.integers(0, len(main_df), len(main_df))
    sample = main_df.iloc[idx].copy().reset_index(drop=True)

    try:
        r = fit_fixed_direction_diag(sample, BULK_L, BULK_B, "bulk_29C")
        r["boot_id"] = i
        boot_rows.append(r)
    except Exception:
        pass

boot = pd.DataFrame(boot_rows)
boot.to_csv(os.path.join(OUTDIR, "teste29C_bootstrap_diagonal.csv"), index=False)

df_boot_summary = pd.DataFrame([{
    "window": MAIN_WINDOW,
    "N_boot": len(boot),
    "A_real_fullcov": main_row["A_bulk"],
    "A_boot_median": float(np.nanmedian(boot["A_bulk"])),
    "A_boot_p025": float(np.nanpercentile(boot["A_bulk"], 2.5)),
    "A_boot_p16": float(np.nanpercentile(boot["A_bulk"], 16)),
    "A_boot_p84": float(np.nanpercentile(boot["A_bulk"], 84)),
    "A_boot_p975": float(np.nanpercentile(boot["A_bulk"], 97.5)),
    "P_boot_A_positive": float(np.mean(boot["A_bulk"] > 0)),
    "P_boot_A_gt_1": float(np.mean(boot["A_bulk"] > 1.0)),
}])

df_boot_summary.to_csv(os.path.join(OUTDIR, "teste29C_bootstrap_summary.csv"), index=False)

print("\n" + "=" * 100)
print("29C - BOOTSTRAP SUMMARY")
print("=" * 100)
print(df_boot_summary)


# Null angular fullcov: permuta direcoes mantendo z, mu e covariancia
null_rows = []

for j in range(N_NULL):
    dnull = main_df.copy()

    perm = rng.permutation(len(dnull))
    dnull[["RA", "DEC", "ux", "uy", "uz"]] = dnull.loc[perm, ["RA", "DEC", "ux", "uy", "uz"]].values

    try:
        r = fit_fixed_direction_fullcov(dnull, Cmain, BULK_L, BULK_B, "bulk_29C")
        r["null_id"] = j
        null_rows.append(r)
    except Exception:
        pass

null = pd.DataFrame(null_rows)
null.to_csv(os.path.join(OUTDIR, "teste29C_null_fullcov.csv"), index=False)

A_real = float(main_row["A_bulk"])
dchi_real = float(main_row["delta_chi2"])

p_A_positive = (np.sum(null["A_bulk"] >= A_real) + 1) / (len(null) + 1)
p_absA = (np.sum(np.abs(null["A_bulk"]) >= abs(A_real)) + 1) / (len(null) + 1)
p_dchi = (np.sum(null["delta_chi2"] >= dchi_real) + 1) / (len(null) + 1)

df_null_summary = pd.DataFrame([{
    "window": MAIN_WINDOW,
    "N_null": len(null),
    "A_real": A_real,
    "delta_chi2_real": dchi_real,
    "A_null_median": float(np.nanmedian(null["A_bulk"])),
    "absA_null_median": float(np.nanmedian(np.abs(null["A_bulk"]))),
    "delta_chi2_null_median": float(np.nanmedian(null["delta_chi2"])),
    "p_A_positive_null": float(p_A_positive),
    "p_absA_null": float(p_absA),
    "p_delta_chi2_null": float(p_dchi),
}])

df_null_summary.to_csv(os.path.join(OUTDIR, "teste29C_null_summary.csv"), index=False)

print("\n" + "=" * 100)
print("29C - NULL SUMMARY")
print("=" * 100)
print(df_null_summary)


# ============================================================
# JACKKNIFE ANGULAR FULLCOV
# ============================================================

def add_jackknife_regions(d):
    out = d.copy()
    ra_bin = pd.cut(out["RA"], bins=[0, 90, 180, 270, 360],
                    labels=[0, 1, 2, 3], include_lowest=True)
    dec_bin = (out["DEC"] >= 0).astype(int)
    out["jk_region"] = ra_bin.astype(int) * 2 + dec_bin.astype(int)
    return out

tmp = add_jackknife_regions(main_df)
jk_rows = []

for reg in sorted(tmp["jk_region"].unique()):
    keep = tmp["jk_region"] != reg
    sub = tmp[keep].copy().reset_index(drop=True)

    original_indices = main_idx[keep.values]
    Csub = COV_FULL[np.ix_(original_indices, original_indices)]

    if len(sub) < 20:
        continue

    r = fit_fixed_direction_fullcov(sub, Csub, BULK_L, BULK_B, "bulk_29C")
    r["left_out_region"] = int(reg)
    r["N_removed"] = int(np.sum(~keep))
    jk_rows.append(r)

df_jk = pd.DataFrame(jk_rows)
df_jk.to_csv(os.path.join(OUTDIR, "teste29C_jackknife_angular.csv"), index=False)

print("\n" + "=" * 100)
print("29C - JACKKNIFE ANGULAR")
print("=" * 100)
print(df_jk)


# ============================================================
# LEAVE-ONE-SURVEY FULLCOV
# ============================================================

survey_rows = []

if "IDSURVEY" in sn.columns:
    main_df["IDSURVEY"] = sn[mask_main]["IDSURVEY"].values
    surveys = sorted(main_df["IDSURVEY"].dropna().unique())

    for sid in surveys:
        keep = main_df["IDSURVEY"] != sid
        sub = main_df[keep].copy().reset_index(drop=True)

        original_indices = main_idx[keep.values]
        Csub = COV_FULL[np.ix_(original_indices, original_indices)]

        if len(sub) < 20:
            continue

        r = fit_fixed_direction_fullcov(sub, Csub, BULK_L, BULK_B, "bulk_29C")
        r["left_out_survey"] = int(sid)
        r["N_removed"] = int(np.sum(~keep))
        survey_rows.append(r)

df_survey = pd.DataFrame(survey_rows)
df_survey.to_csv(os.path.join(OUTDIR, "teste29C_leave_one_survey.csv"), index=False)

print("\n" + "=" * 100)
print("29C - LEAVE ONE SURVEY")
print("=" * 100)
print(df_survey)


# ============================================================
# MASCARA GALACTICA APROXIMADA
# ============================================================

def equatorial_to_galactic(ra_deg, dec_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))

    ra_ngp = np.deg2rad(192.85948)
    dec_ngp = np.deg2rad(27.12825)
    l_omega = np.deg2rad(32.93192)

    b = np.arcsin(
        np.sin(dec) * np.sin(dec_ngp)
        + np.cos(dec) * np.cos(dec_ngp) * np.cos(ra - ra_ngp)
    )

    y = np.cos(dec) * np.sin(ra - ra_ngp)
    x = (
        np.sin(dec) * np.cos(dec_ngp)
        - np.cos(dec) * np.sin(dec_ngp) * np.cos(ra - ra_ngp)
    )

    l = l_omega + np.arctan2(y, x)
    l = np.mod(l, 2.0 * np.pi)

    return np.rad2deg(l), np.rad2deg(b)

glon, glat = equatorial_to_galactic(main_df["RA"].values, main_df["DEC"].values)
main_df["GLON"] = glon
main_df["GLAT"] = glat
main_df["abs_GLAT"] = np.abs(glat)

mask_rows = []

for bcut in [0, 10, 15, 20, 30]:
    keep = main_df["abs_GLAT"] >= bcut
    sub = main_df[keep].copy().reset_index(drop=True)

    original_indices = main_idx[keep.values]
    Csub = COV_FULL[np.ix_(original_indices, original_indices)]

    if len(sub) < 20:
        continue

    r = fit_fixed_direction_fullcov(sub, Csub, BULK_L, BULK_B, "bulk_29C")
    r["mask_abs_GLAT_ge"] = bcut
    mask_rows.append(r)

df_masks = pd.DataFrame(mask_rows)
df_masks.to_csv(os.path.join(OUTDIR, "teste29C_galactic_masks.csv"), index=False)

print("\n" + "=" * 100)
print("29C - MASCARA GALACTICA")
print("=" * 100)
print(df_masks)


# ============================================================
# DIAGNOSTICO GLOBAL
# ============================================================

main_xray = df_results[(df_results["window"] == MAIN_WINDOW) &
                       (df_results["label"] == "xray_control")].iloc[0].to_dict()

main_control = df_results[(df_results["window"] == MAIN_WINDOW) &
                          (df_results["label"] == "opposite_control")].iloc[0].to_dict()

jack_frac_positive = float(np.mean(df_jk["A_bulk"] > 0)) if len(df_jk) else np.nan
jack_frac_p_lt_005 = float(np.mean(df_jk["p_chi2"] < 0.05)) if len(df_jk) else np.nan

survey_frac_positive = float(np.mean(df_survey["A_bulk"] > 0)) if len(df_survey) else np.nan
survey_frac_p_lt_005 = float(np.mean(df_survey["p_chi2"] < 0.05)) if len(df_survey) else np.nan

mask_frac_positive = float(np.mean(df_masks["A_bulk"] > 0)) if len(df_masks) else np.nan
mask_frac_p_lt_005 = float(np.mean(df_masks["p_chi2"] < 0.05)) if len(df_masks) else np.nan

n_robust_controls_ok = int(jack_frac_positive == 1.0) + int(survey_frac_positive == 1.0) + int(mask_frac_positive == 1.0)

if (
    main_row["A_over_sigma"] > 3.0
    and main_row["p_chi2"] < 0.01
    and p_absA < 0.01
    and p_dchi < 0.01
    and n_robust_controls_ok >= 2
):
    final_status = "PASSOU_FORTE_29C_FULLCOV_BULK_ROBUSTO"
elif (
    main_row["A_over_sigma"] > 2.0
    and main_row["p_chi2"] < 0.05
    and n_robust_controls_ok >= 1
):
    final_status = "PASSOU_MODERADO_29C_FULLCOV_BULK"
else:
    final_status = "FALHOU_29C_SEM_EVIDENCIA_ROBUSTA"

df_diag = pd.DataFrame([{
    "MAIN_WINDOW": MAIN_WINDOW,
    "N_main": int(main_row["N"]),
    "bulk_l": BULK_L,
    "bulk_b": BULK_B,
    "A_bulk": float(main_row["A_bulk"]),
    "sigma_A": float(main_row["sigma_A"]),
    "A_over_sigma": float(main_row["A_over_sigma"]),
    "delta_chi2": float(main_row["delta_chi2"]),
    "p_chi2": float(main_row["p_chi2"]),
    "p_A_positive_null": float(p_A_positive),
    "p_absA_null": float(p_absA),
    "p_delta_chi2_null": float(p_dchi),
    "A_boot_median": float(df_boot_summary["A_boot_median"].iloc[0]),
    "A_boot_p16": float(df_boot_summary["A_boot_p16"].iloc[0]),
    "A_boot_p84": float(df_boot_summary["A_boot_p84"].iloc[0]),
    "P_boot_A_positive": float(df_boot_summary["P_boot_A_positive"].iloc[0]),
    "P_boot_A_gt_1": float(df_boot_summary["P_boot_A_gt_1"].iloc[0]),
    "xray_A": float(main_xray["A_bulk"]),
    "xray_p": float(main_xray["p_chi2"]),
    "control_A": float(main_control["A_bulk"]),
    "control_p": float(main_control["p_chi2"]),
    "jack_N": len(df_jk),
    "jack_frac_positive": jack_frac_positive,
    "jack_frac_p_lt_005": jack_frac_p_lt_005,
    "survey_N": len(df_survey),
    "survey_frac_positive": survey_frac_positive,
    "survey_frac_p_lt_005": survey_frac_p_lt_005,
    "mask_N": len(df_masks),
    "mask_frac_positive": mask_frac_positive,
    "mask_frac_p_lt_005": mask_frac_p_lt_005,
    "n_robust_controls_ok": n_robust_controls_ok,
    "FINAL_STATUS_29C": final_status
}])

df_diag.to_csv(os.path.join(OUTDIR, "teste29C_diagnostico_global.csv"), index=False)

print("\n" + "=" * 100)
print("29C - DIAGNOSTICO GLOBAL")
print("=" * 100)
print(df_diag)


# ============================================================
# FIGURAS
# ============================================================

plt.figure(figsize=(10, 5))
for lab in ["bulk_29C", "xray_control", "opposite_control"]:
    sub = df_results[df_results["label"] == lab]
    plt.plot(sub["window"], sub["A_bulk"], marker="o", label=lab)

plt.axhline(0, linestyle="--", linewidth=1)
plt.xticks(rotation=45, ha="right")
plt.ylabel("A_bulk [km/s/Mpc]")
plt.title("TVGD 29C - Bulk amplitude by redshift window")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig29C_A_bulk_by_window.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(boot["A_bulk"], bins=40, alpha=0.75)
plt.axvline(main_row["A_bulk"], linewidth=2, label="fullcov real")
plt.axvline(0, linestyle="--", linewidth=1)
plt.xlabel("A_bulk bootstrap diagonal [km/s/Mpc]")
plt.ylabel("N")
plt.title("TVGD 29C - Bootstrap")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig29C_bootstrap.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(null["delta_chi2"], bins=40, alpha=0.75)
plt.axvline(main_row["delta_chi2"], linewidth=2, label="real")
plt.xlabel("Delta chi2 null")
plt.ylabel("N")
plt.title("TVGD 29C - Null angular permutation")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "fig29C_null_delta_chi2.png"), dpi=220)
plt.show()


# ============================================================
# RELATORIO TXT
# ============================================================

report_path = os.path.join(OUTDIR, "RELATORIO_TESTE29C.txt")

with open(report_path, "w", encoding="utf-8") as f:
    f.write("TVGD - TESTE 29C: Pantheon+ full-covariance bulk anisotropy\n")
    f.write("=" * 100 + "\n\n")

    f.write("Modelo:\n")
    f.write("mu_obs = mu_LCDM + offset + A_bulk * X_bulk\n")
    f.write("X_bulk = -(5/ln10) * (n.d_bulk)/H0_DYN\n\n")

    f.write(f"H0_DYN = {H0_DYN}\n")
    f.write(f"Bulk direction: l={BULK_L}, b={BULK_B}\n\n")

    f.write("Diagnostico global:\n")
    f.write(df_diag.to_string(index=False))
    f.write("\n\nResultados por janela:\n")
    f.write(df_results.to_string(index=False))
    f.write("\n\nBootstrap summary:\n")
    f.write(df_boot_summary.to_string(index=False))
    f.write("\n\nNull summary:\n")
    f.write(df_null_summary.to_string(index=False))
    f.write("\n\nJackknife angular:\n")
    f.write(df_jk.to_string(index=False))
    f.write("\n\nLeave-one-survey:\n")
    f.write(df_survey.to_string(index=False))
    f.write("\n\nGalactic masks:\n")
    f.write(df_masks.to_string(index=False))
    f.write("\n\n")

print("\nRelatorio salvo em:", report_path)

print("\nArquivos principais gerados:")
for fn in os.listdir(OUTDIR):
    if fn.startswith("teste29C") or fn.startswith("fig29C") or fn.startswith("RELATORIO"):
        print(" -", fn)

print("\n" + "=" * 100)
print("TESTE 29C CONCLUIDO")
print("=" * 100)
print("FINAL_STATUS_29C =", final_status)
