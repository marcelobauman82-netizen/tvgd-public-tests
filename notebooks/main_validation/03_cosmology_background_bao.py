# ============================================================
# TVGD PUBLIC TESTS
# TEST 03 — COMPRESSED COSMOLOGY BACKGROUND + BAO
# Versão autônoma para Google Colab
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Configuração
# ============================================================

OUTDIR = "Resultados/test03_cosmology_background_bao"
os.makedirs(OUTDIR, exist_ok=True)

# Constantes
C_KM_S = 299792.458

# Cosmologia fiducial
H0 = 67.4
h = H0 / 100.0

OMEGA_M = 0.31339978339159447
OMEGA_R = 9.0e-5
OMEGA_DE = 1.0 - OMEGA_M - OMEGA_R

# Escala sonora fiducial
RD_MPC = 147.09
RS_STAR_MPC = 144.6
Z_STAR = 1089.92

# TVGD V1.1 efetiva
LAMBDA_TVGD = 7.0
P_TVGD = 0.55

# Grade para comparação visual
Z_MAX_BACKGROUND = 3.0
N_Z = 800

# Critérios automáticos
MAX_MEAN_H_DIFF_PERCENT = 2.0
MAX_MEAN_DISTANCE_DIFF_PERCENT = 1.0
MAX_DELTA_CHI2_BAO_MODERATE = 5.0
MAX_DELTA_CHI2_BAO_STRONG = 2.0
MAX_DELTA_LA_PERCENT = 1.0
MAX_DELTA_R_PERCENT = 1.0

# ============================================================
# Tabela BAO comprimida usada neste teste
# ============================================================
#
# Observáveis:
# - DV/rd em z baixo
# - DM/rd e DH/rd em BOSS-like redshifts
#
# Este é um teste comprimido e reprodutível.
# Em versões futuras, a tabela pode ser substituída por uma likelihood BAO completa.
# ============================================================

bao_data = pd.DataFrame([
    {"name": "6dFGS",      "z": 0.106, "observable": "DV_over_rd", "obs": 3.047,  "sigma": 0.137},
    {"name": "SDSS_MGS",   "z": 0.150, "observable": "DV_over_rd", "obs": 4.480,  "sigma": 0.168},

    {"name": "BOSS_DR12",  "z": 0.380, "observable": "DM_over_rd", "obs": 10.234, "sigma": 0.170},
    {"name": "BOSS_DR12",  "z": 0.380, "observable": "DH_over_rd", "obs": 24.981, "sigma": 0.590},

    {"name": "BOSS_DR12",  "z": 0.510, "observable": "DM_over_rd", "obs": 13.365, "sigma": 0.210},
    {"name": "BOSS_DR12",  "z": 0.510, "observable": "DH_over_rd", "obs": 22.317, "sigma": 0.480},

    {"name": "BOSS_DR12",  "z": 0.610, "observable": "DM_over_rd", "obs": 15.515, "sigma": 0.250},
    {"name": "BOSS_DR12",  "z": 0.610, "observable": "DH_over_rd", "obs": 20.747, "sigma": 0.470},
])

# ============================================================
# Funções cosmológicas
# ============================================================

def omega_psi_tvgd(a, omega_m=OMEGA_M, omega_r=OMEGA_R, lam=LAMBDA_TVGD, p=P_TVGD):
    """
    Densidade escura efetiva TVGD V1.1:

    Omega_psi(a) = Omega_de0 * [1 - exp(-lambda a^p)] / [1 - exp(-lambda)]
    """
    a = np.asarray(a, dtype=float)
    omega_de0 = 1.0 - omega_m - omega_r

    numerator = 1.0 - np.exp(-lam * np.maximum(a, 0.0) ** p)
    denominator = 1.0 - np.exp(-lam)

    return omega_de0 * numerator / denominator


def E_lcdm(z):
    """
    E(z) = H(z)/H0 para Lambda-CDM plano.
    """
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)

    return np.sqrt(
        OMEGA_M * a**(-3)
        + OMEGA_R * a**(-4)
        + OMEGA_DE
    )


def E_tvgd(z):
    """
    E(z) = H(z)/H0 para TVGD V1.1 efetiva.
    """
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)

    return np.sqrt(
        OMEGA_M * a**(-3)
        + OMEGA_R * a**(-4)
        + omega_psi_tvgd(a)
    )


def H_of_z(z, model="lcdm"):
    if model == "lcdm":
        return H0 * E_lcdm(z)
    elif model == "tvgd":
        return H0 * E_tvgd(z)
    else:
        raise ValueError("model must be 'lcdm' or 'tvgd'")


def integrate_trapz(z, inv_E_func, n_grid=6000):
    """
    Integral de 0 até z de dz/E(z).
    """
    z = float(z)

    if z == 0:
        return 0.0

    zz = np.linspace(0.0, z, n_grid)
    yy = inv_E_func(zz)

    return float(np.trapezoid(yy, zz))


def comoving_distance_DM(z, model="lcdm", n_grid=6000):
    """
    Distância comóvel transversal DM(z) em Mpc para universo plano.
    """
    if model == "lcdm":
        integral = integrate_trapz(z, lambda zz: 1.0 / E_lcdm(zz), n_grid=n_grid)
    elif model == "tvgd":
        integral = integrate_trapz(z, lambda zz: 1.0 / E_tvgd(zz), n_grid=n_grid)
    else:
        raise ValueError("model must be 'lcdm' or 'tvgd'")

    return (C_KM_S / H0) * integral


def hubble_distance_DH(z, model="lcdm"):
    """
    DH(z) = c/H(z) em Mpc.
    """
    return C_KM_S / H_of_z(z, model=model)


def volume_distance_DV(z, model="lcdm", n_grid=6000):
    """
    DV(z) = [ z * DM(z)^2 * DH(z) ]^(1/3)
    """
    DM = comoving_distance_DM(z, model=model, n_grid=n_grid)
    DH = hubble_distance_DH(z, model=model)

    return (z * DM**2 * DH) ** (1.0 / 3.0)


def distance_modulus(z, model="lcdm", n_grid=6000):
    """
    mu(z) = 5 log10(DL/Mpc) + 25
    """
    DM = np.array([comoving_distance_DM(zi, model=model, n_grid=n_grid) for zi in np.atleast_1d(z)])
    DL = (1.0 + np.atleast_1d(z)) * DM
    mu = 5.0 * np.log10(np.maximum(DL, 1e-300)) + 25.0
    return mu


def percent_diff(a, b):
    return 100.0 * (np.asarray(a) - np.asarray(b)) / np.maximum(np.asarray(b), 1e-300)


# ============================================================
# Predições BAO
# ============================================================

def predict_bao_row(row, model="lcdm"):
    z = float(row["z"])
    obs_type = row["observable"]

    if obs_type == "DV_over_rd":
        return volume_distance_DV(z, model=model) / RD_MPC

    if obs_type == "DM_over_rd":
        return comoving_distance_DM(z, model=model) / RD_MPC

    if obs_type == "DH_over_rd":
        return hubble_distance_DH(z, model=model) / RD_MPC

    raise ValueError(f"Unknown observable: {obs_type}")


def evaluate_bao_table(model="lcdm"):
    rows = []

    for _, row in bao_data.iterrows():
        pred = predict_bao_row(row, model=model)
        pull = (pred - row["obs"]) / row["sigma"]
        chi2 = pull**2

        rows.append({
            "name": row["name"],
            "z": row["z"],
            "observable": row["observable"],
            "obs": row["obs"],
            "sigma": row["sigma"],
            "pred": pred,
            "pull": pull,
            "chi2": chi2,
            "model": model,
        })

    return pd.DataFrame(rows)


# ============================================================
# CMB compressed proxies
# ============================================================

def cmb_compressed_proxies(model="lcdm"):
    """
    Proxies comprimidos:

    lA = pi * DM(z*) / r_s*
    R  = sqrt(Omega_m) * H0 * DM(z*) / c

    Observação:
    isto não substitui uma likelihood Planck completa.
    É apenas um teste de preservação geométrica de fundo.
    """
    DM_star = comoving_distance_DM(Z_STAR, model=model, n_grid=20000)

    lA = np.pi * DM_star / RS_STAR_MPC
    R = np.sqrt(OMEGA_M) * H0 * DM_star / C_KM_S

    return {
        "DM_star_Mpc": DM_star,
        "lA": lA,
        "R": R,
    }


# ============================================================
# Rodar teste
# ============================================================

print("=" * 72)
print("TESTE 03 — COMPRESSED COSMOLOGY BACKGROUND + BAO")
print("=" * 72)

print("\nParâmetros:")
print("H0:", H0)
print("Omega_m:", OMEGA_M)
print("Omega_r:", OMEGA_R)
print("Omega_de:", OMEGA_DE)
print("TVGD lambda:", LAMBDA_TVGD)
print("TVGD p:", P_TVGD)
print("rd:", RD_MPC)

# Background
z_grid = np.linspace(0.001, Z_MAX_BACKGROUND, N_Z)

H_lcdm = H_of_z(z_grid, model="lcdm")
H_tvgd = H_of_z(z_grid, model="tvgd")
delta_H_percent = percent_diff(H_tvgd, H_lcdm)

DM_lcdm = np.array([comoving_distance_DM(z, model="lcdm", n_grid=4000) for z in z_grid])
DM_tvgd = np.array([comoving_distance_DM(z, model="tvgd", n_grid=4000) for z in z_grid])
delta_DM_percent = percent_diff(DM_tvgd, DM_lcdm)

DH_lcdm = hubble_distance_DH(z_grid, model="lcdm")
DH_tvgd = hubble_distance_DH(z_grid, model="tvgd")
delta_DH_percent = percent_diff(DH_tvgd, DH_lcdm)

mu_lcdm = distance_modulus(z_grid, model="lcdm", n_grid=4000)
mu_tvgd = distance_modulus(z_grid, model="tvgd", n_grid=4000)
delta_mu_mag = mu_tvgd - mu_lcdm

background_summary = {
    "max_abs_delta_H_percent": float(np.max(np.abs(delta_H_percent))),
    "mean_abs_delta_H_percent": float(np.mean(np.abs(delta_H_percent))),
    "max_abs_delta_DM_percent": float(np.max(np.abs(delta_DM_percent))),
    "mean_abs_delta_DM_percent": float(np.mean(np.abs(delta_DM_percent))),
    "max_abs_delta_DH_percent": float(np.max(np.abs(delta_DH_percent))),
    "mean_abs_delta_DH_percent": float(np.mean(np.abs(delta_DH_percent))),
    "max_abs_delta_mu_mag": float(np.max(np.abs(delta_mu_mag))),
    "mean_abs_delta_mu_mag": float(np.mean(np.abs(delta_mu_mag))),
}

print("\nResumo background:")
for k, v in background_summary.items():
    print(f"{k}: {v}")

background_df = pd.DataFrame({
    "z": z_grid,
    "H_LCDM": H_lcdm,
    "H_TVGD": H_tvgd,
    "delta_H_percent": delta_H_percent,
    "DM_LCDM_Mpc": DM_lcdm,
    "DM_TVGD_Mpc": DM_tvgd,
    "delta_DM_percent": delta_DM_percent,
    "DH_LCDM_Mpc": DH_lcdm,
    "DH_TVGD_Mpc": DH_tvgd,
    "delta_DH_percent": delta_DH_percent,
    "mu_LCDM": mu_lcdm,
    "mu_TVGD": mu_tvgd,
    "delta_mu_mag": delta_mu_mag,
})

# BAO
bao_lcdm = evaluate_bao_table(model="lcdm")
bao_tvgd = evaluate_bao_table(model="tvgd")

chi2_bao_lcdm = float(bao_lcdm["chi2"].sum())
chi2_bao_tvgd = float(bao_tvgd["chi2"].sum())
delta_chi2_bao = chi2_bao_tvgd - chi2_bao_lcdm

bao_combined = pd.concat([bao_lcdm, bao_tvgd], ignore_index=True)

print("\nBAO LCDM:")
display(bao_lcdm)

print("\nBAO TVGD:")
display(bao_tvgd)

print("\nChi2 BAO:")
print("chi2_bao_lcdm:", chi2_bao_lcdm)
print("chi2_bao_tvgd:", chi2_bao_tvgd)
print("delta_chi2_bao_TVGD_minus_LCDM:", delta_chi2_bao)

# CMB compressed proxies
cmb_lcdm = cmb_compressed_proxies(model="lcdm")
cmb_tvgd = cmb_compressed_proxies(model="tvgd")

delta_lA_percent = 100.0 * (cmb_tvgd["lA"] - cmb_lcdm["lA"]) / cmb_lcdm["lA"]
delta_R_percent = 100.0 * (cmb_tvgd["R"] - cmb_lcdm["R"]) / cmb_lcdm["R"]

cmb_summary = {
    "LCDM_DM_star_Mpc": cmb_lcdm["DM_star_Mpc"],
    "TVGD_DM_star_Mpc": cmb_tvgd["DM_star_Mpc"],
    "LCDM_lA": cmb_lcdm["lA"],
    "TVGD_lA": cmb_tvgd["lA"],
    "delta_lA_percent": delta_lA_percent,
    "LCDM_R": cmb_lcdm["R"],
    "TVGD_R": cmb_tvgd["R"],
    "delta_R_percent": delta_R_percent,
}

print("\nCMB compressed proxies:")
for k, v in cmb_summary.items():
    print(f"{k}: {v}")

# ============================================================
# Diagnóstico automático
# ============================================================

background_ok = (
    background_summary["mean_abs_delta_H_percent"] < MAX_MEAN_H_DIFF_PERCENT
    and background_summary["mean_abs_delta_DM_percent"] < MAX_MEAN_DISTANCE_DIFF_PERCENT
)

bao_strong = delta_chi2_bao < MAX_DELTA_CHI2_BAO_STRONG
bao_moderate = delta_chi2_bao < MAX_DELTA_CHI2_BAO_MODERATE

cmb_ok = (
    abs(delta_lA_percent) < MAX_DELTA_LA_PERCENT
    and abs(delta_R_percent) < MAX_DELTA_R_PERCENT
)

if background_ok and bao_strong and cmb_ok:
    final_status = "PASSOU_FORTE_COSMOLOGIA_COMPRESSA_BAO"
elif background_ok and bao_moderate and cmb_ok:
    final_status = "PASSOU_MODERADO_COSMOLOGIA_COMPRESSA_BAO"
elif background_ok and cmb_ok:
    final_status = "PASSOU_FRACO_BACKGROUND_CMB_PROXY"
else:
    final_status = "TENSIONADO_COSMOLOGIA_COMPRESSA"

summary = {
    "H0": H0,
    "Omega_m": OMEGA_M,
    "Omega_r": OMEGA_R,
    "Omega_de": OMEGA_DE,
    "lambda_TVGD": LAMBDA_TVGD,
    "p_TVGD": P_TVGD,
    "rd_Mpc": RD_MPC,
    "chi2_bao_lcdm": chi2_bao_lcdm,
    "chi2_bao_tvgd": chi2_bao_tvgd,
    "delta_chi2_bao_TVGD_minus_LCDM": delta_chi2_bao,
    "background_ok": bool(background_ok),
    "bao_strong": bool(bao_strong),
    "bao_moderate": bool(bao_moderate),
    "cmb_ok": bool(cmb_ok),
    "final_status": final_status,
}

summary.update(background_summary)
summary.update(cmb_summary)

summary_df = pd.DataFrame([summary])

print("\nResumo final:")
display(summary_df)

print("\nStatus final:", final_status)

# ============================================================
# Salvar resultados
# ============================================================

background_path = os.path.join(OUTDIR, "test03_background_comparison.csv")
bao_path = os.path.join(OUTDIR, "test03_bao_predictions.csv")
summary_path = os.path.join(OUTDIR, "test03_cosmology_background_bao_summary.csv")

background_df.to_csv(background_path, index=False)
bao_combined.to_csv(bao_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Gráficos
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(z_grid, H_lcdm, label="Lambda-CDM")
plt.plot(z_grid, H_tvgd, "--", label="TVGD V1.1")
plt.xlabel("z")
plt.ylabel("H(z) [km/s/Mpc]")
plt.title("Test 03 — Background expansion")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test03_Hz_comparison.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(z_grid, delta_H_percent)
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("z")
plt.ylabel("Delta H / H [%]")
plt.title("TVGD - Lambda-CDM H(z) difference")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test03_delta_H_percent.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(z_grid, delta_DM_percent, label="DM")
plt.plot(z_grid, delta_DH_percent, label="DH")
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("z")
plt.ylabel("Distance difference [%]")
plt.title("TVGD - Lambda-CDM distance differences")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test03_distance_differences.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(z_grid, delta_mu_mag)
plt.axhline(0, color="black", linewidth=1)
plt.xlabel("z")
plt.ylabel("Delta mu [mag]")
plt.title("TVGD - Lambda-CDM distance modulus difference")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test03_delta_mu.png"), dpi=220)
plt.show()

# BAO pulls
bao_plot = bao_combined.copy()
bao_plot["label"] = bao_plot["name"] + " z=" + bao_plot["z"].astype(str) + " " + bao_plot["observable"]

lcdm_pulls = bao_plot[bao_plot["model"] == "lcdm"]["pull"].values
tvgd_pulls = bao_plot[bao_plot["model"] == "tvgd"]["pull"].values
labels = bao_plot[bao_plot["model"] == "lcdm"]["label"].values

x = np.arange(len(labels))
width = 0.35

plt.figure(figsize=(11, 5))
plt.bar(x - width/2, lcdm_pulls, width, label="Lambda-CDM")
plt.bar(x + width/2, tvgd_pulls, width, label="TVGD")
plt.axhline(0, color="black", linewidth=1)
plt.axhline(2, color="gray", linestyle="--", linewidth=1)
plt.axhline(-2, color="gray", linestyle="--", linewidth=1)
plt.xticks(x, labels, rotation=45, ha="right")
plt.ylabel("Pull")
plt.title("BAO compressed pulls")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test03_bao_pulls.png"), dpi=220)
plt.show()

print("\nArquivos salvos em:")
print(OUTDIR)
print(background_path)
print(bao_path)
print(summary_path)

print("\nInterpretação curta:")
if "PASSOU_FORTE" in final_status:
    print("A TVGD V1.1 preserva a cosmologia comprimida de fundo e fica estatisticamente próxima de Lambda-CDM nos pontos BAO testados.")
elif "PASSOU_MODERADO" in final_status:
    print("A TVGD V1.1 preserva bem o fundo cosmológico e permanece próxima de Lambda-CDM, com tensão BAO moderada.")
elif "PASSOU_FRACO" in final_status:
    print("A TVGD V1.1 preserva o fundo e os proxies CMB, mas a comparação BAO precisa de refinamento.")
else:
    print("A cosmologia comprimida apresentou tensão relevante e deve ser revisada antes de subir ao GitHub.")

print("\nDONE.")
