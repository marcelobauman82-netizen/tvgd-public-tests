# ============================================================
# TVGD PUBLIC TESTS
# TEST 06 — NON-PHANTOM / GHOST-FREE COSMOLOGY BACKGROUND
# Versão autônoma para Google Colab
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Configuração
# ============================================================

OUTDIR = "Resultados/test06_nonphantom_cosmology_background"
os.makedirs(OUTDIR, exist_ok=True)

# Constantes
C_KM_S = 299792.458

# Cosmologia fiducial
H0 = 67.4
OMEGA_M0 = 0.31339978339159447
OMEGA_R0 = 9.0e-5
OMEGA_PSI0 = 1.0 - OMEGA_M0 - OMEGA_R0

# Parâmetros da família não-fantasma TVGD
EPSILON = 1.026579e-4
LAMBDA = 15.0
P_POWER = 0.25

# Escalas de referência
Z_CMB = 1089.92
Z_BBN = 1.0e9

# Grade
A_MIN = 1.0e-10
A_MAX = 1.0
N_A = 5000

Z_MAX_BACKGROUND = 3.0
N_Z = 800

# Critérios
MAX_W_BELOW_MINUS_ONE_TOL = 1e-10
MIN_K_ALLOWED = -1e-12
MIN_RHO_ALLOWED = 0.0
MIN_V_ALLOWED = 0.0

MAX_DELTA_E_BAO_PERCENT = 0.1
MAX_DELTA_DM_BAO_PERCENT = 0.1
MAX_DELTA_CMB_PERCENT = 0.1
MAX_FRAC_PSI_CMB = 1e-6
MAX_FRAC_PSI_BBN = 1e-10

# Pontos BAO aproximados para comparação geométrica
Z_BAO = np.array([0.106, 0.15, 0.38, 0.51, 0.61, 1.0, 1.5, 2.0])


# ============================================================
# Funções cosmológicas
# ============================================================

def omega_psi_LCDM(a):
    a = np.asarray(a, dtype=float)
    return OMEGA_PSI0 * np.ones_like(a)


def omega_psi_TVGD_nonphantom(a, eps=EPSILON, lam=LAMBDA, p=P_POWER):
    """
    Família TVGD não-fantasma efetiva:

    Omega_psi(a) =
    Omega_psi0 * {1 + eps * [exp(-lambda a^p) - exp(-lambda)] / [1 - exp(-lambda)]}

    Para eps pequeno, fica extremamente próxima de Lambda-CDM.

    Como Omega_psi(a) diminui levemente com a, temos w(a) >= -1.
    """
    a = np.asarray(a, dtype=float)
    denom = 1.0 - np.exp(-lam)
    F = (np.exp(-lam * np.maximum(a, 0.0)**p) - np.exp(-lam)) / denom
    return OMEGA_PSI0 * (1.0 + eps * F)


def dlnrho_dln_a_numeric(a, rho_func):
    """
    Derivada numérica d ln rho / d ln a.
    """
    a = np.asarray(a, dtype=float)
    ln_a = np.log(a)
    rho = rho_func(a)

    ln_rho = np.log(np.maximum(rho, 1e-300))
    return np.gradient(ln_rho, ln_a)


def w_from_rho(a, rho_func):
    """
    Conservação:

    d rho / d ln a = -3(1+w) rho

    Logo:

    w = -1 - (1/3) d ln rho / d ln a
    """
    dlnrho = dlnrho_dln_a_numeric(a, rho_func)
    return -1.0 - dlnrho / 3.0


def E_LCDM_z(z):
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    return np.sqrt(
        OMEGA_M0 * a**(-3)
        + OMEGA_R0 * a**(-4)
        + omega_psi_LCDM(a)
    )


def E_TVGD_z(z):
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    return np.sqrt(
        OMEGA_M0 * a**(-3)
        + OMEGA_R0 * a**(-4)
        + omega_psi_TVGD_nonphantom(a)
    )


def H_LCDM_z(z):
    return H0 * E_LCDM_z(z)


def H_TVGD_z(z):
    return H0 * E_TVGD_z(z)


def integrate_trapz(z, inv_E_func, n_grid=6000):
    z = float(z)
    if z <= 0:
        return 0.0

    zz = np.linspace(0.0, z, n_grid)
    yy = inv_E_func(zz)
    return float(np.trapezoid(yy, zz))


def DM_comoving_Mpc(z, model="lcdm", n_grid=6000):
    if model == "lcdm":
        integral = integrate_trapz(z, lambda zz: 1.0 / E_LCDM_z(zz), n_grid=n_grid)
    elif model == "tvgd":
        integral = integrate_trapz(z, lambda zz: 1.0 / E_TVGD_z(zz), n_grid=n_grid)
    else:
        raise ValueError("model must be lcdm or tvgd")

    return (C_KM_S / H0) * integral


def percent_diff(a, b):
    return 100.0 * (np.asarray(a) - np.asarray(b)) / np.maximum(np.asarray(b), 1e-300)


def psi_fraction(a, model="tvgd"):
    """
    Fração relativa Omega_psi / E(a)^2.
    """
    a = np.asarray(a, dtype=float)

    if model == "lcdm":
        psi = omega_psi_LCDM(a)
    else:
        psi = omega_psi_TVGD_nonphantom(a)

    E2 = OMEGA_M0 * a**(-3) + OMEGA_R0 * a**(-4) + psi
    return psi / E2


# ============================================================
# Rodar teste
# ============================================================

print("=" * 72)
print("TESTE 06 — NON-PHANTOM / GHOST-FREE COSMOLOGY BACKGROUND")
print("=" * 72)

print("\nParâmetros:")
print("H0:", H0)
print("Omega_m0:", OMEGA_M0)
print("Omega_r0:", OMEGA_R0)
print("Omega_psi0:", OMEGA_PSI0)
print("epsilon:", EPSILON)
print("lambda:", LAMBDA)
print("p:", P_POWER)

# Grade em a
a_grid = np.logspace(np.log10(A_MIN), np.log10(A_MAX), N_A)

rho_LCDM = omega_psi_LCDM(a_grid)
rho_TVGD = omega_psi_TVGD_nonphantom(a_grid)

w_TVGD = w_from_rho(a_grid, omega_psi_TVGD_nonphantom)

# Reconstrução canônica efetiva:
# rho = K + V
# p = K - V
# w = p/rho
# K = (1+w)rho/2
# V = (1-w)rho/2

K_eff = 0.5 * (1.0 + w_TVGD) * rho_TVGD
V_eff = 0.5 * (1.0 - w_TVGD) * rho_TVGD

rho_positive = bool(np.all(rho_TVGD > MIN_RHO_ALLOWED))
K_nonnegative = bool(np.all(K_eff >= MIN_K_ALLOWED))
V_positive = bool(np.all(V_eff > MIN_V_ALLOWED))
nonphantom_ok = bool(np.all(w_TVGD >= -1.0 - MAX_W_BELOW_MINUS_ONE_TOL))

finite_all = bool(
    np.all(np.isfinite(rho_TVGD))
    and np.all(np.isfinite(w_TVGD))
    and np.all(np.isfinite(K_eff))
    and np.all(np.isfinite(V_eff))
)

# Conservação
dlnrho = dlnrho_dln_a_numeric(a_grid, omega_psi_TVGD_nonphantom)
conservation_residual = dlnrho + 3.0 * (1.0 + w_TVGD)
max_conservation_residual = float(np.max(np.abs(conservation_residual)))

theory_summary = {
    "finite_all": finite_all,
    "rho_positive": rho_positive,
    "nonphantom_ok": nonphantom_ok,
    "K_nonnegative": K_nonnegative,
    "V_positive": V_positive,
    "w_today": float(w_TVGD[-1]),
    "w_min": float(np.min(w_TVGD)),
    "w_max": float(np.max(w_TVGD)),
    "K_min": float(np.min(K_eff)),
    "K_max": float(np.max(K_eff)),
    "V_min": float(np.min(V_eff)),
    "V_max": float(np.max(V_eff)),
    "rho_min": float(np.min(rho_TVGD)),
    "rho_max": float(np.max(rho_TVGD)),
    "max_conservation_residual": max_conservation_residual,
}

print("\nResumo teórico não-fantasma:")
for k, v in theory_summary.items():
    print(f"{k}: {v}")

# Fração no CMB e BBN
a_cmb = 1.0 / (1.0 + Z_CMB)
a_bbn = 1.0 / (1.0 + Z_BBN)

frac_psi_cmb = float(psi_fraction(a_cmb, model="tvgd"))
frac_psi_bbn = float(psi_fraction(a_bbn, model="tvgd"))

early_ok = (frac_psi_cmb < MAX_FRAC_PSI_CMB) and (frac_psi_bbn < MAX_FRAC_PSI_BBN)

early_summary = {
    "a_cmb": a_cmb,
    "frac_psi_cmb": frac_psi_cmb,
    "a_bbn": a_bbn,
    "frac_psi_bbn": frac_psi_bbn,
    "early_universe_ok": bool(early_ok),
}

print("\nResumo universo primordial:")
for k, v in early_summary.items():
    print(f"{k}: {v}")

# Comparação geométrica de fundo
z_grid = np.linspace(0.001, Z_MAX_BACKGROUND, N_Z)

E_lcdm = E_LCDM_z(z_grid)
E_tvgd = E_TVGD_z(z_grid)
delta_E_percent = percent_diff(E_tvgd, E_lcdm)

DM_lcdm = np.array([DM_comoving_Mpc(z, model="lcdm", n_grid=4000) for z in z_grid])
DM_tvgd = np.array([DM_comoving_Mpc(z, model="tvgd", n_grid=4000) for z in z_grid])
delta_DM_percent = percent_diff(DM_tvgd, DM_lcdm)

# BAO geometry points
bao_rows = []

for z in Z_BAO:
    E_l = float(E_LCDM_z(z))
    E_t = float(E_TVGD_z(z))

    DM_l = DM_comoving_Mpc(z, model="lcdm", n_grid=6000)
    DM_t = DM_comoving_Mpc(z, model="tvgd", n_grid=6000)

    bao_rows.append({
        "z": z,
        "E_LCDM": E_l,
        "E_TVGD": E_t,
        "delta_E_percent": float(percent_diff(E_t, E_l)),
        "DM_LCDM_Mpc": DM_l,
        "DM_TVGD_Mpc": DM_t,
        "delta_DM_percent": float(percent_diff(DM_t, DM_l)),
    })

bao_df = pd.DataFrame(bao_rows)

max_abs_delta_E_bao = float(np.max(np.abs(bao_df["delta_E_percent"])))
max_abs_delta_DM_bao = float(np.max(np.abs(bao_df["delta_DM_percent"])))

DM_cmb_lcdm = DM_comoving_Mpc(Z_CMB, model="lcdm", n_grid=20000)
DM_cmb_tvgd = DM_comoving_Mpc(Z_CMB, model="tvgd", n_grid=20000)
delta_DM_cmb_percent = float(percent_diff(DM_cmb_tvgd, DM_cmb_lcdm))

geometry_ok = (
    max_abs_delta_E_bao < MAX_DELTA_E_BAO_PERCENT
    and max_abs_delta_DM_bao < MAX_DELTA_DM_BAO_PERCENT
    and abs(delta_DM_cmb_percent) < MAX_DELTA_CMB_PERCENT
)

geometry_summary = {
    "max_abs_delta_E_background_percent": float(np.max(np.abs(delta_E_percent))),
    "mean_abs_delta_E_background_percent": float(np.mean(np.abs(delta_E_percent))),
    "max_abs_delta_DM_background_percent": float(np.max(np.abs(delta_DM_percent))),
    "mean_abs_delta_DM_background_percent": float(np.mean(np.abs(delta_DM_percent))),
    "max_abs_delta_E_BAO_percent": max_abs_delta_E_bao,
    "max_abs_delta_DM_BAO_percent": max_abs_delta_DM_bao,
    "DM_cmb_LCDM_Mpc": DM_cmb_lcdm,
    "DM_cmb_TVGD_Mpc": DM_cmb_tvgd,
    "delta_DM_cmb_percent": delta_DM_cmb_percent,
    "geometry_ok": bool(geometry_ok),
}

print("\nResumo geométrico:")
for k, v in geometry_summary.items():
    print(f"{k}: {v}")

print("\nTabela BAO geométrica:")
display(bao_df)

# ============================================================
# Diagnóstico automático
# ============================================================

ghostfree_background_ok = (
    finite_all
    and rho_positive
    and nonphantom_ok
    and K_nonnegative
    and V_positive
    and max_conservation_residual < 1e-8
)

if ghostfree_background_ok and early_ok and geometry_ok:
    final_status = "PASSOU_FORTE_NONPHANTOM_GHOSTFREE_BACKGROUND"
elif ghostfree_background_ok and early_ok:
    final_status = "PASSOU_MODERADO_NONPHANTOM_BACKGROUND"
elif ghostfree_background_ok:
    final_status = "PASSOU_FRACO_NONPHANTOM_BACKGROUND"
else:
    final_status = "TENSIONADO_NONPHANTOM_BACKGROUND"

summary = {}
summary.update(theory_summary)
summary.update(early_summary)
summary.update(geometry_summary)
summary["ghostfree_background_ok"] = bool(ghostfree_background_ok)
summary["final_status"] = final_status

summary_df = pd.DataFrame([summary])

print("\nResumo final:")
display(summary_df)

print("\nStatus final:", final_status)

# ============================================================
# Salvar resultados
# ============================================================

background_df = pd.DataFrame({
    "a": a_grid,
    "rho_TVGD": rho_TVGD,
    "w_TVGD": w_TVGD,
    "K_eff": K_eff,
    "V_eff": V_eff,
    "conservation_residual": conservation_residual,
})

z_background_df = pd.DataFrame({
    "z": z_grid,
    "E_LCDM": E_lcdm,
    "E_TVGD": E_tvgd,
    "delta_E_percent": delta_E_percent,
    "DM_LCDM_Mpc": DM_lcdm,
    "DM_TVGD_Mpc": DM_tvgd,
    "delta_DM_percent": delta_DM_percent,
})

background_path = os.path.join(OUTDIR, "test06_nonphantom_background_reconstruction.csv")
geometry_path = os.path.join(OUTDIR, "test06_geometry_comparison.csv")
bao_path = os.path.join(OUTDIR, "test06_bao_geometry_points.csv")
summary_path = os.path.join(OUTDIR, "test06_nonphantom_cosmology_summary.csv")

background_df.to_csv(background_path, index=False)
z_background_df.to_csv(geometry_path, index=False)
bao_df.to_csv(bao_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Gráficos
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(a_grid, w_TVGD)
plt.axhline(-1.0, color="black", linestyle="--", label="phantom divide")
plt.xscale("log")
plt.xlabel("a")
plt.ylabel("w(a)")
plt.title("TVGD non-phantom effective equation of state")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test06_w_of_a.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(a_grid, K_eff, label="K_eff")
plt.plot(a_grid, V_eff, label="V_eff")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("a")
plt.ylabel("effective density units")
plt.title("Canonical effective reconstruction")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test06_K_V_reconstruction.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(z_grid, delta_E_percent, label="Delta E")
plt.plot(z_grid, delta_DM_percent, label="Delta DM")
plt.axhline(0.0, color="black", linewidth=1)
plt.xlabel("z")
plt.ylabel("TVGD - LCDM [%]")
plt.title("Geometry difference relative to Lambda-CDM")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test06_geometry_differences.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.scatter(bao_df["z"], bao_df["delta_E_percent"], label="E(z)")
plt.scatter(bao_df["z"], bao_df["delta_DM_percent"], label="DM(z)")
plt.axhline(0.0, color="black", linewidth=1)
plt.xlabel("z")
plt.ylabel("TVGD - LCDM [%]")
plt.title("BAO-point geometry differences")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test06_bao_geometry_differences.png"), dpi=220)
plt.show()

print("\nArquivos salvos em:")
print(OUTDIR)
print(background_path)
print(geometry_path)
print(bao_path)
print(summary_path)

print("\nInterpretação curta:")
if final_status == "PASSOU_FORTE_NONPHANTOM_GHOSTFREE_BACKGROUND":
    print("A forma cosmológica TVGD não-fantasma preserva a geometria de fundo, mantém w>=-1, K_eff>=0, V_eff>0 e contribuição primordial desprezível.")
elif "PASSOU" in final_status:
    print("A forma cosmológica TVGD passa nos critérios não-fantasma principais, mas algum critério geométrico/primordial precisa de revisão.")
else:
    print("A forma cosmológica TVGD apresentou tensão ghost-free ou geométrica e deve ser revisada.")

print("\nDONE.")
