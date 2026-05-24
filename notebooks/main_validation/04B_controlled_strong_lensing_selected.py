# ============================================================
# TVGD PUBLIC TESTS
# TEST 04B — CONTROLLED STRONG-LENSING SELECTED SAMPLE
# Versão autônoma para Google Colab
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Configuração
# ============================================================

OUTDIR = "Resultados/test04B_controlled_strong_lensing_selected"
os.makedirs(OUTDIR, exist_ok=True)

SEED = 20260518
rng = np.random.default_rng(SEED)

N_LENSES = 300

# Constantes físicas
G = 6.67430e-11
C = 299792458.0
MSUN = 1.98847e30
KPC = 3.0856775814913673e19
MPC = 3.0856775814913673e22
ARCSEC_PER_RAD = 206265.0

# TVGD
A0 = 1.2e-10

# Cosmologia simples para distâncias angulares
H0 = 67.4
OMEGA_M = 0.315
OMEGA_L = 1.0 - OMEGA_M
C_KM_S = 299792.458

# Critérios
MIN_VALID_FRAC = 0.85
MIN_GAIN_FRAC = 0.90
MIN_MEDIAN_GAIN = 1.02
MAX_MEDIAN_GAIN = 1.35
MAX_THETA_P95_ARCSEC = 8.0
MIN_REALISTIC_FRAC = 0.65

# ============================================================
# Funções cosmológicas
# ============================================================

def E_z(z):
    z = np.asarray(z, dtype=float)
    return np.sqrt(OMEGA_M * (1.0 + z)**3 + OMEGA_L)


def comoving_distance_mpc(z, n_grid=3000):
    z = float(z)
    if z <= 0:
        return 0.0

    zz = np.linspace(0.0, z, n_grid)
    yy = 1.0 / E_z(zz)
    integral = np.trapezoid(yy, zz)

    return (C_KM_S / H0) * integral


def angular_diameter_distance_mpc(z):
    return comoving_distance_mpc(z) / (1.0 + z)


def angular_diameter_distance_between_mpc(z1, z2):
    if z2 <= z1:
        return np.nan

    dc1 = comoving_distance_mpc(z1)
    dc2 = comoving_distance_mpc(z2)

    return (dc2 - dc1) / (1.0 + z2)


# ============================================================
# Funções TVGD
# ============================================================

def P_tvgd_from_g(g_bar, a0=A0):
    g_bar = np.asarray(g_bar, dtype=float)
    u = np.sqrt(np.maximum(g_bar, 0.0) / a0)
    return 1.0 - np.exp(-u)


def alpha_tvgd_from_g(g_bar, a0=A0):
    P = P_tvgd_from_g(g_bar, a0=a0)
    return 1.0 / np.maximum(P, 1e-300)


# ============================================================
# Modelo controlado de lente
# ============================================================
#
# Massa bariônica acumulada suave:
#
# Mbar(<R) = Mbar_total * R^2 / (R + Re)^2
#
# O objetivo aqui não é substituir uma modelagem SLACS/HSC real.
# O objetivo é testar estabilidade e escala da amplificação TVGD
# em uma amostra selecionada para produzir lentes fortes.
# ============================================================

def Mbar_enclosed_kg(R_m, Mbar_total_kg, Re_m):
    R_m = np.asarray(R_m, dtype=float)
    return Mbar_total_kg * R_m**2 / np.maximum((R_m + Re_m)**2, 1e-300)


def gbar_at_R(R_m, Mbar_total_kg, Re_m):
    Menc = Mbar_enclosed_kg(R_m, Mbar_total_kg, Re_m)
    return G * Menc / np.maximum(R_m**2, 1e-300)


def Mdyn_tvgd_enclosed_kg(R_m, Mbar_total_kg, Re_m):
    Menc = Mbar_enclosed_kg(R_m, Mbar_total_kg, Re_m)
    gbar = gbar_at_R(R_m, Mbar_total_kg, Re_m)
    alpha = alpha_tvgd_from_g(gbar)
    return alpha * Menc


# ============================================================
# Equação do raio de Einstein
# ============================================================

def einstein_root_arcsec(Mbar_total_msun, Re_kpc, z_l, z_s, model="bar"):
    D_l_mpc = angular_diameter_distance_mpc(z_l)
    D_s_mpc = angular_diameter_distance_mpc(z_s)
    D_ls_mpc = angular_diameter_distance_between_mpc(z_l, z_s)

    if not np.isfinite(D_l_mpc) or not np.isfinite(D_s_mpc) or not np.isfinite(D_ls_mpc):
        return np.nan

    D_l = D_l_mpc * MPC
    D_s = D_s_mpc * MPC
    D_ls = D_ls_mpc * MPC

    Mbar_total_kg = Mbar_total_msun * MSUN
    Re_m = Re_kpc * KPC

    K_lens = (4.0 * G / C**2) * D_ls / (D_l * D_s)

    def f(theta_rad):
        R_m = D_l * theta_rad

        if model == "bar":
            Menc = Mbar_enclosed_kg(R_m, Mbar_total_kg, Re_m)
        elif model == "tvgd":
            Menc = Mdyn_tvgd_enclosed_kg(R_m, Mbar_total_kg, Re_m)
        else:
            raise ValueError("model must be 'bar' or 'tvgd'")

        return theta_rad**2 - K_lens * Menc

    # Intervalo amplo para evitar perder raízes válidas
    th_min = 0.001 / ARCSEC_PER_RAD
    th_max = 30.0 / ARCSEC_PER_RAD

    grid = np.logspace(np.log10(th_min), np.log10(th_max), 700)
    vals = np.array([f(t) for t in grid])

    roots = []

    for i in range(len(grid) - 1):
        if not np.isfinite(vals[i]) or not np.isfinite(vals[i+1]):
            continue

        if vals[i] == 0:
            roots.append(grid[i])

        if vals[i] * vals[i+1] < 0:
            a = grid[i]
            b = grid[i+1]
            fa = vals[i]
            fb = vals[i+1]

            for _ in range(100):
                m = 0.5 * (a + b)
                fm = f(m)

                if not np.isfinite(fm):
                    break

                if abs(fm) < 1e-30:
                    a = b = m
                    break

                if fa * fm <= 0:
                    b = m
                    fb = fm
                else:
                    a = m
                    fa = fm

            roots.append(0.5 * (a + b))

    if len(roots) == 0:
        return np.nan

    # Raiz externa principal
    theta_rad = max(roots)

    return theta_rad * ARCSEC_PER_RAD


# ============================================================
# Gerar amostra selecionada para lentes fortes
# ============================================================

print("=" * 72)
print("TESTE 04B — CONTROLLED STRONG-LENSING SELECTED SAMPLE")
print("=" * 72)

# Amostra mais adequada a lentes fortes:
# - massas maiores
# - raios efetivos menores/moderados
# - fontes suficientemente atrás da lente

logMbar = rng.uniform(11.55, 12.45, N_LENSES)
Mbar_msun = 10.0 ** logMbar

Re_kpc = 10.0 ** rng.uniform(np.log10(1.8), np.log10(8.0), N_LENSES)

z_l = rng.uniform(0.15, 0.55, N_LENSES)
z_s = z_l + rng.uniform(0.60, 2.00, N_LENSES)

rows = []

for i in range(N_LENSES):
    theta_bar = einstein_root_arcsec(
        Mbar_total_msun=Mbar_msun[i],
        Re_kpc=Re_kpc[i],
        z_l=z_l[i],
        z_s=z_s[i],
        model="bar",
    )

    theta_tvgd = einstein_root_arcsec(
        Mbar_total_msun=Mbar_msun[i],
        Re_kpc=Re_kpc[i],
        z_l=z_l[i],
        z_s=z_s[i],
        model="tvgd",
    )

    if np.isfinite(theta_tvgd):
        D_l_mpc = angular_diameter_distance_mpc(z_l[i])
        R_E_m = D_l_mpc * MPC * (theta_tvgd / ARCSEC_PER_RAD)
        Mkg = Mbar_msun[i] * MSUN
        Re_m = Re_kpc[i] * KPC

        g_RE = gbar_at_R(R_E_m, Mkg, Re_m)
        u_RE = np.sqrt(g_RE / A0)
        alpha_RE = alpha_tvgd_from_g(g_RE)
        R_E_kpc = R_E_m / KPC
    else:
        g_RE = np.nan
        u_RE = np.nan
        alpha_RE = np.nan
        R_E_kpc = np.nan

    if np.isfinite(theta_bar) and np.isfinite(theta_tvgd) and theta_bar > 0:
        gain = theta_tvgd / theta_bar
    else:
        gain = np.nan

    rows.append({
        "id": i,
        "logMbar": logMbar[i],
        "Mbar_msun": Mbar_msun[i],
        "Re_kpc": Re_kpc[i],
        "z_l": z_l[i],
        "z_s": z_s[i],
        "theta_bar_arcsec": theta_bar,
        "theta_tvgd_arcsec": theta_tvgd,
        "gain_theta": gain,
        "R_E_tvgd_kpc": R_E_kpc,
        "gbar_RE_m_s2": g_RE,
        "u_RE": u_RE,
        "alpha_RE": alpha_RE,
    })

df = pd.DataFrame(rows)

# ============================================================
# Métricas
# ============================================================

valid_bar = np.isfinite(df["theta_bar_arcsec"])
valid_tvgd = np.isfinite(df["theta_tvgd_arcsec"])
valid_both = valid_bar & valid_tvgd

frac_valid_bar = float(valid_bar.mean())
frac_valid_tvgd = float(valid_tvgd.mean())
frac_valid_both = float(valid_both.mean())

sub = df[valid_both].copy()

frac_tvgd_larger = float((sub["theta_tvgd_arcsec"] > sub["theta_bar_arcsec"]).mean())
median_gain = float(sub["gain_theta"].median())
p16_gain = float(sub["gain_theta"].quantile(0.16))
p84_gain = float(sub["gain_theta"].quantile(0.84))
p95_gain = float(sub["gain_theta"].quantile(0.95))

theta_bar_med = float(sub["theta_bar_arcsec"].median())
theta_tvgd_med = float(sub["theta_tvgd_arcsec"].median())

theta_bar_p16 = float(sub["theta_bar_arcsec"].quantile(0.16))
theta_bar_p84 = float(sub["theta_bar_arcsec"].quantile(0.84))
theta_tvgd_p16 = float(sub["theta_tvgd_arcsec"].quantile(0.16))
theta_tvgd_p84 = float(sub["theta_tvgd_arcsec"].quantile(0.84))
theta_tvgd_p95 = float(sub["theta_tvgd_arcsec"].quantile(0.95))

frac_bar_realistic = float(((sub["theta_bar_arcsec"] >= 0.5) & (sub["theta_bar_arcsec"] <= 3.5)).mean())
frac_tvgd_realistic = float(((sub["theta_tvgd_arcsec"] >= 0.5) & (sub["theta_tvgd_arcsec"] <= 3.5)).mean())

alpha_med = float(sub["alpha_RE"].median())
alpha_p16 = float(sub["alpha_RE"].quantile(0.16))
alpha_p84 = float(sub["alpha_RE"].quantile(0.84))
alpha_p90 = float(sub["alpha_RE"].quantile(0.90))

u_med = float(sub["u_RE"].median())
u_p16 = float(sub["u_RE"].quantile(0.16))
u_p84 = float(sub["u_RE"].quantile(0.84))

# Critérios
valid_ok = frac_valid_both >= MIN_VALID_FRAC
amplification_ok = frac_tvgd_larger >= MIN_GAIN_FRAC
gain_ok = (median_gain >= MIN_MEDIAN_GAIN) and (median_gain <= MAX_MEDIAN_GAIN)
realistic_ok = frac_tvgd_realistic >= MIN_REALISTIC_FRAC
no_explosion_ok = theta_tvgd_p95 < MAX_THETA_P95_ARCSEC

if valid_ok and amplification_ok and gain_ok and realistic_ok and no_explosion_ok:
    final_status = "PASSOU_FORTE_LENTES_CONTROLADAS_SELECIONADAS"
elif valid_ok and amplification_ok and no_explosion_ok:
    final_status = "PASSOU_MODERADO_LENTES_CONTROLADAS_SELECIONADAS"
elif valid_ok:
    final_status = "PASSOU_FRACO_LENTES_CONTROLADAS_SELECIONADAS"
else:
    final_status = "TENSIONADO_LENTES_CONTROLADAS_SELECIONADAS"

summary = {
    "N_lenses": N_LENSES,
    "frac_valid_bar": frac_valid_bar,
    "frac_valid_tvgd": frac_valid_tvgd,
    "frac_valid_both": frac_valid_both,

    "theta_bar_med_arcsec": theta_bar_med,
    "theta_bar_p16_arcsec": theta_bar_p16,
    "theta_bar_p84_arcsec": theta_bar_p84,

    "theta_tvgd_med_arcsec": theta_tvgd_med,
    "theta_tvgd_p16_arcsec": theta_tvgd_p16,
    "theta_tvgd_p84_arcsec": theta_tvgd_p84,
    "theta_tvgd_p95_arcsec": theta_tvgd_p95,

    "median_gain_theta": median_gain,
    "p16_gain_theta": p16_gain,
    "p84_gain_theta": p84_gain,
    "p95_gain_theta": p95_gain,

    "frac_tvgd_larger": frac_tvgd_larger,
    "frac_bar_realistic_0p5_3p5arcsec": frac_bar_realistic,
    "frac_tvgd_realistic_0p5_3p5arcsec": frac_tvgd_realistic,

    "alpha_RE_median": alpha_med,
    "alpha_RE_p16": alpha_p16,
    "alpha_RE_p84": alpha_p84,
    "alpha_RE_p90": alpha_p90,

    "u_RE_median": u_med,
    "u_RE_p16": u_p16,
    "u_RE_p84": u_p84,

    "valid_ok": valid_ok,
    "amplification_ok": amplification_ok,
    "gain_ok": gain_ok,
    "realistic_ok": realistic_ok,
    "no_explosion_ok": no_explosion_ok,

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nResumo:")
display(summary_df)

print("\nStatus final:", final_status)

print("\nAmostra dos resultados:")
display(df.head())

# ============================================================
# Salvar resultados
# ============================================================

results_path = os.path.join(OUTDIR, "test04B_controlled_lensing_selected_results.csv")
summary_path = os.path.join(OUTDIR, "test04B_controlled_lensing_selected_summary.csv")

df.to_csv(results_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Gráficos
# ============================================================

plt.figure(figsize=(7, 6))
plt.scatter(sub["theta_bar_arcsec"], sub["theta_tvgd_arcsec"], s=20, alpha=0.7)
xmax = max(sub["theta_bar_arcsec"].max(), sub["theta_tvgd_arcsec"].max())
plt.plot([0, xmax], [0, xmax], "--", label="equal")
plt.xlabel("Baryonic Einstein radius [arcsec]")
plt.ylabel("TVGD Einstein radius [arcsec]")
plt.title("Test 04B — Controlled selected strong lensing")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test04B_theta_bar_vs_tvgd.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(sub["gain_theta"], bins=35, alpha=0.8)
plt.axvline(1.0, color="black", linewidth=1, label="no gain")
plt.axvline(median_gain, linestyle="--", label="median gain")
plt.xlabel("theta_TVGD / theta_baryonic")
plt.ylabel("N")
plt.title("Einstein-radius gain distribution")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test04B_gain_distribution.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.scatter(sub["u_RE"], sub["gain_theta"], s=20, alpha=0.7)
plt.axhline(1.0, color="black", linewidth=1)
plt.xlabel("u at Einstein radius")
plt.ylabel("theta_TVGD / theta_baryonic")
plt.title("Gain versus TVGD regime")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test04B_gain_vs_u.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.scatter(sub["alpha_RE"], sub["gain_theta"], s=20, alpha=0.7)
plt.axhline(1.0, color="black", linewidth=1)
plt.xlabel("alpha_TVGD at Einstein radius")
plt.ylabel("theta_TVGD / theta_baryonic")
plt.title("Lensing gain versus TVGD amplification")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test04B_gain_vs_alpha.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.hist(sub["theta_tvgd_arcsec"], bins=35, alpha=0.8, label="TVGD")
plt.hist(sub["theta_bar_arcsec"], bins=35, alpha=0.5, label="Baryonic")
plt.xlabel("Einstein radius [arcsec]")
plt.ylabel("N")
plt.title("Einstein-radius distributions")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test04B_theta_distribution.png"), dpi=220)
plt.show()

print("\nArquivos salvos em:")
print(OUTDIR)
print(results_path)
print(summary_path)

print("\nInterpretação curta:")
if final_status == "PASSOU_FORTE_LENTES_CONTROLADAS_SELECIONADAS":
    print("A TVGD aumenta moderadamente o raio de Einstein em uma amostra controlada selecionada para lentes fortes, sem explosão angular e em faixa observacional plausível.")
elif "PASSOU" in final_status:
    print("A TVGD apresenta sinal estável de amplificação de lenteamento, mas algum critério forte ainda precisa ser refinado.")
else:
    print("O teste 04B apresentou tensão e deve ser revisado antes de subir ao GitHub.")

print("\nDONE.")
