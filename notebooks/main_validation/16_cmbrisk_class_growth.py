# ============================================================
# TVGD — 16_CMBRISK_CLASS_GROWTH
# CLASS completo + crescimento linear
#
# Objetivo:
#   Testar se a cosmologia efetiva TVGD V1.1 preserva:
#
#   1. Espectros CMB TT, EE, TE
#   2. Lenteamento CMB phi-phi
#   3. Crescimento linear D(z), f(z), f sigma8(z)
#
# Modelo TVGD V1.1:
#
#   rho_psi(a) ∝ [1 - exp(-lambda a^p)] / [1 - exp(-lambda)]
#
#   com lambda = 7.0, p = 0.55.
#
# O modelo é aproximado em CLASS por um fluido CPL efetivo:
#
#   w(a) = w0 + wa(1-a)
#
# Saída:
#   results/cmbrisk_class_growth/
#
# Uso:
#   python notebooks/main_validation/16_cmbrisk_class_growth.py
# ============================================================

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "16_cmbrisk_class_growth.py"

OUTDIR = "results/cmbrisk_class_growth"
os.makedirs(OUTDIR, exist_ok=True)

# Cosmologia base
H0 = 67.4
h = H0 / 100.0

omega_b = 0.02237
omega_cdm = 0.12000

Omega_b = omega_b / h**2
Omega_cdm = omega_cdm / h**2
Omega_m = Omega_b + Omega_cdm
Omega_r = 9.0e-5
Omega_de = 1.0 - Omega_m - Omega_r

A_s = 2.100e-9
n_s = 0.965
tau_reio = 0.054

# TVGD V1.1
LAMBDA_TVGD = 7.0
P_TVGD = 0.55

# CLASS
L_MAX = 2500

# Crescimento
Z_GROWTH = np.array([0.0, 0.15, 0.25, 0.38, 0.51, 0.61, 0.76, 1.0, 1.5, 2.0])

print("=" * 72)
print("TVGD — 16_CMBRISK_CLASS_GROWTH")
print("CLASS completo + crescimento linear")
print("=" * 72)


# ============================================================
# Instalar CLASS/classy se necessário
# ============================================================

try:
    from classy import Class
    print("classy já disponível.")
except Exception:
    print("classy não encontrado. Tentando instalar...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "classy"])
    from classy import Class
    print("classy instalado com sucesso.")


# ============================================================
# Funções TVGD background
# ============================================================

def A_tvgd(a, lam=LAMBDA_TVGD, p=P_TVGD):
    """
    Fator de densidade estrutural normalizado:
        A(a) = [1 - exp(-lam a^p)] / [1 - exp(-lam)]
    """
    a = np.asarray(a)
    return (1.0 - np.exp(-lam * np.power(a, p))) / (1.0 - np.exp(-lam))


def dlnA_dlna(a, lam=LAMBDA_TVGD, p=P_TVGD):
    """
    d ln A / d ln a.
    """
    a = np.asarray(a)
    x = lam * np.power(a, p)
    numerator = p * x * np.exp(-x)
    denominator = 1.0 - np.exp(-x)
    return numerator / np.maximum(denominator, 1e-300)


def w_tvgd(a, lam=LAMBDA_TVGD, p=P_TVGD):
    """
    Da conservação:
        d ln rho / d ln a = -3(1+w)
    logo:
        w(a) = -1 - (1/3) d ln A / d ln a
    """
    return -1.0 - dlnA_dlna(a, lam, p) / 3.0


def fit_cpl(a_min=0.70, a_max=1.0, n=400):
    """
    Ajusta w(a)=w0+wa(1-a) no intervalo tardio.
    """
    a = np.linspace(a_min, a_max, n)
    w = w_tvgd(a)

    X = np.vstack([np.ones_like(a), 1.0 - a]).T
    coef, *_ = np.linalg.lstsq(X, w, rcond=None)

    w0, wa = coef
    w_fit = w0 + wa * (1.0 - a)

    rms = np.sqrt(np.mean((w - w_fit)**2))
    maxerr = np.max(np.abs(w - w_fit))

    return float(w0), float(wa), float(rms), float(maxerr)


w0, wa, rms_w, max_w = fit_cpl()

print("\nParâmetros:")
print(f"H0 = {H0}")
print(f"h = {h}")
print(f"omega_b = {omega_b}")
print(f"omega_cdm = {omega_cdm}")
print(f"Omega_b = {Omega_b}")
print(f"Omega_cdm = {Omega_cdm}")
print(f"Omega_m = {Omega_m}")
print(f"Omega_de = {Omega_de}")
print(f"lambda = {LAMBDA_TVGD}")
print(f"p = {P_TVGD}")

print("\nAjuste CPL efetivo TVGD V1.1:")
print(f"w0 = {w0}")
print(f"wa = {wa}")
print(f"rms_w = {rms_w}")
print(f"max_w = {max_w}")


# ============================================================
# Rodar CLASS
# ============================================================

def base_class_params():
    return {
        "output": "tCl,pCl,lCl,mPk",
        "lensing": "yes",
        "l_max_scalars": L_MAX,
        "P_k_max_h/Mpc": 10.0,

        "h": h,
        "omega_b": omega_b,
        "omega_cdm": omega_cdm,

        "A_s": A_s,
        "n_s": n_s,
        "tau_reio": tau_reio,

        "N_ur": 3.046,
        "N_ncdm": 0,

        "YHe": 0.2454,
        "reio_parametrization": "reio_camb",

        "gauge": "newtonian",
    }


def run_class_lcdm():
    params = base_class_params()
    params.update({
        # Importante: em LCDM, deixamos CLASS inferir Omega_Lambda.
        # Não passamos Omega_Lambda junto com Omega_fld.
    })

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    return cosmo


def run_class_tvgd_cpl():
    params = base_class_params()
    params.update({
        "Omega_fld": Omega_de,
        "fluid_equation_of_state": "CLP",
        "w0_fld": w0,
        "wa_fld": wa,
        "cs2_fld": 1.0,
        # Não especificar Omega_Lambda quando Omega_fld é usado.
    })

    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    return cosmo


print("\nRodando CLASS LCDM...")
lcdm = run_class_lcdm()
print("LCDM OK.")

print("\nRodando CLASS TVGD-CPL...")
tvgd = run_class_tvgd_cpl()
print("TVGD-CPL OK.")


# ============================================================
# Espectros CMB
# ============================================================

cl_lcdm = lcdm.lensed_cl(L_MAX)
cl_tvgd = tvgd.lensed_cl(L_MAX)

ell = cl_lcdm["ell"]

# Converter para D_ell = ell(ell+1)C_ell/2pi para TT/EE/TE.
fac = ell * (ell + 1.0) / (2.0 * np.pi)

spectra = {}

for key in ["tt", "ee", "te"]:
    y0 = fac * cl_lcdm[key]
    y1 = fac * cl_tvgd[key]
    spectra[key.upper()] = (ell, y0, y1)

# phi-phi vem em pp.
pp0 = cl_lcdm["pp"]
pp1 = cl_tvgd["pp"]
spectra["phiphi"] = (ell, pp0, pp1)


def compare_spectrum(name, x, y0, y1, ell_min=30, ell_max=L_MAX, safe_floor_frac=1e-4):
    mask = (x >= ell_min) & (x <= ell_max) & np.isfinite(y0) & np.isfinite(y1)

    if name == "TE":
        # TE cru cruza zero. Usamos piso relativo ao pico absoluto para não explodir fração.
        scale = np.nanmax(np.abs(y0[mask]))
        denom = np.maximum(np.abs(y0[mask]), safe_floor_frac * scale)
    else:
        denom = np.maximum(np.abs(y0[mask]), 1e-300)

    frac = (y1[mask] - y0[mask]) / denom

    corr = np.corrcoef(y0[mask], y1[mask])[0, 1]

    return {
        f"{name}_N": int(mask.sum()),
        f"{name}_frac_rms_pct": float(100.0 * np.sqrt(np.mean(frac**2))),
        f"{name}_frac_mean_abs_pct": float(100.0 * np.mean(np.abs(frac))),
        f"{name}_frac_max_abs_pct": float(100.0 * np.max(np.abs(frac))),
        f"{name}_corr": float(corr),
    }


metrics = {}

metrics.update(compare_spectrum("TT", *spectra["TT"], ell_min=30, ell_max=L_MAX))
metrics.update(compare_spectrum("EE", *spectra["EE"], ell_min=30, ell_max=L_MAX))
metrics.update(compare_spectrum("TE", *spectra["TE"], ell_min=30, ell_max=L_MAX, safe_floor_frac=1e-3))
metrics.update(compare_spectrum("phiphi", *spectra["phiphi"], ell_min=10, ell_max=1000))

print("\nMétricas CMB:")
for k, v in metrics.items():
    print(f"{k} = {v}")


# ============================================================
# Crescimento linear
# ============================================================

def E2_lcdm_a(a):
    return Omega_m / a**3 + Omega_r / a**4 + (1.0 - Omega_m - Omega_r)


def E2_tvgd_a(a):
    return Omega_m / a**3 + Omega_r / a**4 + Omega_de * A_tvgd(a)


def growth_ode_integrate(E2_func, a_ini=1e-4, a_fin=1.0, n=8000):
    """
    Integra crescimento linear:
        D'' + [3/a + dlnH/da]D' - 3Omega_m/(2a^5E^2)D = 0

    Condição inicial em era de matéria:
        D(a_ini)=a_ini
        D'(a_ini)=1
    """
    a = np.geomspace(a_ini, a_fin, n)
    D = np.zeros_like(a)
    F = np.zeros_like(a)  # dD/da

    D[0] = a[0]
    F[0] = 1.0

    def deriv(ai, Di, Fi):
        # derivada numérica de lnH
        eps = 1e-4
        ap = ai * (1 + eps)
        am = ai * (1 - eps)

        E2p = E2_func(ap)
        E2m = E2_func(am)

        dlnH_da = 0.5 * (np.log(E2p) - np.log(E2m)) / (ap - am)

        E2 = E2_func(ai)

        dD_da = Fi
        dF_da = - (3.0/ai + dlnH_da) * Fi + 1.5 * Omega_m * Di / (ai**5 * E2)

        return dD_da, dF_da

    for i in range(n - 1):
        ai = a[i]
        hstep = a[i+1] - a[i]

        # RK4
        k1D, k1F = deriv(ai, D[i], F[i])
        k2D, k2F = deriv(ai + 0.5*hstep, D[i] + 0.5*hstep*k1D, F[i] + 0.5*hstep*k1F)
        k3D, k3F = deriv(ai + 0.5*hstep, D[i] + 0.5*hstep*k2D, F[i] + 0.5*hstep*k2F)
        k4D, k4F = deriv(ai + hstep, D[i] + hstep*k3D, F[i] + hstep*k3F)

        D[i+1] = D[i] + hstep * (k1D + 2*k2D + 2*k3D + k4D) / 6.0
        F[i+1] = F[i] + hstep * (k1F + 2*k2F + 2*k3F + k4F) / 6.0

    # Normalizar D(a=1)=1.
    D_norm = D / D[-1]

    # f = dlnD/dlna = a D'/D
    F_norm = F / D[-1]
    f = a * F_norm / np.maximum(D_norm, 1e-300)

    return a, D_norm, f


a_lcdm, D_lcdm, f_lcdm = growth_ode_integrate(E2_lcdm_a)
a_tvgd, D_tvgd, f_tvgd = growth_ode_integrate(E2_tvgd_a)

def interp_at_z(a_grid, y_grid, z):
    a = 1.0 / (1.0 + z)
    return float(np.interp(a, a_grid, y_grid))

growth_rows = []

sigma8_ref = 0.811

for z in Z_GROWTH:
    Dl = interp_at_z(a_lcdm, D_lcdm, z)
    Dt = interp_at_z(a_tvgd, D_tvgd, z)

    fl = interp_at_z(a_lcdm, f_lcdm, z)
    ft = interp_at_z(a_tvgd, f_tvgd, z)

    fs8_l = fl * sigma8_ref * Dl
    fs8_t = ft * sigma8_ref * Dt

    growth_rows.append({
        "z": float(z),
        "a": float(1.0 / (1.0 + z)),
        "D_LCDM": Dl,
        "D_TVGD": Dt,
        "frac_D_pct": 100.0 * (Dt - Dl) / Dl,
        "f_LCDM": fl,
        "f_TVGD": ft,
        "frac_f_pct": 100.0 * (ft - fl) / fl,
        "fsigma8_LCDM": fs8_l,
        "fsigma8_TVGD": fs8_t,
        "frac_fsigma8_pct": 100.0 * (fs8_t - fs8_l) / fs8_l,
    })

growth_df = pd.DataFrame(growth_rows)

growth_metrics = {
    "growth_max_abs_frac_D_pct": float(np.max(np.abs(growth_df["frac_D_pct"]))),
    "growth_max_abs_frac_f_pct": float(np.max(np.abs(growth_df["frac_f_pct"]))),
    "growth_max_abs_frac_fsigma8_pct": float(np.max(np.abs(growth_df["frac_fsigma8_pct"]))),
}

print("\nCrescimento linear:")
print(growth_df.to_string(index=False))

print("\nMétricas crescimento:")
for k, v in growth_metrics.items():
    print(f"{k} = {v}")


# ============================================================
# Salvar tabelas
# ============================================================

cmb_rows = []

for name, (x, y0, y1) in spectra.items():
    for i in range(len(x)):
        cmb_rows.append({
            "spectrum": name,
            "ell": int(x[i]),
            "LCDM": float(y0[i]),
            "TVGD_CPL": float(y1[i]),
            "diff": float(y1[i] - y0[i]),
        })

cmb_df = pd.DataFrame(cmb_rows)
cmb_df.to_csv(os.path.join(OUTDIR, "CMBRISK16_cmb_spectra.csv"), index=False)

growth_df.to_csv(os.path.join(OUTDIR, "CMBRISK16_growth_table.csv"), index=False)


# ============================================================
# Diagnóstico final
# ============================================================

summary = {
    "FINAL_STATUS": None,
    "score": None,
    "script": SCRIPT_NAME,

    "H0": H0,
    "Omega_m": float(Omega_m),
    "Omega_de": float(Omega_de),
    "lambda": LAMBDA_TVGD,
    "p": P_TVGD,
    "w0_CPL": w0,
    "wa_CPL": wa,
    "rms_w_CPL": rms_w,
    "max_w_CPL": max_w,
}

summary.update(metrics)
summary.update(growth_metrics)

criteria = {
    "CPL_fit_good": rms_w < 1e-3 and max_w < 1e-2,

    "TT_rms_below_1pct": summary["TT_frac_rms_pct"] < 1.0,
    "TT_max_below_3pct": summary["TT_frac_max_abs_pct"] < 3.0,

    "EE_rms_below_1pct": summary["EE_frac_rms_pct"] < 1.0,
    "EE_max_below_5pct": summary["EE_frac_max_abs_pct"] < 5.0,

    # TE é sensível a cruzamentos por zero. Aceitamos correlação extremamente alta
    # mesmo se fração local explodir.
    "TE_corr_high": summary["TE_corr"] > 0.9999,
    "TE_mean_abs_below_2pct": summary["TE_frac_mean_abs_pct"] < 2.0,

    "phiphi_rms_below_1pct": summary["phiphi_frac_rms_pct"] < 1.0,
    "phiphi_max_below_3pct": summary["phiphi_frac_max_abs_pct"] < 3.0,

    "growth_D_below_1pct": summary["growth_max_abs_frac_D_pct"] < 1.0,
    "growth_f_below_1pct": summary["growth_max_abs_frac_f_pct"] < 1.0,
    "growth_fsigma8_below_1pct": summary["growth_max_abs_frac_fsigma8_pct"] < 1.0,
}

score = sum(bool(v) for v in criteria.values())

if score >= 11:
    status = "PASSOU_FORTE_CMBRISK16_CLASS_GROWTH"
elif score >= 9:
    status = "PASSOU_MODERADO_FORTE_CMBRISK16_COM_ALERTAS"
elif score >= 7:
    status = "PASSOU_MODERADO_CMBRISK16"
else:
    status = "FALHOU_CMBRISK16_CMB_OU_GROWTH"

summary["FINAL_STATUS"] = status
summary["score"] = int(score)
summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "CMBRISK16_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — CMBRISK 16")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

def plot_spectrum(name, filename, yscale="linear"):
    x, y0, y1 = spectra[name]

    plt.figure(figsize=(9, 5))
    plt.plot(x, y0, label="LCDM")
    plt.plot(x, y1, linestyle="--", label="TVGD-CPL")
    plt.xlabel(r"$\ell$")
    plt.ylabel(name)
    plt.title(f"CMBRISK 16 — {name}")
    if yscale == "log":
        plt.yscale("log")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, filename), dpi=170)
    plt.close()

plot_spectrum("TT", "CMBRISK16_TT.png")
plot_spectrum("EE", "CMBRISK16_EE.png")
plot_spectrum("TE", "CMBRISK16_TE.png")
plot_spectrum("phiphi", "CMBRISK16_phiphi.png", yscale="log")

plt.figure(figsize=(8, 5))
plt.plot(growth_df["z"], growth_df["frac_D_pct"], marker="o", label="D")
plt.plot(growth_df["z"], growth_df["frac_f_pct"], marker="s", label="f")
plt.plot(growth_df["z"], growth_df["frac_fsigma8_pct"], marker="^", label=r"$f\sigma_8$")
plt.axhline(0.0, linestyle="--")
plt.xlabel("z")
plt.ylabel("TVGD - LCDM [%]")
plt.title("CMBRISK 16 — Linear growth fractional difference")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "CMBRISK16_growth_fractional_difference.png"), dpi=170)
plt.close()

a_plot = np.linspace(0.05, 1.0, 500)
plt.figure(figsize=(8, 5))
plt.plot(a_plot, w_tvgd(a_plot), label="TVGD exact effective w(a)")
plt.plot(a_plot, w0 + wa*(1-a_plot), linestyle="--", label="CPL fit")
plt.xlabel("a")
plt.ylabel("w(a)")
plt.title("CMBRISK 16 — TVGD effective equation of state")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "CMBRISK16_w_of_a_CPL_fit.png"), dpi=170)
plt.close()


# ============================================================
# Encerrar CLASS
# ============================================================

lcdm.struct_cleanup()
lcdm.empty()

tvgd.struct_cleanup()
tvgd.empty()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste CMBRISK 16.")
