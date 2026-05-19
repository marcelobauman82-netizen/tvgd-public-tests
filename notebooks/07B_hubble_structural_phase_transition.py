# ============================================================
# TVGD PUBLIC TESTS
# TEST 07B — HUBBLE TENSION FROM TVGD STRUCTURAL PHASE TRANSITION
# ============================================================
#
# Goal:
# Test whether the TVGD structural transition
#
#     P(Theta) = 1 - exp(-Theta)
#
# can generate a local H0-like enhancement while suppressing the
# correction before BAO/CMB scales.
#
# This test does not use supernova data. It is a structural/viability
# test of the TVGD phase-transition mechanism for the Hubble sector.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test07B_hubble_structural_phase_transition"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Reference values
# ============================================================

H0_CMB = 67.4
H0_LOCAL = 73.04

A_H = H0_LOCAL / H0_CMB - 1.0

Z_CMB = 1089.92
Z_BAO_MIN = 0.106

# Natural TVGD structural exponent:
#
# rho_m(z) ~ (1+z)^3
# structural amplitude ~ sqrt(rho_m)
# therefore DeltaTheta ~ (1+z)^(3/2)
S_NATURAL = 1.5

# Theta0 scan
THETA0_GRID = np.linspace(0.1, 80.0, 4000)

# Plot grid
Z_GRID = np.linspace(0.0, 0.30, 1000)

# Criteria
MAX_BAO_SHIFT_PERCENT = 0.5
MAX_CMB_SHIFT_ABS = 1e-12
MIN_SHIFT_Z001_PERCENT = 4.0
MIN_SHIFT_Z003_PERCENT = 1.0

# ============================================================
# TVGD structural transition functions
# ============================================================

def P_phase(delta_theta):
    """
    TVGD structural transition:

        P(Theta) = 1 - exp(-Theta)
    """
    delta_theta = np.asarray(delta_theta, dtype=float)
    return 1.0 - np.exp(-np.maximum(delta_theta, 0.0))


def relaxed_fraction(delta_theta):
    """
    Relaxed / non-activated fraction:

        1 - P(Theta) = exp(-Theta)
    """
    delta_theta = np.asarray(delta_theta, dtype=float)
    return np.exp(-np.maximum(delta_theta, 0.0))


def delta_theta_z(z, theta0, s=S_NATURAL):
    """
    Structural transition variable:

        DeltaTheta(z) = Theta0 * [(1+z)^s - 1]

    For s = 3/2 this follows the square-root amplitude of
    the mean matter-density scaling.
    """
    z = np.asarray(z, dtype=float)
    return theta0 * ((1.0 + np.maximum(z, 0.0))**s - 1.0)


def Xi_TVGD(z, theta0, s=S_NATURAL):
    """
    Effective Hubble inference factor:

        Xi_TVGD(z) = 1 + A_H * exp[-DeltaTheta(z)]

    At z = 0:
        DeltaTheta = 0
        Xi = 1 + A_H = H0_LOCAL / H0_CMB

    At high concentration / early universe:
        DeltaTheta >> 1
        Xi -> 1
    """
    dth = delta_theta_z(z, theta0=theta0, s=s)
    return 1.0 + A_H * relaxed_fraction(dth)


def H0_eff(z, theta0, s=S_NATURAL):
    return H0_CMB * Xi_TVGD(z, theta0=theta0, s=s)


def shift_percent(z, theta0, s=S_NATURAL):
    return 100.0 * (Xi_TVGD(z, theta0=theta0, s=s) - 1.0)


# ============================================================
# Run test
# ============================================================

print("=" * 72)
print("TEST 07B — HUBBLE TENSION FROM TVGD STRUCTURAL PHASE TRANSITION")
print("=" * 72)

print("\nReference values:")
print("H0_CMB:", H0_CMB)
print("H0_LOCAL:", H0_LOCAL)
print("A_H = H0_LOCAL/H0_CMB - 1:", A_H)
print("Required local shift [%]:", 100.0 * A_H)
print("Natural structural exponent s:", S_NATURAL)

rows = []

for theta0 in THETA0_GRID:
    H0_z0 = float(H0_eff(0.0, theta0))
    shift_z0 = float(shift_percent(0.0, theta0))
    shift_z001 = float(shift_percent(0.01, theta0))
    shift_z003 = float(shift_percent(0.03, theta0))
    shift_z005 = float(shift_percent(0.05, theta0))
    shift_z010 = float(shift_percent(0.10, theta0))

    # BAO region
    z_bao_grid = np.linspace(Z_BAO_MIN, 2.0, 800)
    bao_shift_max = float(np.max(np.abs(shift_percent(z_bao_grid, theta0))))

    # CMB
    cmb_shift_abs = float(abs(Xi_TVGD(Z_CMB, theta0) - 1.0))

    # Structural phase values
    dtheta_z001 = float(delta_theta_z(0.01, theta0))
    dtheta_z003 = float(delta_theta_z(0.03, theta0))
    dtheta_bao = float(delta_theta_z(Z_BAO_MIN, theta0))
    dtheta_cmb = float(delta_theta_z(Z_CMB, theta0))

    P_z001 = float(P_phase(dtheta_z001))
    P_z003 = float(P_phase(dtheta_z003))
    P_bao = float(P_phase(dtheta_bao))
    P_cmb = float(P_phase(dtheta_cmb))

    bao_ok = bao_shift_max < MAX_BAO_SHIFT_PERCENT
    cmb_ok = cmb_shift_abs < MAX_CMB_SHIFT_ABS

    local_signature_ok = (
        shift_z001 > MIN_SHIFT_Z001_PERCENT
        and shift_z003 > MIN_SHIFT_Z003_PERCENT
    )

    if bao_ok and cmb_ok and local_signature_ok:
        status = "PASSA_TRANSICAO_TVGD_H0_LOCAL_SEM_BAO_CMB"
    elif bao_ok and cmb_ok:
        status = "PASSA_BAO_CMB_MAS_ASSINATURA_LOCAL_FRACA"
    elif cmb_ok:
        status = "PRESERVA_CMB_MAS_TENSIONA_BAO"
    else:
        status = "FALHA_CMB"

    rows.append({
        "theta0": theta0,
        "H0_eff_z0": H0_z0,
        "shift_z0_percent": shift_z0,
        "shift_z0p01_percent": shift_z001,
        "shift_z0p03_percent": shift_z003,
        "shift_z0p05_percent": shift_z005,
        "shift_z0p10_percent": shift_z010,
        "bao_shift_max_z0p106_to_2_percent": bao_shift_max,
        "cmb_shift_abs": cmb_shift_abs,
        "DeltaTheta_z0p01": dtheta_z001,
        "DeltaTheta_z0p03": dtheta_z003,
        "DeltaTheta_BAO_z0p106": dtheta_bao,
        "DeltaTheta_CMB": dtheta_cmb,
        "P_z0p01": P_z001,
        "P_z0p03": P_z003,
        "P_BAO_z0p106": P_bao,
        "P_CMB": P_cmb,
        "bao_ok": bao_ok,
        "cmb_ok": cmb_ok,
        "local_signature_ok": local_signature_ok,
        "status": status,
    })

scan_df = pd.DataFrame(rows)

pass_df = scan_df[
    scan_df["status"] == "PASSA_TRANSICAO_TVGD_H0_LOCAL_SEM_BAO_CMB"
].copy()

N_total = len(scan_df)
N_pass = len(pass_df)

print("\nScan summary:")
print("N_total:", N_total)
print("N_pass:", N_pass)

if N_pass > 0:
    # Best case:
    # preserve as much local signal as possible,
    # while penalizing residual BAO leakage.
    pass_df["score"] = (
        pass_df["shift_z0p01_percent"]
        + 0.5 * pass_df["shift_z0p03_percent"]
        - 10.0 * pass_df["bao_shift_max_z0p106_to_2_percent"]
    )

    best = pass_df.sort_values("score", ascending=False).iloc[0]
    final_status = "PASSOU_TRANSICAO_ESTRUTURAL_TVGD_HUBBLE"
else:
    best = scan_df.sort_values("bao_shift_max_z0p106_to_2_percent").iloc[0]
    final_status = "TENSIONADO_TRANSICAO_ESTRUTURAL_TVGD_HUBBLE"

print("\nBest candidate:")
print(best)

summary = {
    "H0_CMB": H0_CMB,
    "H0_LOCAL": H0_LOCAL,
    "A_H": A_H,
    "required_local_shift_percent": 100.0 * A_H,
    "s_natural": S_NATURAL,
    "N_total": N_total,
    "N_pass": N_pass,
    "best_theta0": float(best["theta0"]),
    "best_H0_eff_z0": float(best["H0_eff_z0"]),
    "best_shift_z0_percent": float(best["shift_z0_percent"]),
    "best_shift_z0p01_percent": float(best["shift_z0p01_percent"]),
    "best_shift_z0p03_percent": float(best["shift_z0p03_percent"]),
    "best_shift_z0p05_percent": float(best["shift_z0p05_percent"]),
    "best_shift_z0p10_percent": float(best["shift_z0p10_percent"]),
    "best_bao_shift_max_percent": float(best["bao_shift_max_z0p106_to_2_percent"]),
    "best_cmb_shift_abs": float(best["cmb_shift_abs"]),
    "best_DeltaTheta_z0p01": float(best["DeltaTheta_z0p01"]),
    "best_DeltaTheta_z0p03": float(best["DeltaTheta_z0p03"]),
    "best_DeltaTheta_BAO_z0p106": float(best["DeltaTheta_BAO_z0p106"]),
    "best_DeltaTheta_CMB": float(best["DeltaTheta_CMB"]),
    "best_P_z0p01": float(best["P_z0p01"]),
    "best_P_z0p03": float(best["P_z0p03"]),
    "best_P_BAO_z0p106": float(best["P_BAO_z0p106"]),
    "best_P_CMB": float(best["P_CMB"]),
    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nFinal summary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

scan_path = os.path.join(OUTDIR, "test07B_hubble_structural_phase_scan.csv")
summary_path = os.path.join(OUTDIR, "test07B_hubble_structural_phase_summary.csv")

scan_df.to_csv(scan_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

theta0_best = float(best["theta0"])

Xi_grid = Xi_TVGD(Z_GRID, theta0_best)
H0_grid = H0_CMB * Xi_grid
shift_grid = 100.0 * (Xi_grid - 1.0)
DeltaTheta_grid = delta_theta_z(Z_GRID, theta0_best)
P_grid = P_phase(DeltaTheta_grid)
relaxed_grid = relaxed_fraction(DeltaTheta_grid)

plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, H0_grid, label="H0_eff(z)")
plt.axhline(H0_CMB, color="black", linestyle="--", label="CMB-like H0")
plt.axhline(H0_LOCAL, color="gray", linestyle=":", label="Local H0 target")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("Effective inferred H0 [km/s/Mpc]")
plt.title("TVGD structural transition: H0 effective")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07B_H0_eff_transition.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, shift_grid, label="local H0 shift")
plt.axhline(0.0, color="black", linewidth=1)
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO tolerance")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("Xi_TVGD(z)-1 [%]")
plt.title("Decay of local TVGD Hubble correction")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07B_shift_decay.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, P_grid, label="P(DeltaTheta)=1-exp(-DeltaTheta)")
plt.plot(Z_GRID, relaxed_grid, label="exp(-DeltaTheta)")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("Structural phase factors")
plt.title("TVGD structural phase transition")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07B_phase_factors.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.scatter(
    scan_df["theta0"],
    scan_df["bao_shift_max_z0p106_to_2_percent"],
    c=(scan_df["status"] == "PASSA_TRANSICAO_TVGD_H0_LOCAL_SEM_BAO_CMB").astype(int),
    s=10
)
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO tolerance")
plt.xlabel("Theta0")
plt.ylabel("Max BAO shift [%]")
plt.title("Allowed structural transition window")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07B_allowed_window.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(scan_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_TRANSICAO_ESTRUTURAL_TVGD_HUBBLE":
    print("The TVGD transition P(Theta)=1-exp(-Theta) can generate a local H0 enhancement and suppress it before BAO/CMB scales.")
else:
    print("The tested TVGD structural transition did not find a safe Hubble window.")

print("\nDONE.")
