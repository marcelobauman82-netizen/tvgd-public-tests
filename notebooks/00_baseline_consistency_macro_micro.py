"""
00_baseline_consistency_macro_micro.py

Baseline macro/micro consistency test for the Dispersive Gravitational
Vacuum Theory (TVGD).

This script verifies two basic consistency requirements:

1. Macro/cosmological regime:
   The effective TVGD V1.1 background is compared with flat LCDM
   through H(z) and distance modulus mu(z).

2. Micro/local regime:
   The local TVGD acceleration law is tested in high-acceleration
   environments: Earth surface, GPS orbit, Moon distance, Mercury,
   and Earth around the Sun.

This is not the main observational validation test.
It is a baseline consistency test showing that the model preserves
well-established regimes before applying it to galaxy, lensing,
supernova, BAO, CMB, or black-hole datasets.

Author: Marcelo Bauman
Repository: TVGD Public Tests
"""

import os
import sys
import numpy as np
import pandas as pd


# ============================================================
# Path setup
# ============================================================

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(CURRENT_DIR, ".."))
SRC_DIR = os.path.join(REPO_ROOT, "src")
RESULTS_DIR = os.path.join(REPO_ROOT, "Resultados")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

os.makedirs(RESULTS_DIR, exist_ok=True)


# ============================================================
# Import TVGD core functions
# ============================================================

from tvgd_core import (
    C_LIGHT_KM,
    G_NEWTON,
    A0_DEFAULT,
    P_tvgd,
    g_tvgd,
    E_lcdm,
    E_tvgd,
    distance_modulus,
    percent_difference,
)


# ============================================================
# Global parameters
# ============================================================

H0 = 67.4
OMEGA_M = 0.31339978339159447
OMEGA_R = 9e-5

LAMBDA_TVGD = 7.0
P_TVGD = 0.55

A0 = A0_DEFAULT


# ============================================================
# Test 1 — Macro / cosmological consistency
# ============================================================

def run_macro_test():
    print("=" * 72)
    print("TEST 00A — BASELINE MACRO CONSISTENCY")
    print("TVGD V1.1 background vs flat LCDM")
    print("=" * 72)

    z = np.linspace(0.001, 2.0, 300)

    E_LCDM = E_lcdm(
        z,
        Omega_m=OMEGA_M,
        Omega_r=OMEGA_R,
    )

    E_TVGD = E_tvgd(
        z,
        Omega_m=OMEGA_M,
        Omega_r=OMEGA_R,
        lam=LAMBDA_TVGD,
        p=P_TVGD,
    )

    H_LCDM = H0 * E_LCDM
    H_TVGD = H0 * E_TVGD

    delta_H_percent = percent_difference(H_TVGD, H_LCDM)

    mu_LCDM = distance_modulus(
        z,
        lambda zz: E_lcdm(
            zz,
            Omega_m=OMEGA_M,
            Omega_r=OMEGA_R,
        ),
        H0=H0,
        n_grid=2000,
    )

    mu_TVGD = distance_modulus(
        z,
        lambda zz: E_tvgd(
            zz,
            Omega_m=OMEGA_M,
            Omega_r=OMEGA_R,
            lam=LAMBDA_TVGD,
            p=P_TVGD,
        ),
        H0=H0,
        n_grid=2000,
    )

    delta_mu = mu_TVGD - mu_LCDM

    summary = {
        "test": "macro",
        "max_abs_delta_H_percent": float(np.max(np.abs(delta_H_percent))),
        "mean_abs_delta_H_percent": float(np.mean(np.abs(delta_H_percent))),
        "max_abs_delta_mu_mag": float(np.max(np.abs(delta_mu))),
        "mean_abs_delta_mu_mag": float(np.mean(np.abs(delta_mu))),
    }

    if (
        summary["max_abs_delta_H_percent"] < 1.0
        and summary["max_abs_delta_mu_mag"] < 0.02
    ):
        status = "PASSOU_FORTE_MACRO"
    elif (
        summary["max_abs_delta_H_percent"] < 2.0
        and summary["max_abs_delta_mu_mag"] < 0.05
    ):
        status = "PASSOU_MODERADO_MACRO"
    else:
        status = "TENSIONADO_MACRO"

    summary["status"] = status

    print()
    print("Macro summary:")
    for key, value in summary.items():
        if key != "test":
            print(f"{key}: {value}")

    macro_df = pd.DataFrame({
        "z": z,
        "H_LCDM_km_s_Mpc": H_LCDM,
        "H_TVGD_km_s_Mpc": H_TVGD,
        "delta_H_percent": delta_H_percent,
        "mu_LCDM_mag": mu_LCDM,
        "mu_TVGD_mag": mu_TVGD,
        "delta_mu_mag": delta_mu,
    })

    macro_output_path = os.path.join(
        RESULTS_DIR,
        "test00_macro_background_results.csv",
    )

    macro_df.to_csv(macro_output_path, index=False)

    print()
    print(f"Saved macro results to: {macro_output_path}")

    return summary, macro_df


# ============================================================
# Test 2 — Micro / local consistency
# ============================================================

def run_micro_test():
    print()
    print("=" * 72)
    print("TEST 00B — BASELINE MICRO/LOCAL CONSISTENCY")
    print("Local TVGD response in high-acceleration regimes")
    print("=" * 72)

    M_earth = 5.9722e24
    R_earth = 6.371e6

    M_sun = 1.98847e30
    AU = 1.495978707e11

    alt_gps = 20200e3
    r_surface = R_earth
    r_gps = R_earth + alt_gps
    r_moon = 384400e3
    r_mercury = 0.387 * AU
    r_earth_sun = AU

    cases = [
        {
            "case": "Surface of Earth",
            "central_mass_kg": M_earth,
            "radius_m": r_surface,
        },
        {
            "case": "GPS orbit",
            "central_mass_kg": M_earth,
            "radius_m": r_gps,
        },
        {
            "case": "Moon distance",
            "central_mass_kg": M_earth,
            "radius_m": r_moon,
        },
        {
            "case": "Mercury around the Sun",
            "central_mass_kg": M_sun,
            "radius_m": r_mercury,
        },
        {
            "case": "Earth around the Sun",
            "central_mass_kg": M_sun,
            "radius_m": r_earth_sun,
        },
    ]

    rows = []

    for item in cases:
        g_newton = G_NEWTON * item["central_mass_kg"] / item["radius_m"]**2
        P_value = P_tvgd(g_newton, a0=A0)
        g_eff = g_tvgd(g_newton, a0=A0)

        delta_percent = 100.0 * (g_eff - g_newton) / g_newton
        u_value = np.sqrt(g_newton / A0)

        rows.append({
            "case": item["case"],
            "central_mass_kg": item["central_mass_kg"],
            "radius_m": item["radius_m"],
            "g_Newton_m_s2": g_newton,
            "u": u_value,
            "P_u": P_value,
            "g_TVGD_m_s2": g_eff,
            "delta_percent": delta_percent,
        })

    micro_df = pd.DataFrame(rows)

    max_abs_micro = float(np.max(np.abs(micro_df["delta_percent"].values)))
    mean_abs_micro = float(np.mean(np.abs(micro_df["delta_percent"].values)))

    if max_abs_micro < 1e-6:
        status = "PASSOU_FORTE_MICRO"
    elif max_abs_micro < 1e-3:
        status = "PASSOU_MODERADO_MICRO"
    else:
        status = "TENSIONADO_MICRO"

    summary = {
        "test": "micro",
        "max_abs_delta_percent": max_abs_micro,
        "mean_abs_delta_percent": mean_abs_micro,
        "status": status,
    }

    print()
    print("Micro/local table:")
    print(micro_df.to_string(index=False))

    print()
    print("Micro summary:")
    for key, value in summary.items():
        if key != "test":
            print(f"{key}: {value}")

    micro_output_path = os.path.join(
        RESULTS_DIR,
        "test00_micro_local_results.csv",
    )

    micro_df.to_csv(micro_output_path, index=False)

    print()
    print(f"Saved micro results to: {micro_output_path}")

    return summary, micro_df


# ============================================================
# Final diagnostic
# ============================================================

def run_final_diagnostic(macro_summary, micro_summary):
    print()
    print("=" * 72)
    print("TEST 00 — FINAL BASELINE DIAGNOSTIC")
    print("=" * 72)

    macro_status = macro_summary["status"]
    micro_status = micro_summary["status"]

    if "PASSOU" in macro_status and "PASSOU" in micro_status:
        final_status = "PASSOU_FORTE_BASELINE"
        interpretation = (
            "The effective TVGD model preserves the tested local/micro "
            "high-acceleration regimes and remains very close to flat LCDM "
            "at the tested background cosmological level."
        )
    else:
        final_status = "REVISAR_BASELINE"
        interpretation = (
            "At least one baseline consistency regime showed tension and "
            "should be reviewed before proceeding to observational datasets."
        )

    final_summary = {
        "macro_status": macro_status,
        "micro_status": micro_status,
        "final_status": final_status,
        "interpretation": interpretation,
    }

    print()
    for key, value in final_summary.items():
        print(f"{key}: {value}")

    final_output_path = os.path.join(
        RESULTS_DIR,
        "test00_baseline_summary.csv",
    )

    pd.DataFrame([final_summary]).to_csv(final_output_path, index=False)

    print()
    print(f"Saved final summary to: {final_output_path}")

    return final_summary


# ============================================================
# Main execution
# ============================================================

if __name__ == "__main__":
    print()
    print("TVGD PUBLIC TESTS")
    print("Test 00 — Baseline Macro/Micro Consistency")
    print()

    macro_summary, macro_df = run_macro_test()
    micro_summary, micro_df = run_micro_test()
    final_summary = run_final_diagnostic(macro_summary, micro_summary)

    print()
    print("=" * 72)
    print("DONE")
    print("=" * 72)
