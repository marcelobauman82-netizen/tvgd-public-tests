# ============================================================
# TVGD PUBLIC TESTS
# TEST 07H — HUBBLE STRUCTURAL PHASE ROBUSTNESS
# ============================================================
#
# Goal:
# Test whether the TVGD Hubble structural transition works only
# at one tuned point or across a robust parameter window.
#
# Structural transition:
#
#     Xi_TVGD(z) = 1 + A_H exp[-DeltaTheta(z)]
#
# with:
#
#     DeltaTheta(z) = theta0 * [(1+z)^s - 1]
#
# Interpretation:
#
# - z ~ 0:
#     relaxed local vacuum -> higher inferred H0.
#
# - BAO/CMB:
#     saturated/early concentrated vacuum -> CMB-like H0 preserved.
#
# This test scans theta0 and s.
#
# It checks:
#
# 1. local H0 enhancement;
# 2. decay before BAO;
# 3. negligible CMB leakage;
# 4. existence of a broad allowed parameter region;
# 5. no extreme fine-tuning around a single theta0.
#
# This is not a Pantheon statistical test.
# It is a structural robustness test of the Hubble transition.
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test07H_hubble_structural_phase_robustness"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Reference values
# ============================================================

H0_CMB = 67.4
H0_LOCAL = 73.04

A_H = H0_LOCAL / H0_CMB - 1.0

Z_CMB = 1089.92
Z_BAO_MIN = 0.106

# Scan ranges
THETA0_GRID = np.linspace(1.0, 80.0, 800)
S_GRID = np.linspace(0.8, 2.2, 281)

# Redshift grid for diagnostics
Z_GRID = np.linspace(0.0, 0.30, 1200)
Z_BAO_GRID = np.linspace(Z_BAO_MIN, 2.0, 1000)

# Criteria
LOCAL_SHIFT_TARGET_PERCENT = 100.0 * A_H

MIN_SHIFT_Z0_PERCENT = 8.0       # must reproduce local H0 enhancement
MIN_SHIFT_Z001_PERCENT = 4.0     # still visible very locally
MIN_SHIFT_Z003_PERCENT = 1.0     # weak but present up to z~0.03

MAX_SHIFT_Z005_PERCENT = 2.0     # should already be decaying
MAX_SHIFT_Z010_PERCENT = 0.5     # should be mostly gone by z~0.10
MAX_BAO_SHIFT_PERCENT = 0.5      # BAO safety
MAX_CMB_SHIFT_ABS = 1e-12        # CMB safety

# Robustness criteria
MIN_ALLOWED_POINTS = 100
MIN_ALLOWED_THETA_WIDTH = 5.0
MIN_ALLOWED_S_WIDTH = 0.10

print("=" * 72)
print("TEST 07H — HUBBLE STRUCTURAL PHASE ROBUSTNESS")
print("=" * 72)

print("\nReference:")
print("H0_CMB =", H0_CMB)
print("H0_LOCAL =", H0_LOCAL)
print("A_H =", A_H)
print("Required local shift [%] =", LOCAL_SHIFT_TARGET_PERCENT)

# ============================================================
# TVGD Hubble transition
# ============================================================

def delta_theta_z(z, theta0, s):
    z = np.asarray(z, dtype=float)
    return theta0 * ((1.0 + np.maximum(z, 0.0))**s - 1.0)


def Xi_TVGD(z, theta0, s):
    dth = delta_theta_z(z, theta0=theta0, s=s)
    return 1.0 + A_H * np.exp(-np.maximum(dth, 0.0))


def H0_eff(z, theta0, s):
    return H0_CMB * Xi_TVGD(z, theta0=theta0, s=s)


def shift_percent(z, theta0, s):
    return 100.0 * (Xi_TVGD(z, theta0=theta0, s=s) - 1.0)


def P_transition(z, theta0, s):
    """
    Activated structural fraction:
        P = 1 - exp[-DeltaTheta].
    """
    dth = delta_theta_z(z, theta0=theta0, s=s)
    return 1.0 - np.exp(-np.maximum(dth, 0.0))


# ============================================================
# Scan
# ============================================================

rows = []

for s in S_GRID:
    for theta0 in THETA0_GRID:
        shift_z0 = float(shift_percent(0.0, theta0, s))
        shift_z001 = float(shift_percent(0.01, theta0, s))
        shift_z003 = float(shift_percent(0.03, theta0, s))
        shift_z005 = float(shift_percent(0.05, theta0, s))
        shift_z010 = float(shift_percent(0.10, theta0, s))

        bao_shift_max = float(np.max(np.abs(shift_percent(Z_BAO_GRID, theta0, s))))
        cmb_shift_abs = float(abs(Xi_TVGD(Z_CMB, theta0, s) - 1.0))

        P_z001 = float(P_transition(0.01, theta0, s))
        P_z003 = float(P_transition(0.03, theta0, s))
        P_bao = float(P_transition(Z_BAO_MIN, theta0, s))
        P_cmb = float(P_transition(Z_CMB, theta0, s))

        local_ok = (
            shift_z0 >= MIN_SHIFT_Z0_PERCENT
            and shift_z001 >= MIN_SHIFT_Z001_PERCENT
            and shift_z003 >= MIN_SHIFT_Z003_PERCENT
        )

        decay_ok = (
            shift_z005 <= MAX_SHIFT_Z005_PERCENT
            and shift_z010 <= MAX_SHIFT_Z010_PERCENT
        )

        bao_ok = bao_shift_max <= MAX_BAO_SHIFT_PERCENT
        cmb_ok = cmb_shift_abs <= MAX_CMB_SHIFT_ABS

        passed = local_ok and decay_ok and bao_ok and cmb_ok

        if passed:
            status = "PASSA_JANELA_HUBBLE_TVGD"
        elif local_ok and bao_ok and cmb_ok:
            status = "PASSA_LOCAL_BAO_CMB_MAS_DECAIMENTO_TENSIONADO"
        elif bao_ok and cmb_ok:
            status = "PASSA_BAO_CMB_MAS_LOCAL_FRACO"
        elif cmb_ok:
            status = "PRESERVA_CMB_MAS_TENSIONA_BAO_LOCAL"
        else:
            status = "FALHA_CMB"

        rows.append({
            "theta0": theta0,
            "s": s,

            "shift_z0_percent": shift_z0,
            "shift_z0p01_percent": shift_z001,
            "shift_z0p03_percent": shift_z003,
            "shift_z0p05_percent": shift_z005,
            "shift_z0p10_percent": shift_z010,

            "bao_shift_max_percent": bao_shift_max,
            "cmb_shift_abs": cmb_shift_abs,

            "P_z0p01": P_z001,
            "P_z0p03": P_z003,
            "P_BAO_z0p106": P_bao,
            "P_CMB": P_cmb,

            "local_ok": bool(local_ok),
            "decay_ok": bool(decay_ok),
            "bao_ok": bool(bao_ok),
            "cmb_ok": bool(cmb_ok),
            "passed": bool(passed),
            "status": status,
        })

scan_df = pd.DataFrame(rows)

pass_df = scan_df[scan_df["passed"]].copy()

N_total = len(scan_df)
N_pass = len(pass_df)
frac_pass = float(N_pass / N_total)

print("\nScan summary:")
print("N_total =", N_total)
print("N_pass =", N_pass)
print("frac_pass =", frac_pass)

# ============================================================
# Robustness diagnostics
# ============================================================

if N_pass > 0:
    theta_min = float(pass_df["theta0"].min())
    theta_max = float(pass_df["theta0"].max())
    theta_width = theta_max - theta_min

    s_min = float(pass_df["s"].min())
    s_max = float(pass_df["s"].max())
    s_width = s_max - s_min

    # Score: prefer local signature while minimizing BAO leakage.
    pass_df["score"] = (
        pass_df["shift_z0p01_percent"]
        + 0.5 * pass_df["shift_z0p03_percent"]
        - 10.0 * pass_df["bao_shift_max_percent"]
        - 0.5 * np.abs(pass_df["s"] - 1.5)
    )

    best = pass_df.sort_values("score", ascending=False).iloc[0]

    robustness_ok = (
        N_pass >= MIN_ALLOWED_POINTS
        and theta_width >= MIN_ALLOWED_THETA_WIDTH
        and s_width >= MIN_ALLOWED_S_WIDTH
    )

else:
    theta_min = np.nan
    theta_max = np.nan
    theta_width = 0.0

    s_min = np.nan
    s_max = np.nan
    s_width = 0.0

    best = scan_df.sort_values("bao_shift_max_percent").iloc[0]
    robustness_ok = False

# Natural point nearest s=1.5 among passed cases
if N_pass > 0:
    pass_df["distance_to_natural_s"] = np.abs(pass_df["s"] - 1.5)
    natural_best = pass_df.sort_values(
        ["distance_to_natural_s", "bao_shift_max_percent"],
        ascending=[True, True]
    ).iloc[0]
else:
    natural_best = best

# ============================================================
# Final status
# ============================================================

if N_pass > 0 and robustness_ok:
    final_status = "PASSOU_FORTE_ROBUSTEZ_HUBBLE_TVGD"
elif N_pass > 0:
    final_status = "PASSOU_MODERADO_JANELA_HUBBLE_TVGD_MAS_ESTREITA"
else:
    final_status = "TENSIONADO_ROBUSTEZ_HUBBLE_TVGD"

summary = {
    "H0_CMB": H0_CMB,
    "H0_LOCAL": H0_LOCAL,
    "A_H": A_H,
    "required_local_shift_percent": LOCAL_SHIFT_TARGET_PERCENT,

    "N_total": N_total,
    "N_pass": N_pass,
    "frac_pass": frac_pass,

    "theta0_min_allowed": theta_min,
    "theta0_max_allowed": theta_max,
    "theta0_width_allowed": theta_width,

    "s_min_allowed": s_min,
    "s_max_allowed": s_max,
    "s_width_allowed": s_width,

    "robustness_ok": bool(robustness_ok),

    "best_theta0": float(best["theta0"]),
    "best_s": float(best["s"]),
    "best_shift_z0_percent": float(best["shift_z0_percent"]),
    "best_shift_z0p01_percent": float(best["shift_z0p01_percent"]),
    "best_shift_z0p03_percent": float(best["shift_z0p03_percent"]),
    "best_shift_z0p05_percent": float(best["shift_z0p05_percent"]),
    "best_shift_z0p10_percent": float(best["shift_z0p10_percent"]),
    "best_bao_shift_max_percent": float(best["bao_shift_max_percent"]),
    "best_cmb_shift_abs": float(best["cmb_shift_abs"]),
    "best_P_BAO_z0p106": float(best["P_BAO_z0p106"]),
    "best_P_CMB": float(best["P_CMB"]),

    "natural_best_theta0": float(natural_best["theta0"]),
    "natural_best_s": float(natural_best["s"]),
    "natural_best_shift_z0p01_percent": float(natural_best["shift_z0p01_percent"]),
    "natural_best_shift_z0p03_percent": float(natural_best["shift_z0p03_percent"]),
    "natural_best_bao_shift_max_percent": float(natural_best["bao_shift_max_percent"]),
    "natural_best_cmb_shift_abs": float(natural_best["cmb_shift_abs"]),

    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nFinal summary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save outputs
# ============================================================

scan_path = os.path.join(OUTDIR, "test07H_hubble_structural_phase_scan.csv")
pass_path = os.path.join(OUTDIR, "test07H_hubble_structural_phase_allowed_window.csv")
summary_path = os.path.join(OUTDIR, "test07H_hubble_structural_phase_summary.csv")

scan_df.to_csv(scan_path, index=False)
pass_df.to_csv(pass_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Plots
# ============================================================

# Plot 1: allowed region in theta0-s plane
plt.figure(figsize=(8, 5))

plt.scatter(
    scan_df["theta0"],
    scan_df["s"],
    c=scan_df["passed"].astype(int),
    s=4,
    alpha=0.6
)

plt.axhline(1.5, color="black", linestyle="--", label="natural s=1.5")
plt.xlabel("theta0")
plt.ylabel("s")
plt.title("Allowed TVGD Hubble structural window")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07H_allowed_theta0_s_window.png"), dpi=220)
plt.show()

# Plot 2: BAO leakage
plt.figure(figsize=(8, 5))

plt.scatter(
    scan_df["theta0"],
    scan_df["bao_shift_max_percent"],
    c=scan_df["passed"].astype(int),
    s=4,
    alpha=0.6
)

plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO tolerance")
plt.xlabel("theta0")
plt.ylabel("max BAO shift [%]")
plt.title("BAO leakage across structural scan")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07H_bao_leakage.png"), dpi=220)
plt.show()

# Plot 3: best transition curve
best_theta0 = float(best["theta0"])
best_s = float(best["s"])

shift_grid = shift_percent(Z_GRID, best_theta0, best_s)
P_grid = P_transition(Z_GRID, best_theta0, best_s)

plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, shift_grid, label="TVGD H0 shift")
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO tolerance")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("Xi(z)-1 [%]")
plt.title("Best TVGD Hubble structural decay")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07H_best_shift_decay.png"), dpi=220)
plt.show()

# Plot 4: structural phase factor
plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, P_grid, label="P=1-exp[-DeltaTheta]")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("P(z)")
plt.title("Best structural phase activation")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07H_best_phase_activation.png"), dpi=220)
plt.show()

# Plot 5: natural s=1.5 curve
natural_theta0 = float(natural_best["theta0"])
natural_s = float(natural_best["s"])

natural_shift_grid = shift_percent(Z_GRID, natural_theta0, natural_s)

plt.figure(figsize=(8, 5))
plt.plot(Z_GRID, natural_shift_grid, label=f"natural-like s={natural_s:.3f}")
plt.axhline(MAX_BAO_SHIFT_PERCENT, color="red", linestyle="--", label="BAO tolerance")
plt.axvline(Z_BAO_MIN, color="red", linestyle="--", label="BAO lower z")
plt.xlabel("z")
plt.ylabel("Xi(z)-1 [%]")
plt.title("Natural structural exponent Hubble decay")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test07H_natural_s_shift_decay.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(scan_path)
print(pass_path)
print(summary_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_ROBUSTEZ_HUBBLE_TVGD":
    print(
        "The TVGD Hubble structural transition is not restricted to a single tuned point. "
        "A finite allowed region in theta0 and s generates local H0 enhancement while "
        "preserving BAO/CMB."
    )
elif "PASSOU_MODERADO" in final_status:
    print(
        "The TVGD Hubble structural transition works, but the allowed region is narrow "
        "and should be treated with caution."
    )
else:
    print(
        "The scanned TVGD Hubble structural transition did not find a robust allowed window."
    )

print("\nDONE.")
