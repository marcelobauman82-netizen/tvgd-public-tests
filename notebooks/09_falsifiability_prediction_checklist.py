# ============================================================
# TVGD PUBLIC TESTS
# TEST 09 — FALSIFIABILITY AND PREDICTION CHECKLIST
# ============================================================
#
# Goal:
# Produce a public falsifiability and prediction checklist for TVGD.
#
# This notebook does not fit new data.
# It defines where the TVGD effective model can be tested,
# confirmed, constrained, or falsified.
#
# It creates:
#
# 1. A falsifiability table.
# 2. A prediction table.
# 3. A public summary CSV.
# 4. A Markdown report for GitHub.
#
# ============================================================

import os
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Output
# ============================================================

OUTDIR = "Resultados/test09_falsifiability_prediction_checklist"
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 72)
print("TEST 09 — FALSIFIABILITY AND PREDICTION CHECKLIST")
print("=" * 72)

# ============================================================
# Central TVGD assumptions
# ============================================================

central_assumptions = [
    {
        "assumption_id": "A1",
        "sector": "Galaxies",
        "statement": "The effective local response is P(u)=1-exp(-u), with u=sqrt(g_bar/a0).",
        "mathematical_form": "g_TVGD = g_bar / [1-exp(-sqrt(g_bar/a0))]",
        "status": "core_effective_assumption",
    },
    {
        "assumption_id": "A2",
        "sector": "High acceleration / Solar System",
        "statement": "For g_bar >> a0, P(u)->1 and the Newtonian/GR local limit is recovered.",
        "mathematical_form": "lim_{u>>1} P(u)=1",
        "status": "screening_requirement",
    },
    {
        "assumption_id": "A3",
        "sector": "Low acceleration",
        "statement": "For g_bar << a0, the model approaches a MOND-like deep-regime behavior.",
        "mathematical_form": "g_TVGD ~ sqrt(g_bar a0)",
        "status": "galactic_prediction",
    },
    {
        "assumption_id": "A4",
        "sector": "Hubble tension",
        "statement": "The vacuum structural phase transition can separate local relaxed and primordial saturated regimes.",
        "mathematical_form": "Xi(z)=1+A_H exp[-DeltaTheta(z)]",
        "status": "cosmological_prediction",
    },
    {
        "assumption_id": "A5",
        "sector": "Black holes",
        "statement": "Strong-field structural saturation regularizes the black-hole core.",
        "mathematical_form": "Q_psi(x)=1-exp[-(x/xc)^n]",
        "status": "strong_field_prediction",
    },
]

assumptions_df = pd.DataFrame(central_assumptions)

# ============================================================
# Falsifiability checklist
# ============================================================

falsifiability_items = [
    {
        "test_id": "F1",
        "sector": "Galaxies / SPARC / RAR",
        "prediction": "Observed galaxy accelerations should follow the TVGD response curve with small systematic residuals.",
        "falsification_condition": "High-quality galaxy data show a statistically significant, repeatable deviation from g_bar/[1-exp(-sqrt(g_bar/a0))] not explainable by baryonic uncertainties.",
        "observable": "g_obs versus g_bar",
        "required_data": "SPARC-like rotation curves, baryonic decomposition, stellar M/L priors",
        "severity_if_failed": "critical",
        "current_public_tests": "01, 01B, 01C",
        "current_status": "tested_preliminarily",
    },
    {
        "test_id": "F2",
        "sector": "Solar System / high acceleration",
        "prediction": "In high acceleration environments, TVGD corrections must be negligible.",
        "falsification_condition": "TVGD predicts measurable deviations in Earth, Moon, GPS, Mercury, or Earth-Sun regimes above observational bounds.",
        "observable": "local gravitational acceleration, orbital dynamics, timing",
        "required_data": "Solar-System bounds, GPS timing, Mercury perihelion constraints",
        "severity_if_failed": "critical",
        "current_public_tests": "00, 05",
        "current_status": "passed_public_sanity",
    },
    {
        "test_id": "F3",
        "sector": "Cosmology / BAO",
        "prediction": "The structural Hubble transition must decay before BAO scales.",
        "falsification_condition": "The required local H0 enhancement necessarily produces BAO leakage above tolerance.",
        "observable": "H(z), D_M/r_d, D_H/r_d, D_V/r_d",
        "required_data": "BAO compilations, compressed likelihoods, full cosmological fits",
        "severity_if_failed": "major",
        "current_public_tests": "03, 06, 07B, 07H",
        "current_status": "passed_structural_tests",
    },
    {
        "test_id": "F4",
        "sector": "CMB",
        "prediction": "The TVGD Hubble transition must be negligible at recombination.",
        "falsification_condition": "The same transition needed for local H0 produces non-negligible CMB leakage or shifts acoustic scales beyond allowed bounds.",
        "observable": "CMB acoustic scale, shift parameter, TT/TE/EE spectra",
        "required_data": "Planck/ACT/SPT likelihoods, CLASS/CAMB implementation",
        "severity_if_failed": "critical",
        "current_public_tests": "03, 07B, 07H",
        "current_status": "passed_compressed_structural_tests",
    },
    {
        "test_id": "F5",
        "sector": "Pantheon+ / Hubble",
        "prediction": "The Hubble structural transition should improve or remain competitive against CMB-like H0 and survive non-calibrator shape tests.",
        "falsification_condition": "Full-covariance Pantheon+ and non-calibrator tests show no improvement beyond offset/placebo controls.",
        "observable": "distance modulus residuals versus redshift",
        "required_data": "Pantheon+SH0ES data and full covariance",
        "severity_if_failed": "major",
        "current_public_tests": "07E, 07G",
        "current_status": "moderate_not_final",
    },
    {
        "test_id": "F6",
        "sector": "Black holes / static",
        "prediction": "Structural saturation should remove the curvature singularity while preserving the external Schwarzschild regime.",
        "falsification_condition": "The saturated metric fails to regularize curvature or destroys horizon/photon-sphere/shadow benchmarks.",
        "observable": "Kretschmann scalar, horizon, photon sphere, shadow impact parameter",
        "required_data": "Analytic/numerical metric audit",
        "severity_if_failed": "major",
        "current_public_tests": "02, 08",
        "current_status": "passed_strong_effective_test",
    },
    {
        "test_id": "F7",
        "sector": "Black holes / rotation",
        "prediction": "A rotating Kerr-like extension should preserve the external Kerr limit and regularize the ring proxy for moderate spin.",
        "falsification_condition": "Moderate rotation destroys the horizon/ergosphere or fails to regularize the ring proxy.",
        "observable": "Delta roots, ergosurface, external Kerr limit, ring curvature proxy",
        "required_data": "Rotating metric proxy audit; future exact Kerr-TVGD solution",
        "severity_if_failed": "major",
        "current_public_tests": "08B",
        "current_status": "passed_moderate_rotation_with_near_extremal_constraint",
    },
    {
        "test_id": "F8",
        "sector": "Lensing",
        "prediction": "TVGD effective amplification should improve lensing consistency in low-acceleration/low-density regimes without breaking strong-lensing constraints.",
        "falsification_condition": "Controlled lensing systems require corrections incompatible with the TVGD response or show no relation to the predicted amplification.",
        "observable": "Einstein radius, lensing mass, weak-lensing excess surface density",
        "required_data": "Strong and weak lensing catalogs",
        "severity_if_failed": "major",
        "current_public_tests": "04B",
        "current_status": "controlled_preliminary",
    },
    {
        "test_id": "F9",
        "sector": "Stability / causality",
        "prediction": "The effective model must avoid obvious ghost-like, acausal, or phantom behavior in tested regimes.",
        "falsification_condition": "The model requires c_s^2<0, c_s^2>1, phantom instability, or uncontrolled local fifth-force behavior.",
        "observable": "sound speed, screening, effective equation of state",
        "required_data": "Analytic stability tests, cosmological reconstruction",
        "severity_if_failed": "critical",
        "current_public_tests": "05, 06",
        "current_status": "passed_public_effective_tests",
    },
]

falsifiability_df = pd.DataFrame(falsifiability_items)

# ============================================================
# Forward predictions
# ============================================================

prediction_items = [
    {
        "prediction_id": "P1",
        "sector": "Galaxies",
        "prediction": "Residuals should be minimized near the TVGD response curve P(u)=1-exp(-u).",
        "future_test": "Bayesian model comparison against NFW halos, MOND/RAR, and flexible empirical models.",
        "expected_TVGD_signature": "Low number of effective parameters with stable residual behavior.",
        "priority": "high",
    },
    {
        "prediction_id": "P2",
        "sector": "Hubble",
        "prediction": "The local H0 enhancement should decay rapidly before BAO redshifts.",
        "future_test": "Use independent SN compilations and BAO likelihoods with consistent priors.",
        "expected_TVGD_signature": "Local enhancement without BAO/CMB leakage.",
        "priority": "high",
    },
    {
        "prediction_id": "P3",
        "sector": "Hubble / environment",
        "prediction": "Low-density relaxed environments may show slightly larger effective expansion than dense environments.",
        "future_test": "Matched-redshift environmental Hubble residuals with void/cluster catalogs.",
        "expected_TVGD_signature": "Small positive H_void - H_dense trend.",
        "priority": "medium",
    },
    {
        "prediction_id": "P4",
        "sector": "Black holes",
        "prediction": "Near-extremal rotation may constrain the structural core scale xc.",
        "future_test": "Exact rotating TVGD solution and comparison with high-spin black-hole observations.",
        "expected_TVGD_signature": "Moderate spin allowed; near-extremal spin restricts large saturated cores.",
        "priority": "high",
    },
    {
        "prediction_id": "P5",
        "sector": "Black-hole shadows",
        "prediction": "External shadow should remain close to Kerr/Schwarzschild for compact cores.",
        "future_test": "Compare with EHT-like constraints on M87* and Sgr A*.",
        "expected_TVGD_signature": "Small deviations unless xc is large or spin is near-extremal.",
        "priority": "medium",
    },
    {
        "prediction_id": "P6",
        "sector": "Structural redshift",
        "prediction": "TVGD may produce sub-km/s structural redshift signatures in galaxies.",
        "future_test": "High-precision spectroscopic comparison across low-acceleration galactic regions.",
        "expected_TVGD_signature": "Correlation with u=sqrt(g_bar/a0), strongest near transitional regimes.",
        "priority": "medium",
    },
    {
        "prediction_id": "P7",
        "sector": "CMB / large-scale structure",
        "prediction": "A full CLASS/CAMB implementation should preserve acoustic scales while possibly modifying late-time growth.",
        "future_test": "Full Boltzmann-code likelihood analysis.",
        "expected_TVGD_signature": "CMB preserved; late-time structure shifts constrained but nonzero.",
        "priority": "high",
    },
]

predictions_df = pd.DataFrame(prediction_items)

# ============================================================
# Scoring
# ============================================================

severity_weight = {
    "critical": 3,
    "major": 2,
    "medium": 1,
    "minor": 0.5,
}

status_score = {
    "passed_public_sanity": 2,
    "passed_structural_tests": 2,
    "passed_compressed_structural_tests": 2,
    "passed_strong_effective_test": 2,
    "passed_public_effective_tests": 2,
    "passed_moderate_rotation_with_near_extremal_constraint": 1.5,
    "moderate_not_final": 1,
    "tested_preliminarily": 1,
    "controlled_preliminary": 1,
}

falsifiability_df["severity_weight"] = falsifiability_df["severity_if_failed"].map(severity_weight)
falsifiability_df["status_score"] = falsifiability_df["current_status"].map(status_score).fillna(0)

falsifiability_df["weighted_public_support"] = (
    falsifiability_df["severity_weight"] * falsifiability_df["status_score"]
)

total_possible_score = float((falsifiability_df["severity_weight"] * 2).sum())
total_public_score = float(falsifiability_df["weighted_public_support"].sum())
public_support_fraction = total_public_score / total_possible_score

N_falsification_tests = len(falsifiability_df)
N_critical = int((falsifiability_df["severity_if_failed"] == "critical").sum())
N_major = int((falsifiability_df["severity_if_failed"] == "major").sum())
N_predictions = len(predictions_df)

# Classification
if public_support_fraction >= 0.75:
    final_status = "PASSOU_FORTE_FALSIFICABILIDADE_PUBLICA_TVGD"
elif public_support_fraction >= 0.50:
    final_status = "PASSOU_MODERADO_FALSIFICABILIDADE_PUBLICA_TVGD"
else:
    final_status = "FALSIFICABILIDADE_PUBLICA_INCOMPLETA_TVGD"

summary = {
    "N_central_assumptions": len(assumptions_df),
    "N_falsification_tests": N_falsification_tests,
    "N_critical_tests": N_critical,
    "N_major_tests": N_major,
    "N_forward_predictions": N_predictions,
    "total_possible_score": total_possible_score,
    "total_public_score": total_public_score,
    "public_support_fraction": public_support_fraction,
    "final_status": final_status,
    "main_interpretation": (
        "TVGD has explicit falsification routes across galaxies, Solar System, "
        "cosmology, Hubble tension, black holes, lensing, and stability. "
        "Several sectors already have public effective tests, while full "
        "covariant and full-likelihood tests remain future work."
    ),
}

summary_df = pd.DataFrame([summary])

print("\nCentral assumptions:")
display(assumptions_df)

print("\nFalsifiability checklist:")
display(falsifiability_df)

print("\nForward predictions:")
display(predictions_df)

print("\nSummary:")
display(summary_df)

print("\nFinal status:", final_status)

# ============================================================
# Save CSV outputs
# ============================================================

assumptions_path = os.path.join(OUTDIR, "test09_central_assumptions.csv")
falsifiability_path = os.path.join(OUTDIR, "test09_falsifiability_checklist.csv")
predictions_path = os.path.join(OUTDIR, "test09_forward_predictions.csv")
summary_path = os.path.join(OUTDIR, "test09_falsifiability_summary.csv")

assumptions_df.to_csv(assumptions_path, index=False)
falsifiability_df.to_csv(falsifiability_path, index=False)
predictions_df.to_csv(predictions_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Create Markdown report
# ============================================================

md_path = os.path.join(OUTDIR, "test09_falsifiability_report.md")

with open(md_path, "w", encoding="utf-8") as f:
    f.write("# TVGD Test 09 — Falsifiability and Prediction Checklist\n\n")

    f.write("This report summarizes explicit falsification routes and forward predictions for the public TVGD effective model.\n\n")

    f.write("## Central assumptions\n\n")
    for _, row in assumptions_df.iterrows():
        f.write(f"### {row['assumption_id']} — {row['sector']}\n\n")
        f.write(f"**Statement:** {row['statement']}\n\n")
        f.write(f"**Mathematical form:** `{row['mathematical_form']}`\n\n")
        f.write(f"**Status:** {row['status']}\n\n")

    f.write("## Falsifiability checklist\n\n")
    for _, row in falsifiability_df.iterrows():
        f.write(f"### {row['test_id']} — {row['sector']}\n\n")
        f.write(f"**Prediction:** {row['prediction']}\n\n")
        f.write(f"**Falsification condition:** {row['falsification_condition']}\n\n")
        f.write(f"**Observable:** {row['observable']}\n\n")
        f.write(f"**Required data:** {row['required_data']}\n\n")
        f.write(f"**Severity if failed:** {row['severity_if_failed']}\n\n")
        f.write(f"**Current public tests:** {row['current_public_tests']}\n\n")
        f.write(f"**Current status:** {row['current_status']}\n\n")

    f.write("## Forward predictions\n\n")
    for _, row in predictions_df.iterrows():
        f.write(f"### {row['prediction_id']} — {row['sector']}\n\n")
        f.write(f"**Prediction:** {row['prediction']}\n\n")
        f.write(f"**Future test:** {row['future_test']}\n\n")
        f.write(f"**Expected TVGD signature:** {row['expected_TVGD_signature']}\n\n")
        f.write(f"**Priority:** {row['priority']}\n\n")

    f.write("## Summary\n\n")
    f.write(f"- Number of falsification tests: **{N_falsification_tests}**\n")
    f.write(f"- Number of critical tests: **{N_critical}**\n")
    f.write(f"- Number of major tests: **{N_major}**\n")
    f.write(f"- Number of forward predictions: **{N_predictions}**\n")
    f.write(f"- Public support fraction: **{public_support_fraction:.3f}**\n")
    f.write(f"- Final status: **{final_status}**\n\n")
    f.write(summary["main_interpretation"] + "\n")

# ============================================================
# Plots
# ============================================================

plt.figure(figsize=(8, 5))
sector_counts = falsifiability_df["sector"].value_counts()
sector_counts.plot(kind="bar")
plt.ylabel("Number of falsification routes")
plt.title("TVGD falsifiability routes by sector")
plt.xticks(rotation=45, ha="right")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test09_falsifiability_by_sector.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
severity_counts = falsifiability_df["severity_if_failed"].value_counts()
severity_counts.plot(kind="bar")
plt.ylabel("Number of tests")
plt.title("Severity of falsification routes")
plt.xticks(rotation=0)
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test09_falsifiability_severity.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.bar(
    falsifiability_df["test_id"],
    falsifiability_df["weighted_public_support"]
)
plt.ylabel("Weighted public support")
plt.title("Current public support by falsification route")
plt.grid(True, axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test09_weighted_public_support.png"), dpi=220)
plt.show()

print("\nSaved files:")
print(OUTDIR)
print(assumptions_path)
print(falsifiability_path)
print(predictions_path)
print(summary_path)
print(md_path)

print("\nShort interpretation:")
if final_status == "PASSOU_FORTE_FALSIFICABILIDADE_PUBLICA_TVGD":
    print(
        "TVGD has a strong public falsifiability structure: multiple sectors contain "
        "explicit predictions, failure conditions, and current public tests."
    )
elif final_status == "PASSOU_MODERADO_FALSIFICABILIDADE_PUBLICA_TVGD":
    print(
        "TVGD has a meaningful falsifiability structure, but several sectors still require "
        "stronger public likelihood or covariant tests."
    )
else:
    print(
        "TVGD still needs a clearer falsifiability structure before public presentation."
    )

print("\nDONE.")
