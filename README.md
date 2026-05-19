# TVGD Public Tests

Public reproducibility tests for the **Dispersive Gravitational Vacuum Theory (TVGD)**.

This repository contains numerical scripts, source code, and output tables designed to allow independent verification, modification, and falsification of the TVGD effective model.

The purpose of this repository is not to hide assumptions, but to expose them clearly.

---

## Central TVGD Response Function

In the local/galactic regime, the basic structural response function tested here is

\[
P(u) = 1 - e^{-u}
\]

with

\[
u = \sqrt{\frac{g_{\rm bar}}{a_0}}
\]

and

\[
g_{\rm TVGD} = \frac{g_{\rm bar}}{P(u)}.
\]

The fiducial acceleration scale used in the public tests is

\[
a_0 = 1.2 \times 10^{-10}\ {\rm m\,s^{-2}}.
\]

In the high-acceleration regime,

\[
P(u) \rightarrow 1,
\]

so the model recovers the Newtonian/General Relativistic local limit.

In the low-acceleration regime,

\[
P(u)<1,
\]

and the effective gravitational acceleration is amplified.

---

## Current Public Tests

### Test 00 — Baseline Macro/Micro Consistency

This test verifies that the effective TVGD model preserves two established regimes before applying it to the main observational datasets:

```text
1. the local high-acceleration Solar-System regime;
2. the background cosmological expansion regime close to flat ΛCDM. ```

Current result:

```text
Macro: PASSOU_FORTE_MACRO
Micro: PASSOU_FORTE_MICRO
Final: PASSOU_FORTE_BASELINE
```

Interpretation:

The tested TVGD background remains very close to flat ΛCDM at the background level, while the local high-acceleration limit satisfies \(P(u)\approx 1\), preserving the micro/local regime.

---

### Test 01 — SPARC/RAR Galactic Regime

This test block evaluates the TVGD effective acceleration law using the public SPARC rotation-curve dataset.

The tested relation is

\[
g_{\rm TVGD}=
\frac{g_{\rm bar}}
{1-\exp\left[-\sqrt{g_{\rm bar}/a_0}\right]}.
\]

The SPARC block contains three public scripts:

```text
01_sparc_rar_tvgd.py
01B_sparc_ml_robustness.py
01C_sparc_rar_comparison.py
```

Current results:

```text
Test 01A — SPARC/RAR:
Status: PASSOU_FORTE_SPARC_RAR_PRELIMINAR
N_points = 3389
N_galaxies = 175
TVGD improvement MAE = 70.61 %
TVGD improvement RMSE = 62.53 %

Test 01B — Mass-to-light robustness:
Status: PASSOU_MODERADO_ROBUSTEZ_ML
N_grid = 66
frac_pass_global_strong = 1.0
median improvement MAE = 68.13 %
minimum improvement MAE = 48.83 %

Test 01C — TVGD vs empirical RAR/MOND:
Status: PASSOU_FORTE_COMO_REALIZACAO_TVGD_DA_RAR
max_relative_TVGD_minus_RAR_fixed = 0.0
TVGD_vs_best_RAR_MAE_percent = 0.766 %
TVGD_vs_best_RAR_RMSE_percent = 1.840 %
```

Interpretation:

The TVGD response function reproduces the empirical RAR behavior in the galactic regime, improves strongly over the baryonic-only model, and remains robust under a reasonable stellar mass-to-light variation.

This result should be interpreted as a strong public reproducibility check of the galactic effective regime, not as a complete Bayesian model-selection analysis.

---

## Repository Structure

Current structure:

```text
tvgd-public-tests/
│
├── README.md
├── requirements.txt
├── LICENSE
│
├── src/
│   └── tvgd_core.py
│
├── cadernos/
│   ├── 00_baseline_consistency_macro_micro.py
│   ├── 01_sparc_rar_tvgd.py
│   ├── 01B_sparc_ml_robustness.py
│   └── 01C_sparc_rar_comparison.py
│
├── Resultados/
│
└── Papel/
```

---

## Reproducibility Principle

All numerical claims should be reproducible from the scripts in this repository.

If a reported trend cannot be reproduced from public data and public code, the corresponding claim should be considered unsupported.

This repository is intended as a falsifiable public record of the TVGD numerical program.

---

## Current Scope and Next Tests

The current public package includes:

```text
1. baseline macro/micro consistency;
2. SPARC/RAR galactic regime;
3. mass-to-light robustness;
4. TVGD vs empirical RAR/MOND control.
```

Future public versions will include:

```text
- black-hole exterior and regular-core checks
- compressed BAO/CMB/Pantheon+ tests
- lensing-sector tests
- environmental Hubble tests
- covariant-limit consistency tests
```

Each new test will be added as a separate reproducibility module.

---

## Scientific Status

TVGD is currently presented as an effective theoretical and phenomenological framework under active development.

The public tests should be interpreted as reproducibility checks of the proposed effective model.

The goal is to make every relevant claim testable, modifiable, and falsifiable by independent readers.

---

## Code and Data Availability

The code and scripts are available directly in this repository.

When external public datasets are used, the corresponding source, citation, and download instructions are provided inside the relevant script.

No hidden or private data should be required to reproduce the public tests.

---

## Citation

If this repository is used or cited, please cite the corresponding manuscript or archived release when available.

A Zenodo DOI will be added after the first stable release.

---

## Author

Marcelo Bauman

---

## License

This project is released under the MIT License.
