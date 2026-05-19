# TVGD Public Tests

Public reproducibility tests for the **Dispersive Gravitational Vacuum Theory (TVGD)**.

This repository contains numerical scripts, source code, and output tables designed to allow independent verification, modification, and falsification of the TVGD effective model.

The purpose of this repository is not to hide assumptions, but to expose them clearly.

---

## 1. Central TVGD Response Function

In the local and galactic regime, the effective TVGD response function tested in this repository is:

```text
P(u) = 1 - exp(-u)
```

with

```text
u = sqrt(g_bar / a0)
```

and

```text
g_TVGD = g_bar / P(u)
```

Equivalently:

```text
g_TVGD = g_bar / [1 - exp(-sqrt(g_bar/a0))]
```

The fiducial acceleration scale used in the public tests is:

```text
a0 = 1.2e-10 m s^-2
```

### High-acceleration limit

For

```text
g_bar >> a0
```

we have:

```text
u >> 1
P(u) -> 1
g_TVGD -> g_bar
```

Therefore, the model recovers the Newtonian / General Relativistic local limit.

### Low-acceleration regime

For

```text
g_bar <= a0
```

we have:

```text
P(u) < 1
g_TVGD > g_bar
```

Therefore, the effective gravitational acceleration is amplified.

This is the regime where the model is intended to reproduce the observed galactic acceleration relation.

---

## 2. Current Public Tests

### Test 00 — Baseline Macro/Micro Consistency

This test verifies that the effective TVGD model preserves two established regimes before applying it to the main observational datasets:

```text
1. the local high-acceleration Solar-System regime;
2. the background cosmological expansion regime close to flat Lambda-CDM.
```

Current result:

```text
Macro: PASSOU_FORTE_MACRO
Micro: PASSOU_FORTE_MICRO
Final: PASSOU_FORTE_BASELINE
```

Interpretation:

The tested TVGD background remains close to flat Lambda-CDM at the background level, while the local high-acceleration limit satisfies:

```text
P(u) approximately 1
```

preserving the micro/local regime.

Script:

```text
cadernos/00_baseline_consistency_macro_micro.py
```

---

### Test 01 — SPARC/RAR Galactic Regime

This test block evaluates the TVGD effective acceleration law using the public SPARC rotation-curve dataset.

The tested relation is:

```text
g_TVGD = g_bar / [1 - exp(-sqrt(g_bar/a0))]
```

The SPARC block contains three public scripts:

```text
cadernos/01_sparc_rar_tvgd.py
cadernos/01B_sparc_ml_robustness.py
cadernos/01C_sparc_rar_comparison.py
```

---

#### Test 01A — SPARC/RAR Main Test

This test compares:

```text
1. baryonic-only acceleration;
2. TVGD effective acceleration.
```

Current result:

```text
Status: PASSOU_FORTE_SPARC_RAR_PRELIMINAR
N_points = 3389
N_galaxies = 175
TVGD improvement MAE = 70.61 %
TVGD improvement RMSE = 62.53 %
```

Interpretation:

TVGD strongly reduces the logarithmic residuals relative to the baryonic-only model in the SPARC/RAR dataset.

---

#### Test 01B — Mass-to-Light Robustness

This test checks whether the SPARC/RAR result depends on a single stellar mass-to-light choice.

The grid tested was:

```text
Upsilon_disk  = 0.30 to 0.80
Upsilon_bulge = 0.50 to 1.00
```

Current result:

```text
Status: PASSOU_MODERADO_ROBUSTEZ_ML
N_grid = 66
frac_pass_global_strong = 1.0
median improvement MAE = 68.13 %
minimum improvement MAE = 48.83 %
```

Interpretation:

The SPARC/RAR result is not a fragile artifact of one specific stellar mass-to-light choice. Across the tested grid, TVGD improves the global logarithmic residuals over the baryonic-only model.

---

#### Test 01C — TVGD vs Empirical RAR/MOND Controls

This test compares TVGD against:

```text
1. baryonic-only acceleration;
2. empirical RAR with fixed g_dagger;
3. empirical RAR with optimized g_dagger;
4. MOND simple interpolating function;
5. MOND standard interpolating function.
```

Current result:

```text
Status: PASSOU_FORTE_COMO_REALIZACAO_TVGD_DA_RAR
max_relative_TVGD_minus_RAR_fixed = 0.0
TVGD_vs_best_RAR_MAE_percent = 0.766 %
TVGD_vs_best_RAR_RMSE_percent = 1.840 %
```

Interpretation:

When

```text
a0 = g_dagger = 1.2e-10 m s^-2
```

the tested TVGD response is mathematically equivalent to the empirical RAR form.

The result should not be interpreted as TVGD "defeating" RAR/MOND in this test. Rather, it shows that TVGD reproduces the empirical RAR behavior as an effective realization of the same acceleration relation, while being embedded in a broader theoretical program.

---

## 3. Summary of Public Results

| Test | Regime | Status | Main meaning |
|---|---|---|---|
| Test 00 | Macro/Micro baseline | PASSOU_FORTE_BASELINE | TVGD preserves the tested local and background limits |
| Test 01A | SPARC/RAR | PASSOU_FORTE_SPARC_RAR_PRELIMINAR | TVGD improves strongly over baryonic-only acceleration |
| Test 01B | Mass-to-light robustness | PASSOU_MODERADO_ROBUSTEZ_ML | Result is robust across a stellar M/L grid |
| Test 01C | RAR/MOND control | PASSOU_FORTE_COMO_REALIZACAO_TVGD_DA_RAR | TVGD reproduces the empirical RAR form |

---

## 4. Repository Structure

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

## 5. How to Run

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/tvgd-public-tests.git
cd tvgd-public-tests
```

Install dependencies:

```bash
pip install -r requirements.txt
```

Run the baseline test:

```bash
python cadernos/00_baseline_consistency_macro_micro.py
```

Run the SPARC/RAR main test:

```bash
python cadernos/01_sparc_rar_tvgd.py
```

Run the mass-to-light robustness test:

```bash
python cadernos/01B_sparc_ml_robustness.py
```

Run the TVGD vs RAR/MOND comparison:

```bash
python cadernos/01C_sparc_rar_comparison.py
```

Each script creates its own output folder inside:

```text
Resultados/
```

---

## 6. Reproducibility Principle

All numerical claims should be reproducible from the scripts in this repository.

If a reported trend cannot be reproduced from public data and public code, the corresponding claim should be considered unsupported.

This repository is intended as a falsifiable public record of the TVGD numerical program.

---

## 7. Scientific Status

TVGD is currently presented as an effective theoretical and phenomenological framework under active development.

The public tests should be interpreted as reproducibility checks of the proposed effective model.

The current results are not presented as a final proof of the theory. They are presented as public, testable evidence that the effective TVGD response function:

```text
1. preserves the tested high-acceleration local regime;
2. remains close to the tested Lambda-CDM background limit;
3. reproduces the empirical galactic acceleration relation;
4. improves strongly over baryonic-only acceleration in SPARC;
5. remains robust under a reasonable stellar mass-to-light variation.
```

The goal is to make every relevant claim testable, modifiable, and falsifiable by independent readers.

---

## 8. Current Scope and Next Tests

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

## 9. Code and Data Availability

The code and scripts are available directly in this repository.

When external public datasets are used, the corresponding source, citation, and download instructions are provided inside the relevant script.

No hidden or private data should be required to reproduce the public tests.

---

## 10. Citation

If this repository is used or cited, please cite the corresponding manuscript or archived release when available.

A Zenodo DOI will be added after the first stable release.

---

## 11. Author

Marcelo Bauman

---

## 12. License

This project is released under the MIT License.

