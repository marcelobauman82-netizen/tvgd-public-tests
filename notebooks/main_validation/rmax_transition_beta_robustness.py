# ============================================================
# TVGD — 12_RMAX_TRANSITION_BETA_ROBUSTNESS
# Robustez da correção beta na janela crítica SPARC
#
# Objetivo:
#   Testar se a preferência por uma correção suave
#
#       P_beta(u) = 1 - exp(-u^beta)
#
#   com beta próximo de 1.05--1.10 é estável contra:
#
#   1. Jackknife por galáxia
#   2. Bootstrap por galáxia
#   3. Comparação com beta=1
#
# Saída:
#   results/rmax_transition_beta_robustness/
#
# Uso:
#   python notebooks/main_validation/12_rmax_transition_beta_robustness.py
# ============================================================

import os
import zipfile
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "12_rmax_transition_beta_robustness.py"

OUTDIR = "results/rmax_transition_beta_robustness"
os.makedirs(OUTDIR, exist_ok=True)

SPARC_URL = "https://astroweb.case.edu/SPARC/Rotmod_LTG.zip"
ZIP_PATH = os.path.join(OUTDIR, "Rotmod_LTG.zip")

A0_REF = 1.2e-10
KPC_TO_M = 3.0856775814913673e19

U_MIN = 0.5
U_MAX = 2.0

EV_FRAC_FLOOR = 0.05

# Grade de beta usada nos testes RMAX.
BETA_GRID = np.round(np.arange(0.90, 1.151, 0.005), 3)

# Grade pequena de a0 para permitir ajuste global simples.
A0_GRID = A0_REF * np.round(np.arange(0.90, 1.151, 0.025), 3)

# Bootstrap.
N_BOOT = 500
RANDOM_SEED = 12345

rng = np.random.default_rng(RANDOM_SEED)


print("=" * 72)
print("TVGD — 12_RMAX_TRANSITION_BETA_ROBUSTNESS")
print("Robustez da correção beta na janela crítica SPARC")
print("=" * 72)


# ============================================================
# Download SPARC
# ============================================================

def download_if_needed(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"Arquivo SPARC já existe: {path}")
        return

    print(f"Tentando baixar: {url}")
    urllib.request.urlretrieve(url, path)

    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        raise RuntimeError("Falha ao baixar Rotmod_LTG.zip.")

    print(f"Download OK: {path}")
    print(f"Tamanho: {os.path.getsize(path)} bytes")


download_if_needed(SPARC_URL, ZIP_PATH)


# ============================================================
# Leitura SPARC
# ============================================================

def read_rotmod_zip(zip_path):
    rows = []

    with zipfile.ZipFile(zip_path, "r") as z:
        dat_files = [f for f in z.namelist() if f.lower().endswith(".dat")]
        print(f"Arquivos .dat encontrados: {len(dat_files)}")

        for fname in dat_files:
            galaxy = os.path.basename(fname).replace("_rotmod.dat", "").replace(".dat", "")

            raw = z.read(fname).decode("utf-8", errors="ignore")
            for line in raw.splitlines():
                s = line.strip()

                if not s or s.startswith("#"):
                    continue

                parts = s.split()
                if len(parts) < 6:
                    continue

                try:
                    rad, vobs, evobs, vgas, vdisk, vbul = [float(x) for x in parts[:6]]
                except Exception:
                    continue

                rows.append({
                    "galaxy": galaxy,
                    "R_kpc": rad,
                    "Vobs": vobs,
                    "eVobs": evobs,
                    "Vgas": vgas,
                    "Vdisk": vdisk,
                    "Vbul": vbul,
                    "source_file": fname,
                })

    df = pd.DataFrame(rows)

    if len(df) == 0:
        raise RuntimeError("Nenhum ponto SPARC foi lido.")

    return df


df = read_rotmod_zip(ZIP_PATH)

print(f"Galáxias lidas: {df['galaxy'].nunique()}")
print(f"Pontos totais: {len(df)}")


# ============================================================
# Preparação física
# ============================================================

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["R_kpc", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul"])

df = df[(df["R_kpc"] > 0) & (df["Vobs"] > 0)].copy()

df["eVobs"] = pd.to_numeric(df["eVobs"], errors="coerce")
df["eVobs"] = df["eVobs"].fillna(EV_FRAC_FLOOR * df["Vobs"])
df["eVobs"] = np.maximum(df["eVobs"], EV_FRAC_FLOOR * df["Vobs"])

df["Vbar2"] = df["Vgas"]**2 + df["Vdisk"]**2 + df["Vbul"]**2
df["Vbar"] = np.sqrt(np.maximum(df["Vbar2"], 0.0))

df = df[df["Vbar"] > 0].copy()

df["gobs"] = (df["Vobs"] * 1000.0)**2 / (df["R_kpc"] * KPC_TO_M)
df["gbar"] = (df["Vbar"] * 1000.0)**2 / (df["R_kpc"] * KPC_TO_M)

df["egobs"] = df["gobs"] * 2.0 * df["eVobs"] / np.maximum(df["Vobs"], 1e-12)
df["egobs"] = np.maximum(df["egobs"], 0.10 * df["gobs"])

df["u_ref"] = np.sqrt(df["gbar"] / A0_REF)

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["gobs", "gbar", "egobs", "u_ref"])
df = df[(df["gobs"] > 0) & (df["gbar"] > 0) & (df["egobs"] > 0)].copy()

print(f"Pontos após sanidade: {len(df)}")
print(f"u_ref min/med/max: {df['u_ref'].min():.6g}, {df['u_ref'].median():.6g}, {df['u_ref'].max():.6g}")


# ============================================================
# Janela de transição
# ============================================================

trans = df[(df["u_ref"] >= U_MIN) & (df["u_ref"] <= U_MAX)].copy()

if len(trans) < 100 or trans["galaxy"].nunique() < 10:
    raise RuntimeError("A janela de transição ficou pequena demais.")

print("\nJanela de transição:")
print(f"U_MIN = {U_MIN}")
print(f"U_MAX = {U_MAX}")
print(f"Pontos: {len(trans)}")
print(f"Galáxias: {trans['galaxy'].nunique()}")

trans.to_csv(os.path.join(OUTDIR, "RMAX12_transition_window.csv"), index=False)


# ============================================================
# Modelos
# ============================================================

def P_beta(u, beta):
    u = np.asarray(u)
    return 1.0 - np.exp(-np.power(u, beta))


def predict_g(gbar, a0, beta):
    u = np.sqrt(np.maximum(gbar, 0.0) / a0)
    P = P_beta(u, beta)
    return gbar / np.maximum(P, 1e-12)


def fit_amp(gobs, gpred, err):
    w = 1.0 / np.maximum(err, 1e-99)**2
    den = np.sum(w * gpred**2)

    if den <= 0 or not np.isfinite(den):
        return 0.0, np.inf

    A = np.sum(w * gobs * gpred) / den
    A = max(float(A), 0.0)

    chi2 = np.sum(((gobs - A * gpred) / err)**2)
    return A, float(chi2)


def evaluate_dataset(data, beta=1.0, a0=A0_REF):
    """
    Ajusta amplitude por galáxia e retorna chi2/AIC/BIC.
    """
    if len(data) == 0:
        return None

    chi2_total = 0.0
    n_used = 0
    n_gal_used = 0
    amps = []

    for gal, sub in data.groupby("galaxy"):
        if len(sub) < 3:
            continue

        gpred = predict_g(sub["gbar"].values, a0=a0, beta=beta)
        A, chi2 = fit_amp(sub["gobs"].values, gpred, sub["egobs"].values)

        if not np.isfinite(chi2):
            continue

        chi2_total += chi2
        n_used += len(sub)
        n_gal_used += 1
        amps.append(A)

    if n_used == 0 or n_gal_used == 0:
        return None

    # Parâmetros: amplitude por galáxia + beta + a0 quando usados como busca global.
    k_amp = n_gal_used
    dof = max(n_used - k_amp, 1)

    return {
        "beta": float(beta),
        "a0": float(a0),
        "a0_over_ref": float(a0 / A0_REF),
        "N": int(n_used),
        "N_gal": int(n_gal_used),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "median_A": float(np.median(amps)),
        "p16_A": float(np.percentile(amps, 16)),
        "p84_A": float(np.percentile(amps, 84)),
    }


def global_grid_search(data, beta_grid=BETA_GRID, a0_grid=A0_GRID):
    """
    Busca beta e a0 em grade.
    AIC/BIC contam amplitude por galáxia + beta + a0 como parâmetros efetivos.
    """
    rows = []

    for a0 in a0_grid:
        for beta in beta_grid:
            row = evaluate_dataset(data, beta=beta, a0=a0)
            if row is None:
                continue

            k_params = row["N_gal"] + 2
            row["k_params"] = int(k_params)
            row["AIC"] = row["chi2"] + 2 * k_params
            row["BIC"] = row["chi2"] + np.log(max(row["N"], 2)) * k_params
            rows.append(row)

    res = pd.DataFrame(rows)

    if len(res) == 0:
        raise RuntimeError("Grid search retornou vazio.")

    res = res.sort_values("AIC").reset_index(drop=True)
    return res


# ============================================================
# Ajuste global
# ============================================================

print("\nRodando ajuste global na janela de transição...")

global_scan = global_grid_search(trans)
global_scan.to_csv(os.path.join(OUTDIR, "RMAX12_global_beta_a0_scan.csv"), index=False)

best_global = global_scan.iloc[0].to_dict()

# Comparação beta=1 no mesmo melhor a0 global
same_a0_beta1 = evaluate_dataset(trans, beta=1.0, a0=best_global["a0"])
k_beta1 = same_a0_beta1["N_gal"] + 1
same_a0_beta1["AIC"] = same_a0_beta1["chi2"] + 2 * k_beta1
same_a0_beta1["BIC"] = same_a0_beta1["chi2"] + np.log(max(same_a0_beta1["N"], 2)) * k_beta1

delta_AIC_beta1_global = float(same_a0_beta1["AIC"] - best_global["AIC"])

print("\nMelhor global:")
for k in ["a0", "a0_over_ref", "beta", "N", "N_gal", "chi2", "chi2red", "AIC", "BIC", "median_A"]:
    print(f"{k}: {best_global[k]}")

print(f"\nDelta AIC beta=1 menos melhor global: {delta_AIC_beta1_global:.6f}")


# ============================================================
# Jackknife por galáxia
# ============================================================

galaxies = sorted(trans["galaxy"].unique())

print(f"\nRodando jackknife por galáxia: {len(galaxies)} galáxias")

jack_rows = []

for i, gal_out in enumerate(galaxies, start=1):
    sub = trans[trans["galaxy"] != gal_out].copy()

    if len(sub) < 100 or sub["galaxy"].nunique() < 10:
        continue

    scan = global_grid_search(sub)
    best = scan.iloc[0].to_dict()

    beta1 = evaluate_dataset(sub, beta=1.0, a0=best["a0"])
    k_beta1 = beta1["N_gal"] + 1
    beta1["AIC"] = beta1["chi2"] + 2 * k_beta1

    best["left_out_galaxy"] = gal_out
    best["delta_AIC_beta1_minus_best"] = float(beta1["AIC"] - best["AIC"])
    jack_rows.append(best)

    if i % 10 == 0:
        print(f"Jackknife: {i} / {len(galaxies)}")

jack = pd.DataFrame(jack_rows)
jack.to_csv(os.path.join(OUTDIR, "RMAX12_jackknife_by_galaxy.csv"), index=False)


# ============================================================
# Bootstrap por galáxia
# ============================================================

print("\nRodando bootstrap por galáxias...")
print(f"N_BOOT = {N_BOOT}")

boot_rows = []

gal_array = np.array(galaxies)

for b in range(1, N_BOOT + 1):
    sample_gals = rng.choice(gal_array, size=len(gal_array), replace=True)

    parts = []
    for j, gal in enumerate(sample_gals):
        sub = trans[trans["galaxy"] == gal].copy()
        sub["galaxy"] = sub["galaxy"].astype(str) + f"__boot{j}"
        parts.append(sub)

    boot_df = pd.concat(parts, ignore_index=True)

    try:
        scan = global_grid_search(boot_df)
        best = scan.iloc[0].to_dict()

        beta1 = evaluate_dataset(boot_df, beta=1.0, a0=best["a0"])
        k_beta1 = beta1["N_gal"] + 1
        beta1["AIC"] = beta1["chi2"] + 2 * k_beta1

        best["boot_id"] = b
        best["delta_AIC_beta1_minus_best"] = float(beta1["AIC"] - best["AIC"])
        boot_rows.append(best)

    except Exception as exc:
        boot_rows.append({
            "boot_id": b,
            "error": str(exc),
        })

    if b % 50 == 0:
        print(f"Bootstrap: {b} / {N_BOOT}")

boot = pd.DataFrame(boot_rows)
boot.to_csv(os.path.join(OUTDIR, "RMAX12_bootstrap_by_galaxy.csv"), index=False)

boot_valid = boot.dropna(subset=["beta", "a0_over_ref", "delta_AIC_beta1_minus_best"]).copy()


# ============================================================
# Diagnóstico
# ============================================================

def safe_percentile(series, q):
    if len(series) == 0:
        return np.nan
    return float(np.percentile(series, q))


jack_valid = jack.dropna(subset=["beta", "a0_over_ref", "delta_AIC_beta1_minus_best"]).copy()

summary = {
    "FINAL_STATUS": None,
    "script": SCRIPT_NAME,
    "N_transition_points": int(len(trans)),
    "N_transition_galaxies": int(trans["galaxy"].nunique()),
    "global_best_a0": float(best_global["a0"]),
    "global_best_a0_over_ref": float(best_global["a0_over_ref"]),
    "global_best_beta": float(best_global["beta"]),
    "global_best_N": int(best_global["N"]),
    "global_best_N_gal": int(best_global["N_gal"]),
    "global_best_chi2": float(best_global["chi2"]),
    "global_best_chi2red": float(best_global["chi2red"]),
    "global_best_AIC": float(best_global["AIC"]),
    "delta_AIC_beta1_global": float(delta_AIC_beta1_global),

    "jack_N": int(len(jack_valid)),
    "jack_beta_median": safe_percentile(jack_valid["beta"], 50),
    "jack_beta_p16": safe_percentile(jack_valid["beta"], 16),
    "jack_beta_p84": safe_percentile(jack_valid["beta"], 84),
    "jack_beta_min": float(jack_valid["beta"].min()) if len(jack_valid) else np.nan,
    "jack_beta_max": float(jack_valid["beta"].max()) if len(jack_valid) else np.nan,
    "jack_frac_beta_gt_1": float(np.mean(jack_valid["beta"] > 1.0)) if len(jack_valid) else np.nan,
    "jack_frac_beta_gt_1p03": float(np.mean(jack_valid["beta"] > 1.03)) if len(jack_valid) else np.nan,
    "jack_frac_AIC_beta_beats_beta1": float(np.mean(jack_valid["delta_AIC_beta1_minus_best"] > 0)) if len(jack_valid) else np.nan,
    "jack_delta_AIC_median": safe_percentile(jack_valid["delta_AIC_beta1_minus_best"], 50),

    "boot_N": int(len(boot_valid)),
    "boot_beta_median": safe_percentile(boot_valid["beta"], 50),
    "boot_beta_p16": safe_percentile(boot_valid["beta"], 16),
    "boot_beta_p84": safe_percentile(boot_valid["beta"], 84),
    "boot_beta_p2p5": safe_percentile(boot_valid["beta"], 2.5),
    "boot_beta_p97p5": safe_percentile(boot_valid["beta"], 97.5),
    "boot_frac_beta_gt_1": float(np.mean(boot_valid["beta"] > 1.0)) if len(boot_valid) else np.nan,
    "boot_frac_beta_gt_1p03": float(np.mean(boot_valid["beta"] > 1.03)) if len(boot_valid) else np.nan,
    "boot_frac_AIC_beta_beats_beta1": float(np.mean(boot_valid["delta_AIC_beta1_minus_best"] > 0)) if len(boot_valid) else np.nan,
    "boot_delta_AIC_median": safe_percentile(boot_valid["delta_AIC_beta1_minus_best"], 50),
}

criteria = {
    "global_beta_gt_1": summary["global_best_beta"] > 1.0,
    "global_beta_near_expected": 1.03 <= summary["global_best_beta"] <= 1.11,
    "global_beta1_disfavored_AIC": summary["delta_AIC_beta1_global"] > 10.0,

    "jack_sufficient": summary["jack_N"] >= 50,
    "jack_beta_gt_1_majority": summary["jack_frac_beta_gt_1"] > 0.80,
    "jack_beta_gt_1p03_majority": summary["jack_frac_beta_gt_1p03"] > 0.60,
    "jack_AIC_majority": summary["jack_frac_AIC_beta_beats_beta1"] > 0.60,

    "boot_sufficient": summary["boot_N"] >= 300,
    "boot_beta_gt_1_majority": summary["boot_frac_beta_gt_1"] > 0.80,
    "boot_beta_gt_1p03_majority": summary["boot_frac_beta_gt_1p03"] > 0.60,
    "boot_AIC_majority": summary["boot_frac_AIC_beta_beats_beta1"] > 0.60,
}

score = sum(bool(v) for v in criteria.values())

if score >= 9:
    status = "PASSOU_FORTE_RMAX12_BETA_ROBUSTO"
elif score >= 7:
    status = "PASSOU_MODERADO_FORTE_RMAX12_BETA_ROBUSTO_COM_ALERTAS"
elif score >= 5:
    status = "PASSOU_MODERADO_RMAX12_BETA_COM_ALERTAS"
else:
    status = "SINAL_FRACO_RMAX12_BETA_NAO_FECHADO"

summary["FINAL_STATUS"] = status
summary["score"] = int(score)
summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "RMAX12_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — RMAX 12")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(global_scan["beta"], global_scan["AIC"], ".", alpha=0.25)
plt.axvline(1.0, linestyle="--", label="beta=1")
plt.axvline(summary["global_best_beta"], linestyle=":", label=f"best beta={summary['global_best_beta']:.3f}")
plt.xlabel(r"$\beta$")
plt.ylabel("AIC")
plt.title("RMAX 12 — Global beta scan")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX12_global_beta_scan_AIC.png"), dpi=170)
plt.close()

if len(jack_valid):
    plt.figure(figsize=(8, 5))
    plt.hist(jack_valid["beta"], bins=25)
    plt.axvline(1.0, linestyle="--", label="beta=1")
    plt.axvline(summary["jack_beta_median"], linestyle=":", label=f"median={summary['jack_beta_median']:.3f}")
    plt.xlabel(r"$\beta_{\rm best}$")
    plt.ylabel("N")
    plt.title("RMAX 12 — Jackknife beta distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "RMAX12_jackknife_beta_hist.png"), dpi=170)
    plt.close()

if len(boot_valid):
    plt.figure(figsize=(8, 5))
    plt.hist(boot_valid["beta"], bins=30)
    plt.axvline(1.0, linestyle="--", label="beta=1")
    plt.axvline(summary["boot_beta_median"], linestyle=":", label=f"median={summary['boot_beta_median']:.3f}")
    plt.xlabel(r"$\beta_{\rm best}$")
    plt.ylabel("N")
    plt.title("RMAX 12 — Bootstrap beta distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "RMAX12_bootstrap_beta_hist.png"), dpi=170)
    plt.close()

if len(boot_valid):
    plt.figure(figsize=(8, 5))
    plt.hist(boot_valid["delta_AIC_beta1_minus_best"], bins=30)
    plt.axvline(0.0, linestyle="--", label=r"$\Delta$AIC=0")
    plt.xlabel(r"$\Delta {\rm AIC}(\beta=1 - \beta_{\rm best})$")
    plt.ylabel("N")
    plt.title("RMAX 12 — Bootstrap AIC gain")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, "RMAX12_bootstrap_delta_AIC.png"), dpi=170)
    plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste RMAX 12.")
