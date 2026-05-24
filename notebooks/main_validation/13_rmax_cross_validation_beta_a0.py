# ============================================================
# TVGD — 13_RMAX_CROSS_VALIDATION_BETA_A0
# Validação cruzada train/test da correção beta e degenerescência beta-a0
#
# Objetivo:
#   Testar se a correção suave
#
#       P_beta(u) = 1 - exp(-u^beta)
#
#   generaliza para galáxias fora do treino, comparando:
#
#   1. beta treinado vs beta=1 com o mesmo a0
#   2. beta fixo 1.055 vs beta=1 com a0 de referência
#   3. beta fixo 1.055 vs beta=1 com a0 alternativo
#
# Interpretação:
#   Este teste não busca o melhor ajuste absoluto do SPARC.
#   Ele mede se a tendência beta>1 é robusta e se não é apenas
#   degenerescência trivial com a0.
#
# Saída:
#   results/rmax_cross_validation_beta_a0/
#
# Uso:
#   python notebooks/main_validation/13_rmax_cross_validation_beta_a0.py
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

SCRIPT_NAME = "13_rmax_cross_validation_beta_a0.py"

OUTDIR = "results/rmax_cross_validation_beta_a0"
os.makedirs(OUTDIR, exist_ok=True)

SPARC_URL = "https://astroweb.case.edu/SPARC/Rotmod_LTG.zip"
ZIP_PATH = os.path.join(OUTDIR, "Rotmod_LTG.zip")

A0_REF = 1.2e-10
KPC_TO_M = 3.0856775814913673e19

U_MIN = 0.5
U_MAX = 2.0

EV_FRAC_FLOOR = 0.05

BETA_GRID = np.round(np.arange(0.90, 1.151, 0.005), 3)
A0_GRID = A0_REF * np.round(np.arange(0.90, 1.151, 0.025), 3)

# Correção beta mínima sugerida pelos testes RMAX.
BETA_FIXED = 1.055

N_SPLITS = 200
TRAIN_FRAC = 0.70
RANDOM_SEED = 20260524

rng = np.random.default_rng(RANDOM_SEED)


print("=" * 72)
print("TVGD — 13_RMAX_CROSS_VALIDATION_BETA_A0")
print("Validação cruzada train/test da correção beta")
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


# ============================================================
# Janela de transição
# ============================================================

trans = df[(df["u_ref"] >= U_MIN) & (df["u_ref"] <= U_MAX)].copy()

if len(trans) < 100 or trans["galaxy"].nunique() < 10:
    raise RuntimeError("A janela de transição ficou pequena demais.")

print("\nJanela de transição:")
print(f"Pontos: {len(trans)}")
print(f"Galáxias: {trans['galaxy'].nunique()}")
print(f"u_ref min/med/max: {trans['u_ref'].min():.6g}, {trans['u_ref'].median():.6g}, {trans['u_ref'].max():.6g}")

trans.to_csv(os.path.join(OUTDIR, "RMAX13_transition_window.csv"), index=False)


# ============================================================
# Modelo TVGD
# ============================================================

def P_beta(u, beta):
    return 1.0 - np.exp(-np.power(np.asarray(u), beta))


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

    chi2 = np.sum(((gobs - A*gpred) / err)**2)
    return A, float(chi2)


def evaluate(data, beta, a0, count_shape_params=0):
    """
    Avalia um beta/a0 fixo com amplitude por galáxia.
    count_shape_params:
      0: beta/a0 tratados como fixos.
      1: conta um parâmetro de forma.
      2: conta beta e a0.
    """
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

    k_params = n_gal_used + count_shape_params
    dof = max(n_used - k_params, 1)

    return {
        "beta": float(beta),
        "a0": float(a0),
        "a0_over_ref": float(a0 / A0_REF),
        "N": int(n_used),
        "N_gal": int(n_gal_used),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "AIC": float(chi2_total + 2*k_params),
        "BIC": float(chi2_total + np.log(max(n_used, 2))*k_params),
        "k_params": int(k_params),
        "median_A": float(np.median(amps)),
    }


def train_grid_search(data):
    rows = []

    for a0 in A0_GRID:
        for beta in BETA_GRID:
            row = evaluate(data, beta=beta, a0=a0, count_shape_params=2)
            if row is not None:
                rows.append(row)

    res = pd.DataFrame(rows)

    if len(res) == 0:
        raise RuntimeError("Grid search vazio.")

    return res.sort_values("AIC").reset_index(drop=True)


# ============================================================
# Validação cruzada
# ============================================================

galaxies = np.array(sorted(trans["galaxy"].unique()))
n_train = int(np.floor(TRAIN_FRAC * len(galaxies)))

print("\nRodando splits train/test...")
print(f"N_SPLITS = {N_SPLITS}")
print(f"TRAIN_FRAC = {TRAIN_FRAC}")
print(f"Galáxias totais na janela = {len(galaxies)}")
print(f"Galáxias treino por split = {n_train}")

rows = []

for split in range(1, N_SPLITS + 1):
    shuffled = rng.permutation(galaxies)
    train_gals = set(shuffled[:n_train])
    test_gals = set(shuffled[n_train:])

    train = trans[trans["galaxy"].isin(train_gals)].copy()
    test = trans[trans["galaxy"].isin(test_gals)].copy()

    if train["galaxy"].nunique() < 10 or test["galaxy"].nunique() < 5:
        continue

    try:
        train_scan = train_grid_search(train)
        best_train = train_scan.iloc[0].to_dict()

        beta_train = float(best_train["beta"])
        a0_train = float(best_train["a0"])

        # Teste A: beta treinado vs beta=1 no mesmo a0 treinado
        test_beta_train = evaluate(test, beta=beta_train, a0=a0_train, count_shape_params=0)
        test_beta1_same_a0 = evaluate(test, beta=1.0, a0=a0_train, count_shape_params=0)

        # Teste B: beta fixo 1.055 vs beta=1 com a0 referência
        test_beta_fixed_ref = evaluate(test, beta=BETA_FIXED, a0=A0_REF, count_shape_params=0)
        test_beta1_ref = evaluate(test, beta=1.0, a0=A0_REF, count_shape_params=0)

        # Teste C: beta fixo 1.055 vs beta=1 com melhor a0 alternativo para beta=1
        beta1_alt_rows = []
        for a0_alt in A0_GRID:
            r_alt = evaluate(test, beta=1.0, a0=a0_alt, count_shape_params=0)
            if r_alt is not None:
                beta1_alt_rows.append(r_alt)

        beta1_alt = pd.DataFrame(beta1_alt_rows).sort_values("AIC").iloc[0].to_dict()

        # Teste D: beta fixo 1.055 com melhor a0 próprio
        beta_fixed_alt_rows = []
        for a0_alt in A0_GRID:
            r_alt = evaluate(test, beta=BETA_FIXED, a0=a0_alt, count_shape_params=0)
            if r_alt is not None:
                beta_fixed_alt_rows.append(r_alt)

        beta_fixed_alt = pd.DataFrame(beta_fixed_alt_rows).sort_values("AIC").iloc[0].to_dict()

        row = {
            "split": split,
            "N_train": int(len(train)),
            "N_test": int(len(test)),
            "Ngal_train": int(train["galaxy"].nunique()),
            "Ngal_test": int(test["galaxy"].nunique()),

            "train_beta": beta_train,
            "train_a0": a0_train,
            "train_a0_over_ref": float(a0_train / A0_REF),
            "train_chi2red": float(best_train["chi2red"]),
            "train_AIC": float(best_train["AIC"]),

            "test_chi2red_beta_train_same_a0": float(test_beta_train["chi2red"]),
            "test_chi2red_beta1_same_a0": float(test_beta1_same_a0["chi2red"]),
            "dAIC_beta_train_minus_beta1_same_a0": float(test_beta_train["AIC"] - test_beta1_same_a0["AIC"]),

            "test_chi2red_beta_fixed_ref": float(test_beta_fixed_ref["chi2red"]),
            "test_chi2red_beta1_ref": float(test_beta1_ref["chi2red"]),
            "dAIC_beta_fixed_minus_beta1_ref": float(test_beta_fixed_ref["AIC"] - test_beta1_ref["AIC"]),

            "beta1_alt_a0_over_ref": float(beta1_alt["a0_over_ref"]),
            "test_chi2red_beta1_alt": float(beta1_alt["chi2red"]),

            "beta_fixed_alt_a0_over_ref": float(beta_fixed_alt["a0_over_ref"]),
            "test_chi2red_beta_fixed_alt": float(beta_fixed_alt["chi2red"]),
            "dAIC_beta_fixed_alt_minus_beta1_alt": float(beta_fixed_alt["AIC"] - beta1_alt["AIC"]),
        }

        rows.append(row)

    except Exception as exc:
        rows.append({
            "split": split,
            "error": str(exc),
        })

    if split % 20 == 0:
        print(f"Split: {split} / {N_SPLITS}")

cv = pd.DataFrame(rows)
cv.to_csv(os.path.join(OUTDIR, "RMAX13_cross_validation_splits.csv"), index=False)

valid = cv.dropna(subset=[
    "train_beta",
    "dAIC_beta_train_minus_beta1_same_a0",
    "dAIC_beta_fixed_minus_beta1_ref",
    "dAIC_beta_fixed_alt_minus_beta1_alt",
]).copy()

if len(valid) == 0:
    raise RuntimeError("Nenhum split válido.")


# ============================================================
# Diagnóstico
# ============================================================

def pct(x, q):
    return float(np.percentile(x, q))


# Convenção:
# dAIC = AIC(beta) - AIC(beta=1)
# Portanto dAIC < 0 favorece beta corrigido.
summary = {
    "FINAL_STATUS": None,
    "script": SCRIPT_NAME,
    "N_splits_requested": int(N_SPLITS),
    "N_splits_valid": int(len(valid)),
    "train_frac": float(TRAIN_FRAC),
    "beta_fixed": float(BETA_FIXED),

    "train_beta_median": pct(valid["train_beta"], 50),
    "train_beta_p16": pct(valid["train_beta"], 16),
    "train_beta_p84": pct(valid["train_beta"], 84),
    "train_beta_p2p5": pct(valid["train_beta"], 2.5),
    "train_beta_p97p5": pct(valid["train_beta"], 97.5),
    "train_a0_over_ref_median": pct(valid["train_a0_over_ref"], 50),
    "train_a0_over_ref_p16": pct(valid["train_a0_over_ref"], 16),
    "train_a0_over_ref_p84": pct(valid["train_a0_over_ref"], 84),

    "frac_train_beta_gt_1": float(np.mean(valid["train_beta"] > 1.0)),
    "frac_train_beta_gt_1p03": float(np.mean(valid["train_beta"] > 1.03)),
    "frac_train_beta_near_1055": float(np.mean(np.abs(valid["train_beta"] - BETA_FIXED) <= 0.04)),

    "dAIC_beta_train_same_a0_median": pct(valid["dAIC_beta_train_minus_beta1_same_a0"], 50),
    "dAIC_beta_train_same_a0_p16": pct(valid["dAIC_beta_train_minus_beta1_same_a0"], 16),
    "dAIC_beta_train_same_a0_p84": pct(valid["dAIC_beta_train_minus_beta1_same_a0"], 84),
    "frac_test_beta_train_beats_beta1_same_a0": float(np.mean(valid["dAIC_beta_train_minus_beta1_same_a0"] < 0.0)),

    "dAIC_beta_fixed_ref_median": pct(valid["dAIC_beta_fixed_minus_beta1_ref"], 50),
    "dAIC_beta_fixed_ref_p16": pct(valid["dAIC_beta_fixed_minus_beta1_ref"], 16),
    "dAIC_beta_fixed_ref_p84": pct(valid["dAIC_beta_fixed_minus_beta1_ref"], 84),
    "frac_test_beta_fixed_beats_beta1_ref": float(np.mean(valid["dAIC_beta_fixed_minus_beta1_ref"] < 0.0)),

    "dAIC_beta_fixed_alt_median": pct(valid["dAIC_beta_fixed_alt_minus_beta1_alt"], 50),
    "dAIC_beta_fixed_alt_p16": pct(valid["dAIC_beta_fixed_alt_minus_beta1_alt"], 16),
    "dAIC_beta_fixed_alt_p84": pct(valid["dAIC_beta_fixed_alt_minus_beta1_alt"], 84),
    "frac_test_beta_fixed_alt_beats_beta1_alt": float(np.mean(valid["dAIC_beta_fixed_alt_minus_beta1_alt"] < 0.0)),

    "test_chi2red_beta_train_median": pct(valid["test_chi2red_beta_train_same_a0"], 50),
    "test_chi2red_beta1_same_a0_median": pct(valid["test_chi2red_beta1_same_a0"], 50),
    "test_chi2red_beta_fixed_ref_median": pct(valid["test_chi2red_beta_fixed_ref"], 50),
    "test_chi2red_beta1_ref_median": pct(valid["test_chi2red_beta1_ref"], 50),
    "test_chi2red_beta_fixed_alt_median": pct(valid["test_chi2red_beta_fixed_alt"], 50),
    "test_chi2red_beta1_alt_median": pct(valid["test_chi2red_beta1_alt"], 50),
}

criteria = {
    "enough_valid_splits": summary["N_splits_valid"] >= 150,
    "train_beta_gt_1_majority": summary["frac_train_beta_gt_1"] > 0.80,
    "train_beta_gt_1p03_majority": summary["frac_train_beta_gt_1p03"] > 0.60,
    "train_beta_near_1055_majority": summary["frac_train_beta_near_1055"] > 0.50,

    "test_same_a0_beats_beta1_majority": summary["frac_test_beta_train_beats_beta1_same_a0"] > 0.55,
    "test_same_a0_median_negative": summary["dAIC_beta_train_same_a0_median"] < 0.0,

    "test_fixed_ref_beats_beta1_majority": summary["frac_test_beta_fixed_beats_beta1_ref"] > 0.55,
    "test_fixed_ref_median_negative": summary["dAIC_beta_fixed_ref_median"] < 0.0,

    "test_fixed_alt_not_destroyed": summary["dAIC_beta_fixed_alt_median"] < 50.0,
    "test_chi2red_reasonable": summary["test_chi2red_beta_fixed_ref_median"] < 12.0,
}

score = sum(bool(v) for v in criteria.values())

if score >= 8:
    status = "PASSOU_FORTE_RMAX13_CROSS_VALIDATION_BETA_A0"
elif score >= 6:
    status = "PASSOU_MODERADO_FORTE_RMAX13_COM_ALERTAS"
elif score >= 4:
    status = "PASSOU_MODERADO_RMAX13_BETA_GENERALIZA_PARCIALMENTE"
else:
    status = "SINAL_FRACO_RMAX13_BETA_NAO_GENERALIZA"

summary["FINAL_STATUS"] = status
summary["score"] = int(score)
summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "RMAX13_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — RMAX 13")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.hist(valid["train_beta"], bins=25)
plt.axvline(1.0, linestyle="--", label="beta=1")
plt.axvline(BETA_FIXED, linestyle=":", label=f"beta fixed={BETA_FIXED}")
plt.axvline(summary["train_beta_median"], linestyle="-.", label=f"median={summary['train_beta_median']:.3f}")
plt.xlabel(r"$\beta_{\rm train}$")
plt.ylabel("N splits")
plt.title("RMAX 13 — Train beta distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX13_train_beta_hist.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(valid["dAIC_beta_train_minus_beta1_same_a0"], bins=30)
plt.axvline(0.0, linestyle="--", label="equal AIC")
plt.xlabel(r"$\Delta{\rm AIC}=\mathrm{AIC}(\beta_{\rm train})-\mathrm{AIC}(\beta=1)$")
plt.ylabel("N splits")
plt.title("RMAX 13 — Test set, same a0")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX13_delta_AIC_same_a0.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(valid["dAIC_beta_fixed_minus_beta1_ref"], bins=30)
plt.axvline(0.0, linestyle="--", label="equal AIC")
plt.xlabel(r"$\Delta{\rm AIC}=\mathrm{AIC}(\beta=1.055)-\mathrm{AIC}(\beta=1)$")
plt.ylabel("N splits")
plt.title("RMAX 13 — Fixed beta, reference a0")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX13_delta_AIC_fixed_beta_ref_a0.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.scatter(valid["train_beta"], valid["dAIC_beta_train_minus_beta1_same_a0"], s=15, alpha=0.6)
plt.axhline(0.0, linestyle="--")
plt.axvline(1.0, linestyle="--")
plt.xlabel(r"$\beta_{\rm train}$")
plt.ylabel(r"$\Delta{\rm AIC}$ on test set")
plt.title("RMAX 13 — Generalization of trained beta")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX13_beta_vs_test_delta_AIC.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste RMAX 13.")
