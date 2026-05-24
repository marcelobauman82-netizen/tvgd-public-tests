# ============================================================
# TVGD — 11_RMAX_SPARC_UNIVERSALITY
# Universalidade de P(u)=1-exp(-u) em SPARC / Rotmod
#
# Objetivo:
#   Testar se a lei estrutural mínima da TVGD,
#
#       P(u) = 1 - exp(-u)
#
#   continua competitiva na janela crítica de transição:
#
#       0.5 <= u <= 2.0
#
#   usando curvas de rotação SPARC públicas.
#
# Saída:
#   results/rmax_sparc_universality/
#
# Uso:
#   python notebooks/main_validation/11_rmax_sparc_universality.py
#
# ============================================================

import os
import io
import zipfile
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "11_rmax_sparc_universality.py"

OUTDIR = "results/rmax_sparc_universality"
os.makedirs(OUTDIR, exist_ok=True)

SPARC_URL = "https://astroweb.case.edu/SPARC/Rotmod_LTG.zip"
ZIP_PATH = os.path.join(OUTDIR, "Rotmod_LTG.zip")

A0_REF = 1.2e-10
KPC_TO_M = 3.0856775814913673e19

U_MIN = 0.5
U_MAX = 2.0

# Grade simples para beta. Mantemos pequena para rodar rápido.
BETA_GRID = np.round(np.arange(0.90, 1.151, 0.005), 3)

# Piso de erro fracional em velocidade.
# SPARC tem incertezas heterogêneas; esse piso evita peso infinito.
EV_FRAC_FLOOR = 0.05


print("=" * 72)
print("TVGD — 11_RMAX_SPARC_UNIVERSALITY")
print("Universalidade de P(u)=1-exp(-u) em SPARC / Rotmod")
print("=" * 72)


# ============================================================
# Download SPARC
# ============================================================

def download_if_needed(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 1000:
        print(f"Arquivo SPARC já existe: {path}")
        return

    print(f"Baixando SPARC Rotmod de:\n{url}")
    urllib.request.urlretrieve(url, path)

    if not os.path.exists(path) or os.path.getsize(path) < 1000:
        raise RuntimeError("Falha no download do Rotmod_LTG.zip.")

    print(f"Download OK: {path}")
    print(f"Tamanho: {os.path.getsize(path)} bytes")


download_if_needed(SPARC_URL, ZIP_PATH)


# ============================================================
# Leitura robusta dos arquivos .dat
# ============================================================

def read_rotmod_zip(zip_path):
    rows = []

    with zipfile.ZipFile(zip_path, "r") as z:
        dat_files = [f for f in z.namelist() if f.lower().endswith(".dat")]

        print(f"Arquivos .dat encontrados: {len(dat_files)}")

        for fname in dat_files:
            galaxy = os.path.basename(fname).replace("_rotmod.dat", "").replace(".dat", "")

            raw = z.read(fname).decode("utf-8", errors="ignore")
            lines = raw.splitlines()

            for line in lines:
                s = line.strip()

                if not s:
                    continue
                if s.startswith("#"):
                    continue

                parts = s.split()

                # Formato típico SPARC Rotmod:
                # Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul
                if len(parts) < 6:
                    continue

                try:
                    nums = [float(x) for x in parts[:6]]
                except Exception:
                    continue

                rad, vobs, evobs, vgas, vdisk, vbul = nums

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
        raise RuntimeError("Nenhum ponto SPARC foi lido. Verifique o arquivo Rotmod.")

    return df


df = read_rotmod_zip(ZIP_PATH)

print(f"Galáxias lidas: {df['galaxy'].nunique()}")
print(f"Pontos totais: {len(df)}")


# ============================================================
# Limpeza e cálculo de acelerações
# ============================================================

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["R_kpc", "Vobs", "eVobs", "Vgas", "Vdisk", "Vbul"])

# Cortes físicos básicos
df = df[
    (df["R_kpc"] > 0) &
    (df["Vobs"] > 0)
].copy()

# Erro mínimo conservador
df["eVobs"] = pd.to_numeric(df["eVobs"], errors="coerce")
df["eVobs"] = df["eVobs"].fillna(EV_FRAC_FLOOR * df["Vobs"])
df["eVobs"] = np.maximum(df["eVobs"], EV_FRAC_FLOOR * df["Vobs"])

# Componentes bariônicas.
# Convenção usual SPARC: Vgas pode ser negativo em alguns pontos por convenção.
# A contribuição em aceleração entra via V^2.
df["Vbar2"] = df["Vgas"]**2 + df["Vdisk"]**2 + df["Vbul"]**2
df["Vbar"] = np.sqrt(np.maximum(df["Vbar2"], 0.0))

df = df[df["Vbar"] > 0].copy()

df["gobs"] = (df["Vobs"] * 1000.0)**2 / (df["R_kpc"] * KPC_TO_M)
df["gbar"] = (df["Vbar"] * 1000.0)**2 / (df["R_kpc"] * KPC_TO_M)

# Propagação aproximada: g = V^2/R -> sigma_g/g ≈ 2 sigma_V/V
df["egobs"] = df["gobs"] * 2.0 * df["eVobs"] / np.maximum(df["Vobs"], 1e-12)
df["egobs"] = np.maximum(df["egobs"], 0.10 * df["gobs"])

df["u"] = np.sqrt(df["gbar"] / A0_REF)

df = df.replace([np.inf, -np.inf], np.nan)
df = df.dropna(subset=["gobs", "gbar", "egobs", "u"])
df = df[(df["gobs"] > 0) & (df["gbar"] > 0) & (df["egobs"] > 0)].copy()

print(f"Pontos após sanidade: {len(df)}")
print(f"u min/med/max: {df['u'].min():.6g}, {df['u'].median():.6g}, {df['u'].max():.6g}")


# ============================================================
# Janela crítica de transição
# ============================================================

trans = df[(df["u"] >= U_MIN) & (df["u"] <= U_MAX)].copy()

print("\nJanela de transição:")
print(f"U_MIN = {U_MIN}")
print(f"U_MAX = {U_MAX}")
print(f"Pontos: {len(trans)}")
print(f"Galáxias: {trans['galaxy'].nunique()}")

if len(trans) < 100 or trans["galaxy"].nunique() < 10:
    raise RuntimeError("A janela de transição ficou pequena demais para o teste RMAX.")


# ============================================================
# Funções de resposta
# ============================================================

def P_exp(u, beta=1.0):
    u = np.asarray(u)
    return 1.0 - np.exp(-np.power(u, beta))


def P_tanh(u):
    return np.tanh(u)


def P_rational(u):
    return u / (1.0 + u)


def P_sqrt(u):
    return u / np.sqrt(1.0 + u*u)


def predict_g(gbar, a0=A0_REF, model="exp", beta=1.0):
    u = np.sqrt(np.maximum(gbar, 0.0) / a0)

    if model == "exp":
        P = P_exp(u, beta=beta)
    elif model == "tanh":
        P = P_tanh(u)
    elif model == "rational":
        P = P_rational(u)
    elif model == "sqrt":
        P = P_sqrt(u)
    elif model == "bar":
        return gbar
    else:
        raise ValueError(f"Modelo desconhecido: {model}")

    P = np.maximum(P, 1e-12)
    return gbar / P


# ============================================================
# Ajuste de amplitude por galáxia
# ============================================================

def fit_amplitude(gobs, gpred, err):
    w = 1.0 / np.maximum(err, 1e-99)**2

    den = np.sum(w * gpred**2)
    if den <= 0 or not np.isfinite(den):
        return 0.0, np.inf

    A = np.sum(w * gobs * gpred) / den
    A = max(float(A), 0.0)

    chi2 = np.sum(((gobs - A*gpred) / err)**2)
    return A, float(chi2)


def evaluate_model(data, model="exp", beta=1.0):
    work = data.copy()
    work["gpred_raw"] = predict_g(
        work["gbar"].values,
        a0=A0_REF,
        model=model,
        beta=beta,
    )

    chi2_total = 0.0
    n_used = 0
    amps = []
    rows = []

    for gal, sub in work.groupby("galaxy"):
        if len(sub) < 3:
            continue

        gobs = sub["gobs"].values
        gpred = sub["gpred_raw"].values
        err = sub["egobs"].values

        A, chi2 = fit_amplitude(gobs, gpred, err)

        if not np.isfinite(chi2):
            continue

        n = len(sub)
        n_used += n
        chi2_total += chi2
        amps.append(A)

        rows.append({
            "galaxy": gal,
            "N": n,
            "A": A,
            "chi2": chi2,
            "chi2red_gal": chi2 / max(n - 1, 1),
        })

    gal_stats = pd.DataFrame(rows)

    n_gal = len(gal_stats)
    if n_used == 0 or n_gal == 0:
        return None, gal_stats

    # Parâmetros:
    #   uma amplitude por galáxia
    #   beta conta como parâmetro de forma se for modelo exp_beta livre
    k_amp = n_gal
    k_shape = 1 if (model == "exp_beta") else 0
    k_params = k_amp + k_shape
    dof = max(n_used - k_params, 1)

    AIC = chi2_total + 2 * k_params
    BIC = chi2_total + np.log(max(n_used, 2)) * k_params

    summary = {
        "model": model,
        "beta": beta if "exp" in model else np.nan,
        "N": int(n_used),
        "N_gal": int(n_gal),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "AIC": float(AIC),
        "BIC": float(BIC),
        "k_params": int(k_params),
        "median_chi2red_gal": float(gal_stats["chi2red_gal"].median()),
        "median_A_gal": float(np.median(amps)),
        "u_min": float(work["u"].min()),
        "u_med": float(work["u"].median()),
        "u_max": float(work["u"].max()),
    }

    return summary, gal_stats


# ============================================================
# Rodar modelos
# ============================================================

all_results = []

# Modelo canônico beta=1
row, gal_stats_exp = evaluate_model(trans, model="exp", beta=1.0)
row["label"] = "exp_beta1"
all_results.append(row)

# Família beta
beta_rows = []
for beta in BETA_GRID:
    row_b, _ = evaluate_model(trans, model="exp_beta", beta=float(beta))
    beta_rows.append(row_b)

beta_df = pd.DataFrame(beta_rows).sort_values("AIC").reset_index(drop=True)
best_beta = float(beta_df.iloc[0]["beta"])

row_best_beta = beta_df.iloc[0].to_dict()
row_best_beta["label"] = f"exp_beta{best_beta:.3f}"
all_results.append(row_best_beta)

# Controles
for model_name in ["tanh", "rational", "sqrt", "bar"]:
    row_c, _ = evaluate_model(trans, model=model_name, beta=1.0)
    row_c["label"] = model_name
    all_results.append(row_c)

res = pd.DataFrame(all_results).sort_values("AIC").reset_index(drop=True)

# Salvar tabelas
df.to_csv(os.path.join(OUTDIR, "RMAX11_sparc_all_points.csv"), index=False)
trans.to_csv(os.path.join(OUTDIR, "RMAX11_transition_window.csv"), index=False)
beta_df.to_csv(os.path.join(OUTDIR, "RMAX11_beta_scan.csv"), index=False)
res.to_csv(os.path.join(OUTDIR, "RMAX11_model_comparison.csv"), index=False)

print("\n" + "=" * 72)
print("RESULTADOS — RMAX 11")
print("=" * 72)
print(res.to_string(index=False))


# ============================================================
# Diagnóstico
# ============================================================

exp1 = res[res["label"] == "exp_beta1"].iloc[0]
best = res.iloc[0]
best_exp = beta_df.iloc[0]

delta_AIC_exp1_minus_best = float(exp1["AIC"] - best["AIC"])
delta_AIC_exp1_minus_best_exp = float(exp1["AIC"] - best_exp["AIC"])

criteria = {
    "has_transition_sample": len(trans) >= 500,
    "has_many_galaxies": trans["galaxy"].nunique() >= 50,
    "exp_available": True,
    "best_beta_near_1": 0.90 <= best_beta <= 1.15,
    "best_beta_close_to_transition_result": 1.03 <= best_beta <= 1.11,
    "exp_beta1_not_failed": exp1["chi2red"] < 15,
    "exp_family_beats_controls_or_competitive": best["label"].startswith("exp"),
    "transition_window_valid": (trans["u"].min() >= U_MIN) and (trans["u"].max() <= U_MAX),
}

score = sum(bool(v) for v in criteria.values())

if score >= 7:
    status = "PASSOU_FORTE_RMAX11_UNIVERSALIDADE_SPARK_TRANSICAO"
elif score >= 5:
    status = "PASSOU_MODERADO_RMAX11_COM_ALERTAS"
elif score >= 3:
    status = "SINAL_FRACO_RMAX11_TRANSICAO_PRESSIONA_MODELO"
else:
    status = "FALHOU_RMAX11_UNIVERSALIDADE"

summary = {
    "FINAL_STATUS": status,
    "score": int(score),
    "script": SCRIPT_NAME,
    "N_total_points": int(len(df)),
    "N_total_galaxies": int(df["galaxy"].nunique()),
    "N_transition_points": int(len(trans)),
    "N_transition_galaxies": int(trans["galaxy"].nunique()),
    "u_min_transition": float(trans["u"].min()),
    "u_median_transition": float(trans["u"].median()),
    "u_max_transition": float(trans["u"].max()),
    "best_model": str(best["label"]),
    "best_beta": float(best_beta),
    "best_chi2red": float(best["chi2red"]),
    "best_AIC": float(best["AIC"]),
    "exp_beta1_chi2red": float(exp1["chi2red"]),
    "exp_beta1_AIC": float(exp1["AIC"]),
    "delta_AIC_exp1_minus_best": delta_AIC_exp1_minus_best,
    "delta_AIC_exp1_minus_best_exp": delta_AIC_exp1_minus_best_exp,
}

summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "RMAX11_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — RMAX 11")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(beta_df["beta"], beta_df["AIC"], marker="o", markersize=3)
plt.axvline(1.0, linestyle="--", label="beta=1")
plt.axvline(best_beta, linestyle=":", label=f"best beta={best_beta:.3f}")
plt.xlabel(r"$\beta$")
plt.ylabel("AIC")
plt.title("RMAX 11 — Beta scan in SPARC transition window")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX11_beta_scan_AIC.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.bar(res["label"], res["chi2red"])
plt.xticks(rotation=35, ha="right")
plt.ylabel(r"$\chi^2_{\rm red}$")
plt.title("RMAX 11 — Model comparison")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX11_model_comparison_chi2red.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(trans["u"], bins=40)
plt.xlabel(r"$u=\sqrt{g_{\rm bar}/a_0}$")
plt.ylabel("N")
plt.title("RMAX 11 — Transition window distribution")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "RMAX11_transition_u_distribution.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste RMAX 11.")
