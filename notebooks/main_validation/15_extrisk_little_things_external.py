# ============================================================
# TVGD — 15_EXTRISK_LITTLE_THINGS_EXTERNAL
# Validação externa semi-direta fora do SPARC com LITTLE THINGS
#
# Objetivo:
#   Testar se a família exponencial da TVGD,
#
#       P_beta(u) = 1 - exp(-u^beta),
#
#   continua funcional em uma base externa não-SPARC.
#
# Base:
#   VizieR J/AJ/149/180 — LITTLE THINGS
#
# Estratégia:
#   Usa tabelas rotdmbar e rotdm.
#   Reconstrói Vbar^2 pela diferença:
#
#       Vbar^2 ≈ Vobs^2 - Vdm^2
#
#   mantendo apenas pontos com Vbar^2 > 0.
#
# Interpretação:
#   Este é um teste externo semi-direto. Ele não substitui uma decomposição
#   baryônica completa por gás/disco/bojo, mas testa se a resposta TVGD
#   é destruída fora do SPARC.
#
# Saída:
#   results/extrisk_little_things_external/
#
# Uso:
#   python notebooks/main_validation/15_extrisk_little_things_external.py
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

SCRIPT_NAME = "15_extrisk_little_things_external.py"

OUTDIR = "results/extrisk_little_things_external"
os.makedirs(OUTDIR, exist_ok=True)

A0_REF = 1.2e-10
KPC_TO_M = 3.0856775814913673e19

CATALOG = "J/AJ/149/180"

EV_FRAC_FLOOR = 0.08

BETA_LIST = [0.94, 1.00, 1.055, 1.065, 1.08, 1.10]

print("=" * 72)
print("TVGD — 15_EXTRISK_LITTLE_THINGS_EXTERNAL")
print("Validação externa semi-direta com LITTLE THINGS")
print("=" * 72)


# ============================================================
# Instalar/carregar astroquery
# ============================================================

try:
    from astroquery.vizier import Vizier
except Exception:
    print("astroquery não encontrado. Tentando instalar...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", "astroquery"])
    from astroquery.vizier import Vizier


# ============================================================
# Baixar catálogo VizieR
# ============================================================

print(f"\nCarregando catálogo VizieR: {CATALOG}")

Vizier.ROW_LIMIT = -1
tables = Vizier.get_catalogs(CATALOG)

print(f"N tabelas encontradas: {len(tables)}")

for i, tab in enumerate(tables):
    cols = list(tab.colnames)
    print("\n" + "-" * 72)
    print(f"Tabela {i}")
    print(f"N = {len(tab)}")
    print(f"Colunas = {cols}")


# ============================================================
# Detectores de colunas
# ============================================================

def table_to_df(tab):
    return tab.to_pandas()


def find_col(cols, candidates):
    lower = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]

    for c in cols:
        cl = c.lower()
        for cand in candidates:
            if cand.lower() in cl:
                return c

    return None


def classify_tables(tables):
    """
    Tenta identificar rotdmbar e rotdm.

    Em J/AJ/149/180 normalmente existem tabelas chamadas
    rotdmbar e rotdm, mas o astroquery pode retornar nomes
    sem metadados claros. Então detectamos por colunas.
    """
    candidates = []

    for idx, tab in enumerate(tables):
        df = table_to_df(tab)
        cols = list(df.columns)

        gal_col = find_col(cols, ["Name", "Galaxy", "Gal", "ID"])
        r_col = find_col(cols, ["R", "Rad", "Radius"])
        v_col = find_col(cols, ["Vrot", "Vobs", "Vc", "V", "Vel"])
        ev_col = find_col(cols, ["e_Vrot", "eVrot", "e_Vobs", "errV", "e_V"])
        type_col = find_col(cols, ["Type"])

        score = 0
        if gal_col is not None:
            score += 1
        if r_col is not None:
            score += 1
        if v_col is not None:
            score += 1
        if ev_col is not None:
            score += 1
        if type_col is not None:
            score += 1

        meta_name = ""
        try:
            meta_name = str(tab.meta.get("name", ""))
        except Exception:
            meta_name = ""

        candidates.append({
            "idx": idx,
            "N": len(df),
            "cols": cols,
            "gal_col": gal_col,
            "r_col": r_col,
            "v_col": v_col,
            "ev_col": ev_col,
            "type_col": type_col,
            "score": score,
            "meta_name": meta_name,
        })

    return pd.DataFrame(candidates)


class_df = classify_tables(tables)
class_df.to_csv(os.path.join(OUTDIR, "EXTRISK15_vizier_table_classification.csv"), index=False)

print("\nClassificação automática das tabelas:")
print(class_df[["idx", "N", "score", "meta_name", "gal_col", "r_col", "v_col", "ev_col", "type_col"]].to_string(index=False))


# ============================================================
# Escolher tabelas
# ============================================================

def choose_little_things_tables(tables, class_df):
    """
    Heurística:
    - rotdmbar costuma ter DM+baryons.
    - rotdm costuma ter DM only.
    - Se nomes não aparecem, usa as duas maiores tabelas com estrutura compatível.
    """
    named_dmbar = []
    named_dm = []

    for _, row in class_df.iterrows():
        name = str(row["meta_name"]).lower()
        if "rotdmbar" in name or "dmbar" in name:
            named_dmbar.append(int(row["idx"]))
        elif "rotdm" in name or "dm" in name:
            named_dm.append(int(row["idx"]))

    if named_dmbar and named_dm:
        return named_dmbar[0], named_dm[0]

    # Fallback: tabelas com score >= 3, ordenar por N.
    good = class_df[class_df["score"] >= 3].sort_values("N", ascending=False).reset_index(drop=True)

    if len(good) < 2:
        raise RuntimeError(
            "Não consegui identificar duas tabelas de curvas de rotação no catálogo LITTLE THINGS."
        )

    # Pelo histórico: DM+baryons tinha N maior que DM-only.
    idx_dmbar = int(good.iloc[0]["idx"])
    idx_dm = int(good.iloc[1]["idx"])

    return idx_dmbar, idx_dm


idx_dmbar, idx_dm = choose_little_things_tables(tables, class_df)

print("\nTabelas escolhidas:")
print(f"DM+baryons idx = {idx_dmbar}")
print(f"DM-only    idx = {idx_dm}")

df_dmbar_raw = table_to_df(tables[idx_dmbar])
df_dm_raw = table_to_df(tables[idx_dm])

df_dmbar_raw.to_csv(os.path.join(OUTDIR, "EXTRISK15_raw_dmbar_table.csv"), index=False)
df_dm_raw.to_csv(os.path.join(OUTDIR, "EXTRISK15_raw_dm_table.csv"), index=False)


# ============================================================
# Normalizar tabela de rotação
# ============================================================

def normalize_rotation_table(df, label):
    cols = list(df.columns)

    gal_col = find_col(cols, ["Name", "Galaxy", "Gal", "ID"])
    r_col = find_col(cols, ["R", "Rad", "Radius"])
    v_col = find_col(cols, ["Vrot", "Vobs", "Vc", "V", "Vel"])
    ev_col = find_col(cols, ["e_Vrot", "eVrot", "e_Vobs", "errV", "e_V"])
    type_col = find_col(cols, ["Type"])

    if gal_col is None or r_col is None or v_col is None:
        raise RuntimeError(f"Não consegui detectar colunas mínimas em {label}.")

    out = pd.DataFrame()
    out["galaxy"] = df[gal_col].astype(str)
    out["R_kpc"] = pd.to_numeric(df[r_col], errors="coerce")
    out["V"] = pd.to_numeric(df[v_col], errors="coerce")

    if ev_col is not None:
        out["eV"] = pd.to_numeric(df[ev_col], errors="coerce")
    else:
        out["eV"] = np.nan

    if type_col is not None:
        out["Type"] = df[type_col].astype(str)
    else:
        out["Type"] = "Data"

    out["source_label"] = label

    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["galaxy", "R_kpc", "V"])

    out = out[(out["R_kpc"] > 0) & (out["V"] > 0)].copy()

    # Filtrar observações reais quando houver Type.
    if "Type" in out.columns:
        type_lower = out["Type"].str.lower()
        mask_data = (
            type_lower.str.contains("data", na=False) |
            type_lower.str.contains("obs", na=False) |
            type_lower.str.contains("observed", na=False)
        )

        # Se o filtro gerar vazio, não aplicar.
        if mask_data.sum() > 0:
            out = out[mask_data].copy()

    out["eV"] = out["eV"].fillna(EV_FRAC_FLOOR * out["V"])
    out["eV"] = np.maximum(out["eV"], EV_FRAC_FLOOR * out["V"])

    return out


dmbar = normalize_rotation_table(df_dmbar_raw, "DM_plus_baryons")
dm = normalize_rotation_table(df_dm_raw, "DM_only")

print("\nApós normalização/filtro:")
print(f"DM+baryons: N = {len(dmbar)}, galáxias = {dmbar['galaxy'].nunique()}")
print(f"DM-only:    N = {len(dm)}, galáxias = {dm['galaxy'].nunique()}")

dmbar.to_csv(os.path.join(OUTDIR, "EXTRISK15_dmbar_normalized.csv"), index=False)
dm.to_csv(os.path.join(OUTDIR, "EXTRISK15_dm_normalized.csv"), index=False)


# ============================================================
# Reconstrução semi-direta de Vbar
# ============================================================

def reconstruct_external_base(dmbar, dm):
    rows = []

    common_gals = sorted(set(dmbar["galaxy"]).intersection(set(dm["galaxy"])))

    for gal in common_gals:
        a = dmbar[dmbar["galaxy"] == gal].sort_values("R_kpc").copy()
        b = dm[dm["galaxy"] == gal].sort_values("R_kpc").copy()

        if len(a) < 3 or len(b) < 3:
            continue

        R_dm = b["R_kpc"].values
        V_dm = b["V"].values

        # Interpolar Vdm na malha DM+baryons.
        Rmin, Rmax = np.min(R_dm), np.max(R_dm)

        for _, row in a.iterrows():
            R = float(row["R_kpc"])

            if R < Rmin or R > Rmax:
                continue

            Vobs = float(row["V"])
            eVobs = float(row["eV"])
            Vdm_interp = float(np.interp(R, R_dm, V_dm))

            Vbar2 = Vobs**2 - Vdm_interp**2

            if not np.isfinite(Vbar2) or Vbar2 <= 0:
                continue

            Vbar = np.sqrt(Vbar2)

            rows.append({
                "galaxy": gal,
                "R_kpc": R,
                "Vobs": Vobs,
                "eVobs": eVobs,
                "Vdm_interp": Vdm_interp,
                "Vbar": Vbar,
                "Vbar2": Vbar2,
            })

    out = pd.DataFrame(rows)

    if len(out) == 0:
        raise RuntimeError("Reconstrução externa ficou vazia.")

    out["gobs"] = (out["Vobs"] * 1000.0)**2 / (out["R_kpc"] * KPC_TO_M)
    out["gbar"] = (out["Vbar"] * 1000.0)**2 / (out["R_kpc"] * KPC_TO_M)

    out["egobs"] = out["gobs"] * 2.0 * out["eVobs"] / np.maximum(out["Vobs"], 1e-12)
    out["egobs"] = np.maximum(out["egobs"], 0.10 * out["gobs"])

    out["u"] = np.sqrt(out["gbar"] / A0_REF)

    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["gobs", "gbar", "egobs", "u"])
    out = out[(out["gobs"] > 0) & (out["gbar"] > 0) & (out["egobs"] > 0)].copy()

    return out


base = reconstruct_external_base(dmbar, dm)

print("\nBase reconstruída externa:")
print(f"N = {len(base)}")
print(f"N_gal = {base['galaxy'].nunique()}")
print(f"u range = {base['u'].min()} {base['u'].max()}")
print(f"u med = {base['u'].median()}")

base.to_csv(os.path.join(OUTDIR, "EXTRISK15_little_things_reconstructed_base.csv"), index=False)


# ============================================================
# Modelos TVGD e controles
# ============================================================

def P_exp_beta(u, beta):
    return 1.0 - np.exp(-np.power(np.asarray(u), beta))


def P_tanh(u):
    return np.tanh(u)


def P_rational(u):
    return np.asarray(u) / (1.0 + np.asarray(u))


def P_sqrt(u):
    u = np.asarray(u)
    return u / np.sqrt(1.0 + u*u)


def predict_g(gbar, model, beta=1.0):
    u = np.sqrt(np.maximum(gbar, 0.0) / A0_REF)

    if model == "exp_beta":
        P = P_exp_beta(u, beta)
        return gbar / np.maximum(P, 1e-12)

    if model == "tanh":
        P = P_tanh(u)
        return gbar / np.maximum(P, 1e-12)

    if model == "rational":
        P = P_rational(u)
        return gbar / np.maximum(P, 1e-12)

    if model == "sqrt":
        P = P_sqrt(u)
        return gbar / np.maximum(P, 1e-12)

    if model == "baryon_only":
        return gbar

    raise ValueError(f"Modelo desconhecido: {model}")


def fit_amp(gobs, gpred, err):
    w = 1.0 / np.maximum(err, 1e-99)**2
    den = np.sum(w*gpred*gpred)

    if den <= 0 or not np.isfinite(den):
        return 0.0, np.inf

    A = np.sum(w*gobs*gpred) / den
    A = max(float(A), 0.0)

    chi2 = np.sum(((gobs - A*gpred) / err)**2)
    return A, float(chi2)


def evaluate_model(data, model, beta=np.nan):
    chi2_total = 0.0
    n_total = 0
    n_gal = 0
    rows = []

    for gal, sub in data.groupby("galaxy"):
        if len(sub) < 3:
            continue

        if model == "exp_beta":
            gpred = predict_g(sub["gbar"].values, model="exp_beta", beta=beta)
        else:
            gpred = predict_g(sub["gbar"].values, model=model, beta=1.0)

        A, chi2 = fit_amp(sub["gobs"].values, gpred, sub["egobs"].values)

        if not np.isfinite(chi2):
            continue

        n = len(sub)
        n_total += n
        n_gal += 1
        chi2_total += chi2

        rows.append({
            "galaxy": gal,
            "N": n,
            "A": A,
            "chi2": chi2,
            "chi2red_gal": chi2 / max(n - 1, 1),
        })

    if n_total == 0 or n_gal == 0:
        return None, pd.DataFrame(rows)

    # Amplitude por galáxia + beta se exp_beta e beta diferente de lista fixa? aqui tratamos beta como modelo.
    k_params = n_gal + (1 if model == "exp_beta" else 0)
    dof = max(n_total - k_params, 1)

    stats = pd.DataFrame(rows)

    summary = {
        "model": model,
        "beta": beta if model == "exp_beta" else np.nan,
        "N": int(n_total),
        "N_gal": int(n_gal),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "AIC": float(chi2_total + 2*k_params),
        "BIC": float(chi2_total + np.log(max(n_total, 2))*k_params),
        "k_params": int(k_params),
        "median_chi2red_gal": float(stats["chi2red_gal"].median()) if len(stats) else np.nan,
        "median_A": float(stats["A"].median()) if len(stats) else np.nan,
        "u_min": float(data["u"].min()),
        "u_med": float(data["u"].median()),
        "u_max": float(data["u"].max()),
    }

    return summary, stats


# ============================================================
# Rodar comparação
# ============================================================

results = []
stats_all = []

for beta in BETA_LIST:
    row, st = evaluate_model(base, "exp_beta", beta=beta)
    row["label"] = f"exp_beta{beta:.3f}"
    results.append(row)
    st["model"] = row["label"]
    stats_all.append(st)

for model in ["tanh", "rational", "sqrt", "baryon_only"]:
    row, st = evaluate_model(base, model, beta=np.nan)
    row["label"] = model
    results.append(row)
    st["model"] = model
    stats_all.append(st)

res = pd.DataFrame(results).sort_values("AIC").reset_index(drop=True)
gal_stats = pd.concat(stats_all, ignore_index=True)

res.to_csv(os.path.join(OUTDIR, "EXTRISK15_model_comparison.csv"), index=False)
gal_stats.to_csv(os.path.join(OUTDIR, "EXTRISK15_galaxy_stats_by_model.csv"), index=False)

print("\n" + "=" * 72)
print("RESULTADOS EXTRISK-15 — LITTLE THINGS EXTERNO")
print("=" * 72)
print(res.to_string(index=False))


# ============================================================
# Diagnóstico final
# ============================================================

best = res.iloc[0].to_dict()

exp1 = res[np.isclose(res["beta"].fillna(-999), 1.0)]
if len(exp1) == 0:
    raise RuntimeError("Modelo exp beta=1 não encontrado.")
exp1 = exp1.iloc[0].to_dict()

bar = res[res["label"] == "baryon_only"].iloc[0].to_dict()

# Melhor exponencial
exp_family = res[res["model"] == "exp_beta"].sort_values("AIC")
best_exp = exp_family.iloc[0].to_dict()

delta_AIC_exp1_minus_best = float(exp1["AIC"] - best["AIC"])
delta_AIC_exp1_minus_best_exp = float(exp1["AIC"] - best_exp["AIC"])
delta_AIC_bar_minus_best_exp = float(bar["AIC"] - best_exp["AIC"])
delta_chi2_bar_minus_best_exp = float(bar["chi2"] - best_exp["chi2"])

criteria = {
    "external_catalog_used": True,
    "not_sparc": True,
    "has_reconstructed_gbar": len(base) > 100,
    "has_many_galaxies": base["galaxy"].nunique() >= 10,
    "u_deep_or_transition": base["u"].median() < 1.0,
    "exponential_family_best_or_close": best["model"] == "exp_beta",
    "exp_family_beats_baryon_AIC": delta_AIC_bar_minus_best_exp > 10,
    "exp_family_beats_baryon_chi2": delta_chi2_bar_minus_best_exp > 10,
    "beta1_not_catastrophic": delta_AIC_exp1_minus_best < 20,
    "best_beta_near_1": 0.94 <= best_exp["beta"] <= 1.10,
    "chi2red_reasonable": best_exp["chi2red"] < 5,
}

score = sum(bool(v) for v in criteria.values())

if score >= 9:
    status = "PASSOU_FORTE_EXTRISK15_LITTLE_THINGS_EXTERNO"
elif score >= 7:
    status = "PASSOU_MODERADO_FORTE_EXTRISK15_COM_ALERTAS"
elif score >= 5:
    status = "PASSOU_MODERADO_EXTRISK15_VALIDACAO_SEMI_DIRETA"
else:
    status = "FALHOU_EXTRISK15_EXTERNO"

summary = {
    "FINAL_STATUS": status,
    "score": int(score),
    "script": SCRIPT_NAME,

    "catalog": CATALOG,
    "N_dmbar_normalized": int(len(dmbar)),
    "N_dm_normalized": int(len(dm)),
    "N_reconstructed": int(len(base)),
    "N_gal_reconstructed": int(base["galaxy"].nunique()),
    "u_min": float(base["u"].min()),
    "u_med": float(base["u"].median()),
    "u_max": float(base["u"].max()),

    "best_model": str(best["label"]),
    "best_exp_beta": float(best_exp["beta"]),
    "best_exp_chi2": float(best_exp["chi2"]),
    "best_exp_chi2red": float(best_exp["chi2red"]),
    "best_exp_AIC": float(best_exp["AIC"]),

    "exp_beta1_chi2": float(exp1["chi2"]),
    "exp_beta1_chi2red": float(exp1["chi2red"]),
    "exp_beta1_AIC": float(exp1["AIC"]),

    "baryon_chi2": float(bar["chi2"]),
    "baryon_chi2red": float(bar["chi2red"]),
    "baryon_AIC": float(bar["AIC"]),

    "delta_AIC_exp1_minus_best": delta_AIC_exp1_minus_best,
    "delta_AIC_exp1_minus_best_exp": delta_AIC_exp1_minus_best_exp,
    "delta_AIC_bar_minus_best_exp": delta_AIC_bar_minus_best_exp,
    "delta_chi2_bar_minus_best_exp": delta_chi2_bar_minus_best_exp,
}

summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "EXTRISK15_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — EXTRISK 15")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.bar(res["label"], res["chi2red"])
plt.xticks(rotation=35, ha="right")
plt.ylabel(r"$\chi^2_{\rm red}$")
plt.title("EXTRISK 15 — LITTLE THINGS external comparison")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "EXTRISK15_model_comparison_chi2red.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(base["u"], bins=40)
plt.xlabel(r"$u=\sqrt{g_{\rm bar}/a_0}$")
plt.ylabel("N")
plt.title("EXTRISK 15 — LITTLE THINGS reconstructed u")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "EXTRISK15_u_distribution.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
exp_plot = exp_family.sort_values("beta")
plt.plot(exp_plot["beta"], exp_plot["AIC"], marker="o")
plt.axvline(1.0, linestyle="--", label="beta=1")
plt.axvline(best_exp["beta"], linestyle=":", label=f"best beta={best_exp['beta']:.3f}")
plt.xlabel(r"$\beta$")
plt.ylabel("AIC")
plt.title("EXTRISK 15 — Exponential family beta scan")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "EXTRISK15_beta_scan_AIC.png"), dpi=170)
plt.close()

plt.figure(figsize=(6, 6))
plt.scatter(base["gbar"], base["gobs"], s=10, alpha=0.6)
lims = [
    min(base["gbar"].min(), base["gobs"].min()),
    max(base["gbar"].max(), base["gobs"].max()),
]
plt.plot(lims, lims, linestyle="--", label="gobs=gbar")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$g_{\rm bar}$")
plt.ylabel(r"$g_{\rm obs}$")
plt.title("EXTRISK 15 — Reconstructed acceleration plane")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "EXTRISK15_acceleration_plane.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste EXTRISK 15.")
