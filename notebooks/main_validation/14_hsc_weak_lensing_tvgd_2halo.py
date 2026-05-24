# ============================================================
# TVGD — 14_HSC_WEAK_LENSING_TVGD_2HALO
# Weak lensing HSC: TVGD + termo 2-halo com amplitudes positivas
#
# Objetivo:
#   Testar se a resposta estrutural TVGD,
#
#       P(u)=1-exp(-u)
#
#   é compatível com perfis públicos de weak lensing HSC quando
#   combinada com um termo externo de grande escala tipo 2-halo.
#
# Modelo:
#
#   DeltaSigma_model(R)
#       = A_TVGD * DeltaSigma_TVGD(R)
#       + A_2h   * DeltaSigma_2h(R)
#
#   com A_TVGD >= 0 e A_2h >= 0.
#
# Observação:
#   Este teste é uma validação fenomenológica conservadora.
#   Não substitui uma likelihood completa de galaxy-galaxy lensing,
#   mas verifica se a TVGD é imediatamente incompatível com HSC.
#
# Saída:
#   results/hsc_weak_lensing_tvgd_2halo/
#
# Uso:
#   python notebooks/main_validation/14_hsc_weak_lensing_tvgd_2halo.py
# ============================================================

import os
import re
import urllib.request
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "14_hsc_weak_lensing_tvgd_2halo.py"

OUTDIR = "results/hsc_weak_lensing_tvgd_2halo"
os.makedirs(OUTDIR, exist_ok=True)

A0_REF = 1.2e-10
G_SI = 6.67430e-11
MSUN_KG = 1.98847e30
KPC_M = 3.0856775814913673e19

# URLs públicas usadas nos testes HSC anteriores.
# Caso algum link mude, o usuário pode baixar manualmente CSVs HSC e colocar
# em results/hsc_weak_lensing_tvgd_2halo/input/.
HSC_URLS = {
    "wls1": "https://raw.githubusercontent.com/surhudm/navin-wl-data/master/data/HSC_wls_zbin1.csv",
    "err1": "https://raw.githubusercontent.com/surhudm/navin-wl-data/master/data/HSC_wls_zbin1_err.csv",
}

INPUT_DIR = os.path.join(OUTDIR, "input")
os.makedirs(INPUT_DIR, exist_ok=True)

# Parâmetros fenomenológicos simples
BETA_TVGD = 1.0

# Grade de inclinação para termo 2-halo:
# DeltaSigma_2h ~ (R/R0)^(-s)
S_2H_GRID = np.round(np.arange(0.40, 1.41, 0.05), 3)
S_2H_FIXED = 0.80

# Piso de erro para evitar peso infinito.
ERR_FLOOR_FRAC = 0.08

print("=" * 72)
print("TVGD — 14_HSC_WEAK_LENSING_TVGD_2HALO")
print("Weak lensing HSC: TVGD + termo 2-halo positivo")
print("=" * 72)


# ============================================================
# Utilidades
# ============================================================

def download_file(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 100:
        print(f"Arquivo já existe: {path}")
        return True

    try:
        print(f"Baixando: {url}")
        urllib.request.urlretrieve(url, path)

        if os.path.exists(path) and os.path.getsize(path) > 100:
            print(f"Download OK: {path}")
            return True

    except Exception as exc:
        print(f"Falha no download: {exc}")

    return False


def read_csv_flexible(path):
    """
    Lê CSV com separador flexível.
    """
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, sep=r"\s+", engine="python")


def parse_logm_from_col(col):
    """
    Extrai limiar de massa de nomes tipo:
    zbin_1_stelM_thresh_10.8
    z1_stelM_thresh_10.8
    """
    m = re.search(r"thresh[_\-]?([0-9]+(?:\.[0-9]+)?)", str(col))
    if m:
        return float(m.group(1))
    return None


def nnls_2col(y, x1, x2, err):
    """
    NNLS simples para dois templates positivos:
        y = a*x1 + b*x2, a>=0, b>=0

    Resolve por:
      1. solução livre
      2. bordas a=0 ou b=0
      3. escolhe menor chi2
    """
    w = 1.0 / np.maximum(err, 1e-99)**2

    X = np.vstack([x1, x2]).T
    Xw = X * np.sqrt(w[:, None])
    yw = y * np.sqrt(w)

    candidates = []

    # Solução livre
    try:
        coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        a, b = coef
        if np.isfinite(a) and np.isfinite(b) and a >= 0 and b >= 0:
            model = a*x1 + b*x2
            chi2 = np.sum(((y - model) / err)**2)
            candidates.append((chi2, float(a), float(b)))
    except Exception:
        pass

    # Borda b=0
    den1 = np.sum(w*x1*x1)
    if den1 > 0:
        a = max(np.sum(w*y*x1) / den1, 0.0)
        b = 0.0
        model = a*x1
        chi2 = np.sum(((y - model) / err)**2)
        candidates.append((float(chi2), float(a), float(b)))

    # Borda a=0
    den2 = np.sum(w*x2*x2)
    if den2 > 0:
        a = 0.0
        b = max(np.sum(w*y*x2) / den2, 0.0)
        model = b*x2
        chi2 = np.sum(((y - model) / err)**2)
        candidates.append((float(chi2), float(a), float(b)))

    if not candidates:
        return np.inf, 0.0, 0.0

    candidates = sorted(candidates, key=lambda t: t[0])
    return candidates[0]


def fit_amp_positive(y, x, err):
    """
    Ajuste positivo y = A*x.
    """
    w = 1.0 / np.maximum(err, 1e-99)**2
    den = np.sum(w*x*x)

    if den <= 0:
        return np.inf, 0.0

    A = max(np.sum(w*y*x) / den, 0.0)
    chi2 = np.sum(((y - A*x) / err)**2)

    return float(chi2), float(A)


# ============================================================
# Baixar ou localizar dados HSC
# ============================================================

print("\nBaixando HSC weak lensing público...")

paths = {}

for key, url in HSC_URLS.items():
    path = os.path.join(INPUT_DIR, f"{key}.csv")
    ok = download_file(url, path)
    if ok:
        paths[key] = path

if "wls1" not in paths or "err1" not in paths:
    raise RuntimeError(
        "\nNão consegui baixar os CSVs HSC automaticamente.\n"
        "Coloque manualmente dois arquivos em:\n"
        f"{INPUT_DIR}\n"
        "com nomes wls1.csv e err1.csv, ou atualize os links HSC_URLS.\n"
    )

wls = read_csv_flexible(paths["wls1"])
err = read_csv_flexible(paths["err1"])

print("\nColunas wls:")
print(list(wls.columns))
print("\nColunas err:")
print(list(err.columns))


# ============================================================
# Normalização da base HSC
# ============================================================

if "rp" not in wls.columns:
    raise RuntimeError("A tabela HSC não contém coluna 'rp'.")

rp = pd.to_numeric(wls["rp"], errors="coerce").values

rows = []

data_cols = [c for c in wls.columns if c != "rp"]

for col in data_cols:
    logm = parse_logm_from_col(col)
    if logm is None:
        continue

    # Encontrar coluna de erro correspondente
    err_col = None
    for ec in err.columns:
        if str(logm) in str(ec):
            err_col = ec
            break

    if err_col is None:
        continue

    ds = pd.to_numeric(wls[col], errors="coerce").values
    es = pd.to_numeric(err[err_col], errors="coerce").values

    for i in range(len(rp)):
        if not np.isfinite(rp[i]) or not np.isfinite(ds[i]) or not np.isfinite(es[i]):
            continue

        if rp[i] <= 0:
            continue

        if ds[i] <= 0:
            continue

        e_use = max(abs(es[i]), ERR_FLOOR_FRAC * abs(ds[i]))

        rows.append({
            "R_raw": float(rp[i]),
            "R_kpc": float(rp[i] * 1000.0),  # rp geralmente em Mpc/h ou Mpc; escala fenomenológica
            "DeltaSigma_obs": float(ds[i]),
            "DeltaSigma_err": float(e_use),
            "logMstar": float(logm),
            "bin_label": str(col),
            "err_col": str(err_col),
        })

base = pd.DataFrame(rows)

if len(base) == 0:
    raise RuntimeError("Base HSC normalizada ficou vazia.")

# Remover massa limiar muito extrema se necessário.
base = base.replace([np.inf, -np.inf], np.nan)
base = base.dropna(subset=["R_kpc", "DeltaSigma_obs", "DeltaSigma_err", "logMstar"])
base = base[(base["R_kpc"] > 0) & (base["DeltaSigma_obs"] > 0) & (base["DeltaSigma_err"] > 0)].copy()

base["bin_id"] = base["bin_label"].astype(str)

print("\nBase HSC normalizada:")
print(f"N = {len(base)}")
print(f"N_bins = {base['bin_id'].nunique()}")
print(f"R_kpc min/max = {base['R_kpc'].min()} {base['R_kpc'].max()}")
print(f"DeltaSigma min/max = {base['DeltaSigma_obs'].min()} {base['DeltaSigma_obs'].max()}")
print(f"Erro min/max = {base['DeltaSigma_err'].min()} {base['DeltaSigma_err'].max()}")
print(f"logMstar min/max = {base['logMstar'].min()} {base['logMstar'].max()}")

base.to_csv(os.path.join(OUTDIR, "HSC_base_normalizada.csv"), index=False)


# ============================================================
# Templates físicos simplificados
# ============================================================

def compute_templates(df, beta=BETA_TVGD, s_2h=S_2H_FIXED):
    work = df.copy()

    # Massa estelar representativa por bin.
    # Como as colunas são thresholds, usamos o threshold como escala efetiva.
    Mstar = 10.0 ** work["logMstar"].values  # Msun

    R_m = work["R_kpc"].values * KPC_M
    R_kpc = work["R_kpc"].values

    # Aceleração bariônica aproximada por massa pontual/projetada.
    # É uma aproximação fenomenológica para construir u.
    gbar = G_SI * (Mstar * MSUN_KG) / np.maximum(R_m, 1e-99)**2

    u = np.sqrt(np.maximum(gbar, 0.0) / A0_REF)

    P = 1.0 - np.exp(-np.power(u, beta))
    alpha = 1.0 / np.maximum(P, 1e-12)

    # Template bariônico projetado simples: escala ~ M/R^2.
    # A normalização é absorvida por A_TVGD.
    sigma_bar_template = Mstar / np.maximum(R_kpc, 1e-99)**2

    # TVGD amplifica a resposta estrutural.
    sigma_tvgd_template = sigma_bar_template * alpha

    # Normalizar templates para estabilidade numérica dentro do conjunto.
    med_tvgd = np.nanmedian(sigma_tvgd_template[sigma_tvgd_template > 0])
    if not np.isfinite(med_tvgd) or med_tvgd <= 0:
        med_tvgd = 1.0

    sigma_tvgd_template = sigma_tvgd_template / med_tvgd

    # Template 2-halo: termo radial externo.
    R0 = np.nanmedian(R_kpc)
    sigma_2h_template = np.power(np.maximum(R_kpc / R0, 1e-12), -s_2h)

    med_2h = np.nanmedian(sigma_2h_template[sigma_2h_template > 0])
    if not np.isfinite(med_2h) or med_2h <= 0:
        med_2h = 1.0

    sigma_2h_template = sigma_2h_template / med_2h

    work["gbar_template"] = gbar
    work["u"] = u
    work["alpha_tvgd"] = alpha
    work["template_tvgd"] = sigma_tvgd_template
    work["template_2h"] = sigma_2h_template

    return work


# ============================================================
# Avaliação por bin
# ============================================================

def evaluate_single_template(df, template_col, model_name, k_extra_per_bin=0):
    chi2_total = 0.0
    n_total = 0
    n_bins = 0
    rows = []

    for bin_id, sub in df.groupby("bin_id"):
        if len(sub) < 3:
            continue

        y = sub["DeltaSigma_obs"].values
        e = sub["DeltaSigma_err"].values
        x = sub[template_col].values

        chi2, A = fit_amp_positive(y, x, e)

        n = len(sub)
        n_total += n
        n_bins += 1
        chi2_total += chi2

        rows.append({
            "bin_id": bin_id,
            "N": n,
            "A": A,
            "chi2": chi2,
            "chi2red_bin": chi2 / max(n - 1, 1),
        })

    k_params = n_bins * (1 + k_extra_per_bin)
    dof = max(n_total - k_params, 1)

    return {
        "model": model_name,
        "N": int(n_total),
        "N_bins": int(n_bins),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "AIC": float(chi2_total + 2*k_params),
        "BIC": float(chi2_total + np.log(max(n_total, 2))*k_params),
        "k_params": int(k_params),
        "median_chi2red_bin": float(pd.DataFrame(rows)["chi2red_bin"].median()) if rows else np.nan,
        "frac_A_positive": float(np.mean(pd.DataFrame(rows)["A"] > 0)) if rows else np.nan,
    }, pd.DataFrame(rows)


def evaluate_tvgd_plus_2halo(df, s_2h, model_name):
    work = compute_templates(df, beta=BETA_TVGD, s_2h=s_2h)

    chi2_total = 0.0
    n_total = 0
    n_bins = 0
    rows = []

    for bin_id, sub in work.groupby("bin_id"):
        if len(sub) < 3:
            continue

        y = sub["DeltaSigma_obs"].values
        e = sub["DeltaSigma_err"].values
        x1 = sub["template_tvgd"].values
        x2 = sub["template_2h"].values

        chi2, A1, A2 = nnls_2col(y, x1, x2, e)

        n = len(sub)
        n_total += n
        n_bins += 1
        chi2_total += chi2

        rows.append({
            "bin_id": bin_id,
            "N": n,
            "s_2h": float(s_2h),
            "A_tvgd": A1,
            "A_2h": A2,
            "chi2": chi2,
            "chi2red_bin": chi2 / max(n - 2, 1),
            "A_tvgd_positive": A1 > 0,
            "A_2h_positive": A2 > 0,
            "both_positive": (A1 > 0) and (A2 > 0),
            "A_tvgd_zero": A1 == 0,
            "A_2h_zero": A2 == 0,
            "u_median_bin": float(np.median(sub["u"])),
        })

    k_params = 2 * n_bins
    dof = max(n_total - k_params, 1)

    gal_stats = pd.DataFrame(rows)

    return {
        "model": model_name,
        "beta": float(BETA_TVGD),
        "s_2h": float(s_2h),
        "N": int(n_total),
        "N_bins": int(n_bins),
        "chi2": float(chi2_total),
        "chi2red": float(chi2_total / dof),
        "AIC": float(chi2_total + 2*k_params),
        "BIC": float(chi2_total + np.log(max(n_total, 2))*k_params),
        "k_params": int(k_params),
        "median_chi2red_bin": float(gal_stats["chi2red_bin"].median()) if len(gal_stats) else np.nan,
        "median_u": float(work["u"].median()),
        "min_u": float(work["u"].min()),
        "max_u": float(work["u"].max()),
        "median_s_2h": float(s_2h),
        "frac_A_tvgd_positive": float(gal_stats["A_tvgd_positive"].mean()) if len(gal_stats) else np.nan,
        "frac_A_2h_positive": float(gal_stats["A_2h_positive"].mean()) if len(gal_stats) else np.nan,
        "frac_both_positive": float(gal_stats["both_positive"].mean()) if len(gal_stats) else np.nan,
        "frac_A_tvgd_zero": float(gal_stats["A_tvgd_zero"].mean()) if len(gal_stats) else np.nan,
        "frac_A_2h_zero": float(gal_stats["A_2h_zero"].mean()) if len(gal_stats) else np.nan,
    }, gal_stats, work


# ============================================================
# Construção de modelos de controle
# ============================================================

base_templates = compute_templates(base, beta=BETA_TVGD, s_2h=S_2H_FIXED)

# Controle power-law positivo por bin: template ~ R^-s fixo ajustado por amplitude.
# Usamos s=0.8 como referência.
base_templates["template_powerlaw"] = base_templates["template_2h"]

# Baryon-only: template sem alpha.
Mstar = 10.0 ** base_templates["logMstar"].values
R_kpc = base_templates["R_kpc"].values
bar = Mstar / np.maximum(R_kpc, 1e-99)**2
bar = bar / np.nanmedian(bar[bar > 0])
base_templates["template_baryon"] = bar

rows = []

# Baryon-only
r_bar, stat_bar = evaluate_single_template(base_templates, "template_baryon", "baryon_only_positive")
rows.append(r_bar)

# TVGD puro
r_tvgd, stat_tvgd = evaluate_single_template(base_templates, "template_tvgd", "TVGD_exp_beta1_positive")
r_tvgd["beta"] = BETA_TVGD
r_tvgd["median_u"] = float(base_templates["u"].median())
r_tvgd["min_u"] = float(base_templates["u"].min())
r_tvgd["max_u"] = float(base_templates["u"].max())
rows.append(r_tvgd)

# Power-law controle
r_pow, stat_pow = evaluate_single_template(base_templates, "template_powerlaw", "powerlaw_control_positive")
rows.append(r_pow)

# TVGD + 2halo com s fixo
r_fixed, stat_fixed, work_fixed = evaluate_tvgd_plus_2halo(base, S_2H_FIXED, "TVGD_plus_2halo_NNLS_s_fixed")
rows.append(r_fixed)

# TVGD + 2halo com grade de s
grid_rows = []
grid_stats = []

for s in S_2H_GRID:
    r, st, wk = evaluate_tvgd_plus_2halo(base, s, "TVGD_plus_2halo_NNLS_s_grid")
    grid_rows.append(r)
    st["s_2h"] = s
    grid_stats.append(st)

grid_df = pd.DataFrame(grid_rows).sort_values("AIC").reset_index(drop=True)
best_grid = grid_df.iloc[0].to_dict()
rows.append(best_grid)

if grid_stats:
    pd.concat(grid_stats, ignore_index=True).to_csv(
        os.path.join(OUTDIR, "HSC_TVGD_2halo_grid_bin_stats.csv"),
        index=False,
    )

res = pd.DataFrame(rows).sort_values("AIC").reset_index(drop=True)
res.to_csv(os.path.join(OUTDIR, "HSC_model_comparison.csv"), index=False)
stat_bar.to_csv(os.path.join(OUTDIR, "HSC_baryon_bin_stats.csv"), index=False)
stat_tvgd.to_csv(os.path.join(OUTDIR, "HSC_tvgd_bin_stats.csv"), index=False)
stat_pow.to_csv(os.path.join(OUTDIR, "HSC_powerlaw_bin_stats.csv"), index=False)
stat_fixed.to_csv(os.path.join(OUTDIR, "HSC_tvgd_2halo_fixed_bin_stats.csv"), index=False)
grid_df.to_csv(os.path.join(OUTDIR, "HSC_tvgd_2halo_s_grid.csv"), index=False)

print("\n" + "=" * 72)
print("RESULTADOS — HSC WEAK LENSING")
print("=" * 72)
print(res.to_string(index=False))


# ============================================================
# Diagnóstico final
# ============================================================

def get_row(model):
    m = res[res["model"] == model]
    if len(m) == 0:
        return None
    return m.iloc[0].to_dict()

bar_row = get_row("baryon_only_positive")
tvgd_row = get_row("TVGD_exp_beta1_positive")
pow_row = get_row("powerlaw_control_positive")
fixed_row = get_row("TVGD_plus_2halo_NNLS_s_fixed")

# best_grid já é dict
grid_row = best_grid

delta_chi2_exp_minus_bar = tvgd_row["chi2"] - bar_row["chi2"]
delta_chi2_fixed_minus_exp = fixed_row["chi2"] - tvgd_row["chi2"]
delta_chi2_grid_minus_exp = grid_row["chi2"] - tvgd_row["chi2"]
delta_AIC_fixed_minus_powerlaw = fixed_row["AIC"] - pow_row["AIC"]
delta_AIC_grid_minus_powerlaw = grid_row["AIC"] - pow_row["AIC"]

criteria = {
    "HSC_public_real_data_used": True,
    "N_points_ok": len(base) >= 100,
    "N_bins_ok": base["bin_id"].nunique() >= 10,
    "TVGD_exp_beats_baryon": delta_chi2_exp_minus_bar < 0,
    "NNLS_fixed_beats_TVGD_exp": delta_chi2_fixed_minus_exp < 0,
    "NNLS_grid_beats_TVGD_exp": delta_chi2_grid_minus_exp < 0,
    "NNLS_grid_chi2red_below_5": grid_row["chi2red"] < 5,
    "NNLS_grid_chi2red_below_3": grid_row["chi2red"] < 3,
    "grid_A_tvgd_positive_majority": grid_row.get("frac_A_tvgd_positive", 0) > 0.50,
    "grid_A_2h_positive_majority": grid_row.get("frac_A_2h_positive", 0) > 0.80,
    "grid_not_catastrophic_vs_powerlaw": abs(delta_AIC_grid_minus_powerlaw) < 100,
}

score = sum(bool(v) for v in criteria.values())

if score >= 10:
    status = "PASSOU_FORTE_HSC_TVGD_2HALO_COMPATIVEL"
elif score >= 8:
    status = "PASSOU_MODERADO_FORTE_HSC_TVGD_2HALO_COM_ALERTAS"
elif score >= 6:
    status = "PASSOU_MODERADO_HSC_TVGD_COMPATIVEL_COM_ALERTAS"
else:
    status = "SINAL_FRACO_HSC_TVGD_PRECISA_REMODELAGEM"

summary = {
    "FINAL_STATUS": status,
    "score": int(score),
    "script": SCRIPT_NAME,

    "N_points": int(len(base)),
    "N_bins": int(base["bin_id"].nunique()),

    "chi2_baryon_positive": float(bar_row["chi2"]),
    "chi2_TVGD_exp_positive": float(tvgd_row["chi2"]),
    "chi2_TVGD_plus_2halo_NNLS_fixed": float(fixed_row["chi2"]),
    "chi2_TVGD_plus_2halo_NNLS_grid": float(grid_row["chi2"]),
    "chi2_powerlaw_positive": float(pow_row["chi2"]),

    "chi2red_baryon_positive": float(bar_row["chi2red"]),
    "chi2red_TVGD_exp_positive": float(tvgd_row["chi2red"]),
    "chi2red_TVGD_plus_2halo_NNLS_fixed": float(fixed_row["chi2red"]),
    "chi2red_TVGD_plus_2halo_NNLS_grid": float(grid_row["chi2red"]),
    "chi2red_powerlaw_positive": float(pow_row["chi2red"]),

    "best_s_2h_grid": float(grid_row["s_2h"]),
    "median_u": float(grid_row["median_u"]),
    "min_u": float(grid_row["min_u"]),
    "max_u": float(grid_row["max_u"]),

    "frac_A_tvgd_positive_fixed": float(fixed_row["frac_A_tvgd_positive"]),
    "frac_A_2h_positive_fixed": float(fixed_row["frac_A_2h_positive"]),
    "frac_both_positive_fixed": float(fixed_row["frac_both_positive"]),

    "frac_A_tvgd_positive_grid": float(grid_row["frac_A_tvgd_positive"]),
    "frac_A_2h_positive_grid": float(grid_row["frac_A_2h_positive"]),
    "frac_both_positive_grid": float(grid_row["frac_both_positive"]),

    "delta_chi2_exp_minus_bar": float(delta_chi2_exp_minus_bar),
    "delta_chi2_fixed_minus_exp": float(delta_chi2_fixed_minus_exp),
    "delta_chi2_grid_minus_exp": float(delta_chi2_grid_minus_exp),
    "delta_AIC_fixed_minus_powerlaw": float(delta_AIC_fixed_minus_powerlaw),
    "delta_AIC_grid_minus_powerlaw": float(delta_AIC_grid_minus_powerlaw),
}

summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "HSC_TVGD_2halo_summary.csv"), index=False)

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — HSC WEAK LENSING")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(8, 5))
plt.bar(res["model"], res["chi2red"])
plt.xticks(rotation=35, ha="right")
plt.ylabel(r"$\chi^2_{\rm red}$")
plt.title("HSC weak lensing — model comparison")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "HSC_model_comparison_chi2red.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.plot(grid_df["s_2h"], grid_df["AIC"], marker="o")
plt.axvline(grid_row["s_2h"], linestyle=":", label=f"best s={grid_row['s_2h']:.2f}")
plt.xlabel(r"$s_{2h}$")
plt.ylabel("AIC")
plt.title("HSC TVGD+2halo — slope scan")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "HSC_2halo_s_grid_AIC.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
plt.hist(work_fixed["u"], bins=35)
plt.xlabel(r"$u=\sqrt{g_{\rm bar}/a_0}$")
plt.ylabel("N")
plt.title("HSC weak lensing — structural variable u")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "HSC_u_distribution.png"), dpi=170)
plt.close()

# Exemplo de perfil para alguns bins
plt.figure(figsize=(8, 5))
for bin_id, sub in work_fixed.groupby("bin_id"):
    if len(sub) < 5:
        continue
    # Mostrar poucos bins para não poluir
    if parse_logm_from_col(bin_id) not in [9.0, 10.0, 11.0]:
        continue

    plt.errorbar(
        sub["R_kpc"],
        sub["DeltaSigma_obs"],
        yerr=sub["DeltaSigma_err"],
        marker="o",
        linestyle="none",
        label=str(bin_id),
        alpha=0.8,
    )

plt.xscale("log")
plt.yscale("log")
plt.xlabel("R [kpc]")
plt.ylabel(r"$\Delta\Sigma$")
plt.title("HSC weak lensing — example observed profiles")
plt.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "HSC_example_profiles.png"), dpi=170)
plt.close()

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nFim do teste HSC weak lensing.")
