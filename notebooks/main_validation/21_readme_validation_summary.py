# ============================================================
# TVGD — 21_README_VALIDATION_SUMMARY
# Gera resumo automático dos testes principais para o README
#
# Objetivo:
#   Ler os resultados públicos salvos em results/ e gerar:
#
#   1. Uma tabela geral dos testes principais
#   2. Um arquivo Markdown pronto para colar no README.md
#   3. Um CSV consolidado de status
#   4. Figuras simples de status/score
#
# Saída:
#   results/readme_validation_summary/
#
# Uso:
#   python notebooks/main_validation/21_readme_validation_summary.py
# ============================================================

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "21_readme_validation_summary.py"

RESULTS_ROOT = "results"
OUTDIR = "results/readme_validation_summary"
os.makedirs(OUTDIR, exist_ok=True)

print("=" * 72)
print("TVGD — 21_README_VALIDATION_SUMMARY")
print("Resumo automático dos testes principais para README")
print("=" * 72)


# ============================================================
# Lista esperada dos testes principais
# ============================================================

EXPECTED_TESTS = [
    {
        "order": 11,
        "script": "11_rmax_sparc_universality.py",
        "sector": "Galaxy dynamics / SPARC",
        "short_name": "RMAX SPARC universality",
        "keywords": ["rmax", "sparc", "universality"],
        "claim": "Tests whether the structural response law P(u) remains viable in SPARC rotation curves.",
        "risk": "high",
    },
    {
        "order": 12,
        "script": "12_rmax_transition_beta_robustness.py",
        "sector": "Galaxy transition window",
        "short_name": "RMAX transition beta robustness",
        "keywords": ["rmax", "transition"],
        "claim": "Tests the critical transition window and robustness of beta near unity.",
        "risk": "high",
    },
    {
        "order": 15,
        "script": "15_extrisk_little_things_external.py",
        "sector": "External galaxy validation",
        "short_name": "EXTRISK LITTLE THINGS",
        "keywords": ["extrisk", "little"],
        "claim": "Checks whether an external non-SPARC dwarf-galaxy dataset breaks the TVGD response family.",
        "risk": "high",
    },
    {
        "order": 16,
        "script": "16_cmbrisk_class_growth.py",
        "sector": "CMB and growth",
        "short_name": "CMBRISK CLASS growth",
        "keywords": ["cmbrisk", "class"],
        "claim": "Tests whether TVGD V1.1 preserves CMB spectra and linear growth.",
        "risk": "maximum",
    },
    {
        "order": 17,
        "script": "17_covrisk_scalar_vector_ppn.py",
        "sector": "Covariant stability / PPN",
        "short_name": "COVRISK scalar-vector-PPN",
        "keywords": ["covrisk", "scalar"],
        "claim": "Audits scalar stability, vector c_T=1 branch, and effective Solar-System PPN safety.",
        "risk": "maximum",
    },
    {
        "order": 18,
        "script": "18_pantheon_full_covariance_sn.py",
        "sector": "Supernovae",
        "short_name": "Pantheon+ full covariance",
        "keywords": ["pantheon"],
        "claim": "Tests TVGD V1.1 against Pantheon+SH0ES with full STAT+SYS covariance.",
        "risk": "high",
    },
    {
        "order": 19,
        "script": "19_black_hole_structural_saturation.py",
        "sector": "Black holes",
        "short_name": "BH structural saturation",
        "keywords": ["black_hole", "bh19"],
        "claim": "Tests whether structural saturation preserves horizon, photon sphere, shadow and ISCO.",
        "risk": "high",
    },
    {
        "order": 20,
        "script": "20_falsifiability_prediction_checklist.py",
        "sector": "Falsifiability",
        "short_name": "Falsifiability checklist",
        "keywords": ["falsifiability"],
        "claim": "Lists explicit ways the model could be falsified or pressured by data.",
        "risk": "meta",
    },
]


# ============================================================
# Funções auxiliares
# ============================================================

def find_summary_files():
    patterns = [
        os.path.join(RESULTS_ROOT, "**", "*summary*.csv"),
        os.path.join(RESULTS_ROOT, "**", "*Summary*.csv"),
        os.path.join(RESULTS_ROOT, "**", "*final_status*.csv"),
        os.path.join(RESULTS_ROOT, "**", "*FINAL_STATUS*.csv"),
    ]

    files = []
    for pat in patterns:
        files.extend(glob.glob(pat, recursive=True))

    return sorted(set(files))


def read_summary(path):
    try:
        df = pd.read_csv(path)
        if len(df) == 0:
            return {}
        return df.iloc[0].to_dict()
    except Exception:
        return {}


def normalize_status(status):
    s = str(status).upper()

    if "FALHOU" in s or "FAILED" in s:
        return "failed"
    if "PASSOU_FORTE" in s:
        return "passed strong"
    if "MODERADO_FORTE" in s:
        return "passed moderate/strong"
    if "PASSOU_MODERADO" in s:
        return "passed moderate"
    if "PASSOU_FRACO" in s or "FRACO" in s:
        return "passed weak"
    if "ALERTA" in s or "ALERT" in s:
        return "passed with alerts"
    if s in ["", "NAN", "NONE", "NA"]:
        return "not found"

    return "unknown"


def status_rank(public_status):
    ranks = {
        "passed strong": 5,
        "passed moderate/strong": 4,
        "passed moderate": 3,
        "passed with alerts": 2.5,
        "passed weak": 2,
        "unknown": 1,
        "not found": 0,
        "failed": -2,
    }
    return ranks.get(public_status, 1)


def find_matching_summary(summary_files, test):
    """
    Procura arquivo de summary que combine com as keywords do teste.
    Usa o caminho do arquivo em minúsculas.
    """
    keywords = [k.lower() for k in test["keywords"]]

    matches = []

    for path in summary_files:
        low = path.lower()
        score = sum(k in low for k in keywords)

        if score > 0:
            matches.append((score, path))

    if len(matches) == 0:
        return None

    matches = sorted(matches, key=lambda x: (-x[0], len(x[1])))
    return matches[0][1]


def extract_key_metrics(data):
    """
    Extrai algumas métricas importantes de diferentes summaries.
    """
    keys_priority = [
        "score",
        "N_points",
        "N_SN",
        "N_models",
        "N_predictions",
        "chi2red_TVGD",
        "chi2red_TVGD_exp",
        "chi2red_TVGD_plus_2halo_NNLS_grid",
        "TT_frac_rms_pct",
        "EE_frac_rms_pct",
        "phiphi_frac_rms_pct",
        "growth_max_abs_frac_fsigma8_pct",
        "delta_chi2_TVGD_minus_LCDM",
        "best_b_shadow_over_rs",
        "best_x_photon",
        "best_x_h_outer",
        "best_cT2",
        "best_ppn_alpha1",
        "weighted_fraction",
    ]

    found = []

    for k in keys_priority:
        if k in data:
            v = data[k]
            try:
                if isinstance(v, float):
                    found.append(f"{k}={v:.6g}")
                else:
                    found.append(f"{k}={v}")
            except Exception:
                found.append(f"{k}={v}")

    if len(found) == 0:
        return ""

    return "; ".join(found[:6])


# ============================================================
# Ler summaries
# ============================================================

summary_files = find_summary_files()

print(f"\nArquivos summary encontrados: {len(summary_files)}")
for p in summary_files:
    print(" -", p)

rows = []

for test in EXPECTED_TESTS:
    path = find_matching_summary(summary_files, test)

    if path is None:
        data = {}
        raw_status = "NOT_FOUND"
        public_status = "not found"
        score = np.nan
        metrics = ""
    else:
        data = read_summary(path)
        raw_status = data.get("FINAL_STATUS", data.get("final_status", data.get("status", "UNKNOWN")))
        public_status = normalize_status(raw_status)
        score = data.get("score", np.nan)
        metrics = extract_key_metrics(data)

    rows.append({
        "order": test["order"],
        "script": test["script"],
        "sector": test["sector"],
        "short_name": test["short_name"],
        "risk": test["risk"],
        "claim": test["claim"],
        "summary_file": path if path is not None else "",
        "raw_status": raw_status,
        "public_status": public_status,
        "status_rank": status_rank(public_status),
        "score": score,
        "key_metrics": metrics,
    })

summary_df = pd.DataFrame(rows).sort_values("order").reset_index(drop=True)
summary_df.to_csv(os.path.join(OUTDIR, "README21_validation_summary.csv"), index=False)

print("\nResumo consolidado:")
print(summary_df[["order", "short_name", "public_status", "score"]].to_string(index=False))


# ============================================================
# Métricas agregadas
# ============================================================

n_total = len(summary_df)
n_found = int((summary_df["public_status"] != "not found").sum())
n_failed = int((summary_df["public_status"] == "failed").sum())
n_strong = int((summary_df["public_status"] == "passed strong").sum())
n_moderate_or_better = int((summary_df["status_rank"] >= 3).sum())

maximum_risk = summary_df[summary_df["risk"].isin(["maximum", "high"])]
n_high_total = len(maximum_risk)
n_high_failed = int((maximum_risk["public_status"] == "failed").sum())
n_high_moderate_or_better = int((maximum_risk["status_rank"] >= 3).sum())

validation_fraction = n_moderate_or_better / max(n_total, 1)
high_risk_fraction = n_high_moderate_or_better / max(n_high_total, 1)

criteria = {
    "all_expected_tests_listed": n_total >= 8,
    "most_summaries_found": n_found >= 6,
    "no_failed_main_tests": n_failed == 0,
    "has_strong_results": n_strong >= 2,
    "majority_moderate_or_better": validation_fraction >= 0.6,
    "high_risk_majority_moderate_or_better": high_risk_fraction >= 0.6,
    "cmb_present": any(summary_df["script"].str.contains("cmbrisk", case=False)),
    "covrisk_present": any(summary_df["script"].str.contains("covrisk", case=False)),
    "pantheon_present": any(summary_df["script"].str.contains("pantheon", case=False)),
    "black_hole_present": any(summary_df["script"].str.contains("black_hole", case=False)),
}

score = sum(bool(v) for v in criteria.values())

if score >= 9:
    final_status = "PASSOU_FORTE_README21_VALIDATION_SUMMARY"
elif score >= 7:
    final_status = "PASSOU_MODERADO_FORTE_README21_COM_ALERTAS"
elif score >= 5:
    final_status = "PASSOU_MODERADO_README21"
else:
    final_status = "FALHOU_README21_SUMMARY_INCOMPLETO"

final_summary = {
    "FINAL_STATUS": final_status,
    "score": int(score),
    "script": SCRIPT_NAME,
    "N_expected_tests": int(n_total),
    "N_summaries_found": int(n_found),
    "N_failed": int(n_failed),
    "N_passed_strong": int(n_strong),
    "N_moderate_or_better": int(n_moderate_or_better),
    "validation_fraction": float(validation_fraction),
    "N_high_or_max_risk": int(n_high_total),
    "N_high_or_max_failed": int(n_high_failed),
    "N_high_or_max_moderate_or_better": int(n_high_moderate_or_better),
    "high_risk_fraction": float(high_risk_fraction),
}

final_summary.update(criteria)

final_summary_df = pd.DataFrame([final_summary])
final_summary_df.to_csv(os.path.join(OUTDIR, "README21_summary.csv"), index=False)


# ============================================================
# Gerar Markdown para README.md
# ============================================================

md = []

md.append("## Main validation chain\n\n")
md.append(
    "The TVGD repository is organized around explicit risk tests. "
    "Each script is designed to check a specific point where an effective modified-gravity model could fail: "
    "galaxy dynamics, external validation, CMB, supernovae, covariant stability, Solar-System/PPN constraints, "
    "black holes, and falsifiability.\n\n"
)

md.append(f"**Automatic summary status:** `{final_status}`  \n")
md.append(f"**Detected validation summaries:** `{n_found}/{n_total}`  \n")
md.append(f"**Failed main tests:** `{n_failed}`  \n")
md.append(f"**Moderate-or-better fraction:** `{validation_fraction:.3f}`  \n\n")

md.append("| # | Script | Sector | Risk | Status | Key metrics |\n")
md.append("|---:|---|---|---|---|---|\n")

for _, row in summary_df.iterrows():
    md.append(
        f"| {int(row['order'])} "
        f"| `{row['script']}` "
        f"| {row['sector']} "
        f"| {row['risk']} "
        f"| `{row['public_status']}` "
        f"| {row['key_metrics']} |\n"
    )

md.append("\n")
md.append("### Interpretation\n\n")
md.append(
    "A `passed strong` or `passed moderate/strong` status does not mean that TVGD is proven as a fundamental theory. "
    "It means that, within the assumptions of the corresponding public test, the model was not ruled out by that risk sector. "
    "The strongest use of this validation chain is reproducibility: an independent reader can run the scripts, inspect the assumptions, "
    "and replace datasets or priors.\n\n"
)

md.append(
    "Exploratory sectors, such as large-scale orientation memory or environment-dependent Hubble modulation, "
    "should be treated as predictions or appendices rather than decisive evidence.\n\n"
)

md_path = os.path.join(OUTDIR, "README21_validation_section.md")
with open(md_path, "w", encoding="utf-8") as f:
    f.write("".join(md))


# ============================================================
# Gerar Markdown em português também
# ============================================================

md_pt = []

md_pt.append("## Cadeia principal de validação\n\n")
md_pt.append(
    "O repositório da TVGD é organizado em torno de testes explícitos de risco. "
    "Cada script verifica um ponto onde uma teoria efetiva de gravidade modificada poderia falhar: "
    "dinâmica de galáxias, validação externa, CMB, supernovas, estabilidade covariante, limites PPN/Sistema Solar, "
    "buracos negros e falsificabilidade.\n\n"
)

md_pt.append(f"**Status automático do resumo:** `{final_status}`  \n")
md_pt.append(f"**Summaries encontrados:** `{n_found}/{n_total}`  \n")
md_pt.append(f"**Testes principais falhos:** `{n_failed}`  \n")
md_pt.append(f"**Fração moderada-ou-melhor:** `{validation_fraction:.3f}`  \n\n")

md_pt.append("| # | Script | Setor | Risco | Status | Métricas-chave |\n")
md_pt.append("|---:|---|---|---|---|---|\n")

for _, row in summary_df.iterrows():
    md_pt.append(
        f"| {int(row['order'])} "
        f"| `{row['script']}` "
        f"| {row['sector']} "
        f"| {row['risk']} "
        f"| `{row['public_status']}` "
        f"| {row['key_metrics']} |\n"
    )

md_pt.append("\n")
md_pt.append("### Interpretação\n\n")
md_pt.append(
    "Um status `passed strong` ou `passed moderate/strong` não significa que a TVGD esteja provada como teoria fundamental. "
    "Significa que, dentro das hipóteses daquele teste público, o modelo não foi derrubado por aquele setor de risco. "
    "O ponto mais forte desta cadeia é a reprodutibilidade: um avaliador independente pode rodar os scripts, inspecionar as hipóteses "
    "e trocar dados ou priors.\n\n"
)

md_pt.append(
    "Setores exploratórios, como memória de orientação em larga escala ou modulação ambiental de H0, "
    "devem ser tratados como previsões/apêndices, não como evidência decisiva.\n\n"
)

md_pt_path = os.path.join(OUTDIR, "README21_validation_section_PT.md")
with open(md_pt_path, "w", encoding="utf-8") as f:
    f.write("".join(md_pt))


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(9, 5))
plt.bar(summary_df["short_name"], summary_df["status_rank"])
plt.xticks(rotation=40, ha="right")
plt.ylabel("status rank")
plt.title("TVGD validation chain — status rank")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "README21_status_rank.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
counts = summary_df["public_status"].value_counts()
plt.bar(counts.index.astype(str), counts.values)
plt.xticks(rotation=30, ha="right")
plt.ylabel("N")
plt.title("TVGD validation chain — status counts")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "README21_status_counts.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
risk_counts = summary_df["risk"].value_counts()
plt.bar(risk_counts.index.astype(str), risk_counts.values)
plt.ylabel("N")
plt.title("TVGD validation chain — risk categories")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "README21_risk_categories.png"), dpi=170)
plt.close()


# ============================================================
# Impressão final
# ============================================================

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — README 21")
print("=" * 72)

for k, v in final_summary.items():
    print(f"{k}: {v}")

print("\nArquivos salvos em:")
print(OUTDIR)

print("\nSeção README em inglês:")
print(md_path)

print("\nSeção README em português:")
print(md_pt_path)

print("\nFim do teste README 21.")
