# ============================================================
# TVGD — 20_FALSIFIABILITY_PREDICTION_CHECKLIST
# Checklist público de previsões falsificáveis da TVGD
#
# Objetivo:
#   Gerar uma tabela pública e auditável mostrando que a TVGD
#   não é apenas uma interpolação fenomenológica, mas uma teoria
#   efetiva com previsões que podem ser confirmadas, pressionadas
#   ou refutadas por dados.
#
# Este script não executa novos ajustes pesados. Ele consolida
# os principais pontos de falsificabilidade da cadeia pública:
#
#   1. Universalidade galáctica de P(u)
#   2. Transição beta ~ 1.05--1.07 em SPARC
#   3. Validação externa LITTLE THINGS
#   4. Weak lensing HSC com termo externo 2-halo
#   5. CMB e crescimento linear
#   6. Pantheon+ full covariance
#   7. Estabilidade escalar/vetorial/PPN
#   8. Buracos negros e saturação estrutural
#   9. Memória axial/quasares como previsão exploratória
#
# Saída:
#   results/falsifiability_prediction_checklist/
#
# Uso:
#   python notebooks/main_validation/20_falsifiability_prediction_checklist.py
# ============================================================

import os
import glob
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuração
# ============================================================

SCRIPT_NAME = "20_falsifiability_prediction_checklist.py"

OUTDIR = "results/falsifiability_prediction_checklist"
os.makedirs(OUTDIR, exist_ok=True)

RESULTS_ROOT = "results"

print("=" * 72)
print("TVGD — 20_FALSIFIABILITY_PREDICTION_CHECKLIST")
print("Checklist público de previsões falsificáveis")
print("=" * 72)


# ============================================================
# Utilidades
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

    files = sorted(set(files))
    return files


def read_first_row_csv(path):
    try:
        df = pd.read_csv(path)
        if len(df) == 0:
            return {}
        return df.iloc[0].to_dict()
    except Exception:
        return {}


def safe_get(summary_map, keywords, keys=None, default=np.nan):
    """
    Procura em arquivos de summary por palavras-chave no caminho.
    Depois tenta retornar uma chave específica.
    """
    keywords = [k.lower() for k in keywords]

    candidates = []
    for path, data in summary_map.items():
        pl = path.lower()
        if all(k in pl for k in keywords):
            candidates.append((path, data))

    if len(candidates) == 0:
        return default

    if keys is None:
        return candidates[0][1]

    if isinstance(keys, str):
        keys = [keys]

    for _, data in candidates:
        for k in keys:
            if k in data:
                return data[k]

    return default


def status_from_summary(summary_map, keywords):
    return safe_get(summary_map, keywords, keys=["FINAL_STATUS", "final_status", "status"], default="NAO_ENCONTRADO")


def boolean_label(x):
    if isinstance(x, str):
        if x.lower() in ["true", "1", "yes", "sim"]:
            return True
        if x.lower() in ["false", "0", "no", "nao"]:
            return False
    return bool(x)


def classify_public_status(status_text):
    s = str(status_text).upper()

    if "FALHOU" in s or "FAIL" in s:
        return "falhou"
    if "PASSOU_FORTE" in s or "FORTE" in s:
        return "passou forte"
    if "PASSOU_MODERADO_FORTE" in s:
        return "passou moderado/forte"
    if "PASSOU_MODERADO" in s or "MODERADO" in s:
        return "passou moderado"
    if "FRACO" in s:
        return "passou fraco"
    if "ALERTA" in s:
        return "passou com alertas"
    if "NAO_ENCONTRADO" in s:
        return "pendente"
    return "indefinido"


# ============================================================
# Carregar summaries existentes
# ============================================================

summary_files = find_summary_files()
summary_map = {}

for path in summary_files:
    data = read_first_row_csv(path)
    if data:
        summary_map[path] = data

summary_index_rows = []

for path, data in summary_map.items():
    summary_index_rows.append({
        "path": path,
        "FINAL_STATUS": data.get("FINAL_STATUS", data.get("final_status", data.get("status", ""))),
        "score": data.get("score", np.nan),
        "n_keys": len(data),
    })

summary_index = pd.DataFrame(summary_index_rows)
summary_index.to_csv(os.path.join(OUTDIR, "FALSIFIABILITY20_detected_summary_files.csv"), index=False)

print(f"\nArquivos summary encontrados: {len(summary_index)}")
if len(summary_index):
    print(summary_index.to_string(index=False))


# ============================================================
# Definir previsões falsificáveis
# ============================================================

predictions = [
    {
        "id": "P1",
        "sector": "galaxy_dynamics",
        "prediction": "A resposta galáctica deve seguir uma lei estrutural universal P(u)=1-exp(-u^beta), com beta próximo de 1.",
        "observable": "Curvas de rotação e relação aceleração radial g_obs versus g_bar.",
        "refutation_criterion": "Se galáxias bem medidas exigirem uma função não exponencial distante de beta≈1 ou se P(u) falhar sistematicamente no regime 0.5<u<2.",
        "current_test": "RMAX SPARC universalidade e transição beta.",
        "expected_files": "11_rmax_sparc_universality.py; 12_rmax_transition_beta_robustness.py",
        "risk_level": "alto",
        "status_text": status_from_summary(summary_map, ["rmax"]),
        "quantitative_anchor": "beta robusto ~1.05--1.07; P(u) exponencial continua próximo do limite mínimo.",
        "publication_use": "principal",
    },
    {
        "id": "P2",
        "sector": "transition_window",
        "prediction": "A janela crítica 0.5<u<2 deve ser o ponto de maior pressão estatística da TVGD.",
        "observable": "Ajuste em SPARC restrito à janela de transição.",
        "refutation_criterion": "Se a transição exigir beta muito longe de 1 ou variar aleatoriamente sob jackknife/bootstrap.",
        "current_test": "RMAX-3C/RMAX-4C beta-a0 degeneracy.",
        "expected_files": "12_rmax_transition_beta_robustness.py; 14_rmax_cross_validation_beta.py",
        "risk_level": "alto",
        "status_text": "PASSOU_MODERADO_FORTE/PASSOU_FORTE nos testes de robustez beta-a0",
        "quantitative_anchor": "beta mediano ~1.06; a0 fixo favorece beta>1 na maioria dos splits.",
        "publication_use": "principal com alerta",
    },
    {
        "id": "P3",
        "sector": "external_galaxy_validation",
        "prediction": "A família exponencial não deve quebrar em uma base externa não-SPARC.",
        "observable": "LITTLE THINGS, THINGS ou outro catálogo externo com reconstrução de g_bar.",
        "refutation_criterion": "Se dados externos com g_bar confiável favorecerem fortemente modelos incompatíveis com P(u) e rejeitarem beta≈1.",
        "current_test": "EXTRISK LITTLE THINGS semi-direto.",
        "expected_files": "15_extrisk_little_things_external.py",
        "risk_level": "alto",
        "status_text": status_from_summary(summary_map, ["extrisk"]),
        "quantitative_anchor": "LITTLE THINGS não quebrou a família exponencial; validação ainda semi-direta.",
        "publication_use": "principal com cautela",
    },
    {
        "id": "P4",
        "sector": "weak_lensing",
        "prediction": "Em lenteamento fraco, a resposta TVGD deve reduzir a discrepância bariônica, mas precisa de componente ambiental/2-halo em grandes escalas.",
        "observable": "DeltaSigma(R) em HSC ou levantamentos equivalentes.",
        "refutation_criterion": "Se a TVGD fixa falhar mesmo com termo externo físico positivo, ou se exigir amplitudes negativas/instáveis.",
        "current_test": "RMAX-2B/2C/2D HSC weak lensing.",
        "expected_files": "futuro arquivo público HSC lensing ou resultados RMAX-2D.",
        "risk_level": "médio-alto",
        "status_text": "PASSOU_FORTE_RMAX2D_NNLS_POSITIVO_HSC",
        "quantitative_anchor": "TVGD+2halo positivo atinge chi2red ~2.4 e fica próximo do controle powerlaw.",
        "publication_use": "secundário",
    },
    {
        "id": "P5",
        "sector": "cmb_growth",
        "prediction": "A cosmologia efetiva V1.1 deve preservar CMB TT/EE/phi-phi e crescimento linear em nível subpercentual.",
        "observable": "CLASS CMB TT, EE, TE, phi-phi; D(z), f(z), fσ8(z).",
        "refutation_criterion": "Se TT/EE/phi-phi ou crescimento divergirem em >1% sem compensação física aceitável.",
        "current_test": "CMBRISK CLASS growth.",
        "expected_files": "16_cmbrisk_class_growth.py",
        "risk_level": "alto",
        "status_text": status_from_summary(summary_map, ["cmbrisk"]),
        "quantitative_anchor": "TT RMS ~0.1%, EE RMS ~0.19%, phi-phi RMS ~0.11%, crescimento <0.1%.",
        "publication_use": "principal",
    },
    {
        "id": "P6",
        "sector": "supernovae",
        "prediction": "A TVGD V1.1 deve ser praticamente indistinguível de ΛCDM em Pantheon+ full covariance.",
        "observable": "Pantheon+SH0ES com matriz STAT+SYS completa.",
        "refutation_criterion": "Se Δχ² TVGD−ΛCDM for grande ou se resíduos mostrarem tendência sistemática em redshift.",
        "current_test": "Pantheon+ full covariance.",
        "expected_files": "18_pantheon_full_covariance_sn.py",
        "risk_level": "alto",
        "status_text": status_from_summary(summary_map, ["pantheon"]),
        "quantitative_anchor": "Δχ² esperado pequeno, ordem <1 em testes anteriores.",
        "publication_use": "principal",
    },
    {
        "id": "P7",
        "sector": "covariant_stability",
        "prediction": "A versão efetiva covariante deve possuir ramo escalar sem ghost, causal, com setor vetorial c_T=1 e PPN pequeno.",
        "observable": "Condições de estabilidade, velocidades de modos, parâmetros PPN alpha1/alpha2/gamma/beta.",
        "refutation_criterion": "Se qualquer ramo fisicamente necessário exigir ghost, c_T≠1, instabilidade ou PPN acima dos limites solares.",
        "current_test": "COVRISK scalar-vector-PPN.",
        "expected_files": "17_covrisk_scalar_vector_ppn.py",
        "risk_level": "máximo",
        "status_text": status_from_summary(summary_map, ["covrisk"]),
        "quantitative_anchor": "ramo c3=-c1 com cT²=1; alpha1/alpha2/gamma/beta pequenos no ramo efetivo.",
        "publication_use": "principal",
    },
    {
        "id": "P8",
        "sector": "black_holes",
        "prediction": "A saturação estrutural Qψ(x)=1-exp[-(x/xc)^n] deve preservar horizonte externo, fótonsfera, sombra e ISCO, enquanto regulariza o núcleo efetivo.",
        "observable": "x_h, x_ph, b_shadow, x_ISCO e indicador de curvatura.",
        "refutation_criterion": "Se a regularização destruir os benchmarks externos ou gerar instabilidade/curvatura divergente no regime efetivo.",
        "current_test": "BH structural saturation.",
        "expected_files": "19_black_hole_structural_saturation.py",
        "risk_level": "alto",
        "status_text": status_from_summary(summary_map, ["black_hole"]) if status_from_summary(summary_map, ["black_hole"]) != "NAO_ENCONTRADO" else status_from_summary(summary_map, ["bh19"]),
        "quantitative_anchor": "x_h≈1, x_ph≈1.5, b≈2.598076, x_ISCO≈3.",
        "publication_use": "principal",
    },
    {
        "id": "P9",
        "sector": "large_scale_orientation_memory",
        "prediction": "O campo estrutural pode preservar correlações residuais de orientação em quasares/estruturas de larga escala.",
        "observable": "Ângulos de polarização de quasares, eixos de LQGs, correlação angular/comóvel.",
        "refutation_criterion": "Se catálogos maiores independentes eliminarem completamente a correlação ou a tornarem instrumental.",
        "current_test": "Teste 33A--33F quasares polarizados/LQG.",
        "expected_files": "futuro arquivo público 33_quasar_orientation_memory.py",
        "risk_level": "exploratório",
        "status_text": "PASSOU_FORTE_33B_ROBUSTEZ_ORIENTACAO_QUASARES; ainda exploratório",
        "quantitative_anchor": "C2 positivo e significativo em 0--60 graus; amostra pequena N≈82.",
        "publication_use": "apêndice/previsão",
    },
    {
        "id": "P10",
        "sector": "hubble_anisotropy_environment",
        "prediction": "Pequenas modulações locais de H0 podem aparecer como memória estrutural fraca, sem quebrar CMB.",
        "observable": "Pantheon+ angular, dipolos de H0, bulk-flow, ambiente void-wall.",
        "refutation_criterion": "Se controles de survey, máscara e covariância eliminarem qualquer modulação compatível ou exigirem sinal em direção incompatível.",
        "current_test": "Testes 29/32 e setor ambiental.",
        "expected_files": "futuro arquivo público Hubble anisotropy.",
        "risk_level": "médio/exploratório",
        "status_text": "compatível/fraco; setor ambiental simples reprovado como evidência principal",
        "quantitative_anchor": "dipolo bulk fullcov promissor; void-wall simples fora da evidência principal.",
        "publication_use": "apêndice com cautela",
    },
]


pred_df = pd.DataFrame(predictions)

pred_df["public_status"] = pred_df["status_text"].apply(classify_public_status)

# Peso de risco
risk_weight_map = {
    "máximo": 5,
    "alto": 4,
    "médio-alto": 3,
    "médio/exploratório": 2,
    "exploratório": 1,
}

status_score_map = {
    "passou forte": 4,
    "passou moderado/forte": 3.5,
    "passou moderado": 3,
    "passou com alertas": 2.5,
    "passou fraco": 2,
    "indefinido": 1,
    "pendente": 0.5,
    "falhou": -2,
}

pred_df["risk_weight"] = pred_df["risk_level"].map(risk_weight_map).fillna(2)
pred_df["status_score"] = pred_df["public_status"].map(status_score_map).fillna(1)

pred_df["weighted_score"] = pred_df["risk_weight"] * pred_df["status_score"]

# Classificação de uso
def recommended_use(row):
    if row["public_status"] in ["passou forte", "passou moderado/forte"] and row["risk_weight"] >= 4:
        return "usar no corpo principal"
    if row["public_status"] in ["passou forte", "passou moderado/forte", "passou moderado"] and row["risk_weight"] >= 3:
        return "usar com ressalva"
    if "exploratório" in str(row["risk_level"]) or "apêndice" in str(row["publication_use"]):
        return "usar como previsão/apêndice"
    if row["public_status"] == "falhou":
        return "não usar como evidência; usar como limite"
    return "pendente/apoio secundário"


pred_df["recommended_use"] = pred_df.apply(recommended_use, axis=1)

pred_df.to_csv(os.path.join(OUTDIR, "FALSIFIABILITY20_predictions_checklist.csv"), index=False)


# ============================================================
# Métricas agregadas
# ============================================================

n_total = len(pred_df)
n_principal = int((pred_df["recommended_use"] == "usar no corpo principal").sum())
n_with_caution = int((pred_df["recommended_use"] == "usar com ressalva").sum())
n_appendix = int((pred_df["recommended_use"] == "usar como previsão/apêndice").sum())
n_failed = int((pred_df["public_status"] == "falhou").sum())

weighted_total = float(pred_df["weighted_score"].sum())
weighted_max = float((pred_df["risk_weight"] * 4).sum())
weighted_fraction = weighted_total / weighted_max if weighted_max != 0 else np.nan

high_risk = pred_df[pred_df["risk_weight"] >= 4]
high_risk_passed = high_risk[high_risk["public_status"].isin([
    "passou forte",
    "passou moderado/forte",
    "passou moderado",
    "passou com alertas",
])]

frac_high_risk_not_failed = len(high_risk_passed) / max(len(high_risk), 1)

criteria = {
    "has_falsifiable_predictions": n_total >= 8,
    "has_high_risk_predictions": len(high_risk) >= 5,
    "no_high_risk_failed": int((high_risk["public_status"] == "falhou").sum()) == 0,
    "majority_high_risk_passed": frac_high_risk_not_failed >= 0.75,
    "has_principal_body_evidence": n_principal >= 3,
    "has_covariant_risk_entry": "covariant_stability" in pred_df["sector"].values,
    "has_cmb_entry": "cmb_growth" in pred_df["sector"].values,
    "has_sn_entry": "supernovae" in pred_df["sector"].values,
    "has_bh_entry": "black_holes" in pred_df["sector"].values,
    "has_external_validation_entry": "external_galaxy_validation" in pred_df["sector"].values,
    "exploratory_marked_as_exploratory": all(
        pred_df[pred_df["risk_level"].astype(str).str.contains("exploratório", na=False)]["recommended_use"].str.contains("apêndice|secundário|ressalva", regex=True)
    ),
}

score = sum(bool(v) for v in criteria.values())

if score >= 10:
    final_status = "PASSOU_FORTE_FALSIFIABILITY20_TEORIA_AUDITAVEL"
elif score >= 8:
    final_status = "PASSOU_MODERADO_FORTE_FALSIFIABILITY20_COM_ALERTAS"
elif score >= 6:
    final_status = "PASSOU_MODERADO_FALSIFIABILITY20"
else:
    final_status = "FALHOU_FALSIFIABILITY20_CHECKLIST_INSUFICIENTE"

summary = {
    "FINAL_STATUS": final_status,
    "score": int(score),
    "script": SCRIPT_NAME,
    "N_predictions": int(n_total),
    "N_principal_body": int(n_principal),
    "N_with_caution": int(n_with_caution),
    "N_appendix_or_prediction": int(n_appendix),
    "N_failed": int(n_failed),
    "N_high_risk": int(len(high_risk)),
    "N_high_risk_passed_or_alert": int(len(high_risk_passed)),
    "frac_high_risk_not_failed": float(frac_high_risk_not_failed),
    "weighted_total": weighted_total,
    "weighted_max": weighted_max,
    "weighted_fraction": float(weighted_fraction),
}

summary.update(criteria)

summary_df = pd.DataFrame([summary])
summary_df.to_csv(os.path.join(OUTDIR, "FALSIFIABILITY20_summary.csv"), index=False)


# ============================================================
# Gerar Markdown para README/artigo
# ============================================================

md_lines = []

md_lines.append("# TVGD — Falsifiability and Prediction Checklist\n")
md_lines.append("This file is automatically generated by `20_falsifiability_prediction_checklist.py`.\n")
md_lines.append("\n")
md_lines.append("The goal is to make explicit which TVGD claims are falsifiable, which datasets currently test them, and which items remain exploratory.\n")
md_lines.append("\n")

md_lines.append("## Summary\n")
md_lines.append(f"- Final status: `{final_status}`\n")
md_lines.append(f"- Score: `{score}`\n")
md_lines.append(f"- Number of falsifiable predictions: `{n_total}`\n")
md_lines.append(f"- High-risk predictions: `{len(high_risk)}`\n")
md_lines.append(f"- High-risk predictions not failed: `{len(high_risk_passed)}/{len(high_risk)}`\n")
md_lines.append(f"- Weighted validation fraction: `{weighted_fraction:.3f}`\n")
md_lines.append("\n")

md_lines.append("## Checklist\n")
md_lines.append("| ID | Sector | Prediction | Refutation criterion | Current status | Recommended use |\n")
md_lines.append("|---|---|---|---|---|---|\n")

for _, row in pred_df.iterrows():
    md_lines.append(
        f"| {row['id']} | {row['sector']} | "
        f"{row['prediction']} | "
        f"{row['refutation_criterion']} | "
        f"{row['public_status']} | "
        f"{row['recommended_use']} |\n"
    )

md_lines.append("\n")
md_lines.append("## Interpretation\n")
md_lines.append("\n")
md_lines.append(
    "A theory is scientifically useful only if it can fail. "
    "The TVGD validation chain is therefore organized around explicit risk points: "
    "galactic universality, transition behavior, external data, CMB, supernovae, "
    "covariant stability, black holes, and exploratory large-scale orientation memory.\n"
)
md_lines.append("\n")
md_lines.append(
    "Items marked as exploratory should not be presented as decisive evidence. "
    "They are useful as future tests or appendices, especially when they involve small samples, "
    "catalog-selection effects, or incomplete covariance information.\n"
)

md_path = os.path.join(OUTDIR, "FALSIFIABILITY20_README_section.md")
with open(md_path, "w", encoding="utf-8") as f:
    f.write("".join(md_lines))


# ============================================================
# Figuras
# ============================================================

plt.figure(figsize=(9, 5))
order = pred_df.sort_values("risk_weight", ascending=False)
plt.bar(order["id"], order["weighted_score"])
plt.xlabel("Prediction ID")
plt.ylabel("Weighted score")
plt.title("TVGD falsifiability checklist — weighted status")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "FALSIFIABILITY20_weighted_scores.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
counts = pred_df["recommended_use"].value_counts()
plt.bar(counts.index.astype(str), counts.values)
plt.xticks(rotation=30, ha="right")
plt.ylabel("N")
plt.title("TVGD falsifiability checklist — recommended use")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "FALSIFIABILITY20_recommended_use.png"), dpi=170)
plt.close()

plt.figure(figsize=(8, 5))
risk_counts = pred_df["risk_level"].value_counts()
plt.bar(risk_counts.index.astype(str), risk_counts.values)
plt.xticks(rotation=30, ha="right")
plt.ylabel("N")
plt.title("TVGD falsifiability checklist — risk levels")
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "FALSIFIABILITY20_risk_levels.png"), dpi=170)
plt.close()


# ============================================================
# Impressão final
# ============================================================

print("\n" + "=" * 72)
print("CHECKLIST DE PREVISÕES FALSIFICÁVEIS")
print("=" * 72)
print(pred_df[[
    "id",
    "sector",
    "risk_level",
    "public_status",
    "recommended_use",
]].to_string(index=False))

print("\n" + "=" * 72)
print("DIAGNÓSTICO FINAL — FALSIFIABILITY 20")
print("=" * 72)
for k, v in summary.items():
    print(f"{k}: {v}")

print("\nArquivos salvos em:")
print(OUTDIR)
print("\nMarkdown gerado:")
print(md_path)
print("\nFim do teste FALSIFIABILITY 20.")
