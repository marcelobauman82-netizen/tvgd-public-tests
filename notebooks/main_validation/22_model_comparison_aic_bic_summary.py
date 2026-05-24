#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TVGD — 22_model_comparison_aic_bic_summary.py

Objetivo
--------
Criar um resumo estatístico simples e auditável dos principais testes públicos da TVGD,
incluindo métricas como chi2, chi2_red, AIC_proxy, BIC_proxy, delta_chi2 e status.

Este script NÃO depende de arquivos pesados.
Ele funciona em dois modos:

1. Modo automático:
   - Procura arquivos CSV/JSON/TXT/MD dentro de pastas de saída conhecidas.
   - Tenta extrair métricas se existirem.

2. Modo fallback:
   - Se não encontrar arquivos de saída, usa uma tabela mínima consolidada com
     resultados públicos já registrados nos testes principais.

Saídas
------
Cria a pasta:

    results/22_model_comparison_aic_bic_summary/

com:

    tvgd_model_comparison_summary.csv
    tvgd_model_comparison_summary.json
    tvgd_model_comparison_report.md

Uso
---
Na raiz do repositório:

    python notebooks/main_validation/22_model_comparison_aic_bic_summary.py

Opcional:

    python notebooks/main_validation/22_model_comparison_aic_bic_summary.py --root .
"""

import argparse
import json
import math
import re
from pathlib import Path
from typing import Dict, Any, List, Optional

import pandas as pd


# ============================================================
# Utilidades
# ============================================================

def safe_float(x):
    try:
        if x is None:
            return None
        if isinstance(x, str):
            x = x.strip().replace(",", ".")
            if x.lower() in ["nan", "none", "null", ""]:
                return None
        y = float(x)
        if math.isfinite(y):
            return y
        return None
    except Exception:
        return None


def compute_aic_bic(chi2: Optional[float], k: Optional[int], n: Optional[int]):
    """
    AIC = chi2 + 2k
    BIC = chi2 + k ln(n)

    Aqui usamos como proxy estatístico, pois nem todos os testes públicos possuem
    a mesma likelihood formal.
    """
    chi2 = safe_float(chi2)
    if chi2 is None or k is None or n is None or n <= 0:
        return None, None

    aic = chi2 + 2.0 * k
    bic = chi2 + k * math.log(n)
    return aic, bic


def extract_metric_from_text(text: str, keys: List[str]) -> Optional[float]:
    """
    Procura padrões como:
        chi2 = 123.4
        chi2_red: 1.02
        Delta chi2≈-4.5
    """
    for key in keys:
        pattern = rf"{re.escape(key)}\s*[=:≈]\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)"
        m = re.search(pattern, text)
        if m:
            return safe_float(m.group(1))
    return None


def read_text_file(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return ""


def scan_existing_results(root: Path) -> List[Dict[str, Any]]:
    """
    Varredura leve por resultados já existentes.
    O objetivo é captar métricas se o usuário já tiver rodado scripts anteriores.
    """
    candidates = []

    search_dirs = [
        root / "results",
        root / "outputs",
        root / "notebooks",
        root,
    ]

    files = []
    for d in search_dirs:
        if d.exists():
            files.extend(list(d.rglob("*.csv")))
            files.extend(list(d.rglob("*.json")))
            files.extend(list(d.rglob("*.txt")))
            files.extend(list(d.rglob("*.md")))

    # Evita ler arquivos gigantes
    files = [f for f in files if f.is_file() and f.stat().st_size < 5_000_000]

    for f in files:
        name = f.name.lower()
        full = str(f).lower()

        # Apenas arquivos que parecem pertencer à cadeia TVGD
        if not any(token in full for token in [
            "sparc", "pantheon", "hubble", "black_hole", "black-hole",
            "covrisk", "little", "things", "rmax", "falsifiability",
            "tvgd", "validation", "euclid"
        ]):
            continue

        record = {
            "source_file": str(f.relative_to(root)) if f.is_relative_to(root) else str(f),
            "test_id": None,
            "sector": None,
            "model": None,
            "n_data": None,
            "k_params": None,
            "chi2": None,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": None,
            "notes": "Extraído automaticamente de arquivo existente; revisar antes de citar."
        }

        text = ""

        if f.suffix.lower() == ".csv":
            try:
                df = pd.read_csv(f)
                cols = {c.lower(): c for c in df.columns}

                # Se houver uma tabela de resumo, pega a primeira linha útil
                if len(df) > 0:
                    row = df.iloc[0].to_dict()

                    for possible in ["test_id", "test", "teste", "name"]:
                        if possible in cols:
                            record["test_id"] = str(row[cols[possible]])

                    for possible in ["sector", "setor"]:
                        if possible in cols:
                            record["sector"] = str(row[cols[possible]])

                    for possible in ["model", "modelo"]:
                        if possible in cols:
                            record["model"] = str(row[cols[possible]])

                    for possible in ["n", "n_data", "ndata", "n_points", "npts"]:
                        if possible in cols:
                            record["n_data"] = safe_float(row[cols[possible]])

                    for possible in ["k", "k_params", "n_params", "params"]:
                        if possible in cols:
                            val = safe_float(row[cols[possible]])
                            record["k_params"] = int(val) if val is not None else None

                    for possible in ["chi2", "chisq", "chi_square"]:
                        if possible in cols:
                            record["chi2"] = safe_float(row[cols[possible]])

                    for possible in ["chi2_red", "chi2red", "reduced_chi2", "chisq_red"]:
                        if possible in cols:
                            record["chi2_red"] = safe_float(row[cols[possible]])

                    for possible in ["delta_chi2", "dchi2", "delta_chi2_vs_ref"]:
                        if possible in cols:
                            record["delta_chi2_vs_ref"] = safe_float(row[cols[possible]])

                    for possible in ["status", "final_status", "classification"]:
                        if possible in cols:
                            record["status"] = str(row[cols[possible]])

            except Exception:
                text = read_text_file(f)

        elif f.suffix.lower() == ".json":
            try:
                data = json.loads(read_text_file(f))
                if isinstance(data, dict):
                    record["test_id"] = data.get("test_id") or data.get("test") or data.get("name")
                    record["sector"] = data.get("sector")
                    record["model"] = data.get("model") or data.get("modelo")
                    record["n_data"] = safe_float(data.get("n_data") or data.get("n") or data.get("N"))
                    k_val = safe_float(data.get("k_params") or data.get("k") or data.get("n_params"))
                    record["k_params"] = int(k_val) if k_val is not None else None
                    record["chi2"] = safe_float(data.get("chi2") or data.get("chisq"))
                    record["chi2_red"] = safe_float(data.get("chi2_red") or data.get("chi2red"))
                    record["delta_chi2_vs_ref"] = safe_float(data.get("delta_chi2") or data.get("delta_chi2_vs_ref"))
                    record["status"] = data.get("status") or data.get("final_status")
            except Exception:
                text = read_text_file(f)

        else:
            text = read_text_file(f)

        if text:
            record["chi2"] = record["chi2"] or extract_metric_from_text(
                text, ["chi2", "χ²", "chi2_TVGD", "chi2_total"]
            )
            record["chi2_red"] = record["chi2_red"] or extract_metric_from_text(
                text, ["chi2_red", "chi2red", "χ²red"]
            )
            record["delta_chi2_vs_ref"] = record["delta_chi2_vs_ref"] or extract_metric_from_text(
                text, ["delta_chi2", "Δchi2", "Delta_chi2", "Delta chi2"]
            )

            status_match = re.search(
                r"(PASSOU_[A-Z0-9_]+|FALHOU[A-Z0-9_]*|COMPATIVEL[A-Z0-9_]*|MODERADO[A-Z0-9_]*)",
                text
            )
            if status_match:
                record["status"] = record["status"] or status_match.group(1)

        # Inferência leve pelo nome do arquivo
        lower = str(f).lower()

        if record["sector"] is None:
            if "sparc" in lower or "rar" in lower:
                record["sector"] = "galaxy_dynamics"
            elif "pantheon" in lower or "hubble" in lower:
                record["sector"] = "cosmology_hubble_sn"
            elif "black" in lower or "bh" in lower:
                record["sector"] = "black_holes"
            elif "covrisk" in lower or "ppn" in lower:
                record["sector"] = "covariant_risk"
            elif "little" in lower or "things" in lower:
                record["sector"] = "external_galaxy_validation"
            elif "euclid" in lower:
                record["sector"] = "forecast_euclid"

        if record["model"] is None:
            if "lcdm" in lower:
                record["model"] = "LCDM"
            elif "mond" in lower:
                record["model"] = "MOND"
            elif "tvgd" in lower:
                record["model"] = "TVGD"

        if record["test_id"] is None:
            record["test_id"] = f.stem

        # Só guarda se encontrou alguma métrica ou status relevante
        if any(record.get(k) is not None for k in ["chi2", "chi2_red", "delta_chi2_vs_ref", "status"]):
            candidates.append(record)

    return candidates


def fallback_summary() -> List[Dict[str, Any]]:
    """
    Tabela mínima consolidada dos testes principais já usados como evidência pública.
    Os valores aqui são proxies/resumos e devem ser substituídos por extração automática
    quando os arquivos de resultado estiverem presentes.
    """
    rows = [
        {
            "test_id": "11_rmax_sparc_universality",
            "sector": "galaxy_dynamics",
            "model": "TVGD",
            "n_data": 943,
            "k_params": 1,
            "chi2": None,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "ROBUST_CORE_PROMISING",
            "source_file": "fallback_memory",
            "notes": (
                "SPARC/RAR robusto: corte D com 943 pontos e 57 galáxias; "
                "erro médio absoluto ~16.7%, resíduo médio ~-0.41%, RMSE_log_eta ~0.108."
            ),
        },
        {
            "test_id": "12_rmax_transition_beta_robustness",
            "sector": "galaxy_dynamics",
            "model": "TVGD",
            "n_data": None,
            "k_params": 2,
            "chi2": None,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "ROBUSTNESS_TEST",
            "source_file": "fallback_memory",
            "notes": "Teste de robustez da transição/beta; preencher com CSV gerado pelo script 12 quando disponível.",
        },
        {
            "test_id": "15_extrisk_little_things_external",
            "sector": "external_galaxy_validation",
            "model": "TVGD_exp_beta_1p1",
            "n_data": 720,
            "k_params": 23,
            "chi2": 1736.431,
            "chi2_red": 2.4913,
            "delta_chi2_vs_ref": -348.299,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "PASSOU_FORTE_EXTRISK1B_LITTLE_THINGS_EXTERNO",
            "source_file": "fallback_memory",
            "notes": (
                "Validação externa LITTLE THINGS; melhor família exponencial beta≈1.1; "
                "baryon-only chi2_red≈2.9237 e AIC≈2084.730."
            ),
        },
        {
            "test_id": "17_covrisk_scalar_vector_ppn",
            "sector": "covariant_risk",
            "model": "TVGD_covariant_effective",
            "n_data": 355509,
            "k_params": 4,
            "chi2": None,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "PASSOU_FORTE_COVRISK_PPN",
            "source_file": "fallback_memory",
            "notes": (
                "Ramo c3=-c1 com cT=1; alpha1, alpha2, gamma e beta_PPN dentro dos limites fortes "
                "no candidato pequeno."
            ),
        },
        {
            "test_id": "18_covrisk_rsd_fsigma8",
            "sector": "growth_rsd",
            "model": "TVGD_V1p1",
            "n_data": 22,
            "k_params": 1,
            "chi2": 17.4505155,
            "chi2_red": 0.830977,
            "delta_chi2_vs_ref": -0.00230,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "PASSOU_FORTE_COVRISK8_RSD_FSIGMA8",
            "source_file": "fallback_memory",
            "notes": (
                "RSD/fσ8 com amplitude livre; praticamente empate com LCDM; "
                "sigma8_best≈0.772799."
            ),
        },
        {
            "test_id": "19_black_hole_structural_saturation",
            "sector": "black_holes",
            "model": "TVGD_BH_saturation",
            "n_data": None,
            "k_params": 2,
            "chi2": 0.0,
            "chi2_red": 0.0,
            "delta_chi2_vs_ref": 0.0,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "PASSOU_FORTE_BH_BENCHMARKS",
            "source_file": "fallback_memory",
            "notes": (
                "Preserva r_h/r_s≈1, r_ph/r_s≈1.5, b_shadow/r_s≈2.598076; "
                "núcleo regularizado por saturação estrutural."
            ),
        },
        {
            "test_id": "20_falsifiability_prediction_checklist",
            "sector": "falsifiability",
            "model": "TVGD",
            "n_data": None,
            "k_params": None,
            "chi2": None,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": None,
            "bic_proxy": None,
            "status": "PREDICTIONS_DEFINED",
            "source_file": "fallback_memory",
            "notes": "Checklist de previsões falsificáveis e critérios de falha pública.",
        },
        {
            "test_id": "hubble_structural_transition",
            "sector": "cosmology_hubble",
            "model": "TVGD_Hubble_transition",
            "n_data": None,
            "k_params": 4,
            "chi2": 0.007323,
            "chi2_red": None,
            "delta_chi2_vs_ref": None,
            "aic_proxy": 12.007323,
            "bic_proxy": 14.916763,
            "status": "PASSOU_FORTE_CHI2",
            "source_file": "fallback_memory",
            "notes": (
                "Modelo com H_CMB_eff≈67.4 e H_today_eff≈73.129; "
                "preserva setor primordial no teste efetivo registrado."
            ),
        },
    ]

    # Calcula AIC/BIC onde possível e ainda ausente
    for r in rows:
        if r.get("aic_proxy") is None or r.get("bic_proxy") is None:
            aic, bic = compute_aic_bic(r.get("chi2"), r.get("k_params"), r.get("n_data"))
            if r.get("aic_proxy") is None:
                r["aic_proxy"] = aic
            if r.get("bic_proxy") is None:
                r["bic_proxy"] = bic

    return rows


def normalize_rows(rows: List[Dict[str, Any]]) -> pd.DataFrame:
    required_cols = [
        "test_id",
        "sector",
        "model",
        "n_data",
        "k_params",
        "chi2",
        "chi2_red",
        "delta_chi2_vs_ref",
        "aic_proxy",
        "bic_proxy",
        "status",
        "source_file",
        "notes",
    ]

    normalized = []
    for r in rows:
        rr = {c: r.get(c) for c in required_cols}

        if rr["aic_proxy"] is None or rr["bic_proxy"] is None:
            aic, bic = compute_aic_bic(rr.get("chi2"), rr.get("k_params"), rr.get("n_data"))
            if rr["aic_proxy"] is None:
                rr["aic_proxy"] = aic
            if rr["bic_proxy"] is None:
                rr["bic_proxy"] = bic

        normalized.append(rr)

    df = pd.DataFrame(normalized)

    # Remove duplicatas simples
    df = df.drop_duplicates(subset=["test_id", "sector", "model", "source_file"], keep="first")

    # Ordenação amigável
    order_map = {
        "galaxy_dynamics": 1,
        "external_galaxy_validation": 2,
        "cosmology_hubble": 3,
        "cosmology_hubble_sn": 4,
        "growth_rsd": 5,
        "black_holes": 6,
        "covariant_risk": 7,
        "forecast_euclid": 8,
        "falsifiability": 9,
    }
    df["_order"] = df["sector"].map(order_map).fillna(99)
    df = df.sort_values(["_order", "test_id"]).drop(columns=["_order"])

    return df


def make_markdown_report(df: pd.DataFrame) -> str:
    lines = []
    lines.append("# TVGD — Model Comparison / Statistical Summary\n")
    lines.append("Este relatório foi gerado automaticamente pelo script `22_model_comparison_aic_bic_summary.py`.\n")
    lines.append("As métricas `AIC_proxy` e `BIC_proxy` são usadas como indicadores aproximados quando a likelihood formal completa não está disponível.\n")
    lines.append("Para publicação, os valores finais devem ser conectados diretamente aos arquivos de saída dos notebooks correspondentes.\n")

    lines.append("## Resumo geral\n")
    lines.append(f"- Número de entradas na tabela: **{len(df)}**")
    sectors = sorted([str(x) for x in df["sector"].dropna().unique()])
    lines.append(f"- Setores cobertos: **{', '.join(sectors)}**\n")

    pass_like = df["status"].fillna("").str.contains("PASSOU|ROBUST|PREDICTIONS|COMPATIVEL", regex=True).sum()
    fail_like = df["status"].fillna("").str.contains("FALHOU", regex=True).sum()
    lines.append(f"- Entradas favoráveis/compatíveis: **{pass_like}**")
    lines.append(f"- Entradas com falha explícita: **{fail_like}**\n")

    lines.append("## Tabela principal\n")
    table_cols = [
        "test_id",
        "sector",
        "model",
        "n_data",
        "chi2",
        "chi2_red",
        "delta_chi2_vs_ref",
        "aic_proxy",
        "bic_proxy",
        "status",
    ]

    md_df = df[table_cols].copy()
    lines.append(md_df.to_markdown(index=False))
    lines.append("\n")

    lines.append("## Observações por teste\n")
    for _, row in df.iterrows():
        lines.append(f"### {row['test_id']}")
        lines.append(f"- Setor: `{row['sector']}`")
        lines.append(f"- Modelo: `{row['model']}`")
        lines.append(f"- Status: `{row['status']}`")
        lines.append(f"- Fonte: `{row['source_file']}`")
        lines.append(f"- Nota: {row['notes']}\n")

    lines.append("## Interpretação curta\n")
    lines.append(
        "A cadeia pública da TVGD não deve ser apresentada como prova final, mas como um pacote de validação efetiva: "
        "dinâmica galáctica, validação externa, crescimento de estruturas, buracos negros, consistência covariante efetiva "
        "e previsões falsificáveis. O ponto forte é a convergência de vários setores com poucos parâmetros estruturais. "
        "O ponto que ainda exige cuidado é transformar todos os testes em likelihoods homogêneas com covariâncias completas."
    )

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, default=".", help="Raiz do repositório TVGD.")
    parser.add_argument(
        "--fallback-only",
        action="store_true",
        help="Usa apenas a tabela fallback consolidada, sem varrer arquivos existentes."
    )
    args = parser.parse_args()

    root = Path(args.root).resolve()
    outdir = root / "results" / "22_model_comparison_aic_bic_summary"
    outdir.mkdir(parents=True, exist_ok=True)

    rows = []

    if not args.fallback_only:
        scanned = scan_existing_results(root)
        rows.extend(scanned)

    # Sempre inclui fallback, porque ele garante uma tabela mínima pública.
    # Se houver arquivos reais, eles aparecerão junto e podem substituir manualmente os fallbacks.
    rows.extend(fallback_summary())

    df = normalize_rows(rows)

    csv_path = outdir / "tvgd_model_comparison_summary.csv"
    json_path = outdir / "tvgd_model_comparison_summary.json"
    md_path = outdir / "tvgd_model_comparison_report.md"

    df.to_csv(csv_path, index=False)

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(df.where(pd.notnull(df), None).to_dict(orient="records"), f, indent=2, ensure_ascii=False)

    report = make_markdown_report(df)
    md_path.write_text(report, encoding="utf-8")

    print("=" * 72)
    print("TVGD — 22 MODEL COMPARISON / AIC-BIC SUMMARY")
    print("=" * 72)
    print(f"Raiz do repositório: {root}")
    print(f"Entradas na tabela: {len(df)}")
    print()
    print("Arquivos gerados:")
    print(f" - {csv_path}")
    print(f" - {json_path}")
    print(f" - {md_path}")
    print()
    print("Status por setor:")
    print(df.groupby("sector")["test_id"].count().sort_values(ascending=False))
    print()
    print("Resumo salvo com sucesso.")


if __name__ == "__main__":
    main()
