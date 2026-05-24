#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TVGD — 23_reproducibility_quickstart.py

Objetivo
--------
Fornecer um teste rápido de reprodutibilidade para o repositório público da TVGD.

Este script verifica:

1. Estrutura básica de pastas.
2. Existência dos principais scripts de validação.
3. Dependências Python mínimas.
4. Sintaxe dos scripts principais.
5. Possibilidade de rodar um subconjunto leve, se ativado com --run-basic.

Ele NÃO baixa dados pesados automaticamente.
Ele NÃO roda CLASS/CAMB/Cobaya.
Ele NÃO exige arquivos grandes como SPARC, Pantheon+, HSC ou COSMOS-Web.

Saídas
------
Cria a pasta:

    results/23_reproducibility_quickstart/

com:

    tvgd_reproducibility_quickstart_report.md
    tvgd_reproducibility_quickstart_summary.json
    tvgd_reproducibility_quickstart_summary.csv

Uso
---
Na raiz do repositório:

    python notebooks/main_validation/23_reproducibility_quickstart.py

Opcional:

    python notebooks/main_validation/23_reproducibility_quickstart.py --run-basic

    python notebooks/main_validation/23_reproducibility_quickstart.py --root .
"""

import argparse
import ast
import importlib.util
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Dict, Any, List

import pandas as pd


# ============================================================
# Scripts principais esperados na pasta notebooks/main_validation
# ============================================================

EXPECTED_MAIN_VALIDATION_SCRIPTS = [
    "00_baseline_consistency_macro_micro.py",
    "01_sparc_rar_tvgd.py",
    "01B_sparc_ml_robustness.py",
    "01C_sparc_rar_comparison.py",
    "02_black_hole_regular_core.py",
    "03_cosmology_background_bao.py",
    "04B_controlled_strong_lensing_selected.py",
    "05_stability_causality_screening.py",
    "06_nonphantom_cosmology_background.py",
    "07B_hubble_structural_phase_transition.py",
    "07E_pantheon_full_covariance_hubble.py",
    "07G_pantheon_noncalibrator_shape.py",
    "07H_hubble_structural_phase_robustness.py",
    "08_black_hole_structural_saturation.py",
    "08B_rotating_black_hole_structural_saturation.py",
    "09_falsifiability_prediction_checklist.py",
    "10_euclid_forecast_structural_response.py",
    "11_rmax_sparc_universality.py",
    "12_rmax_transition_beta_robustness.py",
    "15_extrisk_little_things_external.py",
    "17_covrisk_scalar_vector_ppn.py",
    "18_covrisk_rsd_fsigma8.py",
    "19_black_hole_structural_saturation.py",
    "20_falsifiability_prediction_checklist.py",
    "22_model_comparison_aic_bic_summary.py",
    "23_reproducibility_quickstart.py",
]


# ============================================================
# Dependências mínimas e opcionais
# ============================================================

MINIMAL_DEPENDENCIES = [
    "numpy",
    "pandas",
    "matplotlib",
]

OPTIONAL_DEPENDENCIES = [
    "scipy",
    "astropy",
    "requests",
]


# ============================================================
# Funções auxiliares
# ============================================================

def check_import(module_name: str) -> Dict[str, Any]:
    spec = importlib.util.find_spec(module_name)

    return {
        "item": module_name,
        "category": "dependency",
        "required": True,
        "ok": spec is not None,
        "details": "installed" if spec is not None else "missing",
    }


def check_optional_import(module_name: str) -> Dict[str, Any]:
    spec = importlib.util.find_spec(module_name)

    return {
        "item": module_name,
        "category": "optional_dependency",
        "required": False,
        "ok": spec is not None,
        "details": "installed" if spec is not None else "missing",
    }


def check_path_exists(path: Path, label: str, required: bool = True) -> Dict[str, Any]:
    return {
        "item": label,
        "category": "path",
        "required": required,
        "ok": path.exists(),
        "details": str(path),
    }


def check_python_syntax(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {
            "item": path.name,
            "category": "syntax",
            "required": True,
            "ok": False,
            "details": "file missing",
        }

    try:
        text = path.read_text(encoding="utf-8", errors="ignore")
        ast.parse(text, filename=str(path))

        return {
            "item": path.name,
            "category": "syntax",
            "required": True,
            "ok": True,
            "details": "syntax ok",
        }

    except SyntaxError as e:
        return {
            "item": path.name,
            "category": "syntax",
            "required": True,
            "ok": False,
            "details": f"SyntaxError: {e}",
        }

    except Exception as e:
        return {
            "item": path.name,
            "category": "syntax",
            "required": True,
            "ok": False,
            "details": f"error reading/parsing file: {e}",
        }


def run_script(path: Path, root: Path, timeout: int = 120) -> Dict[str, Any]:
    if not path.exists():
        return {
            "item": path.name,
            "category": "execution",
            "required": False,
            "ok": False,
            "details": "file missing",
        }

    try:
        proc = subprocess.run(
            [sys.executable, str(path)],
            cwd=str(root),
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        ok = proc.returncode == 0

        stdout_tail = proc.stdout[-1500:] if proc.stdout else ""
        stderr_tail = proc.stderr[-1500:] if proc.stderr else ""

        return {
            "item": path.name,
            "category": "execution",
            "required": False,
            "ok": ok,
            "details": (
                f"returncode={proc.returncode}\n\n"
                f"stdout_tail:\n{stdout_tail}\n\n"
                f"stderr_tail:\n{stderr_tail}"
            ),
        }

    except subprocess.TimeoutExpired:
        return {
            "item": path.name,
            "category": "execution",
            "required": False,
            "ok": False,
            "details": f"timeout after {timeout} seconds",
        }

    except Exception as e:
        return {
            "item": path.name,
            "category": "execution",
            "required": False,
            "ok": False,
            "details": f"execution error: {e}",
        }


def make_markdown_report(results: List[Dict[str, Any]], root: Path) -> str:
    df = pd.DataFrame(results)

    required_df = df[df["required"] == True]
    optional_df = df[df["required"] == False]

    n_required = len(required_df)
    n_required_ok = int(required_df["ok"].sum()) if n_required else 0

    n_optional = len(optional_df)
    n_optional_ok = int(optional_df["ok"].sum()) if n_optional else 0

    lines = []

    lines.append("# TVGD — Reproducibility Quickstart Report\n")
    lines.append("Relatório gerado automaticamente por `23_reproducibility_quickstart.py`.\n")

    lines.append("## Ambiente\n")
    lines.append(f"- Python: `{sys.version.split()[0]}`")
    lines.append(f"- Plataforma: `{platform.platform()}`")
    lines.append(f"- Raiz do repositório: `{root}`\n")

    lines.append("## Resultado geral\n")
    lines.append(f"- Itens obrigatórios aprovados: **{n_required_ok}/{n_required}**")
    lines.append(f"- Itens opcionais disponíveis/aprovados: **{n_optional_ok}/{n_optional}**\n")

    if n_required_ok == n_required:
        lines.append("**Status geral:** `REPRODUCIBILITY_QUICKSTART_PASSOU`\n")
    else:
        lines.append("**Status geral:** `REPRODUCIBILITY_QUICKSTART_COM_PENDENCIAS`\n")

    lines.append("## Tabela de verificações\n")

    table_cols = [
        "category",
        "item",
        "required",
        "ok",
        "details",
    ]

    lines.append(df[table_cols].to_markdown(index=False))
    lines.append("\n")

    lines.append("## Interpretação\n")
    lines.append(
        "Este quickstart não valida fisicamente todos os resultados da TVGD. "
        "Ele verifica apenas se a estrutura pública do repositório está coerente, "
        "se as dependências mínimas existem, se os principais scripts estão presentes "
        "e se eles não possuem erro sintático."
    )

    lines.append("\n")
    lines.append(
        "Os testes científicos completos devem ser rodados individualmente, principalmente "
        "os que dependem de SPARC, LITTLE THINGS, Pantheon+, RSD/fσ8, HSC, COSMOS-Web, "
        "CLASS, CAMB ou Cobaya."
    )

    lines.append("\n## Próximos passos recomendados\n")
    lines.append("1. Rodar primeiro os scripts leves.")
    lines.append("2. Depois rodar os scripts com dados externos.")
    lines.append("3. Conferir se cada script salva CSV, JSON ou Markdown em `results/`.")
    lines.append("4. Atualizar o `README.md` com o comando deste quickstart.")
    lines.append("5. Usar o arquivo `22_model_comparison_aic_bic_summary.py` como painel estatístico público.")

    return "\n".join(lines)


# ============================================================
# Programa principal
# ============================================================

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--root",
        type=str,
        default=".",
        help="Raiz do repositório TVGD.",
    )

    parser.add_argument(
        "--run-basic",
        action="store_true",
        help="Roda um subconjunto leve de scripts públicos. Por padrão, apenas checa estrutura e sintaxe.",
    )

    args = parser.parse_args()

    root = Path(args.root).resolve()
    validation_dir = root / "notebooks" / "main_validation"

    outdir = root / "results" / "23_reproducibility_quickstart"
    outdir.mkdir(parents=True, exist_ok=True)

    results: List[Dict[str, Any]] = []

    # ------------------------------------------------------------
    # 1. Estrutura básica do repositório
    # ------------------------------------------------------------

    results.append(check_path_exists(root / "notebooks", "notebooks/", required=True))
    results.append(check_path_exists(validation_dir, "notebooks/main_validation/", required=True))
    results.append(check_path_exists(root / "results", "results/", required=False))
    results.append(check_path_exists(root / "README.md", "README.md", required=False))
    results.append(check_path_exists(root / "requirements.txt", "requirements.txt", required=False))
    results.append(check_path_exists(root / "LICENSE", "LICENSE", required=False))

    # ------------------------------------------------------------
    # 2. Dependências Python
    # ------------------------------------------------------------

    for dep in MINIMAL_DEPENDENCIES:
        results.append(check_import(dep))

    for dep in OPTIONAL_DEPENDENCIES:
        results.append(check_optional_import(dep))

    # ------------------------------------------------------------
    # 3. Existência e sintaxe dos scripts principais
    # ------------------------------------------------------------

    for script_name in EXPECTED_MAIN_VALIDATION_SCRIPTS:
        script_path = validation_dir / script_name

        results.append({
            "item": script_name,
            "category": "script_exists",
            "required": True,
            "ok": script_path.exists(),
            "details": str(script_path),
        })

        if script_path.exists():
            results.append(check_python_syntax(script_path))

    # ------------------------------------------------------------
    # 4. Execução opcional dos scripts leves
    # ------------------------------------------------------------

    if args.run_basic:
        basic_scripts = [
            "00_baseline_consistency_macro_micro.py",
            "02_black_hole_regular_core.py",
            "05_stability_causality_screening.py",
            "06_nonphantom_cosmology_background.py",
            "09_falsifiability_prediction_checklist.py",
            "19_black_hole_structural_saturation.py",
            "20_falsifiability_prediction_checklist.py",
            "22_model_comparison_aic_bic_summary.py",
        ]

        for script_name in basic_scripts:
            script_path = validation_dir / script_name
            results.append(run_script(script_path, root=root, timeout=120))

    # ------------------------------------------------------------
    # 5. Salvar arquivos de saída
    # ------------------------------------------------------------

    df = pd.DataFrame(results)

    csv_path = outdir / "tvgd_reproducibility_quickstart_summary.csv"
    json_path = outdir / "tvgd_reproducibility_quickstart_summary.json"
    md_path = outdir / "tvgd_reproducibility_quickstart_report.md"

    df.to_csv(csv_path, index=False)

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    report = make_markdown_report(results, root)
    md_path.write_text(report, encoding="utf-8")

    # ------------------------------------------------------------
    # 6. Resumo no terminal
    # ------------------------------------------------------------

    required_df = df[df["required"] == True]

    n_required = len(required_df)
    n_required_ok = int(required_df["ok"].sum()) if n_required else 0

    print("=" * 72)
    print("TVGD — 23 REPRODUCIBILITY QUICKSTART")
    print("=" * 72)
    print(f"Raiz do repositório: {root}")
    print(f"Itens obrigatórios OK: {n_required_ok}/{n_required}")

    if n_required_ok == n_required:
        print("STATUS: REPRODUCIBILITY_QUICKSTART_PASSOU")
    else:
        print("STATUS: REPRODUCIBILITY_QUICKSTART_COM_PENDENCIAS")

    print()
    print("Arquivos gerados:")
    print(f" - {csv_path}")
    print(f" - {json_path}")
    print(f" - {md_path}")
    print()

    if not args.run_basic:
        print("Execução de scripts leves não ativada.")
        print("Para rodar o subconjunto básico, use:")
        print()
        print("  python notebooks/main_validation/23_reproducibility_quickstart.py --run-basic")
        print()

    print("Finalizado.")


if __name__ == "__main__":
    main()
