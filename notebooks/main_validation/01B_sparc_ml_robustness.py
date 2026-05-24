#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TVGD — 01B_sparc_ml_robustness.py

Objetivo
--------
Testar a robustez da relação de aceleração radial da TVGD sob variações simples
da razão massa-luminosidade estelar M*/L no setor SPARC.

A ideia é verificar se a resposta

    P(u) = 1 - exp(-u)

com

    u = sqrt(g_bar/a0)

permanece estável quando as componentes estelares do campo bariônico são
reescaladas dentro de uma faixa razoável.

O script é seguro para repositório público:

1. Tenta localizar automaticamente um arquivo SPARC/RAR já processado.
2. Se não encontrar dados reais, roda um modo demonstrativo sintético.
3. Salva CSV, JSON, Markdown e figura em results/01B_sparc_ml_robustness/.
4. Não baixa dados pesados automaticamente.

Entradas aceitas
----------------
O script procura arquivos CSV em locais comuns, por exemplo:

    data/processed/sparc_rar.csv
    data/external/sparc_rar.csv
    results/01_sparc_rar_tvgd/sparc_rar_points.csv
    results/01_sparc_rar_tvgd/tvgd_sparc_rar_points.csv

Colunas aceitas, em ordem de preferência:

Caso A — componentes separadas:
    g_obs, g_gas, g_disk, g_bulge

Caso B — campo bariônico total:
    g_obs, g_bar

Nomes alternativos são reconhecidos automaticamente, como:
    gobs, gtot, g_bar, gbar, ggas, gdisk, gbulge, g_star_disk, g_star_bulge

Unidades
--------
O script assume que g_obs, g_bar e componentes estão na mesma unidade.
A constante a0 é definida na mesma unidade dos dados. Por padrão:

    a0 = 1.2e-10

Se seus dados estiverem normalizados por a0, use:

    --a0 1.0

Uso
---
Na raiz do repositório:

    python cadernos/main_validation/01B_sparc_ml_robustness.py

ou:

    python notebooks/main_validation/01B_sparc_ml_robustness.py

Opções úteis:

    python cadernos/main_validation/01B_sparc_ml_robustness.py --input data/processed/sparc_rar.csv
    python cadernos/main_validation/01B_sparc_ml_robustness.py --a0 1.0
    python cadernos/main_validation/01B_sparc_ml_robustness.py --no-plot

Saídas
------
    results/01B_sparc_ml_robustness/
        ml_robustness_grid.csv
        ml_robustness_summary.json
        ml_robustness_report.md
        ml_robustness_heatmap.png
"""

import argparse
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False


# ============================================================
# Configurações principais
# ============================================================

DEFAULT_A0 = 1.2e-10

DEFAULT_DISK_GRID = np.round(np.linspace(0.30, 0.80, 11), 3)
DEFAULT_BULGE_GRID = np.round(np.linspace(0.30, 0.90, 13), 3)

DEFAULT_GLOBAL_ML_GRID = np.round(np.linspace(0.50, 1.50, 21), 3)

COMMON_INPUTS = [
    "data/processed/sparc_rar.csv",
    "data/processed/tvgd_sparc_rar_points.csv",
    "data/external/sparc_rar.csv",
    "data/raw/sparc_rar.csv",
    "results/01_sparc_rar_tvgd/sparc_rar_points.csv",
    "results/01_sparc_rar_tvgd/tvgd_sparc_rar_points.csv",
    "results/01_sparc_rar_tvgd/rar_points.csv",
]


COLUMN_ALIASES = {
    "g_obs": [
        "g_obs", "gobs", "g_obs_m_s2", "gobs_m_s2", "gtot", "g_tot",
        "gdagger_obs", "g_measured", "grot", "g_rot"
    ],
    "g_bar": [
        "g_bar", "gbar", "g_bar_m_s2", "gbar_m_s2", "gbary", "g_bary",
        "g_baryon", "g_baryonic"
    ],
    "g_gas": [
        "g_gas", "ggas", "gas", "gHI", "g_hi", "g_atomic", "g_gas_m_s2"
    ],
    "g_disk": [
        "g_disk", "gdisk", "disk", "g_star_disk", "gdisk_m_s2",
        "g_stellar_disk", "g_d"
    ],
    "g_bulge": [
        "g_bulge", "gbulge", "bulge", "g_star_bulge", "gbulge_m_s2",
        "g_stellar_bulge", "g_b"
    ],
    "galaxy": [
        "galaxy", "name", "gal", "object", "id", "Galaxy", "Name"
    ],
    "radius": [
        "radius", "r", "R", "rad", "R_kpc", "r_kpc"
    ],
    "err_g_obs": [
        "err_g_obs", "g_obs_err", "e_g_obs", "sigma_g_obs", "dgobs",
        "gobs_err"
    ],
}


# ============================================================
# Utilidades
# ============================================================

def find_column(df: pd.DataFrame, logical_name: str) -> Optional[str]:
    aliases = COLUMN_ALIASES.get(logical_name, [])
    lower_map = {c.lower(): c for c in df.columns}

    for alias in aliases:
        if alias in df.columns:
            return alias
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]

    return None


def safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def tvgd_P(u: np.ndarray) -> np.ndarray:
    """
    Função de resposta TVGD.
    Usa expm1 para estabilidade numérica:
        1 - exp(-u) = -expm1(-u)
    """
    u = np.asarray(u, dtype=float)
    return -np.expm1(-np.clip(u, 0.0, 700.0))


def tvgd_g_from_gbar(gbar: np.ndarray, a0: float) -> np.ndarray:
    gbar = np.asarray(gbar, dtype=float)
    u = np.sqrt(np.maximum(gbar, 0.0) / a0)
    P = tvgd_P(u)
    P = np.maximum(P, 1e-300)
    return gbar / P


def log10_safe(x: np.ndarray) -> np.ndarray:
    return np.log10(np.maximum(np.asarray(x, dtype=float), 1e-300))


def robust_stats(gobs: np.ndarray, gpred: np.ndarray) -> Dict[str, float]:
    """
    Estatísticas em log10 para comparação RAR.
    """
    log_resid = log10_safe(gpred) - log10_safe(gobs)
    frac_resid = (gpred - gobs) / np.maximum(np.abs(gobs), 1e-300)

    rmse_log = float(np.sqrt(np.nanmean(log_resid**2)))
    mae_log = float(np.nanmean(np.abs(log_resid)))
    bias_log = float(np.nanmean(log_resid))
    med_abs_log = float(np.nanmedian(np.abs(log_resid)))

    mean_frac = float(np.nanmean(frac_resid))
    mae_frac = float(np.nanmean(np.abs(frac_resid)))
    med_abs_frac = float(np.nanmedian(np.abs(frac_resid)))

    # Chi2 proxy: usa dispersão log observada efetiva de 0.10 dex se não houver erro formal.
    sigma_log = 0.10
    chi2_proxy = float(np.nansum((log_resid / sigma_log) ** 2))
    n = int(np.sum(np.isfinite(log_resid)))
    k = 1
    chi2_red_proxy = float(chi2_proxy / max(n - k, 1))

    return {
        "n_points": n,
        "rmse_log10": rmse_log,
        "mae_log10": mae_log,
        "median_abs_log10": med_abs_log,
        "bias_log10": bias_log,
        "mean_frac_resid": mean_frac,
        "mean_abs_frac_resid": mae_frac,
        "median_abs_frac_resid": med_abs_frac,
        "chi2_proxy_sigma_0p10dex": chi2_proxy,
        "chi2red_proxy_sigma_0p10dex": chi2_red_proxy,
    }


def locate_input(root: Path, user_input: Optional[str]) -> Optional[Path]:
    if user_input:
        p = Path(user_input)
        if not p.is_absolute():
            p = root / p
        if p.exists():
            return p
        return None

    for rel in COMMON_INPUTS:
        p = root / rel
        if p.exists():
            return p

    # Busca leve por CSVs com nomes prováveis.
    candidates = []
    for folder in ["data", "results", "Resultados"]:
        d = root / folder
        if d.exists():
            for p in d.rglob("*.csv"):
                low = p.name.lower()
                if any(tok in low for tok in ["sparc", "rar", "tvgd"]):
                    candidates.append(p)

    if candidates:
        candidates = sorted(candidates, key=lambda x: len(str(x)))
        return candidates[0]

    return None


def make_synthetic_data(n: int = 600, seed: int = 42, a0: float = DEFAULT_A0) -> pd.DataFrame:
    """
    Cria uma amostra sintética simples para o script ser executável mesmo sem SPARC.
    Essa amostra NÃO deve ser usada como evidência científica.
    """
    rng = np.random.default_rng(seed)

    log_gbar = rng.uniform(-13.0, -9.0, size=n)
    gbar_true = 10.0 ** log_gbar

    # Divide gbar em componentes fictícias
    f_gas = rng.uniform(0.10, 0.70, size=n)
    f_disk = rng.uniform(0.20, 0.80, size=n) * (1.0 - f_gas)
    f_bulge = np.maximum(1.0 - f_gas - f_disk, 0.0)

    g_gas = f_gas * gbar_true
    g_disk = f_disk * gbar_true / 0.50
    g_bulge = f_bulge * gbar_true / 0.70

    gpred = tvgd_g_from_gbar(gbar_true, a0=a0)
    scatter = rng.normal(0.0, 0.08, size=n)
    gobs = gpred * (10.0 ** scatter)

    galaxies = np.array([f"SYN_{i:03d}" for i in rng.integers(0, 60, size=n)])
    radius = rng.uniform(0.5, 30.0, size=n)

    return pd.DataFrame({
        "galaxy": galaxies,
        "radius": radius,
        "g_obs": gobs,
        "g_gas": g_gas,
        "g_disk": g_disk,
        "g_bulge": g_bulge,
        "synthetic": True,
    })


def load_and_standardize(path: Optional[Path], a0: float) -> Tuple[pd.DataFrame, Dict[str, str], bool]:
    """
    Retorna:
        df padronizado
        mapeamento de colunas
        synthetic_mode
    """
    if path is None:
        df = make_synthetic_data(a0=a0)
        mapping = {
            "g_obs": "g_obs",
            "g_gas": "g_gas",
            "g_disk": "g_disk",
            "g_bulge": "g_bulge",
            "galaxy": "galaxy",
            "radius": "radius",
        }
        return df, mapping, True

    df_raw = pd.read_csv(path)
    mapping = {}

    for key in COLUMN_ALIASES:
        col = find_column(df_raw, key)
        if col is not None:
            mapping[key] = col

    if "g_obs" not in mapping:
        raise ValueError(
            "Não encontrei coluna de aceleração observada. "
            "Use uma coluna como g_obs, gobs, g_tot ou g_rot."
        )

    has_components = all(k in mapping for k in ["g_gas", "g_disk", "g_bulge"])
    has_gbar = "g_bar" in mapping

    if not has_components and not has_gbar:
        raise ValueError(
            "Não encontrei componentes bariônicas nem g_bar total. "
            "Use colunas g_obs + g_bar, ou g_obs + g_gas + g_disk + g_bulge."
        )

    out = pd.DataFrame()
    out["g_obs"] = safe_numeric(df_raw[mapping["g_obs"]])

    if "galaxy" in mapping:
        out["galaxy"] = df_raw[mapping["galaxy"]].astype(str)
    else:
        out["galaxy"] = "unknown"

    if "radius" in mapping:
        out["radius"] = safe_numeric(df_raw[mapping["radius"]])
    else:
        out["radius"] = np.nan

    if has_components:
        out["g_gas"] = safe_numeric(df_raw[mapping["g_gas"]])
        out["g_disk"] = safe_numeric(df_raw[mapping["g_disk"]])
        out["g_bulge"] = safe_numeric(df_raw[mapping["g_bulge"]])
    else:
        out["g_bar"] = safe_numeric(df_raw[mapping["g_bar"]])

    if "err_g_obs" in mapping:
        out["err_g_obs"] = safe_numeric(df_raw[mapping["err_g_obs"]])

    out["synthetic"] = False

    return out, mapping, False


def clean_data(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["g_obs"]
    if all(c in df.columns for c in ["g_gas", "g_disk", "g_bulge"]):
        cols += ["g_gas", "g_disk", "g_bulge"]
    else:
        cols += ["g_bar"]

    mask = np.ones(len(df), dtype=bool)
    for c in cols:
        mask &= np.isfinite(df[c].to_numpy(dtype=float))
        mask &= df[c].to_numpy(dtype=float) > 0.0

    out = df.loc[mask].copy()
    out = out.reset_index(drop=True)
    return out


def run_component_grid(df: pd.DataFrame, a0: float, disk_grid: np.ndarray, bulge_grid: np.ndarray) -> pd.DataFrame:
    rows = []

    gobs = df["g_obs"].to_numpy(dtype=float)
    ggas = df["g_gas"].to_numpy(dtype=float)
    gdisk = df["g_disk"].to_numpy(dtype=float)
    gbulge = df["g_bulge"].to_numpy(dtype=float)

    for y_disk in disk_grid:
        for y_bulge in bulge_grid:
            gbar = ggas + y_disk * gdisk + y_bulge * gbulge
            gpred = tvgd_g_from_gbar(gbar, a0=a0)
            st = robust_stats(gobs, gpred)

            rows.append({
                "mode": "component_ml_grid",
                "upsilon_disk": float(y_disk),
                "upsilon_bulge": float(y_bulge),
                "global_ml_scale": np.nan,
                **st,
            })

    return pd.DataFrame(rows)


def run_global_grid(df: pd.DataFrame, a0: float, global_grid: np.ndarray) -> pd.DataFrame:
    rows = []

    gobs = df["g_obs"].to_numpy(dtype=float)
    gbar_base = df["g_bar"].to_numpy(dtype=float)

    for scale in global_grid:
        gbar = scale * gbar_base
        gpred = tvgd_g_from_gbar(gbar, a0=a0)
        st = robust_stats(gobs, gpred)

        rows.append({
            "mode": "global_gbar_scale_grid",
            "upsilon_disk": np.nan,
            "upsilon_bulge": np.nan,
            "global_ml_scale": float(scale),
            **st,
        })

    return pd.DataFrame(rows)


def summarize_grid(grid: pd.DataFrame, synthetic_mode: bool, input_path: Optional[Path], a0: float) -> Dict[str, object]:
    best_idx = grid["rmse_log10"].idxmin()
    best = grid.loc[best_idx].to_dict()

    rmse = grid["rmse_log10"].to_numpy(dtype=float)
    chi2red = grid["chi2red_proxy_sigma_0p10dex"].to_numpy(dtype=float)

    # Robustez: fração da grade dentro de 10%, 20% e 35% do melhor RMSE.
    best_rmse = float(np.nanmin(rmse))
    frac_within_10pct = float(np.mean(rmse <= 1.10 * best_rmse))
    frac_within_20pct = float(np.mean(rmse <= 1.20 * best_rmse))
    frac_within_35pct = float(np.mean(rmse <= 1.35 * best_rmse))

    # Classificação simples e honesta
    if frac_within_20pct >= 0.50 and best_rmse < 0.15:
        status = "PASSOU_FORTE_ML_ROBUSTNESS"
    elif frac_within_20pct >= 0.25 and best_rmse < 0.25:
        status = "PASSOU_MODERADO_ML_ROBUSTNESS"
    elif best_rmse < 0.35:
        status = "COMPATIVEL_FRACO_ML_ROBUSTNESS"
    else:
        status = "SENSIVEL_OU_FALHOU_ML_ROBUSTNESS"

    if synthetic_mode:
        status = "DEMO_SINTETICO_NAO_CIENTIFICO"

    summary = {
        "script": "01B_sparc_ml_robustness.py",
        "input_path": str(input_path) if input_path is not None else None,
        "synthetic_mode": bool(synthetic_mode),
        "a0": float(a0),
        "n_grid": int(len(grid)),
        "n_points": int(best.get("n_points", 0)),
        "best": {
            k: (float(v) if isinstance(v, (np.floating, float, int, np.integer)) and math.isfinite(float(v)) else None)
            for k, v in best.items()
        },
        "rmse_log10_min": float(np.nanmin(rmse)),
        "rmse_log10_median": float(np.nanmedian(rmse)),
        "rmse_log10_p16": float(np.nanpercentile(rmse, 16)),
        "rmse_log10_p84": float(np.nanpercentile(rmse, 84)),
        "chi2red_proxy_min": float(np.nanmin(chi2red)),
        "chi2red_proxy_median": float(np.nanmedian(chi2red)),
        "frac_grid_within_10pct_best_rmse": frac_within_10pct,
        "frac_grid_within_20pct_best_rmse": frac_within_20pct,
        "frac_grid_within_35pct_best_rmse": frac_within_35pct,
        "final_status": status,
    }

    return summary


def write_report(outdir: Path, summary: Dict[str, object], grid: pd.DataFrame, mapping: Dict[str, str]) -> Path:
    path = outdir / "ml_robustness_report.md"

    best = summary["best"]

    lines = []
    lines.append("# TVGD — 01B SPARC M/L Robustness\n")
    lines.append("Este relatório foi gerado automaticamente por `01B_sparc_ml_robustness.py`.\n")

    lines.append("## 1. Objetivo\n")
    lines.append(
        "Testar se a relação de aceleração radial da TVGD permanece estável quando "
        "a contribuição estelar do campo bariônico é reescalada por fatores simples de massa-luminosidade."
    )
    lines.append("\n")

    lines.append("## 2. Equação testada\n")
    lines.append("```math\nP(u)=1-e^{-u}\n```\n")
    lines.append("```math\nu=\\sqrt{g_{\\rm bar}/a_0}\n```\n")
    lines.append("```math\ng_{\\rm TVGD}=\\frac{g_{\\rm bar}}{P(u)}\n```\n")

    lines.append("## 3. Entrada\n")
    lines.append(f"- Arquivo de entrada: `{summary['input_path']}`")
    lines.append(f"- Modo sintético: `{summary['synthetic_mode']}`")
    lines.append(f"- a0 usado: `{summary['a0']}`")
    lines.append(f"- Pontos usados: `{summary['n_points']}`")
    lines.append(f"- Colunas detectadas: `{mapping}`\n")

    lines.append("## 4. Melhor ponto da grade\n")
    lines.append(f"- Modo: `{best.get('mode')}`")
    lines.append(f"- upsilon_disk: `{best.get('upsilon_disk')}`")
    lines.append(f"- upsilon_bulge: `{best.get('upsilon_bulge')}`")
    lines.append(f"- global_ml_scale: `{best.get('global_ml_scale')}`")
    lines.append(f"- RMSE log10: `{best.get('rmse_log10')}`")
    lines.append(f"- MAE log10: `{best.get('mae_log10')}`")
    lines.append(f"- Bias log10: `{best.get('bias_log10')}`")
    lines.append(f"- Chi2red proxy, sigma=0.10 dex: `{best.get('chi2red_proxy_sigma_0p10dex')}`\n")

    lines.append("## 5. Robustez da grade\n")
    lines.append(f"- RMSE mínimo: `{summary['rmse_log10_min']}`")
    lines.append(f"- RMSE mediano: `{summary['rmse_log10_median']}`")
    lines.append(f"- Fração da grade dentro de 10% do melhor RMSE: `{summary['frac_grid_within_10pct_best_rmse']}`")
    lines.append(f"- Fração da grade dentro de 20% do melhor RMSE: `{summary['frac_grid_within_20pct_best_rmse']}`")
    lines.append(f"- Fração da grade dentro de 35% do melhor RMSE: `{summary['frac_grid_within_35pct_best_rmse']}`\n")

    lines.append("## 6. Status\n")
    lines.append(f"```text\n{summary['final_status']}\n```\n")

    lines.append("## 7. Interpretação\n")
    if summary["synthetic_mode"]:
        lines.append(
            "O script rodou em modo sintético porque nenhum arquivo SPARC/RAR real foi encontrado. "
            "Esse modo serve apenas para testar o funcionamento do código e não deve ser usado como evidência científica."
        )
    else:
        lines.append(
            "O resultado mede a sensibilidade da TVGD a variações simples de M*/L. "
            "Se uma região ampla da grade permanece próxima do melhor RMSE, a relação é robusta. "
            "Se apenas um ponto muito fino da grade funciona, o setor é sensível a ajuste fino."
        )

    lines.append("\n## 8. Arquivos gerados\n")
    lines.append("- `ml_robustness_grid.csv`")
    lines.append("- `ml_robustness_summary.json`")
    lines.append("- `ml_robustness_report.md`")
    lines.append("- `ml_robustness_heatmap.png`, se matplotlib estiver disponível\n")

    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def make_plot(outdir: Path, grid: pd.DataFrame) -> Optional[Path]:
    if not HAS_MATPLOTLIB:
        return None

    plot_path = outdir / "ml_robustness_heatmap.png"

    if grid["mode"].iloc[0] == "component_ml_grid":
        pivot = grid.pivot_table(
            index="upsilon_bulge",
            columns="upsilon_disk",
            values="rmse_log10",
            aggfunc="mean",
        )

        plt.figure(figsize=(8, 6))
        plt.imshow(
            pivot.values,
            origin="lower",
            aspect="auto",
            extent=[
                float(pivot.columns.min()),
                float(pivot.columns.max()),
                float(pivot.index.min()),
                float(pivot.index.max()),
            ],
        )
        plt.colorbar(label="RMSE log10")
        plt.xlabel("Upsilon disk")
        plt.ylabel("Upsilon bulge")
        plt.title("TVGD SPARC M/L robustness")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=160)
        plt.close()

    else:
        plt.figure(figsize=(8, 5))
        plt.plot(grid["global_ml_scale"], grid["rmse_log10"], marker="o")
        plt.xlabel("Global M/L scale applied to g_bar")
        plt.ylabel("RMSE log10")
        plt.title("TVGD global g_bar scale robustness")
        plt.tight_layout()
        plt.savefig(plot_path, dpi=160)
        plt.close()

    return plot_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, default=".", help="Raiz do repositório.")
    parser.add_argument("--input", type=str, default=None, help="CSV SPARC/RAR já processado.")
    parser.add_argument("--a0", type=float, default=DEFAULT_A0, help="Escala a0 na mesma unidade dos dados.")
    parser.add_argument("--outdir", type=str, default=None, help="Pasta de saída.")
    parser.add_argument("--no-plot", action="store_true", help="Não gerar figura.")
    args = parser.parse_args()

    root = Path(args.root).resolve()
    input_path = locate_input(root, args.input)

    if args.outdir:
        outdir = Path(args.outdir)
        if not outdir.is_absolute():
            outdir = root / outdir
    else:
        # Compatível com repositórios em português e inglês.
        if (root / "Resultados").exists():
            outdir = root / "Resultados" / "01B_sparc_ml_robustness"
        else:
            outdir = root / "results" / "01B_sparc_ml_robustness"

    outdir.mkdir(parents=True, exist_ok=True)

    df, mapping, synthetic_mode = load_and_standardize(input_path, a0=args.a0)
    df = clean_data(df)

    if len(df) < 5:
        raise ValueError("Poucos pontos válidos após limpeza. Verifique o arquivo de entrada.")

    has_components = all(c in df.columns for c in ["g_gas", "g_disk", "g_bulge"])

    if has_components:
        grid = run_component_grid(
            df,
            a0=args.a0,
            disk_grid=DEFAULT_DISK_GRID,
            bulge_grid=DEFAULT_BULGE_GRID,
        )
    else:
        grid = run_global_grid(
            df,
            a0=args.a0,
            global_grid=DEFAULT_GLOBAL_ML_GRID,
        )

    summary = summarize_grid(grid, synthetic_mode=synthetic_mode, input_path=input_path, a0=args.a0)

    grid_path = outdir / "ml_robustness_grid.csv"
    json_path = outdir / "ml_robustness_summary.json"

    grid.to_csv(grid_path, index=False)

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    report_path = write_report(outdir, summary, grid, mapping)

    plot_path = None
    if not args.no_plot:
        plot_path = make_plot(outdir, grid)

    print("=" * 72)
    print("TVGD — 01B SPARC M/L ROBUSTNESS")
    print("=" * 72)
    print(f"Raiz do repositório: {root}")
    print(f"Arquivo de entrada: {input_path}")
    print(f"Modo sintético: {synthetic_mode}")
    print(f"Pontos válidos: {summary['n_points']}")
    print(f"Grade avaliada: {summary['n_grid']}")
    print(f"Status final: {summary['final_status']}")
    print()
    print("Melhor ponto:")
    for k, v in summary["best"].items():
        if k in [
            "mode",
            "upsilon_disk",
            "upsilon_bulge",
            "global_ml_scale",
            "rmse_log10",
            "mae_log10",
            "bias_log10",
            "chi2red_proxy_sigma_0p10dex",
        ]:
            print(f"  {k}: {v}")

    print()
    print("Arquivos gerados:")
    print(f" - {grid_path}")
    print(f" - {json_path}")
    print(f" - {report_path}")
    if plot_path is not None:
        print(f" - {plot_path}")
    elif args.no_plot:
        print(" - figura desativada por --no-plot")
    else:
        print(" - figura não gerada porque matplotlib não está disponível")

    if synthetic_mode:
        print()
        print("AVISO: nenhum arquivo SPARC/RAR real foi encontrado.")
        print("O script rodou em modo sintético apenas para validar execução.")
        print("Para resultado científico, forneça um CSV real com --input.")


if __name__ == "__main__":
    main()
