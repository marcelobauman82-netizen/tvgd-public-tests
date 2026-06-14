# ============================================================
# TESTE 1 - VALIDACAO DE COERENCIA ESTATICA SEM ROTACAO
#
# Objetivo:
#   Verificar se a bancada mede coerencia/memoria de fase
#   de forma ordenada antes de qualquer teste angular.
#
# Arquivos esperados na pasta /content:
#
#   MEM_00_REFERENCIA.s2p
#   MEM_01_PERT_FRACA.s2p
#   MEM_02_PERT_MEDIA.s2p
#   MEM_03_PERT_FORTE.s2p
#   MEM_04_CONTROLE_CARGA.s2p
#
# Criterio:
#   C_corr(fraca) > C_corr(media) > C_corr(forte)
#   C_corr(controle) aproximadamente 0
#
# Este codigo pode ser usado com dados reais .s2p
# ou com arquivos simulados gerados separadamente.
# ============================================================

import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import skrf as rf
except Exception:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "scikit-rf", "-q"])
    import skrf as rf


# ============================================================
# CONFIGURACAO
# ============================================================

DATA_DIR = Path("/content")

CONTROL_MAX = 0.15
MIN_FREQ_POINTS = 101

FILES = {
    "referencia": "MEM_00_REFERENCIA.s2p",
    "fraca": "MEM_01_PERT_FRACA.s2p",
    "media": "MEM_02_PERT_MEDIA.s2p",
    "forte": "MEM_03_PERT_FORTE.s2p",
    "controle": "MEM_04_CONTROLE_CARGA.s2p",
}


# ============================================================
# FUNCOES
# ============================================================

def load_s2p(path):
    """
    Carrega arquivo .s2p usando scikit-rf.

    Retorna:
        freq_Hz, S11, S21
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Arquivo nao encontrado: {path}")

    ntwk = rf.Network(str(path))

    if ntwk.s.shape[1] < 2 or ntwk.s.shape[2] < 2:
        raise ValueError(f"Arquivo nao parece ser .s2p de 2 portas: {path}")

    freq = np.asarray(ntwk.f, dtype=float)
    S11 = ntwk.s[:, 0, 0]
    S21 = ntwk.s[:, 1, 0]

    return freq, S11, S21


def interp_complex(f_old, y_old, f_new):
    """
    Interpola numero complexo preservando parte real e imaginaria.
    """
    real = np.interp(f_new, f_old, np.real(y_old))
    imag = np.interp(f_new, f_old, np.imag(y_old))
    return real + 1j * imag


def remove_linear_phase(freq, phase):
    """
    Remove tendencia linear de fase associada a atraso de cabo.

    Mantem a estrutura residual de fase, que e usada para medir coerencia.
    """
    x = freq - np.mean(freq)
    X = x.reshape(-1, 1)

    beta, *_ = np.linalg.lstsq(X, phase, rcond=None)
    trend = X @ beta

    residual = phase - trend

    return residual


def coherence_from_phase(theta):
    """
    Calcula:

        C_corr = | < exp(i theta) > |

    Se as fases estao organizadas:
        C_corr ~ 1

    Se as fases estao embaralhadas:
        C_corr ~ 0
    """
    return float(np.abs(np.mean(np.exp(1j * theta))))


def check_missing_files():
    """
    Verifica se todos os arquivos obrigatorios existem.
    """
    missing = []

    for label, filename in FILES.items():
        path = DATA_DIR / filename

        if not path.exists():
            missing.append(path)

    return missing


# ============================================================
# TESTE 1
# ============================================================

def run_test_1_static_memory():
    print("\n===================================================")
    print("TESTE 1 - VALIDACAO DE COERENCIA ESTATICA")
    print("===================================================")

    missing = check_missing_files()

    if missing:
        print("\n[BLOQUEADO] Arquivos faltantes:")

        for path in missing:
            print(" -", path)

        print("\nSTATUS_FINAL: BLOQUEADO_ARQUIVOS_FALTANTES")

        return {
            "status": "BLOQUEADO_ARQUIVOS_FALTANTES",
            "missing": missing,
        }

    # Carregar todos os arquivos
    loaded = {}

    for condition, filename in FILES.items():
        path = DATA_DIR / filename
        freq, S11, S21 = load_s2p(path)

        if len(freq) < MIN_FREQ_POINTS:
            print(f"[AVISO] {filename} tem apenas {len(freq)} pontos de frequencia.")

        loaded[condition] = {
            "freq": freq,
            "S11": S11,
            "S21": S21,
        }

    # Criar grade de frequencia comum
    f_min = max(item["freq"][0] for item in loaded.values())
    f_max = min(item["freq"][-1] for item in loaded.values())
    n_common = min(len(item["freq"]) for item in loaded.values())

    if f_max <= f_min:
        print("\n[ERRO] Os arquivos nao possuem faixa de frequencia comum.")
        print("STATUS_FINAL: BLOQUEADO_FAIXA_FREQUENCIA_INCOMPATIVEL")

        return {
            "status": "BLOQUEADO_FAIXA_FREQUENCIA_INCOMPATIVEL",
        }

    freq_common = np.linspace(f_min, f_max, n_common)

    # Fase da referencia
    ref = loaded["referencia"]
    S21_ref = interp_complex(ref["freq"], ref["S21"], freq_common)
    phi_ref = np.unwrap(np.angle(S21_ref))

    rows = []

    for condition, item in loaded.items():
        S21 = interp_complex(item["freq"], item["S21"], freq_common)

        phi = np.unwrap(np.angle(S21))

        # Fase relativa a referencia
        dphi = np.unwrap(phi - phi_ref)

        # Remove fase linear de cabo
        residual = remove_linear_phase(freq_common, dphi)

        # Coerencia corrigida
        C_corr = coherence_from_phase(residual)

        rows.append({
            "condition": condition,
            "C_corr": C_corr,
            "N_freq": len(freq_common),
        })

    df = pd.DataFrame(rows)

    print("\nResultado por condicao:")
    print(df)

    values = dict(zip(df["condition"], df["C_corr"]))

    monotonic_ok = (
        values["fraca"] > values["media"] > values["forte"]
    )

    control_ok = (
        values["controle"] < CONTROL_MAX
    )

    print("\nVerificacao:")
    print("C_fraca  =", values["fraca"])
    print("C_media  =", values["media"])
    print("C_forte  =", values["forte"])
    print("C_controle =", values["controle"])

    print("\nmonotonic_ok:", monotonic_ok)
    print("control_ok:", control_ok)

    if monotonic_ok and control_ok:
        status = "MEMORIA_DE_FASE_COMPATIVEL_EM_BANCADA"
    else:
        status = "SEM_MEMORIA_DE_FASE_REPRODUTIVEL"

    print("\nSTATUS_FINAL:", status)

    # Salvar resultado
    df.to_csv(DATA_DIR / "resultado_teste_1_coerencia_estatica.csv", index=False)

    return {
        "status": status,
        "summary": df,
        "monotonic_ok": monotonic_ok,
        "control_ok": control_ok,
    }


# ============================================================
# EXECUTAR
# ============================================================

result = run_test_1_static_memory()
