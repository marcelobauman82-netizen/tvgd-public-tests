# ============================================================
# TESTE 1 - GERAR ARQUIVOS DEMO DE MEMORIA DE FASE
#
# Objetivo:
#   Criar os 5 arquivos .s2p necessarios para validar
#   a coerencia estatica sem rotacao.
#
# Este codigo gera:
#
#   MEM_00_REFERENCIA.s2p
#   MEM_01_PERT_FRACA.s2p
#   MEM_02_PERT_MEDIA.s2p
#   MEM_03_PERT_FORTE.s2p
#   MEM_04_CONTROLE_CARGA.s2p
#
# IMPORTANTE:
#   Estes arquivos sao DEMO/SIMULADOS.
#   Eles servem apenas para testar o pipeline.
#   Nao sao dados reais de VNA.
# ============================================================

from pathlib import Path
import numpy as np


# ============================================================
# CONFIGURACAO
# ============================================================

DATA_DIR = Path("/content")

N_FREQ = 201
FREQ_MIN_HZ = 1.0e9
FREQ_MAX_HZ = 2.0e9

SEED = 20260614
rng = np.random.default_rng(SEED)


# ============================================================
# FUNCOES AUXILIARES
# ============================================================

def db20(z):
    """
    Converte magnitude complexa para dB.
    """
    return 20.0 * np.log10(np.abs(z) + 1e-15)


def phase_deg(z):
    """
    Retorna fase em graus.
    """
    return np.rad2deg(np.angle(z))


def write_s2p(path, freq_hz, S11, S21, S12=None, S22=None):
    """
    Escreve arquivo Touchstone .s2p no formato:

        # Hz S DB R 50

    Colunas:
        freq
        S11_dB S11_deg
        S21_dB S21_deg
        S12_dB S12_deg
        S22_dB S22_deg
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if S12 is None:
        S12 = S21.copy()

    if S22 is None:
        S22 = S11.copy()

    with open(path, "w", encoding="utf-8") as f:
        f.write("# Hz S DB R 50\n")
        f.write("! Arquivo DEMO/SIMULADO - nao e dado real de VNA\n")
        f.write("! freq S11_dB S11_deg S21_dB S21_deg S12_dB S12_deg S22_dB S22_deg\n")

        for i in range(len(freq_hz)):
            f.write(
                f"{freq_hz[i]:.8e} "
                f"{db20(S11[i]):.8f} {phase_deg(S11[i]):.8f} "
                f"{db20(S21[i]):.8f} {phase_deg(S21[i]):.8f} "
                f"{db20(S12[i]):.8f} {phase_deg(S12[i]):.8f} "
                f"{db20(S22[i]):.8f} {phase_deg(S22[i]):.8f}\n"
            )


def make_demo_s2p(filename, sigma_phase, control_random=False):
    """
    Cria arquivo .s2p sintetico.

    sigma_phase controla o embaralhamento da fase:

        sigma pequeno  -> coerencia alta
        sigma grande   -> coerencia baixa

    control_random=True:
        cria fase quase aleatoria, simulando carga/controle nulo.
    """
    freq_hz = np.linspace(FREQ_MIN_HZ, FREQ_MAX_HZ, N_FREQ)

    x = (freq_hz - freq_hz.mean()) / (freq_hz.max() - freq_hz.min())

    # Fase comum de cabo
    cable_phase = 2.0 * np.pi * freq_hz * 1.2e-9

    # Pequena ondulacao instrumental comum
    ripple = 0.015 * np.sin(4.0 * np.pi * x + 0.2)

    if control_random:
        phase_noise = rng.uniform(-np.pi, np.pi, size=N_FREQ)
    else:
        phase_noise = rng.normal(0.0, sigma_phase, size=N_FREQ)

    phase = cable_phase + ripple + phase_noise

    # S21 principal
    S21_amp = 10.0 ** (-3.0 / 20.0)
    S21 = S21_amp * np.exp(1j * phase)

    # S11 pequeno, apenas para manter o arquivo completo
    S11_amp = 10.0 ** (-18.0 / 20.0)
    S11_phase = 0.1 * np.sin(2.0 * np.pi * x)
    S11 = S11_amp * np.exp(1j * S11_phase)

    write_s2p(DATA_DIR / filename, freq_hz, S11, S21)


# ============================================================
# GERAR OS ARQUIVOS MEM
# ============================================================

def gerar_arquivos_demo_memoria():
    print("\n===================================================")
    print("TESTE 1 - GERAR ARQUIVOS DEMO DE MEMORIA")
    print("===================================================")

    print("\n[INFO] Pasta de saida:", DATA_DIR)
    print("[INFO] Estes arquivos sao SIMULADOS, nao reais.\n")

    make_demo_s2p(
        filename="MEM_00_REFERENCIA.s2p",
        sigma_phase=0.05,
        control_random=False,
    )

    make_demo_s2p(
        filename="MEM_01_PERT_FRACA.s2p",
        sigma_phase=0.30,
        control_random=False,
    )

    make_demo_s2p(
        filename="MEM_02_PERT_MEDIA.s2p",
        sigma_phase=0.75,
        control_random=False,
    )

    make_demo_s2p(
        filename="MEM_03_PERT_FORTE.s2p",
        sigma_phase=1.45,
        control_random=False,
    )

    make_demo_s2p(
        filename="MEM_04_CONTROLE_CARGA.s2p",
        sigma_phase=3.0,
        control_random=True,
    )

    arquivos = [
        "MEM_00_REFERENCIA.s2p",
        "MEM_01_PERT_FRACA.s2p",
        "MEM_02_PERT_MEDIA.s2p",
        "MEM_03_PERT_FORTE.s2p",
        "MEM_04_CONTROLE_CARGA.s2p",
    ]

    print("[OK] Arquivos gerados:\n")

    for nome in arquivos:
        path = DATA_DIR / nome
        print(" -", path)

    print("\nSTATUS_FINAL: ARQUIVOS_DEMO_MEMORIA_GERADOS")


# ============================================================
# EXECUTAR
# ============================================================

if __name__ == "__main__":
    gerar_arquivos_demo_memoria()
