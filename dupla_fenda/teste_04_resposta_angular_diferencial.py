# ============================================================
# TESTE 4 - RESPOSTA ANGULAR DIFERENCIAL
#
# Objetivo:
#   Verificar se o protocolo angular esta montado corretamente.
#
# Este teste usa:
#   - DUT_NORMAL
#   - DUT_INVERTED
#   - LOAD_NORMAL
#   - LOAD_INVERTED
#   - CW
#   - CCW
#
# Observavel diferencial:
#
#   R_DUT_diff(chi) =
#       1/2 [R_DUT_normal(chi) - R_DUT_inverted(chi)]
#
#   R_LOAD_diff(chi) =
#       1/2 [R_LOAD_normal(chi) - R_LOAD_inverted(chi)]
#
# Modelo angular:
#
#   R_diff(chi) = A cos(chi) + B sin(chi) + C
#
# IMPORTANTE:
#   Este codigo gera dados DEMO/SIMULADOS.
#   Ele serve para validar a logica do protocolo.
#   Nao e dado real de VNA.
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# CONFIGURACAO
# ============================================================

DATA_DIR = Path("/content")

N_FREQ = 201
FREQ_MIN_HZ = 1.0e9
FREQ_MAX_HZ = 2.0e9

ANGLES_DEG = [0, 45, 90, 135, 180, 225, 270, 315]
DIRECTIONS = ["CW", "CCW"]
DEVICES = ["DUT", "LOAD"]
ORIENTATIONS = ["NORMAL", "INVERTED"]

SEED = 20260614
rng = np.random.default_rng(SEED)

# Criterios estatisticos
AIC_LIMIT = -10.0
P_PERM_LIMIT = 0.05
SIGMA_LIMIT = 3.0
B_LEAK_LIMIT = 0.15


# ============================================================
# FUNCOES PARA ARQUIVO .s2p
# ============================================================

def db20(z):
    return 20.0 * np.log10(np.abs(z) + 1e-15)


def phase_deg(z):
    return np.rad2deg(np.angle(z))


def write_s2p(path, freq_hz, S11, S21, S12=None, S22=None):
    """
    Escreve arquivo Touchstone .s2p simples.
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


def read_s2p_simple(path):
    """
    Le o .s2p gerado por este proprio codigo.
    Retorna:
        freq_hz, S11_complexo, S21_complexo
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Arquivo nao encontrado: {path}")

    rows = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("#") or line.startswith("!"):
                continue

            parts = line.split()

            if len(parts) < 9:
                continue

            rows.append([float(x) for x in parts[:9]])

    if not rows:
        raise ValueError(f"Nenhuma linha numerica encontrada em {path}")

    arr = np.array(rows, dtype=float)

    freq_hz = arr[:, 0]

    S11_db = arr[:, 1]
    S11_deg = arr[:, 2]

    S21_db = arr[:, 3]
    S21_deg = arr[:, 4]

    S11 = 10.0 ** (S11_db / 20.0) * np.exp(1j * np.deg2rad(S11_deg))
    S21 = 10.0 ** (S21_db / 20.0) * np.exp(1j * np.deg2rad(S21_deg))

    return freq_hz, S11, S21


# ============================================================
# NOMENCLATURA DOS ARQUIVOS ANGULARES
# ============================================================

def angular_filename(device, orientation, direction, angle_deg):
    return f"{device}_{orientation}_{direction}_{angle_deg:03d}.s2p"


# ============================================================
# GERADOR DEMO ANGULAR
# ============================================================

def make_angular_demo_s2p(device, orientation, direction, angle_deg, A_true=0.030):
    """
    Gera um arquivo angular sintetico.

    Regras:
      DUT_NORMAL   recebe +A cos(chi)
      DUT_INVERTED recebe -A cos(chi)

      LOAD_NORMAL e LOAD_INVERTED nao recebem sinal fisico.

    Artefatos comuns de bancada sao adicionados,
    mas devem cancelar no diferencial Normal/Invertido.
    """
    freq_hz = np.linspace(FREQ_MIN_HZ, FREQ_MAX_HZ, N_FREQ)

    x = (freq_hz - freq_hz.mean()) / (freq_hz.max() - freq_hz.min())
    chi = np.deg2rad(angle_deg)

    # Fase de cabo comum
    cable_phase = 2.0 * np.pi * freq_hz * 1.1e-9

    # Artefatos comuns de bancada
    bench_common = (
        0.08 * np.sin(2.0 * chi + 0.3)
        + 0.04 * np.sin(3.0 * chi - 0.5)
    )

    # Pequena histerese comum por sentido de giro
    direction_sign = +1.0 if direction == "CW" else -1.0
    hysteresis = direction_sign * 0.002 * np.sin(chi + 0.2)

    # Sinal diferencial esperado
    if device == "DUT":
        orientation_sign = +1.0 if orientation == "NORMAL" else -1.0
        angular_signal = orientation_sign * A_true * np.cos(chi)
    else:
        angular_signal = 0.0

    # Ruido pequeno
    noise = rng.normal(0.0, 0.002, size=N_FREQ)

    phase = cable_phase + bench_common + hysteresis + angular_signal + noise

    # S21
    S21_amp = 10.0 ** (-3.0 / 20.0)
    S21 = S21_amp * np.exp(1j * phase)

    # S11 pequeno
    S11_amp = 10.0 ** (-18.0 / 20.0)
    S11_phase = 0.1 * np.sin(2.0 * np.pi * x)
    S11 = S11_amp * np.exp(1j * S11_phase)

    filename = angular_filename(device, orientation, direction, angle_deg)
    write_s2p(DATA_DIR / filename, freq_hz, S11, S21)


def gerar_arquivos_demo_angulares():
    print("\n[GERANDO] Arquivos DEMO angulares...")

    for angle in ANGLES_DEG:
        for direction in DIRECTIONS:
            for device in DEVICES:
                for orientation in ORIENTATIONS:
                    make_angular_demo_s2p(
                        device=device,
                        orientation=orientation,
                        direction=direction,
                        angle_deg=angle,
                        A_true=0.030,
                    )

    print("[OK] Arquivos DEMO angulares gerados.")


# ============================================================
# ANALISE DE FASE
# ============================================================

def remove_linear_phase(freq_hz, phase):
    """
    Remove apenas a inclinacao linear da fase.
    Mantem o offset medio angular.
    """
    x = freq_hz - np.mean(freq_hz)
    X = x.reshape(-1, 1)

    beta, *_ = np.linalg.lstsq(X, phase, rcond=None)
    trend = X @ beta

    residual = phase - trend

    return residual


def phase_feature_from_file(path):
    """
    Extrai uma grandeza escalar R(chi) de cada arquivo.
    """
    freq_hz, S11, S21 = read_s2p_simple(path)

    phase = np.unwrap(np.angle(S21))

    residual = remove_linear_phase(freq_hz, phase)

    R = float(np.mean(residual))

    return R


# ============================================================
# AJUSTE ANGULAR
# ============================================================

def aic_bic(y, yhat, k):
    n = len(y)

    rss = float(np.sum((y - yhat) ** 2))
    rss = max(rss, 1e-30)

    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)

    return aic, bic, rss


def fit_angular(theta, y, n_perm=2000):
    """
    Ajusta:

        y = A cos(theta) + B sin(theta) + C

    e compara com modelo nulo:

        y = C
    """
    theta = np.array(theta, dtype=float)
    y = np.array(y, dtype=float)

    # Modelo nulo
    X0 = np.ones((len(theta), 1))
    beta0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    y0 = X0 @ beta0

    aic0, bic0, rss0 = aic_bic(y, y0, k=1)

    # Modelo angular
    X = np.column_stack([
        np.cos(theta),
        np.sin(theta),
        np.ones_like(theta),
    ])

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta

    aic, bic, rss = aic_bic(y, yhat, k=3)

    A, B, C = beta

    # Incertezas
    n, k = X.shape
    dof = max(n - k, 1)
    sigma2 = rss / dof

    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
        err = np.sqrt(np.diag(cov))
    except Exception:
        err = np.array([np.nan, np.nan, np.nan])

    dA, dB, dC = err

    # Teste de permutacao
    rng_perm = np.random.default_rng(123)
    A_perm = []

    for _ in range(n_perm):
        y_perm = rng_perm.permutation(y)
        beta_perm, *_ = np.linalg.lstsq(X, y_perm, rcond=None)
        A_perm.append(abs(beta_perm[0]))

    A_perm = np.array(A_perm)
    p_perm = float(np.mean(A_perm >= abs(A)))

    return {
        "A": float(A),
        "B": float(B),
        "C": float(C),
        "dA": float(dA),
        "dB": float(dB),
        "dC": float(dC),
        "delta_AIC": float(aic - aic0),
        "delta_BIC": float(bic - bic0),
        "p_perm": p_perm,
        "rss": float(rss),
        "rss_null": float(rss0),
        "leak_B_over_A": float(abs(B) / (abs(A) + 1e-30)),
    }


# ============================================================
# TESTE 4
# ============================================================

def run_teste_4_resposta_angular():
    print("\n===================================================")
    print("TESTE 4 - RESPOSTA ANGULAR DIFERENCIAL")
    print("===================================================")

    # 1. Gerar arquivos DEMO
    gerar_arquivos_demo_angulares()

    # 2. Verificar arquivos
    missing = []

    for angle in ANGLES_DEG:
        for direction in DIRECTIONS:
            for device in DEVICES:
                for orientation in ORIENTATIONS:
                    filename = angular_filename(device, orientation, direction, angle)
                    path = DATA_DIR / filename

                    if not path.exists():
                        missing.append(path)

    if missing:
        print("\n[BLOQUEADO] Arquivos faltantes:")

        for path in missing:
            print(" -", path)

        print("\nSTATUS_FINAL: BLOQUEADO_ARQUIVOS_ANGULARES_FALTANTES")

        return {
            "status": "BLOQUEADO_ARQUIVOS_ANGULARES_FALTANTES",
            "missing": missing,
        }

    # 3. Extrair R de cada arquivo
    rows = []

    for angle in ANGLES_DEG:
        theta = np.deg2rad(angle)

        for direction in DIRECTIONS:
            for device in DEVICES:
                for orientation in ORIENTATIONS:
                    filename = angular_filename(device, orientation, direction, angle)
                    path = DATA_DIR / filename

                    R = phase_feature_from_file(path)

                    rows.append({
                        "angle_deg": angle,
                        "theta_rad": theta,
                        "direction": direction,
                        "device": device,
                        "orientation": orientation,
                        "R": R,
                        "file": filename,
                    })

    raw = pd.DataFrame(rows)

    # 4. Montar diferencial Normal/Invertido
    diff_rows = []

    for (angle, theta, direction, device), group in raw.groupby(
        ["angle_deg", "theta_rad", "direction", "device"]
    ):
        R_normal = group[group["orientation"] == "NORMAL"]["R"].mean()
        R_inverted = group[group["orientation"] == "INVERTED"]["R"].mean()

        R_diff = 0.5 * (R_normal - R_inverted)
        R_common = 0.5 * (R_normal + R_inverted)

        diff_rows.append({
            "angle_deg": angle,
            "theta_rad": theta,
            "direction": direction,
            "device": device,
            "R_normal": R_normal,
            "R_inverted": R_inverted,
            "R_diff": R_diff,
            "R_common": R_common,
        })

    diff = pd.DataFrame(diff_rows)

    # 5. Ajustar por device/direction e tambem ALL
    report_rows = []

    # Separado CW e CCW
    for (device, direction), group in diff.groupby(["device", "direction"]):
        fit = fit_angular(group["theta_rad"].values, group["R_diff"].values)

        A = fit["A"]
        B = fit["B"]
        dA = fit["dA"]
        dB = fit["dB"]

        A_sig = abs(A) > SIGMA_LIMIT * dA
        B_ok = abs(B) < max(SIGMA_LIMIT * dB, B_LEAK_LIMIT * abs(A))
        AIC_ok = fit["delta_AIC"] < AIC_LIMIT
        perm_ok = fit["p_perm"] < P_PERM_LIMIT

        if device == "LOAD":
            if A_sig and AIC_ok and perm_ok:
                status = "FALSO_POSITIVO_INSTRUMENTAL"
            else:
                status = "LOAD_PASSOU"
        else:
            if A_sig and B_ok and AIC_ok and perm_ok:
                status = "CANDIDATO_ANGULAR_DIFERENCIAL"
            else:
                status = "SEM_CANDIDATO_ANGULAR_FORTE"

        report_rows.append({
            "device": device,
            "direction": direction,
            **fit,
            "A_sig": A_sig,
            "B_ok": B_ok,
            "AIC_ok": AIC_ok,
            "perm_ok": perm_ok,
            "status": status,
        })

    # Agregado ALL
    for device, group in diff.groupby("device"):
        fit = fit_angular(group["theta_rad"].values, group["R_diff"].values)

        A = fit["A"]
        B = fit["B"]
        dA = fit["dA"]
        dB = fit["dB"]

        A_sig = abs(A) > SIGMA_LIMIT * dA
        B_ok = abs(B) < max(SIGMA_LIMIT * dB, B_LEAK_LIMIT * abs(A))
        AIC_ok = fit["delta_AIC"] < AIC_LIMIT
        perm_ok = fit["p_perm"] < P_PERM_LIMIT

        if device == "LOAD":
            if A_sig and AIC_ok and perm_ok:
                status = "FALSO_POSITIVO_INSTRUMENTAL"
            else:
                status = "LOAD_PASSOU"
        else:
            if A_sig and B_ok and AIC_ok and perm_ok:
                status = "CANDIDATO_ANGULAR_DIFERENCIAL"
            else:
                status = "SEM_CANDIDATO_ANGULAR_FORTE"

        report_rows.append({
            "device": device,
            "direction": "ALL",
            **fit,
            "A_sig": A_sig,
            "B_ok": B_ok,
            "AIC_ok": AIC_ok,
            "perm_ok": perm_ok,
            "status": status,
        })

    report = pd.DataFrame(report_rows)

    # 6. Checar CW/CCW
    direction_rows = []

    for device in DEVICES:
        cw = report[(report["device"] == device) & (report["direction"] == "CW")]
        ccw = report[(report["device"] == device) & (report["direction"] == "CCW")]

        if len(cw) > 0 and len(ccw) > 0:
            A_CW = float(cw["A"].iloc[0])
            A_CCW = float(ccw["A"].iloc[0])
            dA_CW = float(cw["dA"].iloc[0])
            dA_CCW = float(ccw["dA"].iloc[0])

            sigma = np.sqrt(dA_CW**2 + dA_CCW**2)

            if sigma > 0:
                z_diff = (A_CW - A_CCW) / sigma
            else:
                z_diff = np.nan

            direction_consistent = abs(z_diff) < 3

            direction_rows.append({
                "device": device,
                "A_CW": A_CW,
                "A_CCW": A_CCW,
                "z_diff": z_diff,
                "direction_consistent": direction_consistent,
            })

    direction_check = pd.DataFrame(direction_rows)

    # 7. Veredito global
    dut_all = report[(report["device"] == "DUT") & (report["direction"] == "ALL")]
    load_all = report[(report["device"] == "LOAD") & (report["direction"] == "ALL")]

    if len(dut_all) == 0 or len(load_all) == 0:
        global_status = "BLOQUEADO_RELATORIO_INCOMPLETO"
    else:
        dut_status = dut_all["status"].iloc[0]
        load_status = load_all["status"].iloc[0]

        if len(direction_check) > 0:
            direction_ok = bool(direction_check["direction_consistent"].all())
        else:
            direction_ok = False

        if load_status != "LOAD_PASSOU":
            global_status = "FALSO_POSITIVO_INSTRUMENTAL"
        elif dut_status == "CANDIDATO_ANGULAR_DIFERENCIAL" and direction_ok:
            global_status = "PASSOU_FORTE_RESPOSTA_ANGULAR_DIFERENCIAL"
        elif not direction_ok:
            global_status = "FALHOU_CW_CCW_INCONSISTENTE"
        else:
            global_status = "SEM_CANDIDATO_ANGULAR_FORTE"

    # 8. Imprimir resultados
    print("\nRelatorio angular:")
    print(report[[
        "device",
        "direction",
        "A",
        "dA",
        "B",
        "dB",
        "C",
        "delta_AIC",
        "delta_BIC",
        "p_perm",
        "leak_B_over_A",
        "status",
    ]])

    print("\nCheck CW/CCW:")
    print(direction_check)

    print("\nSTATUS_FINAL:", global_status)

    # 9. Salvar CSVs
    raw.to_csv(DATA_DIR / "teste_04_raw_angular.csv", index=False)
    diff.to_csv(DATA_DIR / "teste_04_diferencial_angular.csv", index=False)
    report.to_csv(DATA_DIR / "teste_04_relatorio_angular.csv", index=False)
    direction_check.to_csv(DATA_DIR / "teste_04_check_cw_ccw.csv", index=False)

    return {
        "status": global_status,
        "raw": raw,
        "diff": diff,
        "report": report,
        "direction_check": direction_check,
    }


# ============================================================
# EXECUTAR
# ============================================================

if __name__ == "__main__":
    result = run_teste_4_resposta_angular()
