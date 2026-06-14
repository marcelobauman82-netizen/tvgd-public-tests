# ============================================================
# TESTE 5 - PIPELINE COMPLETO DO PROTOCOLO
#
# Objetivo:
#   Executar o protocolo inteiro na ordem correta:
#
#   1) Gerar arquivos DEMO de memoria
#   2) Validar coerencia estatica
#   3) Testar robustez temporal
#   4) Liberar ou bloquear o angular
#   5) Rodar resposta angular diferencial
#   6) Emitir um veredito final unico
#
# IMPORTANTE:
#   Este codigo usa dados DEMO/SIMULADOS.
#   Ele valida a logica do protocolo.
#   Nao e evidencia experimental real.
# ============================================================

from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# CONFIGURACAO GERAL
# ============================================================

DATA_DIR = Path("/content")

N_FREQ = 201
FREQ_MIN_HZ = 1.0e9
FREQ_MAX_HZ = 2.0e9

SEED = 20260614
rng = np.random.default_rng(SEED)

CONTROL_MAX = 0.15

ANGLES_DEG = [0, 45, 90, 135, 180, 225, 270, 315]
DIRECTIONS = ["CW", "CCW"]
DEVICES = ["DUT", "LOAD"]
ORIENTATIONS = ["NORMAL", "INVERTED"]

AIC_LIMIT = -10.0
P_PERM_LIMIT = 0.05
SIGMA_LIMIT = 3.0
B_LEAK_LIMIT = 0.15


# ============================================================
# FUNCOES GERAIS .s2p
# ============================================================

def db20(z):
    return 20.0 * np.log10(np.abs(z) + 1e-15)


def phase_deg(z):
    return np.rad2deg(np.angle(z))


def write_s2p(path, freq_hz, S11, S21, S12=None, S22=None):
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


def remove_linear_phase(freq_hz, phase):
    x = freq_hz - np.mean(freq_hz)
    X = x.reshape(-1, 1)

    beta, *_ = np.linalg.lstsq(X, phase, rcond=None)
    trend = X @ beta

    residual = phase - trend

    return residual


def coherence_from_phase(theta):
    return float(np.abs(np.mean(np.exp(1j * theta))))


# ============================================================
# PARTE 1 - GERAR ARQUIVOS DE MEMORIA
# ============================================================

def make_memory_s2p(filename, sigma_phase, seed, control_random=False):
    local_rng = np.random.default_rng(seed)

    freq_hz = np.linspace(FREQ_MIN_HZ, FREQ_MAX_HZ, N_FREQ)

    x = (freq_hz - freq_hz.mean()) / (freq_hz.max() - freq_hz.min())

    cable_phase = 2.0 * np.pi * freq_hz * 1.2e-9
    ripple = 0.015 * np.sin(4.0 * np.pi * x + 0.2)

    if control_random:
        phase_noise = local_rng.uniform(-np.pi, np.pi, size=N_FREQ)
    else:
        phase_noise = local_rng.normal(0.0, sigma_phase, size=N_FREQ)

    phase = cable_phase + ripple + phase_noise

    S21_amp = 10.0 ** (-3.0 / 20.0)
    S21 = S21_amp * np.exp(1j * phase)

    S11_amp = 10.0 ** (-18.0 / 20.0)
    S11_phase = 0.1 * np.sin(2.0 * np.pi * x)
    S11 = S11_amp * np.exp(1j * S11_phase)

    write_s2p(DATA_DIR / filename, freq_hz, S11, S21)


def gerar_arquivos_memoria(prefix="", seed_base=1000):
    make_memory_s2p(
        filename=prefix + "MEM_00_REFERENCIA.s2p",
        sigma_phase=0.05,
        seed=seed_base + 0,
        control_random=False,
    )

    make_memory_s2p(
        filename=prefix + "MEM_01_PERT_FRACA.s2p",
        sigma_phase=0.30,
        seed=seed_base + 1,
        control_random=False,
    )

    make_memory_s2p(
        filename=prefix + "MEM_02_PERT_MEDIA.s2p",
        sigma_phase=0.75,
        seed=seed_base + 2,
        control_random=False,
    )

    make_memory_s2p(
        filename=prefix + "MEM_03_PERT_FORTE.s2p",
        sigma_phase=1.45,
        seed=seed_base + 3,
        control_random=False,
    )

    make_memory_s2p(
        filename=prefix + "MEM_04_CONTROLE_CARGA.s2p",
        sigma_phase=3.0,
        seed=seed_base + 4,
        control_random=True,
    )


# ============================================================
# PARTE 2 - VALIDAR COERENCIA ESTATICA
# ============================================================

MEM_FILES = {
    "referencia": "MEM_00_REFERENCIA.s2p",
    "fraca": "MEM_01_PERT_FRACA.s2p",
    "media": "MEM_02_PERT_MEDIA.s2p",
    "forte": "MEM_03_PERT_FORTE.s2p",
    "controle": "MEM_04_CONTROLE_CARGA.s2p",
}


def validar_coerencia_estatica(file_map=MEM_FILES, label="STATIC"):
    missing = []

    for condition, filename in file_map.items():
        path = DATA_DIR / filename

        if not path.exists():
            missing.append(path)

    if missing:
        return {
            "status": "BLOQUEADO_ARQUIVOS_FALTANTES",
            "missing": missing,
            "summary": pd.DataFrame(),
        }

    loaded = {}

    for condition, filename in file_map.items():
        freq, S11, S21 = read_s2p_simple(DATA_DIR / filename)
        loaded[condition] = {
            "freq": freq,
            "S11": S11,
            "S21": S21,
        }

    f_min = max(item["freq"][0] for item in loaded.values())
    f_max = min(item["freq"][-1] for item in loaded.values())
    n_common = min(len(item["freq"]) for item in loaded.values())

    if f_max <= f_min:
        return {
            "status": "BLOQUEADO_FAIXA_FREQUENCIA_INCOMPATIVEL",
            "summary": pd.DataFrame(),
        }

    freq_common = np.linspace(f_min, f_max, n_common)

    ref = loaded["referencia"]
    S21_ref = np.interp(freq_common, ref["freq"], np.real(ref["S21"])) + 1j * np.interp(
        freq_common, ref["freq"], np.imag(ref["S21"])
    )
    phi_ref = np.unwrap(np.angle(S21_ref))

    rows = []

    for condition, item in loaded.items():
        S21 = np.interp(freq_common, item["freq"], np.real(item["S21"])) + 1j * np.interp(
            freq_common, item["freq"], np.imag(item["S21"])
        )

        phi = np.unwrap(np.angle(S21))
        dphi = np.unwrap(phi - phi_ref)

        residual = remove_linear_phase(freq_common, dphi)

        C_corr = coherence_from_phase(residual)

        rows.append({
            "condition": condition,
            "C_corr": C_corr,
            "N_freq": len(freq_common),
        })

    df = pd.DataFrame(rows)

    values = dict(zip(df["condition"], df["C_corr"]))

    monotonic_ok = values["fraca"] > values["media"] > values["forte"]
    control_ok = values["controle"] < CONTROL_MAX

    if monotonic_ok and control_ok:
        status = "MEMORIA_DE_FASE_COMPATIVEL_EM_BANCADA"
    else:
        status = "SEM_MEMORIA_DE_FASE_REPRODUTIVEL"

    df.to_csv(DATA_DIR / f"{label}_coerencia_estatica.csv", index=False)

    return {
        "status": status,
        "summary": df,
        "monotonic_ok": monotonic_ok,
        "control_ok": control_ok,
    }


# ============================================================
# PARTE 3 - ROBUSTEZ TEMPORAL
# ============================================================

TEMP_MAP = {
    "T0": {
        "referencia": "T0_MEM_00_REFERENCIA.s2p",
        "fraca": "T0_MEM_01_PERT_FRACA.s2p",
        "media": "T0_MEM_02_PERT_MEDIA.s2p",
        "forte": "T0_MEM_03_PERT_FORTE.s2p",
        "controle": "T0_MEM_04_CONTROLE_CARGA.s2p",
    },
    "T1": {
        "referencia": "T1_MEM_00_REFERENCIA.s2p",
        "fraca": "T1_MEM_01_PERT_FRACA.s2p",
        "media": "T1_MEM_02_PERT_MEDIA.s2p",
        "forte": "T1_MEM_03_PERT_FORTE.s2p",
        "controle": "T1_MEM_04_CONTROLE_CARGA.s2p",
    },
    "T2": {
        "referencia": "T2_MEM_00_REFERENCIA.s2p",
        "fraca": "T2_MEM_01_PERT_FRACA.s2p",
        "media": "T2_MEM_02_PERT_MEDIA.s2p",
        "forte": "T2_MEM_03_PERT_FORTE.s2p",
        "controle": "T2_MEM_04_CONTROLE_CARGA.s2p",
    },
}


def gerar_arquivos_temporais():
    # tres tempos com sementes diferentes, mas mesma ordem de perturbacao
    gerar_arquivos_memoria(prefix="T0_", seed_base=2000)
    gerar_arquivos_memoria(prefix="T1_", seed_base=3000)
    gerar_arquivos_memoria(prefix="T2_", seed_base=4000)


def validar_robustez_temporal():
    time_results = []
    summaries = []

    for time_label, fmap in TEMP_MAP.items():
        result = validar_coerencia_estatica(
            file_map=fmap,
            label=f"{time_label}_temporal",
        )

        if result["status"] != "MEMORIA_DE_FASE_COMPATIVEL_EM_BANCADA":
            return {
                "status": "DRIFT_INSTRUMENTAL_DOMINANTE_OU_MEMORIA_NAO_REPRODUTIVEL",
                "time_results": time_results,
                "summary": pd.DataFrame(),
            }

        df = result["summary"].copy()
        df["time"] = time_label
        summaries.append(df)
        time_results.append(result)

    df_all = pd.concat(summaries, ignore_index=True)

    grouped = df_all.groupby("condition").agg(
        C_mean=("C_corr", "mean"),
        C_std=("C_corr", "std"),
        C_min=("C_corr", "min"),
        C_max=("C_corr", "max"),
        N_reps=("C_corr", "size"),
    ).reset_index()

    values = dict(zip(grouped["condition"], grouped["C_mean"]))

    monotonic_ok = values["fraca"] > values["media"] > values["forte"]
    control_ok = values["controle"] < CONTROL_MAX

    separation_1 = values["fraca"] - values["media"]
    separation_2 = values["media"] - values["forte"]
    separation_ok = separation_1 > 0.05 and separation_2 > 0.05

    stability_ok = grouped["C_std"].max() < 0.08

    if monotonic_ok and control_ok and separation_ok and stability_ok:
        status = "PASSOU_FORTE_MEMORIA_DE_FASE_REPRODUTIVEL_TEMPORALMENTE"
    else:
        status = "DRIFT_INSTRUMENTAL_DOMINANTE_OU_MEMORIA_NAO_REPRODUTIVEL"

    grouped.to_csv(DATA_DIR / "teste_05_resumo_temporal.csv", index=False)

    return {
        "status": status,
        "summary": grouped,
        "monotonic_ok": monotonic_ok,
        "control_ok": control_ok,
        "separation_ok": separation_ok,
        "stability_ok": stability_ok,
    }


# ============================================================
# PARTE 4 - RESPOSTA ANGULAR DIFERENCIAL
# ============================================================

def angular_filename(device, orientation, direction, angle_deg):
    return f"{device}_{orientation}_{direction}_{angle_deg:03d}.s2p"


def make_angular_demo_s2p(device, orientation, direction, angle_deg, A_true=0.030):
    freq_hz = np.linspace(FREQ_MIN_HZ, FREQ_MAX_HZ, N_FREQ)

    x = (freq_hz - freq_hz.mean()) / (freq_hz.max() - freq_hz.min())
    chi = np.deg2rad(angle_deg)

    cable_phase = 2.0 * np.pi * freq_hz * 1.1e-9

    bench_common = (
        0.08 * np.sin(2.0 * chi + 0.3)
        + 0.04 * np.sin(3.0 * chi - 0.5)
    )

    direction_sign = +1.0 if direction == "CW" else -1.0
    hysteresis = direction_sign * 0.002 * np.sin(chi + 0.2)

    if device == "DUT":
        orientation_sign = +1.0 if orientation == "NORMAL" else -1.0
        angular_signal = orientation_sign * A_true * np.cos(chi)
    else:
        angular_signal = 0.0

    noise = rng.normal(0.0, 0.002, size=N_FREQ)

    phase = cable_phase + bench_common + hysteresis + angular_signal + noise

    S21_amp = 10.0 ** (-3.0 / 20.0)
    S21 = S21_amp * np.exp(1j * phase)

    S11_amp = 10.0 ** (-18.0 / 20.0)
    S11_phase = 0.1 * np.sin(2.0 * np.pi * x)
    S11 = S11_amp * np.exp(1j * S11_phase)

    filename = angular_filename(device, orientation, direction, angle_deg)
    write_s2p(DATA_DIR / filename, freq_hz, S11, S21)


def gerar_arquivos_angulares():
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


def phase_feature_from_file(path):
    freq_hz, S11, S21 = read_s2p_simple(path)

    phase = np.unwrap(np.angle(S21))

    residual = remove_linear_phase(freq_hz, phase)

    R = float(np.mean(residual))

    return R


def aic_bic(y, yhat, k):
    n = len(y)
    rss = float(np.sum((y - yhat) ** 2))
    rss = max(rss, 1e-30)

    aic = n * np.log(rss / n) + 2 * k
    bic = n * np.log(rss / n) + k * np.log(n)

    return aic, bic, rss


def fit_angular(theta, y, n_perm=2000):
    theta = np.array(theta, dtype=float)
    y = np.array(y, dtype=float)

    X0 = np.ones((len(theta), 1))
    beta0, *_ = np.linalg.lstsq(X0, y, rcond=None)
    y0 = X0 @ beta0

    aic0, bic0, rss0 = aic_bic(y, y0, k=1)

    X = np.column_stack([
        np.cos(theta),
        np.sin(theta),
        np.ones_like(theta),
    ])

    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ beta

    aic, bic, rss = aic_bic(y, yhat, k=3)

    A, B, C = beta

    n, k = X.shape
    dof = max(n - k, 1)
    sigma2 = rss / dof

    try:
        cov = sigma2 * np.linalg.inv(X.T @ X)
        err = np.sqrt(np.diag(cov))
    except Exception:
        err = np.array([np.nan, np.nan, np.nan])

    dA, dB, dC = err

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


def validar_resposta_angular():
    rows = []

    for angle in ANGLES_DEG:
        theta = np.deg2rad(angle)

        for direction in DIRECTIONS:
            for device in DEVICES:
                for orientation in ORIENTATIONS:
                    filename = angular_filename(device, orientation, direction, angle)
                    path = DATA_DIR / filename

                    if not path.exists():
                        return {
                            "status": "BLOQUEADO_ARQUIVOS_ANGULARES_FALTANTES",
                            "report": pd.DataFrame(),
                        }

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

    report_rows = []

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
            z_diff = (A_CW - A_CCW) / sigma if sigma > 0 else np.nan

            direction_rows.append({
                "device": device,
                "A_CW": A_CW,
                "A_CCW": A_CCW,
                "z_diff": z_diff,
                "direction_consistent": abs(z_diff) < 3,
            })

    direction_check = pd.DataFrame(direction_rows)

    dut_all = report[(report["device"] == "DUT") & (report["direction"] == "ALL")]
    load_all = report[(report["device"] == "LOAD") & (report["direction"] == "ALL")]

    if len(dut_all) == 0 or len(load_all) == 0:
        global_status = "BLOQUEADO_RELATORIO_INCOMPLETO"
    else:
        dut_status = dut_all["status"].iloc[0]
        load_status = load_all["status"].iloc[0]

        direction_ok = bool(direction_check["direction_consistent"].all())

        if load_status != "LOAD_PASSOU":
            global_status = "FALSO_POSITIVO_INSTRUMENTAL"
        elif dut_status == "CANDIDATO_ANGULAR_DIFERENCIAL" and direction_ok:
            global_status = "PASSOU_FORTE_RESPOSTA_ANGULAR_DIFERENCIAL"
        elif not direction_ok:
            global_status = "FALHOU_CW_CCW_INCONSISTENTE"
        else:
            global_status = "SEM_CANDIDATO_ANGULAR_FORTE"

    raw.to_csv(DATA_DIR / "teste_05_raw_angular.csv", index=False)
    diff.to_csv(DATA_DIR / "teste_05_diferencial_angular.csv", index=False)
    report.to_csv(DATA_DIR / "teste_05_relatorio_angular.csv", index=False)
    direction_check.to_csv(DATA_DIR / "teste_05_check_cw_ccw.csv", index=False)

    return {
        "status": global_status,
        "report": report,
        "direction_check": direction_check,
    }


# ============================================================
# EXECUCAO DO PIPELINE COMPLETO
# ============================================================

def run_pipeline_completo():
    print("\n===================================================")
    print("TESTE 5 - PIPELINE COMPLETO DO PROTOCOLO")
    print("===================================================")

    print("\nETAPA 1 - Gerar arquivos DEMO de memoria")
    gerar_arquivos_memoria(prefix="", seed_base=1000)
    print("STATUS_ETAPA_1: ARQUIVOS_DEMO_MEMORIA_GERADOS")

    print("\nETAPA 2 - Validar coerencia estatica")
    static_result = validar_coerencia_estatica(label="teste_05_static")

    print(static_result["summary"])
    print("STATUS_ETAPA_2:", static_result["status"])

    if static_result["status"] != "MEMORIA_DE_FASE_COMPATIVEL_EM_BANCADA":
        print("\n[BLOQUEADO] Robustez temporal e angular nao liberados.")
        print("STATUS_FINAL: BLOQUEADO_NA_COHERENCIA_ESTATICA")
        return

    print("\nETAPA 3 - Gerar e validar robustez temporal")
    gerar_arquivos_temporais()
    temporal_result = validar_robustez_temporal()

    print(temporal_result["summary"])
    print("STATUS_ETAPA_3:", temporal_result["status"])

    if temporal_result["status"] != "PASSOU_FORTE_MEMORIA_DE_FASE_REPRODUTIVEL_TEMPORALMENTE":
        print("\n[BLOQUEADO] Angular nao liberado.")
        print("STATUS_FINAL: BLOQUEADO_NA_ROBUSTEZ_TEMPORAL")
        return

    print("\nETAPA 4 - Gerar e validar resposta angular diferencial")
    gerar_arquivos_angulares()
    angular_result = validar_resposta_angular()

    print("\nRelatorio angular:")
    print(angular_result["report"][[
        "device",
        "direction",
        "A",
        "dA",
        "B",
        "dB",
        "delta_AIC",
        "p_perm",
        "leak_B_over_A",
        "status",
    ]])

    print("\nCheck CW/CCW:")
    print(angular_result["direction_check"])

    print("STATUS_ETAPA_4:", angular_result["status"])

    if angular_result["status"] != "PASSOU_FORTE_RESPOSTA_ANGULAR_DIFERENCIAL":
        print("\nSTATUS_FINAL: BLOQUEADO_NO_TESTE_ANGULAR")
        return

    print("\n===================================================")
    print("VEREDITO FINAL")
    print("===================================================")
    print("Teste 2 - Coerencia estatica:", static_result["status"])
    print("Teste 3 - Robustez temporal:", temporal_result["status"])
    print("Teste 4 - Angular diferencial:", angular_result["status"])
    print("\nSTATUS_FINAL: PROTOCOLO_COMPLETO_PASSOU_NA_ORDEM_CORRETA")
    print("\nOBS: Resultado DEMO/SIMULADO. Nao e evidencia experimental real.")


# ============================================================
# EXECUTAR
# ============================================================

if __name__ == "__main__":
    run_pipeline_completo()
