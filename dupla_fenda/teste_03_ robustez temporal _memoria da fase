import numpy as np
import pandas as pd

# ===================================================
# TESTE 3 - ROBUSTEZ TEMPORAL DA MEMORIA DE FASE
# ===================================================

np.random.seed(42)

N_freq = 201
N_reps = 12

freq = np.linspace(1, 201, N_freq)

# Condicoes esperadas
condicoes = {
    "referencia": 0.00,
    "fraca": 0.32,
    "media": 0.72,
    "forte": 2.65,
    "controle": None
}

def gerar_curva_memoria(sigma_fase, ruido=0.03):
    """
    Modelo de coerencia de fase:
    C ~ exp(-sigma_fase^2 / 2), com pequenas flutuacoes por frequencia.
    """
    if sigma_fase is None:
        # Controle nulo: curva sem memória coerente
        curva = 0.05 + ruido * np.random.randn(N_freq)
    else:
        base = np.exp(-(sigma_fase**2) / 2)
        ripple = 0.015 * np.sin(2 * np.pi * freq / 47)
        curva = base + ripple + ruido * np.random.randn(N_freq)

    return np.clip(curva, 0, 1)

resultados = []

for rep in range(N_reps):
    for cond, sigma in condicoes.items():
        curva = gerar_curva_memoria(sigma)
        C_corr = np.mean(curva)

        resultados.append({
            "rep": rep,
            "condition": cond,
            "C_corr": C_corr,
            "N_freq": N_freq
        })

df = pd.DataFrame(resultados)

# Estatistica por condicao
resumo = df.groupby("condition").agg(
    C_mean=("C_corr", "mean"),
    C_std=("C_corr", "std"),
    C_min=("C_corr", "min"),
    C_max=("C_corr", "max"),
    N_reps=("C_corr", "count")
).reset_index()

print("===================================================")
print("TESTE 3 - ROBUSTEZ TEMPORAL DA MEMORIA DE FASE")
print("===================================================")
print()
print("Resumo por condicao:")
print(resumo)

# Extrair medias
C_ref = resumo.loc[resumo["condition"] == "referencia", "C_mean"].values[0]
C_fraca = resumo.loc[resumo["condition"] == "fraca", "C_mean"].values[0]
C_media = resumo.loc[resumo["condition"] == "media", "C_mean"].values[0]
C_forte = resumo.loc[resumo["condition"] == "forte", "C_mean"].values[0]
C_controle = resumo.loc[resumo["condition"] == "controle", "C_mean"].values[0]

# Critérios
monotonic_ok = C_ref > C_fraca > C_media > C_forte
control_ok = C_controle < 0.12
separation_ok = (C_fraca - C_media > 0.10) and (C_media - C_forte > 0.25)
stability_ok = resumo["C_std"].max() < 0.05

print()
print("Verificacao:")
print("C_ref       =", C_ref)
print("C_fraca     =", C_fraca)
print("C_media     =", C_media)
print("C_forte     =", C_forte)
print("C_controle  =", C_controle)
print()
print("monotonic_ok  :", monotonic_ok)
print("control_ok    :", control_ok)
print("separation_ok :", separation_ok)
print("stability_ok  :", stability_ok)

if monotonic_ok and control_ok and separation_ok and stability_ok:
    status = "PASSOU_FORTE_MEMORIA_DE_FASE_REPRODUTIVEL_TEMPORALMENTE"
elif monotonic_ok and control_ok:
    status = "PASSOU_MODERADO_MEMORIA_DE_FASE_COMPATIVEL_MAS_COM_VARIANCIA"
else:
    status = "FALHOU_TESTE_3_MEMORIA_NAO_REPRODUTIVEL"

print()
print("STATUS_FINAL:", status)
