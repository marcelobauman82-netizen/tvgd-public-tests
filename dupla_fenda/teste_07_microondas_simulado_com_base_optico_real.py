# ============================================================
# TESTE SIMULADO MICRO-ONDAS USANDO DADO ÓPTICO REAL/DIGITALIZADO
# TVGD-MEM | Ponte: M_optico_real -> sigma_med -> C_corr_micro
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

np.random.seed(2026)

print("=====================================================================")
print("TESTE SIMULADO MICRO-ONDAS COM BASE NO DADO ÓPTICO REAL/DIGITALIZADO")
print("TVGD-MEM | M_optico_real -> sigma_med -> C_corr_micro")
print("=====================================================================\n")

# ============================================================
# 1. DADO ÓPTICO REAL/DIGITALIZADO JÁ EXTRAÍDO
# ============================================================

M_optico_real = 0.6276459525749891

# Pela hipótese TVGD-MEM:
# M_struct ≈ exp(-sigma_med^2 / 2)
sigma_equiv = np.sqrt(-2 * np.log(M_optico_real))

print("DADO DE ENTRADA")
print("M_optico_real/dig =", M_optico_real)
print("sigma_med equivalente =", sigma_equiv)

# ============================================================
# 2. SIMULAÇÃO DE MICRO-ONDAS
# ============================================================

N_freq = 201
N_trials = 1000

freq = np.linspace(8e9, 12e9, N_freq)  # faixa simulada: 8 GHz a 12 GHz

def calcula_C_corr(delta_phi):
    """
    Métrica fasorial:
        C_corr = | < exp(i delta_phi) > |
    """
    return np.abs(np.mean(np.exp(1j * delta_phi)))

def simular_microondas_por_sigma(sigma, N_freq=201):
    """
    Simula uma varredura de fase relativa em micro-ondas.
    delta_phi ~ Normal(0, sigma)
    """
    delta_phi = np.random.normal(0, sigma, N_freq)
    C_corr = calcula_C_corr(delta_phi)
    return C_corr, delta_phi

# Rodar muitas repetições para estimar estabilidade
resultados = []

for i in range(N_trials):
    C_corr, delta_phi = simular_microondas_por_sigma(sigma_equiv, N_freq=N_freq)
    
    resultados.append({
        "trial": i,
        "C_corr_micro": C_corr,
        "M_optico_real": M_optico_real,
        "erro_abs": abs(C_corr - M_optico_real)
    })

df = pd.DataFrame(resultados)

C_mean = df["C_corr_micro"].mean()
C_std = df["C_corr_micro"].std()
C_median = df["C_corr_micro"].median()
erro_mean = df["erro_abs"].mean()
erro_median = df["erro_abs"].median()

# Uma execução representativa para plot
C_rep, delta_phi_rep = simular_microondas_por_sigma(sigma_equiv, N_freq=N_freq)

print("\n============================================================")
print("RESULTADO ESTATÍSTICO MICRO-ONDAS SIMULADO")
print("============================================================")
print("N_trials              =", N_trials)
print("N_freq                =", N_freq)
print("M_optico_real/dig     =", M_optico_real)
print("sigma_equiv           =", sigma_equiv)
print("C_corr_micro_mean     =", C_mean)
print("C_corr_micro_std      =", C_std)
print("C_corr_micro_median   =", C_median)
print("erro_abs_mean         =", erro_mean)
print("erro_abs_median       =", erro_median)
print("C_corr_representativo =", C_rep)

# ============================================================
# 3. CRITÉRIOS
# ============================================================

erro_tolerancia_forte = 0.05
erro_tolerancia_moderado = 0.10

erro_final = abs(C_mean - M_optico_real)

forte_ok = erro_final < erro_tolerancia_forte
moderado_ok = erro_final < erro_tolerancia_moderado

print("\n============================================================")
print("VERIFICAÇÕES")
print("============================================================")
print("erro_final_mean_vs_optico =", erro_final)
print("forte_ok                  =", forte_ok)
print("moderado_ok               =", moderado_ok)

if forte_ok:
    status = "PASSOU_FORTE_MICROONDAS_SIMULADO_COM_BASE_EM_DADO_OPTICO_REAL"
elif moderado_ok:
    status = "PASSOU_MODERADO_MICROONDAS_SIMULADO_COM_BASE_EM_DADO_OPTICO_REAL"
else:
    status = "FALHOU_MICROONDAS_SIMULADO_COM_BASE_EM_DADO_OPTICO_REAL"

print("\nSTATUS_FINAL:", status)

# ============================================================
# 4. GRÁFICOS
# ============================================================

plt.figure(figsize=(8,4))
plt.hist(df["C_corr_micro"], bins=40, alpha=0.7)
plt.axvline(M_optico_real, linestyle="--", linewidth=2, label=f"M óptico real = {M_optico_real:.4f}")
plt.axvline(C_mean, linestyle="-", linewidth=2, label=f"Média micro = {C_mean:.4f}")
plt.xlabel("C_corr micro-ondas simulado")
plt.ylabel("contagem")
plt.title("Distribuição de C_corr em micro-ondas simulado")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

plt.figure(figsize=(8,4))
plt.bar(["Óptico real/dig", "Micro simulado médio"], [M_optico_real, C_mean])
plt.ylim(0, 1.05)
plt.ylabel("M_struct / C_corr")
plt.title(f"Comparação M_optico real vs C_corr micro simulado | erro={erro_final:.4f}")
plt.grid(True, alpha=0.3)
plt.show()

plt.figure(figsize=(8,4))
plt.plot(freq/1e9, np.unwrap(delta_phi_rep), label="fase relativa simulada")
plt.xlabel("frequência [GHz]")
plt.ylabel("fase relativa [rad]")
plt.title(f"Exemplo de fase relativa micro-ondas | C_corr={C_rep:.4f}")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# Círculo fasorial
phasors = np.exp(1j * delta_phi_rep)
mean_phasor = np.mean(phasors)

plt.figure(figsize=(5,5))
plt.scatter(phasors.real, phasors.imag, s=15, alpha=0.5, label="fasores")
plt.arrow(0, 0, mean_phasor.real, mean_phasor.imag,
          head_width=0.04, length_includes_head=True,
          label="vetor médio")
theta = np.linspace(0, 2*np.pi, 400)
plt.plot(np.cos(theta), np.sin(theta), linestyle="--", alpha=0.5)
plt.axis("equal")
plt.xlabel("Re")
plt.ylabel("Im")
plt.title(f"Círculo fasorial micro-ondas | C_corr={C_rep:.4f}")
plt.grid(True, alpha=0.3)
plt.legend()
plt.show()

# ============================================================
# 5. EXPORTAR RESULTADOS
# ============================================================

df.to_csv("resultado_microondas_simulado_com_base_optico_real.csv", index=False)

df_summary = pd.DataFrame([{
    "M_optico_real_dig": M_optico_real,
    "sigma_equiv": sigma_equiv,
    "C_corr_micro_mean": C_mean,
    "C_corr_micro_std": C_std,
    "C_corr_micro_median": C_median,
    "erro_final_mean_vs_optico": erro_final,
    "N_trials": N_trials,
    "N_freq": N_freq,
    "STATUS_FINAL": status
}])

df_summary.to_csv("summary_microondas_simulado_com_base_optico_real.csv", index=False)

print("\nArquivos salvos:")
print("resultado_microondas_simulado_com_base_optico_real.csv")
print("summary_microondas_simulado_com_base_optico_real.csv")
