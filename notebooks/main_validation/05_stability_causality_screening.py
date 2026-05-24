# ============================================================
# TVGD PUBLIC TESTS
# TEST 05 — STABILITY, CAUSALITY AND LOCAL SCREENING
# Versão autônoma para Google Colab
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Configuração
# ============================================================

OUTDIR = "Resultados/test05_stability_causality_screening"
os.makedirs(OUTDIR, exist_ok=True)

# Constantes
A0 = 1.2e-10  # m/s^2
G = 6.67430e-11
C = 299792458.0

# Massa e raios principais
M_EARTH = 5.9722e24
R_EARTH = 6.371e6

M_SUN = 1.98847e30
R_SUN = 6.957e8

AU = 1.495978707e11

# Raios orbitais aproximados
R_GPS = R_EARTH + 20200e3
R_MOON = 384400e3
R_MERCURY = 0.387098 * AU
R_VENUS = 0.723332 * AU
R_EARTH_ORBIT = AU
R_MARS = 1.523679 * AU
R_JUPITER = 5.2044 * AU
R_SATURN = 9.5826 * AU
R_NEPTUNE = 30.07 * AU

# Grade de u
U_MIN = 1e-8
U_MAX = 1e8
N_U = 20000

# Critérios
MAX_LOCAL_ALPHA_MINUS_1 = 1e-12
CS2_MIN_ALLOWED = 0.0
CS2_MAX_ALLOWED = 1.0

# ============================================================
# Funções TVGD
# ============================================================

def P_tvgd_u(u):
    u = np.asarray(u, dtype=float)
    return 1.0 - np.exp(-u)


def D_tvgd_u(u):
    """
    D(u) aparece como denominador efetivo no teste de estabilidade:

    D(u) = P(u) + (u/2) exp(-u)
    """
    u = np.asarray(u, dtype=float)
    return P_tvgd_u(u) + 0.5 * u * np.exp(-u)


def cs2_tvgd_u(u):
    """
    Velocidade de som efetiva ao quadrado:

    c_s^2 = P(u) / D(u)
    """
    P = P_tvgd_u(u)
    D = D_tvgd_u(u)
    return P / np.maximum(D, 1e-300)


def alpha_tvgd_u(u):
    """
    Fator de amplificação estrutural:

    alpha = 1/P(u)
    """
    P = P_tvgd_u(u)
    return 1.0 / np.maximum(P, 1e-300)


def u_from_g(g):
    return np.sqrt(np.asarray(g, dtype=float) / A0)


def g_newton(M, r):
    return G * M / r**2


# ============================================================
# Teste matemático global em u
# ============================================================

print("=" * 72)
print("TESTE 05 — STABILITY, CAUSALITY AND LOCAL SCREENING")
print("=" * 72)

u_grid = np.logspace(np.log10(U_MIN), np.log10(U_MAX), N_U)

P_grid = P_tvgd_u(u_grid)
D_grid = D_tvgd_u(u_grid)
cs2_grid = cs2_tvgd_u(u_grid)
alpha_grid = alpha_tvgd_u(u_grid)

finite_all = (
    np.all(np.isfinite(P_grid))
    and np.all(np.isfinite(D_grid))
    and np.all(np.isfinite(cs2_grid))
    and np.all(np.isfinite(alpha_grid))
)

P_positive = np.all(P_grid > 0)
D_positive = np.all(D_grid > 0)
cs2_positive = np.all(cs2_grid > CS2_MIN_ALLOWED)
cs2_subluminal = np.all(cs2_grid <= CS2_MAX_ALLOWED + 1e-14)

math_summary = {
    "finite_all": bool(finite_all),
    "P_positive": bool(P_positive),
    "D_positive": bool(D_positive),
    "cs2_positive": bool(cs2_positive),
    "cs2_subluminal": bool(cs2_subluminal),
    "P_min": float(np.min(P_grid)),
    "P_max": float(np.max(P_grid)),
    "D_min": float(np.min(D_grid)),
    "D_max": float(np.max(D_grid)),
    "cs2_min": float(np.min(cs2_grid)),
    "cs2_max": float(np.max(cs2_grid)),
    "alpha_min": float(np.min(alpha_grid)),
    "alpha_max": float(np.max(alpha_grid)),
}

print("\nResumo matemático global:")
for k, v in math_summary.items():
    print(f"{k}: {v}")

# ============================================================
# Pontos físicos representativos
# ============================================================

physical_points = []

def add_point(name, g_value, source):
    u = float(u_from_g(g_value))
    P = float(P_tvgd_u(u))
    D = float(D_tvgd_u(u))
    cs2 = float(cs2_tvgd_u(u))
    alpha = float(alpha_tvgd_u(u))
    alpha_minus_1 = alpha - 1.0

    physical_points.append({
        "name": name,
        "source": source,
        "g_m_s2": float(g_value),
        "u": u,
        "P": P,
        "D": D,
        "cs2": cs2,
        "alpha": alpha,
        "alpha_minus_1": alpha_minus_1,
    })

# Regimes abstratos úteis
for u0 in [1e-4, 0.1, 0.5, 1.0, 2.94, 10.0, 100.0]:
    g0 = A0 * u0**2
    add_point(f"u={u0}", g0, "dimensionless_reference")

# Terra superfície
g_earth_surface = g_newton(M_EARTH, R_EARTH)
add_point("Earth surface", g_earth_surface, "Earth gravity")

# GPS
g_gps = g_newton(M_EARTH, R_GPS)
add_point("GPS orbit", g_gps, "Earth gravity")

# Lua
g_moon_orbit = g_newton(M_EARTH, R_MOON)
add_point("Moon orbit around Earth", g_moon_orbit, "Earth gravity")

# Sol superfície
g_sun_surface = g_newton(M_SUN, R_SUN)
add_point("Solar surface", g_sun_surface, "Solar gravity")

# Órbitas planetárias ao redor do Sol
planet_orbits = [
    ("Mercury orbit", R_MERCURY),
    ("Venus orbit", R_VENUS),
    ("Earth orbit", R_EARTH_ORBIT),
    ("Mars orbit", R_MARS),
    ("Jupiter orbit", R_JUPITER),
    ("Saturn orbit", R_SATURN),
    ("Neptune orbit", R_NEPTUNE),
]

for name, r in planet_orbits:
    g = g_newton(M_SUN, r)
    add_point(name, g, "Solar gravity")

points_df = pd.DataFrame(physical_points)

print("\nPontos físicos:")
display(points_df)

# ============================================================
# Blindagem local
# ============================================================

local_names = [
    "Earth surface",
    "GPS orbit",
    "Moon orbit around Earth",
    "Solar surface",
    "Mercury orbit",
    "Venus orbit",
    "Earth orbit",
    "Mars orbit",
    "Jupiter orbit",
    "Saturn orbit",
    "Neptune orbit",
]

local_df = points_df[points_df["name"].isin(local_names)].copy()

max_local_alpha_minus_1 = float(np.max(np.abs(local_df["alpha_minus_1"])))
max_local_delta_g_over_g = max_local_alpha_minus_1

local_screening_ok = max_local_alpha_minus_1 <= MAX_LOCAL_ALPHA_MINUS_1

screening_summary = {
    "N_local_points": int(len(local_df)),
    "max_local_alpha_minus_1": max_local_alpha_minus_1,
    "max_local_delta_g_over_g": max_local_delta_g_over_g,
    "local_screening_ok": bool(local_screening_ok),
}

print("\nResumo blindagem local:")
for k, v in screening_summary.items():
    print(f"{k}: {v}")

# ============================================================
# Diagnóstico automático
# ============================================================

stability_ok = (
    finite_all
    and P_positive
    and D_positive
    and cs2_positive
    and cs2_subluminal
)

if stability_ok and local_screening_ok:
    final_status = "PASSOU_FORTE_ESTABILIDADE_CAUSALIDADE_BLINDAGEM"
elif stability_ok:
    final_status = "PASSOU_MODERADO_ESTABILIDADE_CAUSALIDADE"
else:
    final_status = "TENSIONADO_ESTABILIDADE_CAUSALIDADE"

summary = {}
summary.update(math_summary)
summary.update(screening_summary)
summary["stability_ok"] = bool(stability_ok)
summary["final_status"] = final_status

summary_df = pd.DataFrame([summary])

print("\nResumo final:")
display(summary_df)

print("\nStatus final:", final_status)

# ============================================================
# Salvar resultados
# ============================================================

grid_df = pd.DataFrame({
    "u": u_grid,
    "P": P_grid,
    "D": D_grid,
    "cs2": cs2_grid,
    "alpha": alpha_grid,
})

grid_path = os.path.join(OUTDIR, "test05_stability_grid.csv")
points_path = os.path.join(OUTDIR, "test05_physical_points.csv")
summary_path = os.path.join(OUTDIR, "test05_stability_causality_screening_summary.csv")

grid_df.to_csv(grid_path, index=False)
points_df.to_csv(points_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Gráficos
# ============================================================

plt.figure(figsize=(8, 5))
plt.plot(u_grid, P_grid, label="P(u)")
plt.xscale("log")
plt.xlabel("u")
plt.ylabel("P(u)")
plt.title("TVGD response function")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test05_P_u.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(u_grid, cs2_grid, label="c_s^2")
plt.axhline(1.0, color="black", linestyle="--", linewidth=1, label="luminal bound")
plt.xscale("log")
plt.xlabel("u")
plt.ylabel("c_s^2")
plt.title("Effective sound-speed stability")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test05_cs2_u.png"), dpi=220)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(u_grid, alpha_grid, label="alpha = 1/P(u)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("u")
plt.ylabel("alpha")
plt.title("TVGD amplification factor")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test05_alpha_u.png"), dpi=220)
plt.show()

plt.figure(figsize=(9, 5))
plot_local = local_df.copy()
plot_local["abs_alpha_minus_1"] = np.abs(plot_local["alpha_minus_1"])

plt.bar(plot_local["name"], plot_local["abs_alpha_minus_1"])
plt.axhline(MAX_LOCAL_ALPHA_MINUS_1, color="black", linestyle="--", label="criterion")
plt.yscale("log")
plt.ylabel("|alpha - 1|")
plt.title("Local screening in high-acceleration regimes")
plt.xticks(rotation=45, ha="right")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test05_local_screening.png"), dpi=220)
plt.show()

print("\nArquivos salvos em:")
print(OUTDIR)
print(grid_path)
print(points_path)
print(summary_path)

print("\nInterpretação curta:")
if final_status == "PASSOU_FORTE_ESTABILIDADE_CAUSALIDADE_BLINDAGEM":
    print("O núcleo P(u)=1-exp(-u) passa nos critérios de positividade, estabilidade, causalidade subluminal e blindagem local.")
elif "PASSOU" in final_status:
    print("O núcleo passa em estabilidade/causalidade, mas a blindagem local precisa de revisão.")
else:
    print("O núcleo apresentou tensão matemática e deve ser revisado antes de uso público.")

print("\nDONE.")
