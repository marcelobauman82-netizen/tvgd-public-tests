# ============================================================
# TVGD PUBLIC TESTS
# TEST 02 — BLACK-HOLE EXTERIOR AND REGULAR-CORE CHECK
# Versão autônoma para Google Colab
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# Configuração
# ============================================================

OUTDIR = "/content/tvgd_test02_black_hole_regular_core"
os.makedirs(OUTDIR, exist_ok=True)

# Parâmetros testados
N_VALUES = [5, 6, 8]
XC_VALUES = [0.001, 0.003, 0.01, 0.03, 0.1]

# Referências Schwarzschild em unidades x = r/r_s
X_H_REF = 1.0
X_PH_REF = 1.5
B_SHADOW_REF = 3.0 * np.sqrt(3.0) / 2.0   # 2.598076...
X_ISCO_REF = 3.0

# Critérios
COMPACT_XC_MAX = 0.01
N_STRONG_MIN = 6

TOL_HORIZON_PERCENT = 1e-3
TOL_PH_PERCENT = 1e-3
TOL_SHADOW_PERCENT = 1e-3
TOL_ISCO_PERCENT = 1e-2

EPS = 1e-12

# ============================================================
# Funções TVGD BH
# ============================================================

def S_tvgd(x, n, xc):
    """
    Saturação estrutural TVGD:

    S(x) = 1 - exp[-(x/xc)^n]
    """
    x = np.asarray(x, dtype=float)
    y = np.maximum(x / xc, 0.0)
    return 1.0 - np.exp(-(y ** n))


def A_metric(x, n, xc):
    """
    Métrica efetiva:

    A(x) = 1 - S(x)/x
    """
    x = np.asarray(x, dtype=float)
    return 1.0 - S_tvgd(x, n, xc) / np.maximum(x, EPS)


def numerical_derivative(func, x, n, xc, h=1e-5):
    x = np.asarray(x, dtype=float)
    xp = x + h
    xm = np.maximum(x - h, EPS)
    return (func(xp, n, xc) - func(xm, n, xc)) / (xp - xm)


def numerical_second_derivative(func, x, n, xc, h=1e-4):
    x = np.asarray(x, dtype=float)
    xp = x + h
    xm = np.maximum(x - h, EPS)
    return (func(xp, n, xc) - 2.0 * func(x, n, xc) + func(xm, n, xc)) / (h ** 2)


# ============================================================
# Horizonte
# ============================================================

def find_roots_by_sign(func_values, x_grid):
    roots = []
    y = func_values

    for i in range(len(x_grid) - 1):
        if not np.isfinite(y[i]) or not np.isfinite(y[i + 1]):
            continue

        if y[i] == 0:
            roots.append(x_grid[i])

        if y[i] * y[i + 1] < 0:
            x1, x2 = x_grid[i], x_grid[i + 1]
            y1, y2 = y[i], y[i + 1]

            # interpolação linear simples
            xr = x1 - y1 * (x2 - x1) / (y2 - y1)
            roots.append(xr)

    return roots


def find_horizons(n, xc):
    """
    Horizontes são raízes de A(x)=0.
    """
    x_grid = np.logspace(-8, 2, 30000)
    A = A_metric(x_grid, n, xc)
    roots = find_roots_by_sign(A, x_grid)

    # Remover duplicados próximos
    clean = []
    for r in roots:
        if r > 0 and all(abs(r - q) > 1e-4 for q in clean):
            clean.append(r)

    return clean


# ============================================================
# Esfera de fótons e sombra
# ============================================================

def photon_sphere_function(x, n, xc):
    """
    Condição para órbita circular nula:

    d/dx [ A(x)/x^2 ] = 0

    Equivalente:

    x A'(x) - 2 A(x) = 0
    """
    A = A_metric(x, n, xc)
    Ap = numerical_derivative(A_metric, x, n, xc, h=1e-5)
    return x * Ap - 2.0 * A


def find_photon_sphere(n, xc):
    x_grid = np.linspace(1.0001, 10.0, 40000)
    F = photon_sphere_function(x_grid, n, xc)
    roots = find_roots_by_sign(F, x_grid)

    if len(roots) == 0:
        return np.nan

    # Queremos a raiz externa mais próxima de 1.5
    roots = np.array(roots)
    roots = roots[np.isfinite(roots)]
    if len(roots) == 0:
        return np.nan

    return float(roots[np.argmin(np.abs(roots - X_PH_REF))])


def shadow_impact_parameter(x_ph, n, xc):
    """
    b = x / sqrt(A(x))
    em unidades de r_s.
    """
    if not np.isfinite(x_ph):
        return np.nan

    A = A_metric(x_ph, n, xc)
    if A <= 0:
        return np.nan

    return float(x_ph / np.sqrt(A))


# ============================================================
# ISCO timelike
# ============================================================

def L2_circular(x, n, xc):
    """
    Momento angular específico ao quadrado para órbita circular timelike:

    L^2 = x^3 A'(x) / [2A(x) - xA'(x)]

    em unidades r_s^2.
    """
    A = A_metric(x, n, xc)
    Ap = numerical_derivative(A_metric, x, n, xc, h=1e-5)
    denom = 2.0 * A - x * Ap

    out = np.full_like(np.asarray(x, dtype=float), np.nan, dtype=float)
    mask = np.isfinite(A) & np.isfinite(Ap) & np.isfinite(denom) & (denom > 0)
    out[mask] = x[mask] ** 3 * Ap[mask] / denom[mask]
    return out


def find_isco(n, xc):
    """
    ISCO: mínimo de L^2(x) para órbitas circulares estáveis externas.

    Para Schwarzschild em x=r/r_s, x_ISCO=3.
    """
    x_grid = np.linspace(1.5001, 20.0, 60000)
    L2 = L2_circular(x_grid, n, xc)

    mask = np.isfinite(L2) & (L2 > 0)

    if np.sum(mask) < 10:
        return np.nan

    x_valid = x_grid[mask]
    L2_valid = L2[mask]

    # mínimo global externo
    idx = np.argmin(L2_valid)
    return float(x_valid[idx])


# ============================================================
# Regularidade do núcleo
# ============================================================

def core_regular_diagnostics(n, xc):
    """
    Verificação simples de regularidade no núcleo.

    Para x -> 0:

    S(x) ~ (x/xc)^n
    A(x) = 1 - S(x)/x ~ 1 - x^(n-1)/xc^n

    Para n >= 3, o núcleo fica regular no sentido de não divergir
    como Schwarzschild em A ~ 1 - 1/x.
    """
    x_core = np.logspace(-8, -3, 2000)

    A = A_metric(x_core, n, xc)
    Ap = numerical_derivative(A_metric, x_core, n, xc, h=1e-8)
    App = numerical_second_derivative(A_metric, x_core, n, xc, h=1e-7)

    finite_A = np.all(np.isfinite(A))
    finite_Ap = np.all(np.isfinite(Ap))
    finite_App = np.all(np.isfinite(App))

    max_abs_A = float(np.nanmax(np.abs(A)))
    max_abs_Ap = float(np.nanmax(np.abs(Ap)))
    max_abs_App = float(np.nanmax(np.abs(App)))

    A0_like = float(A[0])

    regular_core = finite_A and finite_Ap and finite_App and (n >= 3)

    return {
        "finite_A": bool(finite_A),
        "finite_Ap": bool(finite_Ap),
        "finite_App": bool(finite_App),
        "max_abs_A_core": max_abs_A,
        "max_abs_Ap_core": max_abs_Ap,
        "max_abs_App_core": max_abs_App,
        "A_near_zero": A0_like,
        "regular_core": bool(regular_core),
    }


# ============================================================
# Rodar varredura
# ============================================================

print("=" * 72)
print("TESTE 02 — BLACK-HOLE EXTERIOR AND REGULAR-CORE CHECK")
print("=" * 72)

rows = []

for n in N_VALUES:
    for xc in XC_VALUES:

        horizons = find_horizons(n, xc)
        x_h_outer = max(horizons) if len(horizons) > 0 else np.nan

        x_ph = find_photon_sphere(n, xc)
        b_shadow = shadow_impact_parameter(x_ph, n, xc)
        x_isco = find_isco(n, xc)

        core = core_regular_diagnostics(n, xc)

        err_h_percent = 100.0 * abs(x_h_outer - X_H_REF) / X_H_REF if np.isfinite(x_h_outer) else np.nan
        err_ph_percent = 100.0 * abs(x_ph - X_PH_REF) / X_PH_REF if np.isfinite(x_ph) else np.nan
        err_shadow_percent = 100.0 * abs(b_shadow - B_SHADOW_REF) / B_SHADOW_REF if np.isfinite(b_shadow) else np.nan
        err_isco_percent = 100.0 * abs(x_isco - X_ISCO_REF) / X_ISCO_REF if np.isfinite(x_isco) else np.nan

        compact_core = (xc <= COMPACT_XC_MAX) and (n >= N_STRONG_MIN)

        exterior_ok = (
            np.isfinite(x_h_outer)
            and np.isfinite(x_ph)
            and np.isfinite(b_shadow)
            and np.isfinite(x_isco)
            and err_h_percent < TOL_HORIZON_PERCENT
            and err_ph_percent < TOL_PH_PERCENT
            and err_shadow_percent < TOL_SHADOW_PERCENT
            and err_isco_percent < TOL_ISCO_PERCENT
        )

        strong_ok = compact_core and exterior_ok and core["regular_core"]

        if strong_ok:
            status = "PASSA_FORTE"
        elif core["regular_core"] and np.isfinite(x_h_outer):
            status = "PASSA_REGULAR_MAS_NAO_FORTE"
        elif core["regular_core"]:
            status = "REGULAR_SEM_HORIZONTE_EXTERNO_CLARO"
        else:
            status = "FALHA_REGULARIDADE"

        row = {
            "n": n,
            "xc": xc,
            "compact_core": compact_core,
            "N_horizons": len(horizons),
            "horizons": str([float(h) for h in horizons]),
            "x_h_outer": x_h_outer,
            "x_ph": x_ph,
            "b_shadow": b_shadow,
            "x_isco": x_isco,
            "err_h_percent": err_h_percent,
            "err_ph_percent": err_ph_percent,
            "err_shadow_percent": err_shadow_percent,
            "err_isco_percent": err_isco_percent,
            "exterior_ok": exterior_ok,
            "regular_core": core["regular_core"],
            "finite_A": core["finite_A"],
            "finite_Ap": core["finite_Ap"],
            "finite_App": core["finite_App"],
            "A_near_zero": core["A_near_zero"],
            "max_abs_A_core": core["max_abs_A_core"],
            "max_abs_Ap_core": core["max_abs_Ap_core"],
            "max_abs_App_core": core["max_abs_App_core"],
            "status": status,
        }

        rows.append(row)

df = pd.DataFrame(rows)

print("\nTabela principal:")
display(df)

# ============================================================
# Resumo
# ============================================================

N_total = len(df)
N_regular = int(df["regular_core"].sum())
N_exterior_ok = int(df["exterior_ok"].sum())
N_strong = int((df["status"] == "PASSA_FORTE").sum())

compact_df = df[df["compact_core"] == True].copy()
N_compact = len(compact_df)
N_compact_strong = int((compact_df["status"] == "PASSA_FORTE").sum())

if N_compact > 0 and N_compact_strong == N_compact:
    final_status = "PASSOU_FORTE_BH_EXTERIOR_REGULAR_CORE"
elif N_regular == N_total and N_exterior_ok > 0:
    final_status = "PASSOU_MODERADO_BH_REGULARIDADE_COM_EXTERIOR"
elif N_regular == N_total:
    final_status = "PASSOU_FRACO_BH_REGULARIDADE"
else:
    final_status = "TENSIONADO_BH_REGULAR_CORE"

summary = {
    "N_total": N_total,
    "N_regular": N_regular,
    "N_exterior_ok": N_exterior_ok,
    "N_strong": N_strong,
    "N_compact": N_compact,
    "N_compact_strong": N_compact_strong,
    "final_status": final_status,
}

summary_df = pd.DataFrame([summary])

print("\nResumo:")
display(summary_df)

print("\nStatus final:", final_status)

# ============================================================
# Salvar resultados
# ============================================================

main_path = os.path.join(OUTDIR, "test02_black_hole_regular_core_results.csv")
summary_path = os.path.join(OUTDIR, "test02_black_hole_regular_core_summary.csv")

df.to_csv(main_path, index=False)
summary_df.to_csv(summary_path, index=False)

# ============================================================
# Gráficos
# ============================================================

# Gráfico A(x) para casos representativos
plt.figure(figsize=(8, 5))

for n, xc in [(6, 0.001), (6, 0.01), (8, 0.001), (8, 0.01), (6, 0.1)]:
    x = np.linspace(0.001, 5.0, 3000)
    A = A_metric(x, n, xc)
    plt.plot(x, A, label=f"n={n}, xc={xc}")

plt.axhline(0, color="black", linewidth=1)
plt.axvline(1.0, linestyle="--", label="Schwarzschild horizon x=1")
plt.xlabel("x = r / r_s")
plt.ylabel("A(x)")
plt.title("TVGD regular-core metric function")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test02_A_metric_profiles.png"), dpi=220)
plt.show()

# Erros externos
plot_df = df.copy()
plot_df["case"] = plot_df.apply(lambda r: f"n={int(r['n'])}, xc={r['xc']}", axis=1)

plt.figure(figsize=(10, 5))
plt.bar(plot_df["case"], plot_df["err_ph_percent"])
plt.axhline(TOL_PH_PERCENT, color="black", linestyle="--", label="criterion")
plt.ylabel("Photon-sphere error [%]")
plt.title("Photon-sphere preservation")
plt.xticks(rotation=45, ha="right")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test02_photon_sphere_errors.png"), dpi=220)
plt.show()

plt.figure(figsize=(10, 5))
plt.bar(plot_df["case"], plot_df["err_shadow_percent"])
plt.axhline(TOL_SHADOW_PERCENT, color="black", linestyle="--", label="criterion")
plt.ylabel("Shadow-impact error [%]")
plt.title("Shadow preservation")
plt.xticks(rotation=45, ha="right")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test02_shadow_errors.png"), dpi=220)
plt.show()

plt.figure(figsize=(10, 5))
plt.bar(plot_df["case"], plot_df["err_isco_percent"])
plt.axhline(TOL_ISCO_PERCENT, color="black", linestyle="--", label="criterion")
plt.ylabel("ISCO error [%]")
plt.title("ISCO preservation")
plt.xticks(rotation=45, ha="right")
plt.grid(True, axis="y", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(OUTDIR, "test02_isco_errors.png"), dpi=220)
plt.show()

# ============================================================
# Impressão final
# ============================================================

print("\nArquivos salvos em:")
print(OUTDIR)
print(main_path)
print(summary_path)

print("\nInterpretação curta:")
if final_status == "PASSOU_FORTE_BH_EXTERIOR_REGULAR_CORE":
    print("A saturação TVGD regulariza o núcleo e preserva os observáveis externos clássicos para núcleos compactos n>=6 e xc<=0.01.")
elif "PASSOU" in final_status:
    print("A regularidade do núcleo foi preservada, mas alguns critérios externos fortes exigem revisão ou parametrização mais restrita.")
else:
    print("O setor de buracos negros apresentou tensão e precisa ser revisado.")

print("\nDONE.")
