"""
tvgd_core.py

Core numerical functions for the Dispersive Gravitational Vacuum Theory (TVGD).

This module contains the basic effective functions used in the public
reproducibility tests.

Author: Marcelo Bauman
Repository: TVGD Public Tests
"""

import numpy as np


# ============================================================
# Constants
# ============================================================

C_LIGHT = 299792458.0          # speed of light in m/s
C_LIGHT_KM = 299792.458        # speed of light in km/s
G_NEWTON = 6.67430e-11         # Newtonian gravitational constant in SI
A0_DEFAULT = 1.2e-10           # TVGD acceleration scale in m/s^2


# ============================================================
# Local / Galactic TVGD Response
# ============================================================

def tvgd_u(g_bar, a0=A0_DEFAULT):
    """
    Dimensionless TVGD structural variable.

    u = sqrt(g_bar / a0)

    Parameters
    ----------
    g_bar : float or array-like
        Baryonic/Newtonian acceleration in m/s^2.
    a0 : float
        TVGD acceleration scale in m/s^2.

    Returns
    -------
    u : float or ndarray
        Dimensionless structural variable.
    """
    g_bar = np.asarray(g_bar, dtype=float)
    return np.sqrt(np.maximum(g_bar, 0.0) / a0)


def P_tvgd(g_bar, a0=A0_DEFAULT):
    """
    TVGD structural response function.

    P(u) = 1 - exp(-u)

    with

    u = sqrt(g_bar / a0)

    Parameters
    ----------
    g_bar : float or array-like
        Baryonic/Newtonian acceleration in m/s^2.
    a0 : float
        TVGD acceleration scale in m/s^2.

    Returns
    -------
    P : float or ndarray
        TVGD response function.
    """
    u = tvgd_u(g_bar, a0=a0)
    return 1.0 - np.exp(-u)


def g_tvgd(g_bar, a0=A0_DEFAULT):
    """
    Effective TVGD acceleration.

    g_TVGD = g_bar / P(u)

    Parameters
    ----------
    g_bar : float or array-like
        Baryonic/Newtonian acceleration in m/s^2.
    a0 : float
        TVGD acceleration scale in m/s^2.

    Returns
    -------
    g_eff : float or ndarray
        Effective TVGD acceleration in m/s^2.
    """
    g_bar = np.asarray(g_bar, dtype=float)
    P = P_tvgd(g_bar, a0=a0)

    return g_bar / np.maximum(P, 1e-300)


def g_tvgd_approx(g_bar, a0=A0_DEFAULT):
    """
    Approximate RAR-like TVGD form.

    g_TVGD approx = g_bar + sqrt(g_bar * a0)

    This is useful for simple comparisons, but the primary tested
    structural form is g_bar / P(u).
    """
    g_bar = np.asarray(g_bar, dtype=float)
    return g_bar + np.sqrt(np.maximum(g_bar, 0.0) * a0)


def fractional_difference(model, reference):
    """
    Fractional difference between a model and a reference.

    delta = (model - reference) / reference
    """
    model = np.asarray(model, dtype=float)
    reference = np.asarray(reference, dtype=float)

    return (model - reference) / np.maximum(reference, 1e-300)


def percent_difference(model, reference):
    """
    Percent difference between a model and a reference.
    """
    return 100.0 * fractional_difference(model, reference)


# ============================================================
# Cosmological Background: LCDM and Effective TVGD V1.1
# ============================================================

def E_lcdm(z, Omega_m=0.31339978339159447, Omega_r=9e-5):
    """
    Dimensionless Hubble function for flat LCDM.

    E(z) = H(z) / H0

    Parameters
    ----------
    z : float or array-like
        Redshift.
    Omega_m : float
        Matter density parameter today.
    Omega_r : float
        Radiation density parameter today.

    Returns
    -------
    E : float or ndarray
        Dimensionless Hubble parameter.
    """
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)
    Omega_de = 1.0 - Omega_m - Omega_r

    return np.sqrt(
        Omega_m * a**(-3)
        + Omega_r * a**(-4)
        + Omega_de
    )


def Omega_psi_tvgd(
    a,
    Omega_m=0.31339978339159447,
    Omega_r=9e-5,
    lam=7.0,
    p=0.55,
):
    """
    Effective TVGD cosmological density component.

    Omega_psi(a) = Omega_de * [1 - exp(-lambda a^p)] / [1 - exp(-lambda)]

    Parameters
    ----------
    a : float or array-like
        Scale factor.
    Omega_m : float
        Matter density parameter today.
    Omega_r : float
        Radiation density parameter today.
    lam : float
        TVGD effective transition parameter.
    p : float
        TVGD effective exponent.

    Returns
    -------
    Omega_psi : float or ndarray
        Effective TVGD dark-sector-like density contribution.
    """
    a = np.asarray(a, dtype=float)
    Omega_de = 1.0 - Omega_m - Omega_r

    numerator = 1.0 - np.exp(-lam * a**p)
    denominator = 1.0 - np.exp(-lam)

    return Omega_de * numerator / denominator


def E_tvgd(
    z,
    Omega_m=0.31339978339159447,
    Omega_r=9e-5,
    lam=7.0,
    p=0.55,
):
    """
    Dimensionless Hubble function for the effective TVGD V1.1 background.

    E(z) = H(z) / H0
    """
    z = np.asarray(z, dtype=float)
    a = 1.0 / (1.0 + z)

    return np.sqrt(
        Omega_m * a**(-3)
        + Omega_r * a**(-4)
        + Omega_psi_tvgd(
            a,
            Omega_m=Omega_m,
            Omega_r=Omega_r,
            lam=lam,
            p=p,
        )
    )


def H_of_z(z, H0=67.4, model="lcdm", **kwargs):
    """
    Hubble function H(z) in km/s/Mpc.

    Parameters
    ----------
    z : float or array-like
        Redshift.
    H0 : float
        Hubble constant in km/s/Mpc.
    model : str
        Either 'lcdm' or 'tvgd'.

    Returns
    -------
    H : float or ndarray
        Hubble parameter in km/s/Mpc.
    """
    if model.lower() == "lcdm":
        return H0 * E_lcdm(z, **kwargs)

    if model.lower() == "tvgd":
        return H0 * E_tvgd(z, **kwargs)

    raise ValueError("model must be either 'lcdm' or 'tvgd'")


def comoving_distance(z_values, E_function, H0=67.4, n_grid=2000):
    """
    Comoving distance in Mpc using simple trapezoidal integration.

    D_C(z) = c/H0 * integral_0^z dz'/E(z')
    """
    z_values = np.asarray(z_values, dtype=float)
    distances = np.zeros_like(z_values, dtype=float)

    for i, z in enumerate(z_values):
        zz = np.linspace(0.0, z, n_grid)
        integrand = 1.0 / E_function(zz)
        distances[i] = (C_LIGHT_KM / H0) * np.trapezoid(integrand, zz)

    return distances


def luminosity_distance(z_values, E_function, H0=67.4, n_grid=2000):
    """
    Luminosity distance in Mpc.

    D_L = (1 + z) D_C
    """
    z_values = np.asarray(z_values, dtype=float)
    return (1.0 + z_values) * comoving_distance(
        z_values,
        E_function,
        H0=H0,
        n_grid=n_grid,
    )


def distance_modulus(z_values, E_function, H0=67.4, n_grid=2000):
    """
    Distance modulus.

    mu = 5 log10(D_L / Mpc) + 25
    """
    D_L = luminosity_distance(
        z_values,
        E_function,
        H0=H0,
        n_grid=n_grid,
    )

    return 5.0 * np.log10(np.maximum(D_L, 1e-300)) + 25.0


# ============================================================
# Summary metrics
# ============================================================

def summary_absolute(values):
    """
    Simple summary for absolute values.

    Returns maximum and mean absolute value.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return {
            "max_abs": np.nan,
            "mean_abs": np.nan,
        }

    return {
        "max_abs": float(np.max(np.abs(values))),
        "mean_abs": float(np.mean(np.abs(values))),
    }


def log10_residual(observed, model):
    """
    Logarithmic residual:

    log10(model) - log10(observed)
    """
    observed = np.asarray(observed, dtype=float)
    model = np.asarray(model, dtype=float)

    result = np.full_like(observed, np.nan, dtype=float)

    mask = (
        np.isfinite(observed)
        & np.isfinite(model)
        & (observed > 0)
        & (model > 0)
    )

    result[mask] = np.log10(model[mask]) - np.log10(observed[mask])

    return result


def log_metrics(observed, model):
    """
    Basic log-space metrics.
    """
    residual = log10_residual(observed, model)
    valid = np.isfinite(residual)

    if valid.sum() == 0:
        return {
            "N": 0,
            "mean_log_residual": np.nan,
            "median_log_residual": np.nan,
            "mean_abs_log_residual": np.nan,
            "rmse_log": np.nan,
        }

    r = residual[valid]

    return {
        "N": int(valid.sum()),
        "mean_log_residual": float(np.mean(r)),
        "median_log_residual": float(np.median(r)),
        "mean_abs_log_residual": float(np.mean(np.abs(r))),
        "rmse_log": float(np.sqrt(np.mean(r**2))),
    }


# ============================================================
# Self-test
# ============================================================

if __name__ == "__main__":
    print("TVGD core module loaded successfully.")

    g_test = np.array([9.8, 1e-10, 1e-12])
    print("g_bar:", g_test)
    print("u:", tvgd_u(g_test))
    print("P(u):", P_tvgd(g_test))
    print("g_TVGD:", g_tvgd(g_test))

    z_test = np.array([0.1, 0.5, 1.0])
    print("E_LCDM:", E_lcdm(z_test))
    print("E_TVGD:", E_tvgd(z_test))
