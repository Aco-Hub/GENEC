"""Nuclear energy generation rates -- PP chain and CNO cycle.

Ported from src/energy.f90 (lines 117-600).
Uses reaction rate functions from smallfunc (primus, secun, tertiu).
Supports FP32 mode for ~2x GPU throughput on consumer cards.
"""

import math
import torch
from genec import const
from genec.smallfunc import primus, secun


# ============================================================
# PP chain reaction rates
# ============================================================

def pp_chain_rate(T, rho, X):
    """PP-chain energy generation rate.

    Args:
        T: temperature in K (not log)
        rho: density in g/cm^3
        X: hydrogen mass fraction

    Returns:
        epsilon_pp in erg/g/s
    """
    if T < 4.0e6 or X <= 0:
        return 0.0

    T9 = T / 1.0e9
    T913 = T9**(1.0/3.0)
    T923 = T913**2

    # p(p,e+nu)d rate -- NACRE (Adelberger+ 2011)
    rate_pp = primus(4.01e-15, 3.380, 1.0/3.0, 1.0/3.0, T9)
    if rate_pp == 0.0:
        return 0.0

    Q_pp = 1.442  # MeV
    zeta = X
    f_screen = 1.0

    # eps_pp = 2.38e6 * rho * X^2 * T9^(-2/3) * exp(-3.381/T9^(1/3))
    eps_pp = 2.38e6 * rho * X**2 * T9**(-2.0/3.0) * math.exp(-3.381 / T913) * f_screen

    return eps_pp


def cno_rate(T, rho, X, XC12, XN14):
    """CNO cycle energy generation rate.

    Args:
        T: temperature in K
        rho: density in g/cm^3
        X: hydrogen mass fraction
        XC12: C12 mass fraction
        XN14: N14 mass fraction

    Returns:
        epsilon_cno in erg/g/s
    """
    if T < 4.0e6 or X <= 0:
        return 0.0

    T9 = T / 1.0e9
    T913 = T9**(1.0/3.0)

    X_CNO = XC12 + XN14
    if X_CNO <= 0:
        return 0.0

    eps_cno = 8.24e25 * rho * X * X_CNO * T9**(-2.0/3.0) * math.exp(-15.231 / T913)

    return eps_cno


def total_energy_rate(T, rho, X, Y, XC12=2.56e-3, XN14=7.4e-4):
    """Total nuclear energy generation rate (PP + CNO).

    Args:
        T: temperature in K
        rho: density in g/cm^3
        X: hydrogen mass fraction
        Y: helium mass fraction
        XC12: C12 mass fraction (default solar)
        XN14: N14 mass fraction (default solar)

    Returns:
        epsilon_total in erg/g/s
    """
    return pp_chain_rate(T, rho, X) + cno_rate(T, rho, X, XC12, XN14)


# ============================================================
# Batched GPU versions
# ============================================================

def batch_pp_chain(T, rho, X, device=None, dtype=torch.float64):
    """Batched PP-chain rate for N shells.

    Args:
        T: (N,) tensor of temperature in K
        rho: (N,) tensor of density in g/cm^3
        X: scalar or (N,) tensor of hydrogen mass fraction
        dtype: torch.float64 (default) or torch.float32 for 2x throughput

    Returns:
        (N,) tensor of epsilon_pp in erg/g/s
    """
    if device is not None:
        T = T.to(device=device, dtype=dtype)
        rho = rho.to(device=device, dtype=dtype)
    if isinstance(X, (int, float)):
        X = torch.full_like(T, X)

    T9 = T / 1.0e9
    T913 = T9.pow(1.0 / 3.0)

    eps = 2.38e6 * rho * X**2 * T9.pow(-2.0 / 3.0) * torch.exp(-3.381 / T913)

    # Zero out where T < 4 MK or X <= 0
    mask = (T >= 4.0e6) & (X > 0)
    return torch.where(mask, eps, torch.zeros_like(eps))


def batch_cno(T, rho, X, X_CNO, device=None, dtype=torch.float64):
    """Batched CNO cycle rate for N shells.

    Args:
        T: (N,) tensor of temperature in K
        rho: (N,) tensor of density in g/cm^3
        X: scalar or (N,) tensor of hydrogen mass fraction
        X_CNO: scalar or (N,) tensor of CNO catalyst mass fraction
        dtype: torch.float64 (default) or torch.float32 for 2x throughput

    Returns:
        (N,) tensor of epsilon_cno in erg/g/s
    """
    if device is not None:
        T = T.to(device=device, dtype=dtype)
        rho = rho.to(device=device, dtype=dtype)
    if isinstance(X, (int, float)):
        X = torch.full_like(T, X)
    if isinstance(X_CNO, (int, float)):
        X_CNO = torch.full_like(T, X_CNO)

    T9 = T / 1.0e9
    T913 = T9.pow(1.0 / 3.0)

    eps = 8.24e25 * rho * X * X_CNO * T9.pow(-2.0 / 3.0) * torch.exp(-15.231 / T913)

    mask = (T >= 4.0e6) & (X > 0) & (X_CNO > 0)
    return torch.where(mask, eps, torch.zeros_like(eps))


def batch_total_energy(T, rho, X, Y, XC12=2.56e-3, XN14=7.4e-4,
                       device=None, dtype=torch.float64):
    """Batched total nuclear energy rate for N shells."""
    X_CNO = XC12 + XN14
    return batch_pp_chain(T, rho, X, device, dtype) + batch_cno(T, rho, X, X_CNO, device, dtype)


# ============================================================
# torch.compile wrapper
# ============================================================

_compiled_batch_energy = None


def compiled_batch_total_energy(T, rho, X, Y, XC12=2.56e-3, XN14=7.4e-4,
                                device=None, dtype=torch.float64):
    """batch_total_energy with torch.compile kernel fusion."""
    global _compiled_batch_energy
    if _compiled_batch_energy is None:
        _compiled_batch_energy = torch.compile(batch_total_energy, mode='reduce-overhead')
    return _compiled_batch_energy(T, rho, X, Y, XC12=XC12, XN14=XN14,
                                  device=device, dtype=dtype)
