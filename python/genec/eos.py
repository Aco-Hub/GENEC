"""Equation of State -- ideal gas + radiation (ieos=0 path from src/EOS.f90).

Computes thermodynamic properties from (P, T, composition).
Supports FP32 mode for ~2x GPU throughput on consumer cards.
"""

import math
import torch
from genec import const


def eos_ideal(log_p, log_t, mu):
    """Ideal gas + radiation EOS.

    Args:
        log_p: log10(total pressure) in cgs (dyne/cm^2)
        log_t: log10(temperature) in K
        mu: mean molecular weight

    Returns:
        dict with: log_rho, beta, Cv, dlnP_dlnrho_T, dlnP_dlnT_rho, nabla_ad
    """
    P = 10.0**log_p
    T = 10.0**log_t

    # Radiation pressure: P_rad = (a/3) * T^4
    P_rad = (const.cst_a / 3.0) * T**4

    # Gas pressure
    P_gas = P - P_rad
    P_gas = max(P_gas, P * 1.0e-10)  # Numerical safety

    # Beta = P_gas / P_total
    beta = P_gas / P

    # Density from ideal gas: P_gas = rho * k_B * T / (mu * m_H)
    # rho = P_gas * mu * m_H / (k_B * T)
    rho = P_gas * mu * const.cst_mh / (const.cst_k * T)
    log_rho = math.log10(rho) if rho > 0 else -30.0

    # Thermodynamic derivatives
    # (dlnP/dlnrho)_T = P_gas/P = beta
    dlnP_dlnrho_T = beta

    # (dlnP/dlnT)_rho = (P_gas + 4*P_rad) / P = beta + 4*(1-beta)= 4 - 3*beta
    dlnP_dlnT_rho = 4.0 - 3.0 * beta

    # Adiabatic gradient
    delta = dlnP_dlnT_rho / dlnP_dlnrho_T if dlnP_dlnrho_T > 0 else 1.0
    alpha = 1.0 / dlnP_dlnrho_T if dlnP_dlnrho_T > 0 else 1.0

    # Cv per unit mass (ideal gas + radiation)
    Cv = 1.5 * const.cst_k / (mu * const.cst_mh)
    if rho > 0:
        Cv += 4.0 * const.cst_a * T**3 / rho

    # nabla_ad for ideal gas + radiation
    denom = 24.0 * (1.0 - beta)**2 + 3.0 * (2.0 - beta) * beta**2
    gamma1 = (32.0 - 24.0 * beta - 3.0 * beta**2) * (1.0 - beta) + beta**2 * 5.0 / 2.0
    gamma1 = gamma1 / max(denom / beta + 2.5 * beta, 1e-30) if beta > 0 else 5.0/3.0
    nabla_ad = 2.0 / 5.0 * beta / max(beta, 1e-10)

    return {
        'log_rho': log_rho,
        'beta': beta,
        'Cv': Cv,
        'dlnP_dlnrho_T': dlnP_dlnrho_T,
        'dlnP_dlnT_rho': dlnP_dlnT_rho,
        'nabla_ad': nabla_ad,
        'gamma1': gamma1,
    }


def batch_eos_ideal(log_p, log_t, mu, device=None, dtype=torch.float64):
    """Batched GPU EOS for N shells.

    Args:
        log_p: (N,) tensor of log10(P)
        log_t: (N,) tensor of log10(T)
        mu: (N,) tensor or scalar of mean molecular weight
        dtype: torch.float64 (default) or torch.float32 for 2x throughput

    Returns:
        dict of (N,) tensors
    """
    if device is not None:
        log_p = log_p.to(device=device, dtype=dtype)
        log_t = log_t.to(device=device, dtype=dtype)
        if isinstance(mu, torch.Tensor):
            mu = mu.to(device=device, dtype=dtype)

    P = torch.pow(10.0, log_p)
    T = torch.pow(10.0, log_t)

    P_rad = (const.cst_a / 3.0) * T.pow(4)
    P_gas = (P - P_rad).clamp(min=P * 1.0e-10)
    beta = P_gas / P

    rho = P_gas * mu * const.cst_mh / (const.cst_k * T)
    log_rho = torch.log10(rho.clamp(min=1e-30))

    dlnP_dlnrho_T = beta
    dlnP_dlnT_rho = 4.0 - 3.0 * beta

    Cv = 1.5 * const.cst_k / (mu * const.cst_mh)
    Cv = Cv + 4.0 * const.cst_a * T.pow(3) / rho.clamp(min=1e-30)

    gamma1 = beta + (4.0 - 3.0 * beta)**2 * beta / \
             (8.0 * (4.0 - 3.0 * beta) * (1.0 - beta) + beta * (32.0 - 24.0 * beta - 3.0 * beta**2)).clamp(min=1e-30)
    nabla_ad = 1.0 - 1.0 / gamma1.clamp(min=1e-30)

    return {
        'log_rho': log_rho,
        'beta': beta,
        'Cv': Cv,
        'dlnP_dlnrho_T': dlnP_dlnrho_T,
        'dlnP_dlnT_rho': dlnP_dlnT_rho,
        'nabla_ad': nabla_ad,
        'gamma1': gamma1,
    }


# ============================================================
# torch.compile wrapper
# ============================================================

_compiled_batch_eos = None


def compiled_batch_eos_ideal(log_p, log_t, mu, device=None, dtype=torch.float64):
    """batch_eos_ideal with torch.compile kernel fusion."""
    global _compiled_batch_eos
    if _compiled_batch_eos is None:
        _compiled_batch_eos = torch.compile(batch_eos_ideal, mode='reduce-overhead')
    return _compiled_batch_eos(log_p, log_t, mu, device=device, dtype=dtype)
