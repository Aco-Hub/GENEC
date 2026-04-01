"""Ionisation module -- full port of src/ionisation.f90.

Implements the Saha equation solver for partial ionization of
H, He, C, O, Ne, Mg (6 elements, up to 12 ionization states).

Provides both scalar (matching Fortran) and batched GPU versions.
"""

import math
import torch
from genec import const

# ============================================================
# Atomic data (from ionisation.f90 lines 18-40)
# ============================================================

iatoms = 6
ionstates = 12

# Atomic weights (g/mol)
a_ion = [1.00794, 4.002602, 12.0107, 15.9994, 20.1797, 24.305]

# Atomic numbers
iz = [1, 2, 6, 8, 10, 12]

# Ionization potentials in eV -- chi1 flat array reshaped to (iatoms, ionstates)
# Fortran stores as chi(iatoms, 0:ionstates-1) column-major
_chi1 = [
    13.59844, 24.58741, 11.26030, 13.61806, 21.5646, 7.64624,
    0.0, 54.41778, 24.38332, 35.1173, 40.96328, 15.03528,
    0.0, 0.0, 47.8878, 54.9355, 63.45, 80.1437,
    0.0, 0.0, 64.4939, 77.41353, 97.12, 109.2655,
    0.0, 0.0, 392.0870, 113.8990, 126.21, 141.27,
    0.0, 0.0, 489.99334, 138.1197, 157.93, 186.76,
    0.0, 0.0, 0.0, 739.29, 207.2759, 225.02,
    0.0, 0.0, 0.0, 871.4101, 239.0989, 265.96,
    0.0, 0.0, 0.0, 0.0, 1195.8286, 328.06,
    0.0, 0.0, 0.0, 0.0, 1362.1995, 367.5,
    0.0, 0.0, 0.0, 0.0, 0.0, 1761.805,
    0.0, 0.0, 0.0, 0.0, 0.0, 1962.665,
]
# Reshape: Fortran (iatoms, ionstates) column-major → chi[atom][state]
chi = [[_chi1[j * iatoms + i] for j in range(ionstates)] for i in range(iatoms)]

# Log partition function ratios
_vlu1 = [
    0.0, 0.602059991, 0.124938737, -0.051152522, 1.079181246, 0.602059991,
    0.0, 0.0, -0.477121255, 0.653212514, 0.477121255, 0.0,
    0.0, 0.0, 0.602059991, 0.124938737, -0.051152522, 1.079181246,
    0.0, 0.0, 0.0, -0.477121255, 0.653212514, 0.477121255,
    0.0, 0.0, 0.602059991, 0.602059991, 0.124938737, -0.051152522,
    0.0, 0.0, 0.0, 0.0, -0.477121255, 0.653212514,
    0.0, 0.0, 0.0, 0.602059991, 0.602059991, 0.124938737,
    0.0, 0.0, 0.0, 0.0, 0.0, -0.477121255,
    0.0, 0.0, 0.0, 0.0, 0.602059991, 0.602059991,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.602059991,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
]
vlu = [[_vlu1[j * iatoms + i] for j in range(ionstates)] for i in range(iatoms)]


# ============================================================
# Composition setup
# ============================================================

def setup_composition(abond_in):
    """Initialize composition state from mass fractions.

    Args:
        abond_in: list of 6 mass fractions [X_H, X_He, X_C, X_O, X_Ne, X_Mg]

    Returns:
        dict with keys: abond, vnu, vmol, vmyion, active_list
    """
    abond = list(abond_in)
    vmol = sum(ab / ai for ab, ai in zip(abond, a_ion))
    # Add contribution from "other" elements (Z metals not in the 6)
    vmol += 0.5 * (1.0 - sum(abond))
    vmol = 1.0 / vmol if vmol > 0 else 0.0

    # Active elements list (Fortran: elements with abond > 0.02)
    active_list = []
    for i in range(iatoms):
        if abond[i] > 0.02:
            active_list.append(i)

    vnu = [0.0] * iatoms
    for i in active_list:
        vnu[i] = vmol * abond[i] / a_ion[i]

    # Mean molecular weight if fully ionized
    e_full = sum(iz[i] * vnu[i] for i in active_list)
    vmyion = vmol / (1.0 + e_full) if (1.0 + e_full) > 0 else 0.0

    return {
        'abond': abond,
        'vnu': vnu,
        'vmol': vmol,
        'vmyion': vmyion,
        'active_list': active_list,
    }


# ============================================================
# Scalar Saha solver (matching Fortran exactly)
# ============================================================

def saha(p, t, beta, e_init, vnu_arr, active_list):
    """Iterative Saha equation solver.

    Args:
        p: log10(Pressure) in cgs
        t: log10(Temperature) in K
        beta: Pgas/Ptotal
        e_init: initial estimate of free electrons per atom
        vnu_arr: relative number of atoms per element
        active_list: indices of active elements

    Returns:
        (xion, e) where xion[i][j] is ionization fraction of element i in state j,
        and e is the converged electron count.
    """
    thet = const.cst_e / (1.0e-7 * const.cst_k * const.um * 10.0**t)

    h = (2.5 * t - p - math.log10(beta)
         + 1.5 * (math.log10(2.0) + const.lgpi + const.cstlg_me)
         + 2.5 * const.cstlg_k - 3.0 * const.cstlg_h)

    xion = [[0.0] * (ionstates + 1) for _ in range(iatoms)]
    e = max(e_init, 1.0e-10)

    for iteration in range(41):
        diff = 0.0
        enew = 0.0
        e = max(e, 1.0e-10)

        for idx in active_list:
            i = idx
            max_state = iz[i]

            # Compute Saha ratios
            vlk = [0.0] * ionstates
            for j in range(max_state):
                vlk[j] = h - chi[i][j] * thet + vlu[i][j] + math.log10((1.0 + e) / e)

            # Compute ionization fractions
            for j in range(max_state + 1):
                xi_j = 1.0
                skip = False

                # Sum down: lower ionization states
                temp = 0.0
                for k in range(j - 1, -1, -1):
                    temp -= vlk[k]
                    if temp > 12.0:
                        xi_j = 0.0
                        skip = True
                        break
                    else:
                        xi_j += 10.0**temp

                if skip:
                    diff += abs(xion[i][j] - xi_j)
                    xion[i][j] = xi_j
                    continue

                # Sum up: higher ionization states
                temp = 0.0
                for k in range(j, max_state):
                    temp += vlk[k]
                    if temp > 12.0:
                        xi_j = 0.0
                        skip = True
                        break
                    else:
                        xi_j += 10.0**temp

                if skip:
                    diff += abs(xion[i][j] - xi_j)
                    xion[i][j] = xi_j
                    continue

                xi_j = 1.0 / xi_j
                diff += abs(xion[i][j] - xi_j)
                xion[i][j] = xi_j

            # Electron contribution from this element
            temp = sum(k * xion[i][k] for k in range(1, max_state + 1))
            enew += vnu_arr[i] * temp

        e = enew
        if diff <= 5.0e-3:
            break

    return xion, e


def ionpart(p, t, comp, vmion_prev=0.0, ionized=0, chem=0.0, ychem=0.0):
    """Full ionization calculation.

    Args:
        p: log10(Pressure) in cgs
        t: log10(Temperature) in K
        comp: composition dict from setup_composition()
        vmion_prev: previous mean molecular weight (0 for first call)
        ionized: 0=partial, 1=fully ionized
        chem: sum of heavy element mass fractions
        ychem: helium mass fraction

    Returns:
        dict with: xion, vmion, beta_env, cp, vna, vmionp, vmiont,
                   HII, HeII, HeIII, mu
    """
    vmol = comp['vmol']
    vnu_arr = comp['vnu']
    active_list = comp['active_list']

    # Initial electron estimate
    if vmion_prev != 0.0 and vmion_prev != vmol:
        e = vmol / vmion_prev - 1.0
    else:
        e = 0.0
        for i in active_list:
            if i >= 1:  # He or heavier
                e = vnu_arr[i] * comp['abond'][i] / a_ion[i]
            if e != 0.0:
                break

    # Beta_env = Pgas/Ptotal
    beta_env = 1.0 - 10.0**(const.cstlg_a - math.log10(3.0) + 4.0 * t - p)
    beta_env = max(beta_env, 1.0e-4)

    # Initialize xion
    xion = [[0.0] * (ionstates + 1) for _ in range(iatoms)]

    # Saha calculation if not fully ionized
    if ionized == 0:
        xion, e = saha(p, t, beta_env, e, vnu_arr, active_list)

    # Check pressure ionization thresholds
    fully_ionized = False
    if chem + ychem < 0.95:
        if (p - t + math.log10(beta_env) - math.log10(1.0 + e) > 5.3447) or ionized == 1:
            fully_ionized = True
    else:
        if (t - 0.123 * p > 4.09) or ionized == 1:
            fully_ionized = True

    if fully_ionized:
        e = 0.0
        for i in active_list:
            for j in range(iz[i]):
                xion[i][j] = 0.0
            xion[i][iz[i]] = 1.0
            e += iz[i] * vnu_arr[i]

    # Thermodynamic derivatives
    d = 1.0 + e
    h3 = 4.0 * (1.0 - beta_env) / beta_env + 2.5
    h4 = (1.0 - beta_env) * (4.0 + beta_env) * d / (beta_env * beta_env)

    vng = 0.0
    vngp = 0.0
    vngpp = 0.0

    if not fully_ionized and ionized != 1:
        kB_SI = const.cst_k * 1.0e-7 / const.cst_e  # k_B in eV/K

        for idx in active_list:
            i = idx
            max_z = iz[i]

            # Find dominant ionization transition
            sums = [xion[i][j] + xion[i][j + 1] for j in range(max_z)]
            k = 0
            temp_max = 0.0
            for j in range(max_z):
                if sums[j] > temp_max:
                    temp_max = sums[j]
                    k = j

            # Compute cumulative y values
            y = [0.0] * ionstates
            l = 0 if k == 0 else 1
            temp = 0.0
            for j in range(max_z - 1, k - l - 1, -1):
                temp += xion[i][j + 1]
                y[j] = temp

            # Compute gi and phi for dominant term
            if k == 0:
                temp = y[0] * (1.0 - y[0])
                denom = d * e + vnu_arr[i] * temp
                gi = d * e * temp / denom if denom != 0.0 else 0.0
                phi = h3 + chi[i][0] / (10.0**t) / kB_SI
            else:
                temp = y[k] * (y[k - 1] - y[k])
                denom = d * e * y[k - 1] + vnu_arr[i] * temp
                gi = d * e * temp / denom if denom != 0.0 else 0.0
                phi = h3 + chi[i][k] / (10.0**t) / kB_SI

            h5 = vnu_arr[i] * gi
            vng += h5
            vngp += h5 * phi
            vngpp += h5 * phi * phi

    cp = 2.5 * d + 4.0 * h4 + vngpp
    vna = (d + h4 + vngp / beta_env) / cp if cp != 0.0 else 0.0
    vmionp = vng / (beta_env * d) if beta_env * d != 0.0 else 0.0
    vmiont = -vngp / d if d != 0.0 else 0.0
    vmion = vmol / d if d != 0.0 else 0.0

    return {
        'xion': xion,
        'vmion': vmion,
        'beta_env': beta_env,
        'cp': cp,
        'vna': vna,
        'vmionp': vmionp,
        'vmiont': vmiont,
        'HII': xion[0][1] if len(xion[0]) > 1 else 0.0,    # H+ fraction
        'HeII': xion[1][1] if len(xion[1]) > 1 else 0.0,   # He+ fraction
        'HeIII': xion[1][2] if len(xion[1]) > 2 else 0.0,  # He++ fraction
        'mu': vmion,
    }


# ============================================================
# Batched GPU version
# ============================================================

# Pre-compute tensors for GPU
_chi_tensor = None
_vlu_tensor = None
_iz_tensor = None
_a_ion_tensor = None


def _ensure_tensors(device, dtype=torch.float64):
    global _chi_tensor, _vlu_tensor, _iz_tensor, _a_ion_tensor
    if _chi_tensor is None or _chi_tensor.device != device:
        _chi_tensor = torch.tensor(chi, dtype=dtype, device=device)  # (6, 12)
        _vlu_tensor = torch.tensor(vlu, dtype=dtype, device=device)  # (6, 12)
        _iz_tensor = torch.tensor(iz, dtype=torch.long, device=device)  # (6,)
        _a_ion_tensor = torch.tensor(a_ion, dtype=dtype, device=device)  # (6,)


def batch_ionpart(p_batch, t_batch, abond, device=None, dtype=torch.float64):
    """Batched GPU ionization calculation for N shells.

    Args:
        p_batch: (N,) tensor of log10(Pressure)
        t_batch: (N,) tensor of log10(Temperature)
        abond: list of 6 mass fractions [X_H, X_He, X_C, X_O, X_Ne, X_Mg]
        device: torch device
        dtype: torch dtype

    Returns:
        dict with tensors: HII(N,), HeII(N,), HeIII(N,), mu(N,), beta_env(N,)
    """
    if device is None:
        device = p_batch.device
    _ensure_tensors(device, dtype)

    N = p_batch.shape[0]
    p = p_batch.to(device=device, dtype=dtype)
    t = t_batch.to(device=device, dtype=dtype)

    # Setup composition (scalar, same for all shells)
    comp = setup_composition(abond)
    vmol = comp['vmol']
    vnu_t = torch.tensor(comp['vnu'], dtype=dtype, device=device)  # (6,)
    active = torch.tensor(comp['active_list'], dtype=torch.long, device=device)

    # Beta_env
    beta_env = 1.0 - torch.pow(10.0, const.cstlg_a - math.log10(3.0) + 4.0 * t - p)
    beta_env = beta_env.clamp(min=1.0e-4)

    # Saha constants
    thet = const.cst_e / (1.0e-7 * const.cst_k * const.um * torch.pow(10.0, t))  # (N,)
    h_saha = (2.5 * t - p - torch.log10(beta_env)
              + 1.5 * (math.log10(2.0) + const.lgpi + const.cstlg_me)
              + 2.5 * const.cstlg_k - 3.0 * const.cstlg_h)  # (N,)

    # Initialize xion: (N, iatoms, ionstates+1)
    xion = torch.zeros(N, iatoms, ionstates + 1, dtype=dtype, device=device)

    # Initial electron estimate
    e = torch.full((N,), 1.0e-10, dtype=dtype, device=device)
    for i_idx in range(len(active)):
        i = active[i_idx].item()
        if i >= 1:  # He or heavier
            e_val = vnu_t[i] * abond[i] / a_ion[i]
            mask = e < 1.0e-9
            e = torch.where(mask, torch.full_like(e, max(e_val, 1e-10)), e)

    # Iterative Saha solver (all shells in parallel)
    for iteration in range(41):
        e = e.clamp(min=1.0e-10)
        diff = torch.zeros(N, dtype=dtype, device=device)
        enew = torch.zeros(N, dtype=dtype, device=device)

        for i_idx in range(len(active)):
            i = active[i_idx].item()
            max_state = iz[i]

            # Saha ratios: vlk[j] = h - chi*thet + vlu + log10((1+e)/e)
            # Shape: (N, max_state)
            log_ratio = torch.log10((1.0 + e) / e)  # (N,)
            vlk = torch.zeros(N, max_state, dtype=dtype, device=device)
            for j in range(max_state):
                vlk[:, j] = h_saha - _chi_tensor[i, j] * thet + _vlu_tensor[i, j] + log_ratio

            # Compute ionization fractions for each state
            for j in range(max_state + 1):
                xi_j = torch.ones(N, dtype=dtype, device=device)

                # Sum down
                temp_down = torch.zeros(N, dtype=dtype, device=device)
                overflow_down = torch.zeros(N, dtype=torch.bool, device=device)
                for k in range(j - 1, -1, -1):
                    temp_down = temp_down - vlk[:, k]
                    overflow_down = overflow_down | (temp_down > 12.0)
                    xi_j = torch.where(overflow_down, xi_j, xi_j + torch.pow(10.0, temp_down.clamp(max=12.0)))

                # Sum up
                temp_up = torch.zeros(N, dtype=dtype, device=device)
                overflow_up = torch.zeros(N, dtype=torch.bool, device=device)
                for k in range(j, max_state):
                    temp_up = temp_up + vlk[:, k]
                    overflow_up = overflow_up | (temp_up > 12.0)
                    xi_j = torch.where(overflow_up, xi_j, xi_j + torch.pow(10.0, temp_up.clamp(max=12.0)))

                # If any overflow occurred, set to 0
                any_overflow = overflow_down | overflow_up
                xi_j = torch.where(any_overflow, torch.zeros_like(xi_j), 1.0 / xi_j)

                diff += torch.abs(xion[:, i, j] - xi_j)
                xion[:, i, j] = xi_j

            # Electron contribution
            for k in range(1, max_state + 1):
                enew += vnu_t[i] * k * xion[:, i, k]

        e = enew

        # Check convergence (all shells)
        if (diff <= 5.0e-3).all():
            break

    # Pressure ionization check
    fully_ionized = (p - t + torch.log10(beta_env) - torch.log10(1.0 + e) > 5.3447)

    # Apply full ionization where needed
    if fully_ionized.any():
        e_full = torch.zeros(N, dtype=dtype, device=device)
        for i_idx in range(len(active)):
            i = active[i_idx].item()
            # Set all states to 0 except fully ionized
            for j in range(iz[i]):
                xion[:, i, j] = torch.where(fully_ionized, torch.zeros_like(xion[:, i, j]), xion[:, i, j])
            xion[:, i, iz[i]] = torch.where(fully_ionized, torch.ones_like(xion[:, i, iz[i]]), xion[:, i, iz[i]])
            e_full += iz[i] * vnu_t[i]
        e = torch.where(fully_ionized, e_full, e)

    # Mean molecular weight
    d = 1.0 + e
    vmion = vmol / d

    return {
        'HII': xion[:, 0, 1],
        'HeII': xion[:, 1, 1],
        'HeIII': xion[:, 1, 2],
        'mu': vmion,
        'beta_env': beta_env,
        'xion': xion,
    }
