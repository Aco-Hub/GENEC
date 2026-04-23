"""Triton-fused Saha ionization kernel.

Replaces batch_ionpart's Python loop (6 elements x 12 states x 40 iterations
= ~2880 separate kernel launches) with a single fused GPU kernel.

Each Triton program handles one stellar shell, running the complete 40-iteration
Saha solver in registers.
"""

import math

import torch

try:
    import triton
    import triton.language as tl

    HAS_TRITON = True
except ImportError:
    HAS_TRITON = False

from genec import const
from genec.ionisation import setup_composition

# ============================================================
# Atomic data (duplicated from ionisation.py for kernel access)
# ============================================================

_IZ = [1, 2, 6, 8, 10, 12]  # max ionization state per element
_A_ION = [1.00794, 4.002602, 12.0107, 15.9994, 20.1797, 24.305]

# chi[i][j] ionization potentials (eV) -- (6, 12)
_CHI_FLAT = [
    13.59844,
    24.58741,
    11.26030,
    13.61806,
    21.5646,
    7.64624,
    0.0,
    54.41778,
    24.38332,
    35.1173,
    40.96328,
    15.03528,
    0.0,
    0.0,
    47.8878,
    54.9355,
    63.45,
    80.1437,
    0.0,
    0.0,
    64.4939,
    77.41353,
    97.12,
    109.2655,
    0.0,
    0.0,
    392.0870,
    113.8990,
    126.21,
    141.27,
    0.0,
    0.0,
    489.99334,
    138.1197,
    157.93,
    186.76,
    0.0,
    0.0,
    0.0,
    739.29,
    207.2759,
    225.02,
    0.0,
    0.0,
    0.0,
    871.4101,
    239.0989,
    265.96,
    0.0,
    0.0,
    0.0,
    0.0,
    1195.8286,
    328.06,
    0.0,
    0.0,
    0.0,
    0.0,
    1362.1995,
    367.5,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1761.805,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    1962.665,
]
# Reshape: Fortran col-major → chi[atom][state]
_CHI = [[_CHI_FLAT[j * 6 + i] for j in range(12)] for i in range(6)]

# vlu[i][j] log partition function ratios -- (6, 12)
_VLU_FLAT = [
    0.0,
    0.602059991,
    0.124938737,
    -0.051152522,
    1.079181246,
    0.602059991,
    0.0,
    0.0,
    -0.477121255,
    0.653212514,
    0.477121255,
    0.0,
    0.0,
    0.0,
    0.602059991,
    0.124938737,
    -0.051152522,
    1.079181246,
    0.0,
    0.0,
    0.0,
    -0.477121255,
    0.653212514,
    0.477121255,
    0.0,
    0.0,
    0.602059991,
    0.602059991,
    0.124938737,
    -0.051152522,
    0.0,
    0.0,
    0.0,
    0.0,
    -0.477121255,
    0.653212514,
    0.0,
    0.0,
    0.0,
    0.602059991,
    0.602059991,
    0.124938737,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    -0.477121255,
    0.0,
    0.0,
    0.0,
    0.0,
    0.602059991,
    0.602059991,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.602059991,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
]
_VLU = [[_VLU_FLAT[j * 6 + i] for j in range(12)] for i in range(6)]

# Constants needed inside kernel
_LOG10_E = 0.4342944819032518  # 1/ln(10)
_LN_10 = 2.302585092994046
_H_CONST = (
    1.5 * (math.log10(2.0) + const.lgpi + const.cstlg_me)
    + 2.5 * const.cstlg_k
    - 3.0 * const.cstlg_h
)
_THET_COEFF = const.cst_e / (1.0e-7 * const.cst_k * const.um)


if HAS_TRITON:

    @triton.jit
    def _log10(x):
        """log10(x) via log2."""
        return tl.log2(x) * 0.30102999566398114

    @triton.jit
    def _pow10(x):
        """10^x via exp2."""
        return tl.exp2(x * 3.321928094887362)

    @triton.jit
    def _saha_kernel(
        # Inputs
        p_ptr,
        t_ptr,
        # Composition arrays (all flattened)
        chi_ptr,  # (72,) = 6*12 flattened chi[i][j]  row-major
        vlu_ptr,  # (72,) = 6*12 flattened vlu[i][j]  row-major
        iz_ptr,  # (6,) max ionization states
        vnu_ptr,  # (6,) number fractions
        active_ptr,  # (6,) active flags (1 or 0)
        # Scalar params
        vmol: tl.constexpr,
        n_active: tl.constexpr,
        # Outputs
        hii_ptr,
        heii_ptr,
        heiii_ptr,
        mu_ptr,
        beta_ptr,
        # Grid
        N,
        BLOCK: tl.constexpr,
    ):
        """One program per shell. Full Saha iteration in registers."""
        pid = tl.program_id(0)
        idx = pid * BLOCK + tl.arange(0, BLOCK)
        mask = idx < N

        p = tl.load(p_ptr + idx, mask=mask, other=0.0)
        t = tl.load(t_ptr + idx, mask=mask, other=0.0)

        # Beta_env = 1 - 10^(cstlg_a - log10(3) + 4t - p)
        beta_exp = -14.1210378 - 0.4771212547196624 + 4.0 * t - p
        beta = 1.0 - _pow10(beta_exp)
        beta = tl.maximum(beta, 1.0e-4)

        # Saha constants
        t_lin = _pow10(t)  # 10^t
        thet = 1.602176487e-19 / (1.0e-7 * 1.3806504e-16 * 2.3025850929940457 * t_lin)
        h_saha = (
            2.5 * t
            - p
            - _log10(beta)
            + 1.5 * (0.30102999566398114 + 0.497149873 + (-27.04051106))
            + 2.5 * (-15.85991628)
            - 3.0 * (-26.17874405)
        )

        # Initialize electron count
        e = tl.zeros_like(p) + 1.0e-10

        # Load vnu for initial e estimate (use He if available)
        vnu_1 = tl.load(vnu_ptr + 1)  # He
        active_1 = tl.load(active_ptr + 1)
        e = tl.where(active_1 > 0.5, tl.zeros_like(e) + tl.maximum(vnu_1, 1.0e-10), e)

        # ---- Saha iteration (40 steps) ----
        # We handle H (iz=1) and He (iz=2) explicitly since they dominate
        # and produce HII, HeII, HeIII. For heavier elements (C,O,Ne,Mg)
        # we approximate as fully ionized in hot regions, neutral in cool.

        # xion storage: H has states 0,1; He has states 0,1,2
        xion_h0 = tl.zeros_like(p)
        xion_h1 = tl.zeros_like(p)
        xion_he0 = tl.zeros_like(p)
        xion_he1 = tl.zeros_like(p)
        xion_he2 = tl.zeros_like(p)

        vnu_0 = tl.load(vnu_ptr + 0)  # H
        active_0 = tl.load(active_ptr + 0)

        # Load chi and vlu for H and He
        chi_h0 = tl.load(chi_ptr + 0)  # chi[0][0]
        vlu_h0 = tl.load(vlu_ptr + 0)  # vlu[0][0]
        chi_he0 = tl.load(chi_ptr + 12)  # chi[1][0]
        chi_he1 = tl.load(chi_ptr + 13)  # chi[1][1]
        vlu_he0 = tl.load(vlu_ptr + 12)  # vlu[1][0]
        vlu_he1 = tl.load(vlu_ptr + 13)  # vlu[1][1]

        # Heavy element vnu and iz
        vnu_2 = tl.load(vnu_ptr + 2)  # C
        vnu_3 = tl.load(vnu_ptr + 3)  # O
        vnu_4 = tl.load(vnu_ptr + 4)  # Ne
        vnu_5 = tl.load(vnu_ptr + 5)  # Mg
        active_2 = tl.load(active_ptr + 2)
        active_3 = tl.load(active_ptr + 3)
        active_4 = tl.load(active_ptr + 4)
        active_5 = tl.load(active_ptr + 5)

        for _iter in range(41):
            e = tl.maximum(e, 1.0e-10)
            log_ratio = _log10((1.0 + e) / e)

            # ---- Hydrogen (iz=1, states 0 and 1) ----
            vlk_h0 = h_saha - chi_h0 * thet + vlu_h0 + log_ratio

            # state 0: xi = 1 / (1 + 10^vlk_h0)
            pow_h0 = _pow10(tl.minimum(vlk_h0, 12.0))
            overflow_h = vlk_h0 > 12.0
            denom_h0 = 1.0 + pow_h0
            xion_h0_new = tl.where(overflow_h, tl.zeros_like(p), 1.0 / denom_h0)

            # state 1: xi = 1 / (1 + 10^(-vlk_h0))
            pow_h0_neg = _pow10(tl.minimum(-vlk_h0, 12.0))
            overflow_h_neg = -vlk_h0 > 12.0
            denom_h1 = 1.0 + pow_h0_neg
            xion_h1_new = tl.where(overflow_h_neg, tl.zeros_like(p), 1.0 / denom_h1)

            xion_h0 = tl.where(active_0 > 0.5, xion_h0_new, xion_h0)
            xion_h1 = tl.where(active_0 > 0.5, xion_h1_new, xion_h1)

            # ---- Helium (iz=2, states 0, 1, 2) ----
            vlk_he0 = h_saha - chi_he0 * thet + vlu_he0 + log_ratio
            vlk_he1 = h_saha - chi_he1 * thet + vlu_he1 + log_ratio

            # state 0: 1 / (1 + 10^vlk_he0 + 10^(vlk_he0+vlk_he1))
            p10_he0 = _pow10(tl.minimum(vlk_he0, 12.0))
            p10_he01 = _pow10(tl.minimum(vlk_he0 + vlk_he1, 12.0))
            of_he0 = vlk_he0 > 12.0
            of_he01 = (vlk_he0 + vlk_he1) > 12.0
            denom_he0 = (
                1.0
                + tl.where(of_he0, tl.zeros_like(p), p10_he0)
                + tl.where(of_he01, tl.zeros_like(p), p10_he01)
            )
            xion_he0_new = tl.where(of_he0 | of_he01, tl.zeros_like(p), 1.0 / denom_he0)

            # state 1: 1 / (10^(-vlk_he0) + 1 + 10^vlk_he1)
            p10_mhe0 = _pow10(tl.minimum(-vlk_he0, 12.0))
            p10_he1 = _pow10(tl.minimum(vlk_he1, 12.0))
            of_mhe0 = -vlk_he0 > 12.0
            of_he1 = vlk_he1 > 12.0
            denom_he1 = (
                tl.where(of_mhe0, tl.zeros_like(p), p10_mhe0)
                + 1.0
                + tl.where(of_he1, tl.zeros_like(p), p10_he1)
            )
            xion_he1_new = tl.where(of_mhe0 | of_he1, tl.zeros_like(p), 1.0 / denom_he1)

            # state 2: 1 / (10^(-vlk_he0-vlk_he1) + 10^(-vlk_he1) + 1)
            p10_mhe01 = _pow10(tl.minimum(-vlk_he0 - vlk_he1, 12.0))
            p10_mhe1 = _pow10(tl.minimum(-vlk_he1, 12.0))
            of_mhe01 = (-vlk_he0 - vlk_he1) > 12.0
            of_mhe1 = -vlk_he1 > 12.0
            denom_he2 = (
                tl.where(of_mhe01, tl.zeros_like(p), p10_mhe01)
                + tl.where(of_mhe1, tl.zeros_like(p), p10_mhe1)
                + 1.0
            )
            xion_he2_new = tl.where(
                of_mhe01 | of_mhe1, tl.zeros_like(p), 1.0 / denom_he2
            )

            xion_he0 = tl.where(active_1 > 0.5, xion_he0_new, xion_he0)
            xion_he1 = tl.where(active_1 > 0.5, xion_he1_new, xion_he1)
            xion_he2 = tl.where(active_1 > 0.5, xion_he2_new, xion_he2)

            # Electron count: e = sum_i vnu_i * sum_j (j * xion_i_j)
            enew = tl.zeros_like(p)
            # H contribution
            enew = enew + tl.where(active_0 > 0.5, vnu_0 * xion_h1, tl.zeros_like(p))
            # He contribution
            enew = enew + tl.where(
                active_1 > 0.5, vnu_1 * (xion_he1 + 2.0 * xion_he2), tl.zeros_like(p)
            )

            # Heavy elements: approximate as fully ionized at high T
            # (they contribute iz[i] * vnu[i] electrons when fully ionized)
            hot = t > 5.0  # log(T) > 5 => T > 100,000 K
            enew = enew + tl.where(
                (active_2 > 0.5) & hot, vnu_2 * 6.0, tl.zeros_like(p)
            )
            enew = enew + tl.where(
                (active_3 > 0.5) & hot, vnu_3 * 8.0, tl.zeros_like(p)
            )
            enew = enew + tl.where(
                (active_4 > 0.5) & hot, vnu_4 * 10.0, tl.zeros_like(p)
            )
            enew = enew + tl.where(
                (active_5 > 0.5) & hot, vnu_5 * 12.0, tl.zeros_like(p)
            )

            e = enew

        # Pressure ionization check
        fully_ionized = (p - t + _log10(beta) - _log10(1.0 + e)) > 5.3447
        e_full = tl.zeros_like(p)
        e_full = e_full + tl.where(active_0 > 0.5, vnu_0 * 1.0, tl.zeros_like(p))
        e_full = e_full + tl.where(active_1 > 0.5, vnu_1 * 2.0, tl.zeros_like(p))
        e_full = e_full + tl.where(active_2 > 0.5, vnu_2 * 6.0, tl.zeros_like(p))
        e_full = e_full + tl.where(active_3 > 0.5, vnu_3 * 8.0, tl.zeros_like(p))
        e_full = e_full + tl.where(active_4 > 0.5, vnu_4 * 10.0, tl.zeros_like(p))
        e_full = e_full + tl.where(active_5 > 0.5, vnu_5 * 12.0, tl.zeros_like(p))

        xion_h1 = tl.where(fully_ionized, tl.zeros_like(p) + 1.0, xion_h1)
        xion_he1 = tl.where(fully_ionized, tl.zeros_like(p), xion_he1)
        xion_he2 = tl.where(fully_ionized, tl.zeros_like(p) + 1.0, xion_he2)
        e = tl.where(fully_ionized, e_full, e)

        # Mean molecular weight
        d = 1.0 + e
        mu_val = vmol / d

        # Store outputs
        tl.store(hii_ptr + idx, xion_h1, mask=mask)
        tl.store(heii_ptr + idx, xion_he1, mask=mask)
        tl.store(heiii_ptr + idx, xion_he2, mask=mask)
        tl.store(mu_ptr + idx, mu_val, mask=mask)
        tl.store(beta_ptr + idx, beta, mask=mask)


def triton_batch_ionpart(p_batch, t_batch, abond, device=None, dtype=torch.float64):
    """Triton-fused batch ionization -- drop-in replacement for batch_ionpart.

    Args:
        p_batch: (N,) tensor of log10(Pressure)
        t_batch: (N,) tensor of log10(Temperature)
        abond: list of 6 mass fractions [X_H, X_He, X_C, X_O, X_Ne, X_Mg]
        device: torch device
        dtype: torch dtype

    Returns:
        dict with tensors: HII(N,), HeII(N,), HeIII(N,), mu(N,), beta_env(N,)
    """
    if not HAS_TRITON:
        raise RuntimeError("Triton not available. Install with: pip install triton")

    if device is None:
        device = p_batch.device

    N = p_batch.shape[0]
    p = p_batch.to(device=device, dtype=dtype).contiguous()
    t = t_batch.to(device=device, dtype=dtype).contiguous()

    # Setup composition
    comp = setup_composition(abond)
    vmol = comp["vmol"]
    vnu = comp["vnu"]
    active_list = comp["active_list"]

    # Flatten chi and vlu to (72,) row-major: chi[i*12+j]
    chi_flat = []
    vlu_flat = []
    for i in range(6):
        for j in range(12):
            chi_flat.append(_CHI[i][j])
            vlu_flat.append(_VLU[i][j])

    chi_t = torch.tensor(chi_flat, dtype=dtype, device=device)
    vlu_t = torch.tensor(vlu_flat, dtype=dtype, device=device)
    iz_t = torch.tensor(_IZ, dtype=torch.int64, device=device)
    vnu_t = torch.tensor(vnu, dtype=dtype, device=device)

    # Active flags: 1.0 if element is active, 0.0 otherwise
    active_flags = [0.0] * 6
    for i in active_list:
        active_flags[i] = 1.0
    active_t = torch.tensor(active_flags, dtype=dtype, device=device)

    n_active = len(active_list)

    # Allocate outputs
    hii = torch.empty(N, dtype=dtype, device=device)
    heii = torch.empty(N, dtype=dtype, device=device)
    heiii = torch.empty(N, dtype=dtype, device=device)
    mu = torch.empty(N, dtype=dtype, device=device)
    beta = torch.empty(N, dtype=dtype, device=device)

    # Launch kernel
    BLOCK = 1024
    grid = ((N + BLOCK - 1) // BLOCK,)

    _saha_kernel[grid](
        p,
        t,
        chi_t,
        vlu_t,
        iz_t,
        vnu_t,
        active_t,
        vmol,
        n_active,
        hii,
        heii,
        heiii,
        mu,
        beta,
        N,
        BLOCK=BLOCK,
    )

    return {
        "HII": hii,
        "HeII": heii,
        "HeIII": heiii,
        "mu": mu,
        "beta_env": beta,
    }
