"""Small mathematical and physics functions -- ported from src/SmallFunc.f90.

Each function has a scalar version (matching Fortran exactly)
and a batched GPU-accelerated version using PyTorch tensors.
"""

import torch
import math
from genec.interpolation import indic as _indic, batch_indic as _batch_indic


# ============================================================
# Scalar versions
# ============================================================

def neg_root(a, b):
    """Signed power function: sign(a) * |a|^b."""
    c = abs(a) ** b
    return -c if a < 0.0 else c


def expf10(x):
    """Compute 10^x."""
    vloge = math.log10(math.e)
    return math.exp(x / vloge)


def exphi(x):
    """Numerically stable 1 - exp(x). Uses Taylor series for small |x|."""
    if abs(x) <= 3.0e-2:
        return -((x * (1.0 / 3.0) + 1.0) * x * 0.5 + 1.0) * x
    else:
        return 1.0 - math.exp(x)


def primus(a, b, e1, e2, t9):
    """Reaction rate: a / t9^e1 * exp(-b / t9^e2) with underflow protection."""
    exponent = -b / t9**e2
    if exponent > -706.0:
        return a / t9**e1 * math.exp(exponent)
    else:
        return 0.0


def drimus(b, e1, e2, t9):
    """Derivative of log reaction rate: -e1 + e2 * b / t9^e2."""
    return -e1 + e2 * b / t9**e2


def secun(a, b, c, e1, e2, e3, t9):
    """Extended reaction rate: a / t9^e1 * exp(-b/t9^e2 - (t9/c)^e3)."""
    return a / t9**e1 * math.exp(-b / t9**e2 - (t9 / c)**e3)


def dsecun(b, c, e1, e2, e3, t9):
    """Derivative of extended reaction rate."""
    return -e1 + e2 * b / t9**e2 - e3 * (t9 / c)**e3


def tertiu(a, b, c, d, e, f, e1, e2, e3, e4, e5, t9):
    """Polynomial reaction rate: a + b*t9^e1 + c*t9^e2 + d*t9^e3 + e*t9^e4 + f*t9^e5."""
    return a + b * t9**e1 + c * t9**e2 + d * t9**e3 + e * t9**e4 + f * t9**e5


def dterti(b, c, d, e, f, e1, e2, e3, e4, e5, t9):
    """Derivative of polynomial reaction rate."""
    return e1 * b * t9**e1 + e2 * c * t9**e2 + e3 * d * t9**e3 + e4 * e * t9**e4 + e5 * f * t9**e5


def pos(x0, x, m):
    """Binary search (same as indic)."""
    return _indic(x0, x, m)


def ValInterp(x1, x2, y1, y2, y3):
    """2D interpolation: (x1*(y2-y3) + x2*(y3-y1)) / (y2-y1)."""
    return (x1 * (y2 - y3) + x2 * (y3 - y1)) / (y2 - y1)


def girl(a_flat, n, m):
    """Matrix inversion by Gauss elimination.

    Args:
        a_flat: flat list/array of size n*(n+m) representing [A|I] in column-major order
        n: matrix dimension
        m: number of right-hand side columns

    Returns:
        (b, flag) where b is flat list of size n*m, flag is 0 on success
    """
    npm = n + m
    a = list(a_flat)  # work copy
    flag = 0

    for j in range(1, n + 1):
        nj = (j - 1) * n
        jj = nj + j - 1  # 0-indexed
        j1 = j + 1

        amax = abs(a[jj])
        jm = j

        if j1 - 1 <= n - 1:
            for i in range(j1, n + 1):
                ij = nj + i - 1  # 0-indexed
                if abs(a[ij]) > amax:
                    amax = abs(a[ij])
                    jm = i

            if jm < j:
                flag = 1
                return [0.0] * (n * m), flag

            if jm > j:
                # Partial pivoting
                i1 = jm - 1 + nj  # 0-indexed
                i2 = jj
                for ii in range(j, npm + 1):
                    a[i1], a[i2] = a[i2], a[i1]
                    i1 += n
                    i2 += n

        if a[jj] == 0.0:
            flag = 2
            return [0.0] * (n * m), flag

        # Elimination
        for i in range(1, n + 1):
            if i != j:
                ij = nj + i - 1  # 0-indexed
                faktor = -a[ij] / a[jj]
                jk = jj
                ik = ij
                for k in range(j1, npm + 1):
                    jk += n
                    ik += n
                    a[ik] = a[ik] + faktor * a[jk]

        # Division of pivot row
        jk = jj
        faktor = 1.0 / a[jj]
        for k in range(j1, npm + 1):
            jk += n
            a[jk] = a[jk] * faktor

    # Extract result
    b = [0.0] * (n * m)
    i2 = n * n
    for i in range(n * m):
        b[i] = a[i2]
        i2 += 1

    return b, flag


# ============================================================
# Batched GPU-accelerated versions
# ============================================================

def batch_neg_root(a, b, device=None):
    """Batched signed power: sign(a) * |a|^b."""
    if device is not None:
        a = a.to(device)
        b = b.to(device)
    return torch.sign(a) * torch.abs(a).pow(b)


def batch_expf10(x, device=None):
    """Batched 10^x."""
    if device is not None:
        x = x.to(device)
    return torch.pow(10.0, x)


def batch_exphi(x, device=None):
    """Batched numerically stable 1 - exp(x)."""
    if device is not None:
        x = x.to(device)
    taylor = -((x * (1.0 / 3.0) + 1.0) * x * 0.5 + 1.0) * x
    standard = 1.0 - torch.exp(x)
    return torch.where(torch.abs(x) <= 3.0e-2, taylor, standard)


def batch_primus(a, b, e1, e2, t9, device=None):
    """Batched reaction rate with underflow protection."""
    if device is not None:
        t9 = t9.to(device)
    exponent = -b / t9.pow(e2)
    result = a / t9.pow(e1) * torch.exp(exponent)
    return torch.where(exponent > -706.0, result, torch.zeros_like(result))


def batch_drimus(b, e1, e2, t9, device=None):
    """Batched derivative of log reaction rate."""
    if device is not None:
        t9 = t9.to(device)
    return -e1 + e2 * b / t9.pow(e2)


def batch_secun(a, b, c, e1, e2, e3, t9, device=None):
    """Batched extended reaction rate."""
    if device is not None:
        t9 = t9.to(device)
    return a / t9.pow(e1) * torch.exp(-b / t9.pow(e2) - (t9 / c).pow(e3))


def batch_dsecun(b, c, e1, e2, e3, t9, device=None):
    """Batched derivative of extended reaction rate."""
    if device is not None:
        t9 = t9.to(device)
    return -e1 + e2 * b / t9.pow(e2) - e3 * (t9 / c).pow(e3)


def batch_tertiu(a, b, c, d, e, f, e1, e2, e3, e4, e5, t9, device=None):
    """Batched polynomial reaction rate."""
    if device is not None:
        t9 = t9.to(device)
    return a + b * t9.pow(e1) + c * t9.pow(e2) + d * t9.pow(e3) + e * t9.pow(e4) + f * t9.pow(e5)


def batch_dterti(b, c, d, e, f, e1, e2, e3, e4, e5, t9, device=None):
    """Batched derivative of polynomial reaction rate."""
    if device is not None:
        t9 = t9.to(device)
    return e1 * b * t9.pow(e1) + e2 * c * t9.pow(e2) + e3 * d * t9.pow(e3) + e4 * e * t9.pow(e4) + e5 * f * t9.pow(e5)


def batch_ValInterp(x1, x2, y1, y2, y3, device=None):
    """Batched 2D interpolation."""
    return (x1 * (y2 - y3) + x2 * (y3 - y1)) / (y2 - y1)


def batch_girl(matrices, device=None):
    """Batched matrix inversion using torch.linalg.inv.

    Args:
        matrices: (N, k, k) tensor of matrices to invert

    Returns:
        (N, k, k) tensor of inverses
    """
    if device is not None:
        matrices = matrices.to(device)
    return torch.linalg.inv(matrices)
