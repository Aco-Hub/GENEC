"""Interpolation functions -- ported from src/interpolation.f90.

Each function has a scalar version (matching Fortran behavior exactly)
and a batched GPU-accelerated version using PyTorch tensors.
"""

import torch
import math


# ============================================================
# Scalar versions (matching Fortran exactly)
# ============================================================

def indic(x0, x, m):
    """Binary search in a monotonic table (ascending or descending).

    Returns k such that x0 is between x[k] and x[k+1] (1-indexed).
    """
    n = m
    k = 1
    while n - k - 1 != 0:
        i = (k + n) // 2
        if (x[i - 1] - x0) * (x[n - 1] - x[0]) <= 0:
            k = i
        else:
            n = i
    return k


def fipoi(x, n, a, b):
    """Linear interpolation (function form)."""
    k = indic(x, a, n)
    return b[k - 1] + (b[k] - b[k - 1]) * (x - a[k - 1]) / (a[k] - a[k - 1])


def fipoi1(x, n, a, b):
    """Linear interpolation (returns value)."""
    k = indic(x, a, n)
    r1 = b[k - 1] + (b[k] - b[k - 1]) * (x - a[k - 1]) / (a[k] - a[k - 1])
    return r1


def fipoi2(x, n, a, b):
    """Linear interpolation with derivative. Returns (value, derivative)."""
    k = indic(x, a, n)
    r1 = b[k - 1] + (b[k] - b[k - 1]) * (x - a[k - 1]) / (a[k] - a[k - 1])
    r2 = (b[k] - b[k - 1]) / (a[k] - a[k - 1])
    return r1, r2


def quint(x, x0, h, y0, y1, y2):
    """Quadratic interpolation for equidistant points."""
    d1 = y1 - y0
    d2 = y2 - 2.0 * y1 + y0
    t = (x - x0) / h
    return y0 + t * d1 + 0.5 * t * (t - 1) * d2


def qua(x1, x2, x3, y1, y2, y3, x0):
    """Lagrange quadratic interpolation."""
    a23 = x2 - x3
    a12 = x1 - x2
    a13 = x1 - x3
    a01 = x0 - x1
    a02 = x0 - x2
    a03 = x0 - x3
    return a02 * a03 / a12 / a13 * y1 - a01 * a03 / a12 / a23 * y2 + a01 * a02 / a13 / a23 * y3


def quad_gg(x1, x2, x3, y1, y2, y3, x0):
    """Derivative of Lagrange quadratic interpolation."""
    a23 = x2 - x3
    a12 = x1 - x2
    a13 = x1 - x3
    a01 = x0 - x1
    a02 = x0 - x2
    a03 = x0 - x3
    return (a02 + a03) / a12 / a13 * y1 - (a01 + a03) / a12 / a23 * y2 + (a01 + a02) / a13 / a23 * y3


def flin(x1, x2, y1, y2, x0):
    """Linear interpolation / extrapolation."""
    a12 = x1 - x2
    a01 = x0 - x1
    a02 = x0 - x2
    return (a02 * y1 - a01 * y2) / a12


def spline(x, y, n):
    """Compute cubic spline second derivatives. Returns y2 array."""
    nmax = 100
    u = [0.0] * nmax
    y2 = [0.0] * n

    # Endpoint derivatives via cubic fit
    yp1 = ((y[2] - y[0]) * (x[1] - x[0])**2 - (y[1] - y[0]) * (x[2] - x[0])**2) / \
          ((x[2] - x[0]) * (x[1] - x[0]) * (x[1] - x[2]))
    ypn = ((y[n - 3] - y[n - 1]) * (x[n - 2] - x[n - 1])**2 - (y[n - 2] - y[n - 1]) * (x[n - 3] - x[n - 1])**2) / \
          ((x[n - 3] - x[n - 1]) * (x[n - 2] - x[n - 1]) * (x[n - 2] - x[n - 3]))

    y2[0] = -0.5
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1)

    for i in range(1, n - 1):
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1])
        p = sig * y2[i - 1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (6.0 * ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1])) / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p

    qn = 0.5
    un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]))
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0)

    for k in range(n - 2, -1, -1):
        y2[k] = y2[k] * y2[k + 1] + u[k]

    return y2


def splint(xa, ya, n, y2a, x):
    """Evaluate cubic spline. Returns (y, yp)."""
    klo = 0
    khi = n - 1
    while khi - klo > 1:
        k = (khi + klo) // 2
        if xa[k] > x:
            khi = k
        else:
            klo = k
    h = xa[khi] - xa[klo]
    a = (xa[khi] - x) / h
    b = (x - xa[klo]) / h
    yval = a * ya[klo] + b * ya[khi] + ((a**3 - a) * y2a[klo] + (b**3 - b) * y2a[khi]) * h**2 / 6.0
    yp = 0.05 * ((-ya[klo] + ya[khi]) / h + (-(3.0 * a**2 - 1.0) * y2a[klo] + (3.0 * b**2 - 1.0) * y2a[khi]) * h / 6.0)
    return yval, yp


def getd(f, n):
    """Spline coefficients for unit intervals. Returns (d, fp1, fpn)."""
    t = [0.0] * max(n, 85)
    d = [0.0] * n

    fp1 = (-11.0 * f[0] + 18.0 * f[1] - 9.0 * f[2] + 2.0 * f[3]) / 6.0
    fpn = (11.0 * f[n - 1] - 18.0 * f[n - 2] + 9.0 * f[n - 3] - 2.0 * f[n - 4]) / 6.0

    d[0] = -0.5
    t[0] = 0.5 * (-f[0] + f[1] - fp1)

    for j in range(1, n - 1):
        d[j] = -1.0 / (4.0 + d[j - 1])
        t[j] = -d[j] * (f[j - 1] - 2.0 * f[j] + f[j + 1] - t[j - 1])

    d[n - 1] = (fpn + f[n - 2] - f[n - 1] - t[n - 2]) / (2.0 + d[n - 2])

    for j in range(n - 2, -1, -1):
        d[j] = d[j] * d[j + 1] + t[j]

    return d, fp1, fpn


# ============================================================
# Batched GPU-accelerated versions
# ============================================================

def batch_indic(x0, x, device=None):
    """Batched binary search using torch.searchsorted.

    Args:
        x0: (N,) tensor of query points
        x: (M,) tensor of table values (ascending or descending)
        device: torch device

    Returns:
        (N,) tensor of 1-indexed positions (matching Fortran convention)
    """
    if device is not None:
        x0 = x0.to(device)
        x = x.to(device)
    m = x.shape[0]
    ascending = x[-1] > x[0]
    if ascending:
        k = torch.searchsorted(x, x0, right=False)
        # Clamp: searchsorted returns [0, m], we need [1, m-1] (1-indexed)
        k = k.clamp(1, m - 1)
    else:
        x_flip = x.flip(0)
        k_flip = torch.searchsorted(x_flip, x0, right=False)
        k = m - k_flip
        k = k.clamp(1, m - 1)
    return k


def batch_fipoi(x0, a, b, device=None):
    """Batched linear interpolation.

    Args:
        x0: (N,) tensor of query points
        a: (M,) tensor of x-values (table)
        b: (M,) tensor of y-values (table)

    Returns:
        (N,) tensor of interpolated values
    """
    if device is not None:
        x0 = x0.to(device)
        a = a.to(device)
        b = b.to(device)
    k = batch_indic(x0, a, device) - 1  # Convert to 0-indexed
    bk = b[k]
    bk1 = b[k + 1]
    ak = a[k]
    ak1 = a[k + 1]
    return bk + (bk1 - bk) * (x0 - ak) / (ak1 - ak)


def batch_fipoi2(x0, a, b, device=None):
    """Batched linear interpolation with derivative.

    Returns:
        (values, derivatives) both (N,) tensors
    """
    if device is not None:
        x0 = x0.to(device)
        a = a.to(device)
        b = b.to(device)
    k = batch_indic(x0, a, device) - 1
    bk = b[k]
    bk1 = b[k + 1]
    ak = a[k]
    ak1 = a[k + 1]
    slope = (bk1 - bk) / (ak1 - ak)
    vals = bk + slope * (x0 - ak)
    return vals, slope


def batch_quint(x, x0, h, y0, y1, y2, device=None):
    """Batched quadratic interpolation for equidistant points."""
    d1 = y1 - y0
    d2 = y2 - 2.0 * y1 + y0
    t = (x - x0) / h
    return y0 + t * d1 + 0.5 * t * (t - 1) * d2


def batch_qua(x1, x2, x3, y1, y2, y3, x0, device=None):
    """Batched Lagrange quadratic interpolation."""
    a23 = x2 - x3
    a12 = x1 - x2
    a13 = x1 - x3
    a01 = x0 - x1
    a02 = x0 - x2
    a03 = x0 - x3
    return a02 * a03 / a12 / a13 * y1 - a01 * a03 / a12 / a23 * y2 + a01 * a02 / a13 / a23 * y3


def batch_quad_gg(x1, x2, x3, y1, y2, y3, x0, device=None):
    """Batched derivative of Lagrange quadratic interpolation."""
    a23 = x2 - x3
    a12 = x1 - x2
    a13 = x1 - x3
    a01 = x0 - x1
    a02 = x0 - x2
    a03 = x0 - x3
    return (a02 + a03) / a12 / a13 * y1 - (a01 + a03) / a12 / a23 * y2 + (a01 + a02) / a13 / a23 * y3


def batch_flin(x1, x2, y1, y2, x0, device=None):
    """Batched linear interpolation / extrapolation."""
    a12 = x1 - x2
    a01 = x0 - x1
    a02 = x0 - x2
    return (a02 * y1 - a01 * y2) / a12


def batch_spline(x, y, device=None):
    """Batched cubic spline coefficient computation.

    Args:
        x: (B, N) tensor of x-values for B independent splines
        y: (B, N) tensor of y-values

    Returns:
        (B, N) tensor of second derivatives
    """
    if device is not None:
        x = x.to(device)
        y = y.to(device)
    B, N = x.shape

    # Endpoint derivatives via cubic fit
    yp1 = ((y[:, 2] - y[:, 0]) * (x[:, 1] - x[:, 0])**2 - (y[:, 1] - y[:, 0]) * (x[:, 2] - x[:, 0])**2) / \
          ((x[:, 2] - x[:, 0]) * (x[:, 1] - x[:, 0]) * (x[:, 1] - x[:, 2]))
    ypn = ((y[:, -3] - y[:, -1]) * (x[:, -2] - x[:, -1])**2 - (y[:, -2] - y[:, -1]) * (x[:, -3] - x[:, -1])**2) / \
          ((x[:, -3] - x[:, -1]) * (x[:, -2] - x[:, -1]) * (x[:, -2] - x[:, -3]))

    y2 = torch.zeros(B, N, dtype=x.dtype, device=x.device)
    u = torch.zeros(B, N, dtype=x.dtype, device=x.device)

    y2[:, 0] = -0.5
    u[:, 0] = (3.0 / (x[:, 1] - x[:, 0])) * ((y[:, 1] - y[:, 0]) / (x[:, 1] - x[:, 0]) - yp1)

    for i in range(1, N - 1):
        sig = (x[:, i] - x[:, i - 1]) / (x[:, i + 1] - x[:, i - 1])
        p = sig * y2[:, i - 1] + 2.0
        y2[:, i] = (sig - 1.0) / p
        u[:, i] = (6.0 * ((y[:, i + 1] - y[:, i]) / (x[:, i + 1] - x[:, i]) - (y[:, i] - y[:, i - 1]) / (x[:, i] - x[:, i - 1])) / (x[:, i + 1] - x[:, i - 1]) - sig * u[:, i - 1]) / p

    qn = 0.5
    un = (3.0 / (x[:, -1] - x[:, -2])) * (ypn - (y[:, -1] - y[:, -2]) / (x[:, -1] - x[:, -2]))
    y2[:, -1] = (un - qn * u[:, -2]) / (qn * y2[:, -2] + 1.0)

    for k in range(N - 2, -1, -1):
        y2[:, k] = y2[:, k] * y2[:, k + 1] + u[:, k]

    return y2


def batch_splint(xa, ya, y2a, x_query, device=None):
    """Batched cubic spline evaluation.

    Args:
        xa: (B, N) or (N,) tensor of x-values
        ya: (B, N) or (N,) tensor of y-values
        y2a: (B, N) or (N,) tensor of second derivatives
        x_query: (B,) tensor of query points

    Returns:
        (y, yp) tensors of shape (B,)
    """
    if device is not None:
        xa = xa.to(device)
        ya = ya.to(device)
        y2a = y2a.to(device)
        x_query = x_query.to(device)

    if xa.dim() == 1:
        # Single table, multiple queries: use searchsorted
        khi = torch.searchsorted(xa, x_query, right=False).clamp(1, xa.shape[0] - 1)
        klo = khi - 1
        h = xa[khi] - xa[klo]
        a = (xa[khi] - x_query) / h
        b = (x_query - xa[klo]) / h
        yval = a * ya[klo] + b * ya[khi] + ((a**3 - a) * y2a[klo] + (b**3 - b) * y2a[khi]) * h**2 / 6.0
        yp = 0.05 * ((-ya[klo] + ya[khi]) / h + (-(3.0 * a**2 - 1.0) * y2a[klo] + (3.0 * b**2 - 1.0) * y2a[khi]) * h / 6.0)
        return yval, yp
    else:
        # Batched tables: need per-batch search
        B, N = xa.shape
        khi = torch.searchsorted(xa, x_query.unsqueeze(1), right=False).squeeze(1).clamp(1, N - 1)
        klo = khi - 1
        batch_idx = torch.arange(B, device=xa.device)
        h = xa[batch_idx, khi] - xa[batch_idx, klo]
        a = (xa[batch_idx, khi] - x_query) / h
        b_val = (x_query - xa[batch_idx, klo]) / h
        yval = a * ya[batch_idx, klo] + b_val * ya[batch_idx, khi] + \
               ((a**3 - a) * y2a[batch_idx, klo] + (b_val**3 - b_val) * y2a[batch_idx, khi]) * h**2 / 6.0
        yp = 0.05 * ((-ya[batch_idx, klo] + ya[batch_idx, khi]) / h +
                      (-(3.0 * a**2 - 1.0) * y2a[batch_idx, klo] + (3.0 * b_val**2 - 1.0) * y2a[batch_idx, khi]) * h / 6.0)
        return yval, yp


def batch_getd(f, device=None):
    """Batched spline coefficients for unit intervals.

    Args:
        f: (B, N) tensor

    Returns:
        (d, fp1, fpn) where d is (B, N), fp1 and fpn are (B,)
    """
    if device is not None:
        f = f.to(device)
    B, N = f.shape

    fp1 = (-11.0 * f[:, 0] + 18.0 * f[:, 1] - 9.0 * f[:, 2] + 2.0 * f[:, 3]) / 6.0
    fpn = (11.0 * f[:, -1] - 18.0 * f[:, -2] + 9.0 * f[:, -3] - 2.0 * f[:, -4]) / 6.0

    d = torch.zeros(B, N, dtype=f.dtype, device=f.device)
    t = torch.zeros(B, N, dtype=f.dtype, device=f.device)

    d[:, 0] = -0.5
    t[:, 0] = 0.5 * (-f[:, 0] + f[:, 1] - fp1)

    for j in range(1, N - 1):
        d[:, j] = -1.0 / (4.0 + d[:, j - 1])
        t[:, j] = -d[:, j] * (f[:, j - 1] - 2.0 * f[:, j] + f[:, j + 1] - t[:, j - 1])

    d[:, -1] = (fpn + f[:, -2] - f[:, -1] - t[:, -2]) / (2.0 + d[:, -2])

    for j in range(N - 2, -1, -1):
        d[:, j] = d[:, j] * d[:, j + 1] + t[:, j]

    return d, fp1, fpn
