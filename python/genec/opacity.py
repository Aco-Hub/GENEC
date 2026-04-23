"""GPU-accelerated OPAL opacity table interpolation.

Ported from the Fortran ``opacgn93`` routine in ``src/opacity.f90``.
Reads GN93-format opacity tables (e.g. ``opaSol_GN93.dat``) and provides
both scalar and batched GPU interpolation of log10(kappa) as a function
of composition (Z, X) and thermodynamic state (logT, logR).

The batch path uses bilinear interpolation in (logT, logR) at each
(X, Z) bracket corner, followed by linear interpolation in X and then Z.
Derivatives dlog(kappa)/dlogT and dlog(kappa)/dlogR are computed
analytically from the bilinear weights.
"""

from __future__ import annotations

import os
import re
from typing import Dict, Optional, Tuple

import torch

# ---------------------------------------------------------------------------
# Table dimensions (matching Fortran parameters)
# ---------------------------------------------------------------------------
_MX: int = 10       # number of X slots (9 fixed + 1 dynamic = 1-Z)
_MZ: int = 13       # number of Z grid points
_NTM: int = 85      # number of logT rows per table
_NRM: int = 19      # number of logR columns per table
_HEADER_LINES: int = 240  # lines to skip at top of .dat file
_LINES_PER_TABLE: int = 92  # 7 formatting + 85 data rows

# Number of Z sub-tables available for each X slot
_N_TABLES_PER_X: Tuple[int, ...] = (13, 13, 13, 13, 13, 13, 13, 13, 10, 12)

# Canonical X grid (slot 10 is set to 1-Z at runtime)
_XA_DEFAULT: Tuple[float, ...] = (
    0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9, 0.95, 0.0,
)

# Canonical Z grid
_ZA: Tuple[float, ...] = (
    0.0, 1.0e-4, 3.0e-4, 1.0e-3, 2.0e-3, 4.0e-3,
    1.0e-2, 2.0e-2, 3.0e-2, 4.0e-2, 6.0e-2, 8.0e-2, 1.0e-1,
)

# Sentinel value in the tables indicating out-of-range data
_SENTINEL: float = 9.999

# Default path to the opacity data file
_DEFAULT_PATH: str = os.path.join(
    os.path.dirname(__file__), '..', '..', 'src', 'inputs', 'opaSol_GN93.dat',
)


# =====================================================================
# Fixed-width line parser
# =====================================================================

def _parse_data_row(line: str) -> Tuple[float, list[float]]:
    """Parse one 85-row data line in Fortran f4.2,19f7.3 format.

    Args:
        line: Raw text line from the opacity file.

    Returns:
        (logT, values) where *values* has exactly ``_NRM`` entries.
        Missing columns (line too short) are filled with ``_SENTINEL``.
    """
    # Pad to full width so slicing never goes out of bounds
    padded = line.ljust(4 + _NRM * 7)
    log_t = float(padded[0:4])
    values: list[float] = []
    for i in range(_NRM):
        start = 4 + 7 * i
        field = padded[start:start + 7].strip()
        if field == '' or field == '.':
            values.append(_SENTINEL)
        else:
            try:
                values.append(float(field))
            except ValueError:
                values.append(_SENTINEL)
    return log_t, values


# =====================================================================
# OpacityTable
# =====================================================================

class OpacityTable:
    """Load and store OPAL GN93 opacity tables for GPU interpolation.

    After construction the table data lives on CPU.  Call :meth:`to` to
    move everything to a CUDA device before using :meth:`batch_interpolate`.

    Attributes:
        data: (mx, mz, nt, nr) tensor of log10(kappa) values.
        log_T_grid: (nt,) tensor of logT values read from the file.
        log_R_grid: (nr,) tensor of logR values (-8.0 .. +1.0).
        X_grid: (mx,) tensor of hydrogen fractions (slot 10 = 0 initially).
        Z_grid: (mz,) tensor of metallicities.
        n_per_X: (mx,) tensor giving number of valid Z sub-tables per X.
    """

    # -----------------------------------------------------------------
    # Construction
    # -----------------------------------------------------------------

    def __init__(self, filepath: Optional[str] = None) -> None:
        """Load an OPAL opacity table from *filepath*.

        Args:
            filepath: Path to the ``.dat`` file.  When *None* the default
                ``opaSol_GN93.dat`` shipped with GENEC is used.
        """
        if filepath is None:
            filepath = _DEFAULT_PATH
        filepath = os.path.abspath(filepath)
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"Opacity table not found: {filepath}")

        # Allocate storage -- fill with sentinel so missing sub-tables
        # are automatically marked out-of-range.
        data_np = [[[[_SENTINEL] * _NRM for _ in range(_NTM)]
                     for _ in range(_MZ)]
                    for _ in range(_MX)]

        log_T_grid: list[float] = []
        log_R_grid: Optional[list[float]] = None
        _first_table = True  # flag: build logT grid from first table only

        with open(filepath, 'r') as fh:
            lines = fh.readlines()

        # Skip the header block
        pos = _HEADER_LINES

        # Read tables in Fortran order: for each X slot, for each Z sub-table
        for m in range(_MX):
            for i in range(_N_TABLES_PER_X[m]):
                # --- 7 formatting / header lines ---
                # line 1: blank
                pos += 1
                # line 2: TABLE header  ->  extract X, Y, Z
                header = lines[pos]
                pos += 1
                match = re.search(
                    r'X=([0-9.]+)\s+Y=([0-9.]+)\s+Z=([0-9.]+)', header,
                )
                if match is None:
                    raise ValueError(
                        f"Cannot parse TABLE header at line {pos}: {header!r}"
                    )
                _x_val = float(match.group(1))
                _z_val = float(match.group(3))

                # Determine Z sub-index: the i-th sub-table for this X
                # maps directly to mz index i (0-based)
                mz_idx = i

                # line 3: blank
                pos += 1
                # line 4: "log R" text
                pos += 1
                # line 5: blank
                pos += 1
                # line 6: logR header values
                if log_R_grid is None:
                    rline = lines[pos]
                    # skip "logT" prefix (4 chars), parse f6.1 + 18*f7.1
                    rpadded = rline.ljust(4 + 6 + 18 * 7)
                    log_R_grid = []
                    # first value: chars 4:10 (f6.1)
                    log_R_grid.append(float(rpadded[4:10].strip()))
                    # remaining 18 values: chars 10+7*j : 10+7*(j+1)
                    for j in range(18):
                        start = 10 + 7 * j
                        log_R_grid.append(float(rpadded[start:start + 7].strip()))
                pos += 1
                # line 7: blank
                pos += 1

                # --- 85 data rows ---
                for k in range(_NTM):
                    log_t_k, vals = _parse_data_row(lines[pos])
                    pos += 1

                    if _first_table:
                        log_T_grid.append(log_t_k)

                    for ll in range(_NRM):
                        data_np[m][mz_idx][k][ll] = vals[ll]

                _first_table = False

        assert len(log_T_grid) == _NTM
        assert log_R_grid is not None and len(log_R_grid) == _NRM

        # Convert to tensors (CPU, float64)
        self.data = torch.tensor(data_np, dtype=torch.float64)       # (mx,mz,nt,nr)
        self.log_T_grid = torch.tensor(log_T_grid, dtype=torch.float64)  # (nt,)
        self.log_R_grid = torch.tensor(log_R_grid, dtype=torch.float64)  # (nr,)
        self.X_grid = torch.tensor(list(_XA_DEFAULT), dtype=torch.float64)  # (mx,)
        self.Z_grid = torch.tensor(list(_ZA), dtype=torch.float64)        # (mz,)
        self.n_per_X = torch.tensor(list(_N_TABLES_PER_X), dtype=torch.int64)

        # Pre-process: replace sentinel (9.999) values by propagating the
        # nearest valid value downward along the logT axis.  This avoids
        # special-case branching during interpolation.
        self._fill_sentinels()

    # -----------------------------------------------------------------

    def _fill_sentinels(self) -> None:
        """Replace 9.999 sentinel entries by propagating valid neighbours.

        For each (X, Z, logR) column the method performs two passes:

        1. **Low-T pass** -- scan from low logT upward; fill any leading
           sentinels with the first valid value encountered.
        2. **High-T pass** -- scan from low logT upward; once a valid
           value has been seen, any subsequent sentinels (from short
           data lines at the jagged table edge) are replaced with
           the most recent valid value.

        Columns that are entirely sentinel are left unchanged (they
        correspond to X/Z combinations outside the table coverage).
        """
        sentinel_thresh = _SENTINEL - 0.01  # 9.989
        for m in range(_MX):
            for j in range(int(self.n_per_X[m].item())):
                for lr in range(_NRM):
                    col = self.data[m, j, :, lr]

                    # --- Pass 1: fill sentinels at the low-T end ---
                    first_valid: Optional[int] = None
                    for k in range(_NTM):
                        if col[k].item() < sentinel_thresh:
                            first_valid = k
                            break
                    if first_valid is not None and first_valid > 0:
                        fill_val = col[first_valid].item()
                        for k in range(first_valid):
                            self.data[m, j, k, lr] = fill_val

                    # --- Pass 2: fill sentinels at the high-T end ---
                    # Forward scan: once we have seen a valid value,
                    # propagate it into any subsequent sentinel slots.
                    last_good: Optional[float] = None
                    for k in range(_NTM):
                        v = self.data[m, j, k, lr].item()
                        if v < sentinel_thresh:
                            last_good = v
                        elif last_good is not None:
                            self.data[m, j, k, lr] = last_good

    # -----------------------------------------------------------------
    # Device management
    # -----------------------------------------------------------------

    def to(self, device: torch.device | str, dtype: torch.dtype = torch.float64) -> "OpacityTable":
        """Move all internal tensors to *device* and *dtype*.

        Args:
            device: Target device (e.g. ``"cuda"`` or ``torch.device("cpu")``).
            dtype: Floating-point dtype.  Default ``torch.float64``.

        Returns:
            ``self`` for chaining.
        """
        self.data = self.data.to(device=device, dtype=dtype)
        self.log_T_grid = self.log_T_grid.to(device=device, dtype=dtype)
        self.log_R_grid = self.log_R_grid.to(device=device, dtype=dtype)
        self.X_grid = self.X_grid.to(device=device, dtype=dtype)
        self.Z_grid = self.Z_grid.to(device=device, dtype=dtype)
        self.n_per_X = self.n_per_X.to(device=device)
        return self

    # -----------------------------------------------------------------
    # Scalar interpolation
    # -----------------------------------------------------------------

    def interpolate(
        self,
        z: float,
        xh: float,
        log_t: float,
        log_r: float,
    ) -> Dict[str, float]:
        """Scalar interpolation for a single shell.

        Uses bilinear interpolation in (logT, logR) at each (X, Z)
        bracket corner, then linear interpolation in X and Z.

        Args:
            z: Metallicity Z.
            xh: Hydrogen mass fraction X.
            log_t: log10(T) in Kelvin.
            log_r: log10(R) where R = rho / T6**3.

        Returns:
            Dictionary with keys ``log_kappa``, ``dlogk_dlogT``,
            ``dlogk_dlogR``.
        """
        log_t_t = torch.tensor([log_t], dtype=torch.float64, device=self.data.device)
        log_r_t = torch.tensor([log_r], dtype=torch.float64, device=self.data.device)
        result = self.batch_interpolate(z, xh, log_t_t, log_r_t)
        return {
            'log_kappa': result['log_kappa'].item(),
            'dlogk_dlogT': result['dlogk_dlogT'].item(),
            'dlogk_dlogR': result['dlogk_dlogR'].item(),
        }

    # -----------------------------------------------------------------
    # Batched GPU interpolation
    # -----------------------------------------------------------------

    def batch_interpolate(
        self,
        z: float,
        xh: float,
        log_t_batch: torch.Tensor,
        log_r_batch: torch.Tensor,
    ) -> Dict[str, torch.Tensor]:
        """GPU-batched opacity interpolation for N shells.

        All shells share the same composition (z, xh).  Temperature and
        density vary per shell.

        Args:
            z: Scalar metallicity.
            xh: Scalar hydrogen mass fraction.
            log_t_batch: (N,) tensor of log10(T).
            log_r_batch: (N,) tensor of log10(R), R = rho/T6**3.

        Returns:
            Dictionary with (N,) tensors ``log_kappa``, ``dlogk_dlogT``,
            ``dlogk_dlogR``.
        """
        dev = self.data.device
        dtype = self.data.dtype
        log_t_batch = log_t_batch.to(device=dev, dtype=dtype)
        log_r_batch = log_r_batch.to(device=dev, dtype=dtype)

        # --- 1. Build the X grid with dynamic 10th slot ---------------
        X_grid = self.X_grid.clone()
        X_grid[_MX - 1] = 1.0 - z
        # If X(10) < X(9), collapse the last slot onto X(9)
        if X_grid[_MX - 1] < X_grid[_MX - 2]:
            X_grid[_MX - 2] = X_grid[_MX - 1]

        # Transform to log10(0.005 + X) for smoother interpolation
        # (matches Fortran variable xxx / xx)
        xx_grid = torch.log10(0.005 + X_grid)
        xx_h = torch.log10(torch.tensor(0.005 + xh, device=dev, dtype=dtype))

        # --- 2. Bracket in Z ------------------------------------------
        j0, j1, wz = self._bracket_z(z)

        # --- 3. Bracket in X (using transformed variable) -------------
        i0, i1, wx = self._bracket_x(xx_h, xx_grid)

        # --- 4. For each (X_i, Z_j) corner, bilinear in (logT, logR) -
        corners_i = [i0, i1]
        corners_j = [j0, j1]
        corner_vals = {}   # (ci, cj) -> (N,) tensor of log_kappa
        corner_dT = {}     # (ci, cj) -> (N,) tensor of dlogk/dlogT
        corner_dR = {}     # (ci, cj) -> (N,) tensor of dlogk/dlogR

        for ci in corners_i:
            for cj in corners_j:
                # Check that the table actually has data at this (X,Z)
                if cj >= int(self.n_per_X[ci].item()):
                    # No data -- use sentinel-like values
                    n = log_t_batch.shape[0]
                    corner_vals[(ci, cj)] = torch.full((n,), _SENTINEL, device=dev, dtype=dtype)
                    corner_dT[(ci, cj)] = torch.zeros(n, device=dev, dtype=dtype)
                    corner_dR[(ci, cj)] = torch.zeros(n, device=dev, dtype=dtype)
                    continue
                table = self.data[ci, cj]  # (nt, nr)
                v, dt, dr = self._bilinear_TR(table, log_t_batch, log_r_batch)
                corner_vals[(ci, cj)] = v
                corner_dT[(ci, cj)] = dt
                corner_dR[(ci, cj)] = dr

        # --- 5. Linear interpolation in X -----------------------------
        def _lerp_x(field: dict) -> dict:
            """Linearly interpolate across X brackets."""
            out = {}
            for cj in corners_j:
                v0 = field[(i0, cj)]
                v1 = field[(i1, cj)]
                out[cj] = v0 * (1.0 - wx) + v1 * wx
            return out

        val_by_z = _lerp_x(corner_vals)
        dT_by_z = _lerp_x(corner_dT)
        dR_by_z = _lerp_x(corner_dR)

        # --- 6. Linear interpolation in Z -----------------------------
        log_kappa = val_by_z[j0] * (1.0 - wz) + val_by_z[j1] * wz
        dlogk_dlogT = dT_by_z[j0] * (1.0 - wz) + dT_by_z[j1] * wz
        dlogk_dlogR = dR_by_z[j0] * (1.0 - wz) + dR_by_z[j1] * wz

        return {
            'log_kappa': log_kappa,
            'dlogk_dlogT': dlogk_dlogT,
            'dlogk_dlogR': dlogk_dlogR,
        }

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    def _bracket_z(self, z: float) -> Tuple[int, int, float]:
        """Find bracketing Z indices and interpolation weight.

        Returns:
            (j0, j1, wz) such that Z_grid[j0] <= z <= Z_grid[j1] and
            wz = (z - Z[j0]) / (Z[j1] - Z[j0]).
        """
        Z = self.Z_grid
        mz = Z.shape[0]

        # Exact match check
        for j in range(mz):
            if abs(z - Z[j].item()) < 1.0e-10:
                return j, j, 0.0

        # Binary search
        j1 = int(torch.searchsorted(Z, torch.tensor(z, device=Z.device, dtype=Z.dtype)).item())
        j1 = max(1, min(j1, mz - 1))
        j0 = j1 - 1

        dz = Z[j1].item() - Z[j0].item()
        if dz < 1.0e-30:
            wz = 0.0
        else:
            wz = (z - Z[j0].item()) / dz
        wz = max(0.0, min(1.0, wz))
        return j0, j1, wz

    def _bracket_x(
        self,
        xx_h: torch.Tensor,
        xx_grid: torch.Tensor,
    ) -> Tuple[int, int, float]:
        """Find bracketing X indices and weight in log(0.005+X) space.

        Args:
            xx_h: Scalar tensor log10(0.005 + xh).
            xx_grid: (mx,) tensor of log10(0.005 + X_grid).

        Returns:
            (i0, i1, wx) with linear weight in transformed X space.
        """
        mx = xx_grid.shape[0]
        val = xx_h.item()

        # Determine the effective upper limit of the grid
        # (may be mx-1 or mx depending on the 10th slot adjustment)
        mx_end = mx
        if xx_grid[mx - 1].item() < xx_grid[mx - 2].item():
            mx_end = mx - 1

        # Binary search in the first mx_end entries
        grid_slice = xx_grid[:mx_end]
        j = int(torch.searchsorted(grid_slice, xx_h.unsqueeze(0)).item())
        j = max(1, min(j, mx_end - 1))
        i0 = j - 1
        i1 = j

        dx = xx_grid[i1].item() - xx_grid[i0].item()
        if abs(dx) < 1.0e-30:
            wx = 0.0
        else:
            wx = (val - xx_grid[i0].item()) / dx
        wx = max(0.0, min(1.0, wx))
        return i0, i1, wx

    def _bilinear_TR(
        self,
        table: torch.Tensor,
        log_t: torch.Tensor,
        log_r: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """Bilinear interpolation in (logT, logR) on a single 2-D table.

        Args:
            table: (nt, nr) tensor of log10(kappa).
            log_t: (N,) tensor of log10(T) query points.
            log_r: (N,) tensor of log10(R) query points.

        Returns:
            (values, dval_dlogT, dval_dlogR), each (N,) tensors.
        """
        T_grid = self.log_T_grid  # (nt,)
        R_grid = self.log_R_grid  # (nr,)
        nt = T_grid.shape[0]
        nr = R_grid.shape[0]

        # Clamp inputs to the valid table range
        t_clamped = log_t.clamp(T_grid[0], T_grid[-1])
        r_clamped = log_r.clamp(R_grid[0], R_grid[-1])

        # Find bracketing indices in T and R
        kt = torch.searchsorted(T_grid, t_clamped, right=False).clamp(1, nt - 1)
        kr = torch.searchsorted(R_grid, r_clamped, right=False).clamp(1, nr - 1)

        kt0 = kt - 1  # (N,)
        kt1 = kt       # (N,)
        kr0 = kr - 1
        kr1 = kr

        # Grid spacings at the bracket points
        dt = T_grid[kt1] - T_grid[kt0]  # (N,)
        dr = R_grid[kr1] - R_grid[kr0]  # (N,)

        # Fractional weights
        wt = (t_clamped - T_grid[kt0]) / dt.clamp(min=1.0e-30)  # (N,)
        wr = (r_clamped - R_grid[kr0]) / dr.clamp(min=1.0e-30)  # (N,)
        wt = wt.clamp(0.0, 1.0)
        wr = wr.clamp(0.0, 1.0)

        # Gather the four corner values from the table
        f00 = table[kt0, kr0]  # (N,)
        f10 = table[kt1, kr0]
        f01 = table[kt0, kr1]
        f11 = table[kt1, kr1]

        # Bilinear interpolation
        val = (
            f00 * (1.0 - wt) * (1.0 - wr)
            + f10 * wt * (1.0 - wr)
            + f01 * (1.0 - wt) * wr
            + f11 * wt * wr
        )

        # Analytical derivatives of the bilinear form
        dval_dt = ((f10 - f00) * (1.0 - wr) + (f11 - f01) * wr) / dt.clamp(min=1.0e-30)
        dval_dr = ((f01 - f00) * (1.0 - wt) + (f11 - f10) * wt) / dr.clamp(min=1.0e-30)

        return val, dval_dt, dval_dr


# =====================================================================
# Convenience constructor
# =====================================================================

def load_opacity_table(
    filepath: Optional[str] = None,
    device: Optional[str] = None,
    dtype: torch.dtype = torch.float64,
) -> OpacityTable:
    """Load an OPAL opacity table and optionally move it to a device.

    Args:
        filepath: Path to the ``.dat`` file, or *None* for the default.
        device: ``"cpu"``, ``"cuda"``, etc.  *None* keeps CPU.
        dtype: Floating-point dtype.

    Returns:
        An :class:`OpacityTable` ready for interpolation.
    """
    table = OpacityTable(filepath)
    if device is not None:
        table.to(device, dtype)
    return table
