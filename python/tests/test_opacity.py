"""Tests for opacity table loading and GPU interpolation."""

import pytest
import os
import numpy as np
import torch

TABLE_PATH = os.path.join(
    os.path.dirname(__file__), '..', '..', 'src', 'inputs', 'opaSol_GN93.dat'
)


@pytest.fixture(scope='session')
def opacity_table():
    if not os.path.isfile(TABLE_PATH):
        pytest.skip("Opacity table not found (data files not distributed with git)")
    from genec.opacity import OpacityTable
    return OpacityTable(TABLE_PATH)


class TestOpacityTable:
    """Test opacity table loading and interpolation."""

    def test_table_dimensions(self, opacity_table):
        """Table should be (10, 13, 85, 19)."""
        assert opacity_table.data.shape == (10, 13, 85, 19)

    def test_grids_loaded(self, opacity_table):
        """logT and logR grids should be correct."""
        assert opacity_table.log_T_grid.shape == (85,)
        assert opacity_table.log_R_grid.shape == (19,)
        assert abs(opacity_table.log_T_grid[0].item() - 3.00) < 0.01
        assert abs(opacity_table.log_T_grid[-1].item() - 8.70) < 0.01
        assert abs(opacity_table.log_R_grid[0].item() - (-8.0)) < 0.1
        assert abs(opacity_table.log_R_grid[-1].item() - 1.0) < 0.1

    def test_sentinels_filled(self, opacity_table):
        """After _fill_sentinels, 9.999 values should be replaced in valid sub-tables."""
        # Table 1 (X=0, Z=0): first two logR columns at logT=3.0 were 9.999
        # After fill, they should match the first valid value in that column
        val = opacity_table.data[0, 0, 0, 0].item()
        assert val < 9.0, f"Sentinel should be filled, got {val}"

    def test_scalar_interpolation(self, opacity_table):
        """Scalar interpolation should return reasonable log(kappa)."""
        result = opacity_table.interpolate(z=0.014, xh=0.72, log_t=6.0, log_r=-3.0)
        # At solar composition, T=1e6K, R=1e-3, log(kappa) ~ 0.3-1.0
        assert -5.0 < result['log_kappa'] < 5.0, f"log_kappa={result['log_kappa']}"

    def test_batch_interpolation(self, opacity_table, device, dtype):
        """GPU batch should match scalar results."""
        table = opacity_table
        table.to(device, dtype)

        log_t = torch.tensor([5.0, 6.0, 7.0, 4.0], dtype=dtype, device=device)
        log_r = torch.tensor([-3.0, -2.0, -1.0, -4.0], dtype=dtype, device=device)

        batch = table.batch_interpolate(0.014, 0.72, log_t, log_r)
        assert batch['log_kappa'].shape == (4,)

        for i in range(4):
            scalar = table.interpolate(0.014, 0.72, log_t[i].item(), log_r[i].item())
            assert abs(batch['log_kappa'][i].item() - scalar['log_kappa']) < 0.01, \
                f"i={i}: batch={batch['log_kappa'][i].item():.4f} vs scalar={scalar['log_kappa']:.4f}"

    def test_kappa_reference(self, opacity_table, reference_data, device, dtype):
        """Compare batch kappa against 690-zone reference data."""
        rd = reference_data
        table = opacity_table
        table.to(device, dtype)

        log_T = torch.tensor(rd['log_T'], dtype=dtype, device=device)
        log_rho = torch.tensor(rd['log_rho'], dtype=dtype, device=device)
        T6 = torch.pow(10.0, log_T - 6.0)
        log_R = log_rho - 3.0 * torch.log10(T6.clamp(min=1e-30))

        result = table.batch_interpolate(0.014, 0.72, log_T, log_R)
        kappa_gpu = result['log_kappa'].cpu().numpy()
        kappa_ref = rd['log_kappa']

        # Only compare zones with valid reference data
        valid = np.abs(kappa_ref) < 90
        if valid.sum() < 10:
            pytest.skip("Too few valid kappa reference values")

        errors = np.abs(kappa_gpu[valid] - kappa_ref[valid])
        mean_err = np.mean(errors)
        # Table interpolation precision: bilinear vs Fortran quadratic
        assert mean_err < 0.5, f"Mean kappa error too large: {mean_err:.4f}"
        print(f"  kappa: mean abs error={mean_err:.4f}, max={np.max(errors):.4f} ({valid.sum()} zones)")
