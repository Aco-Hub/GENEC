"""Tests for nuclear energy generation rates against reference model."""

import pytest
import numpy as np
import torch
from genec.energy import pp_chain_rate, cno_rate, total_energy_rate
from genec.energy import batch_pp_chain, batch_cno, batch_total_energy


class TestEnergyScalar:
    """Test scalar energy rates against reference data."""

    def test_core_has_energy(self, reference_data):
        """Core zones should produce nuclear energy."""
        rd = reference_data
        T = 10.0**rd['log_T'][0]
        rho = 10.0**rd['log_rho'][0]
        X = rd['X_H1'][0]
        eps = total_energy_rate(T, rho, X, rd['X_He4'][0])
        assert eps > 0, f"Core epsilon should be > 0, got {eps}"

    def test_surface_no_energy(self, reference_data):
        """Surface zones (cold) should produce no nuclear energy."""
        rd = reference_data
        T = 10.0**rd['log_T'][-1]
        rho = 10.0**rd['log_rho'][-1]
        X = rd['X_H1'][-1]
        eps = total_energy_rate(T, rho, X, rd['X_He4'][-1])
        assert eps == 0.0 or eps < 1e-20, f"Surface epsilon should be ~0, got {eps}"

    def test_pp_dominates_at_low_T(self, reference_data):
        """PP chain should dominate at T < 15 MK."""
        # Use a zone where T ~ 10 MK
        rd = reference_data
        for idx in range(rd['n_zones']):
            T = 10.0**rd['log_T'][idx]
            if 8e6 < T < 15e6:
                rho = 10.0**rd['log_rho'][idx]
                X = rd['X_H1'][idx]
                eps_pp = pp_chain_rate(T, rho, X)
                eps_cno = cno_rate(T, rho, X, 2.56e-3, 7.4e-4)
                if eps_pp > 0:
                    assert eps_pp >= eps_cno, \
                        f"PP should dominate at T={T/1e6:.1f}MK: pp={eps_pp:.2e}, cno={eps_cno:.2e}"
                break

    def test_energy_order_of_magnitude(self, reference_data):
        """Computed energy should be within ~2 orders of magnitude of reference."""
        rd = reference_data
        # Find zones with significant energy production
        for idx in range(rd['n_zones']):
            ref_eps = rd['epsilon'][idx]
            if ref_eps > 1.0:  # Significant energy production
                T = 10.0**rd['log_T'][idx]
                rho = 10.0**rd['log_rho'][idx]
                X = rd['X_H1'][idx]
                eps = total_energy_rate(T, rho, X, rd['X_He4'][idx])
                if eps > 0 and ref_eps > 0:
                    ratio = eps / ref_eps
                    # Within factor of 100 (we implement simplified PP+CNO only)
                    assert 0.01 < ratio < 100, \
                        f"Zone {idx}: eps={eps:.2e} vs ref={ref_eps:.2e} (ratio={ratio:.2f})"
                break

    def test_energy_increases_with_temperature(self):
        """Energy rate should increase with temperature."""
        rho = 100.0  # typical core density
        X = 0.72
        Y = 0.266
        temps = [5e6, 10e6, 15e6, 20e6]
        rates = [total_energy_rate(T, rho, X, Y) for T in temps]
        for i in range(len(rates) - 1):
            if rates[i] > 0:
                assert rates[i + 1] >= rates[i], \
                    f"Rate should increase: eps({temps[i]/1e6:.0f}MK)={rates[i]:.2e} > eps({temps[i+1]/1e6:.0f}MK)={rates[i+1]:.2e}"


class TestEnergyBatch:
    """Test batched GPU energy rates."""

    def test_batch_matches_scalar(self, reference_data, device, dtype):
        """GPU batch should match scalar for core zones."""
        rd = reference_data
        indices = np.arange(0, min(50, rd['n_zones']), 5)
        T = torch.tensor(10.0**rd['log_T'][indices], dtype=dtype)
        rho = torch.tensor(10.0**rd['log_rho'][indices], dtype=dtype)
        X = rd['X_H1'][0]  # Same composition

        gpu_eps = batch_total_energy(T, rho, X, rd['X_He4'][0],
                                     device=device, dtype=dtype)

        for i, idx in enumerate(indices):
            cpu_eps = total_energy_rate(
                10.0**rd['log_T'][idx], 10.0**rd['log_rho'][idx],
                rd['X_H1'][idx], rd['X_He4'][idx])
            gpu_val = gpu_eps[i].item()
            if cpu_eps > 0:
                rel_err = abs(gpu_val - cpu_eps) / cpu_eps
                assert rel_err < 0.01, \
                    f"Zone {idx}: GPU={gpu_val:.4e} vs CPU={cpu_eps:.4e}"
