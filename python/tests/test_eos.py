"""Tests for ideal gas + radiation EOS against reference model."""

import pytest
import numpy as np
import torch
from genec.eos import eos_ideal, batch_eos_ideal


class TestEOSScalar:
    """Test scalar EOS against reference data."""

    def test_core_beta_near_one(self, reference_data):
        """In the core, gas pressure dominates (beta ~ 1)."""
        rd = reference_data
        result = eos_ideal(rd['log_P'][0], rd['log_T'][0], rd['mu'][0])
        # For 1 Msol, beta should be very close to 1 (radiation negligible)
        assert result['beta'] > 0.99, f"Core beta={result['beta']:.6f}, expected ~1"

    def test_core_dlnP_dlnrho(self, reference_data):
        """(dlnP/dlnrho)_T should match reference for core zones."""
        rd = reference_data
        # Test first 50 core zones where EOS is simplest
        for idx in range(0, 50, 5):
            result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
            ref_val = rd['dlnP_dlnrho_T'][idx]
            rel_err = abs(result['dlnP_dlnrho_T'] - ref_val) / max(abs(ref_val), 1e-10)
            assert rel_err < 0.02, \
                f"Zone {idx}: dlnP_dlnrho={result['dlnP_dlnrho_T']:.6f} vs ref={ref_val:.6f}"

    def test_core_dlnP_dlnT(self, reference_data):
        """(dlnP/dlnT)_rho should be close to reference for core zones.
        Note: reference includes ionization corrections we don't model, so tolerance is wider."""
        rd = reference_data
        for idx in range(0, 50, 5):
            result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
            ref_val = rd['dlnP_dlnT_rho'][idx]
            rel_err = abs(result['dlnP_dlnT_rho'] - ref_val) / max(abs(ref_val), 1e-10)
            assert rel_err < 0.05, \
                f"Zone {idx}: dlnP_dlnT={result['dlnP_dlnT_rho']:.6f} vs ref={ref_val:.6f}"

    def test_nabla_ad_range(self, reference_data):
        """nabla_ad should be in a physically reasonable range.
        Note: reference nabla_ad includes ionization zone effects (partial ionization
        reduces nabla_ad significantly). Our pure ideal gas+radiation EOS gives
        different values in ionization zones, so we test only the core."""
        rd = reference_data
        # Test deep core zones where ionization is complete (simple ideal gas)
        for idx in [0, 5, 10]:
            result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
            # nabla_ad should be physically reasonable (0.1 to 0.5)
            assert 0.1 < result['nabla_ad'] < 0.5, \
                f"Zone {idx}: nabla_ad={result['nabla_ad']:.4f} outside range"
        # For a zone far from ionization boundaries, nabla_ad should be closer to ref
        idx = 300  # mid-zone, likely fully ionized interior
        result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
        ref_val = rd['nabla_ad'][idx]
        assert abs(result['nabla_ad'] - ref_val) < 0.2, \
            f"Zone {idx}: nabla_ad={result['nabla_ad']:.4f} vs ref={ref_val:.4f}"

    def test_density_reasonable(self, reference_data):
        """Computed density should match reference log(rho)."""
        rd = reference_data
        for idx in range(0, 50, 5):
            result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
            ref_val = rd['log_rho'][idx]
            assert abs(result['log_rho'] - ref_val) < 0.1, \
                f"Zone {idx}: log_rho={result['log_rho']:.4f} vs ref={ref_val:.4f}"


class TestEOSBatch:
    """Test batched GPU EOS."""

    def test_batch_matches_scalar(self, reference_data, device, dtype):
        """GPU batch results should match scalar loop."""
        rd = reference_data
        indices = np.arange(0, 50, 5)
        log_P = torch.tensor(rd['log_P'][indices], dtype=dtype)
        log_T = torch.tensor(rd['log_T'][indices], dtype=dtype)
        mu = torch.tensor(rd['mu'][indices], dtype=dtype)

        batch_result = batch_eos_ideal(log_P, log_T, mu, device=device, dtype=dtype)

        for i, idx in enumerate(indices):
            scalar_result = eos_ideal(rd['log_P'][idx], rd['log_T'][idx], rd['mu'][idx])
            for key in ['log_rho', 'beta', 'dlnP_dlnrho_T', 'dlnP_dlnT_rho']:
                gpu_val = batch_result[key][i].item()
                cpu_val = scalar_result[key]
                assert abs(gpu_val - cpu_val) < 1e-10, \
                    f"Zone {idx} {key}: GPU={gpu_val:.10f}, scalar={cpu_val:.10f}"
