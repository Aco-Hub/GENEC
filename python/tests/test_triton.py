"""Tests for Triton-fused Saha kernel -- validates against standard batch_ionpart."""

import pytest
import numpy as np
import torch

try:
    from genec.triton_saha import triton_batch_ionpart, HAS_TRITON
except ImportError:
    HAS_TRITON = False

from genec.ionisation import batch_ionpart

ABOND = [0.72, 0.266, 2.56e-3, 6.42e-3, 1.65e-3, 5.13e-4]


def _rel_err(got, ref):
    return abs(got - ref) / max(abs(ref), 1e-30)


@pytest.mark.skipif(not HAS_TRITON, reason="Triton not available")
class TestTritonSaha:
    """Validate Triton kernel against standard PyTorch batch_ionpart."""

    def test_core_fully_ionized(self, reference_data, device, dtype):
        """Core zones: Triton should produce HII~1, HeIII~1."""
        rd = reference_data
        log_P = torch.tensor(rd['log_P'][:10], dtype=dtype, device=device)
        log_T = torch.tensor(rd['log_T'][:10], dtype=dtype, device=device)

        result = triton_batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        hii = result['HII'].cpu().numpy()
        heiii = result['HeIII'].cpu().numpy()

        assert np.all(hii > 0.95), f"Core HII: {hii}"
        assert np.all(heiii > 0.95), f"Core HeIII: {heiii}"

    def test_triton_matches_pytorch(self, reference_data, device, dtype):
        """Triton results should match standard batch_ionpart within 5%."""
        rd = reference_data
        indices = np.linspace(0, rd['n_zones'] - 1, 50, dtype=int)
        log_P = torch.tensor(rd['log_P'][indices], dtype=dtype, device=device)
        log_T = torch.tensor(rd['log_T'][indices], dtype=dtype, device=device)

        ref = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        tri = triton_batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)

        for key in ['HII', 'HeII', 'HeIII', 'mu']:
            ref_vals = ref[key].cpu().numpy()
            tri_vals = tri[key].cpu().numpy()
            for i in range(len(indices)):
                r, t = ref_vals[i], tri_vals[i]
                # Triton approximates heavy elements, so transition zones differ
                if abs(r) < 0.03:
                    assert abs(t - r) < 0.03, \
                        f"Zone {indices[i]} {key}: triton={t:.6f} ref={r:.6f}"
                else:
                    assert _rel_err(t, r) < 0.05, \
                        f"Zone {indices[i]} {key}: triton={t:.6f} ref={r:.6f} rel={_rel_err(t,r)*100:.2f}%"

    def test_triton_mu_all_zones(self, reference_data, device, dtype):
        """mu from Triton should match reference within 5% average."""
        rd = reference_data
        log_P = torch.tensor(rd['log_P'], dtype=dtype, device=device)
        log_T = torch.tensor(rd['log_T'], dtype=dtype, device=device)

        result = triton_batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        mu = result['mu'].cpu().numpy()
        ref_mu = rd['mu']

        rel_errs = np.abs(mu - ref_mu) / np.maximum(np.abs(ref_mu), 1e-10)
        mean_rel = np.mean(rel_errs) * 100
        print(f"  Triton mu: mean_rel_err={mean_rel:.4f}%")
        assert mean_rel < 10.0, f"Mean mu relative error too large: {mean_rel:.4f}%"

    def test_triton_large_batch(self, device, dtype):
        """Triton should handle large batches (100k shells)."""
        N = 100_000
        rng = np.random.default_rng(42)
        log_T = torch.tensor(rng.uniform(3.7, 7.15, N), dtype=dtype, device=device)
        log_P = torch.tensor(rng.uniform(0.5, 17.3, N), dtype=dtype, device=device)

        result = triton_batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        assert result['HII'].shape == (N,)
        assert result['mu'].shape == (N,)
        # All mu should be positive and reasonable
        mu = result['mu'].cpu().numpy()
        assert np.all(mu > 0.1), f"Some mu too small: min={mu.min():.6f}"
        assert np.all(mu < 10.0), f"Some mu too large: max={mu.max():.6f}"
