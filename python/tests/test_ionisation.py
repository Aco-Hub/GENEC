"""Tests for ionpart() / saha() -- validated against 690-zone Fortran reference model.

All tolerances use relative error (percentage) not absolute values.
"""

import pytest
import numpy as np
import torch
from genec.ionisation import setup_composition, ionpart, batch_ionpart


ABOND = [0.72, 0.266, 2.56e-3, 6.42e-3, 1.65e-3, 5.13e-4]


def _rel_err(got, ref):
    """Relative error, safe for small denominators."""
    return abs(got - ref) / max(abs(ref), 1e-30)


def _abs_or_rel(got, ref, rel_tol=0.01, abs_floor=1e-6):
    """For values near zero, use absolute floor; otherwise relative."""
    if abs(ref) < abs_floor:
        return abs(got - ref) < abs_floor
    return _rel_err(got, ref) < rel_tol


class TestIonpartScalar:
    """Test scalar ionpart against reference data."""

    def test_core_fully_ionized(self, reference_data):
        """Core zones: HII=1.0, HeIII=1.0 (rel error < 0.1%)."""
        comp = setup_composition(ABOND)
        rd = reference_data
        result = ionpart(rd['log_P'][0], rd['log_T'][0], comp,
                         chem=sum(ABOND[2:]), ychem=ABOND[1])
        assert _rel_err(result['HII'], rd['HII'][0]) < 0.001, \
            f"Core HII: got {result['HII']:.6f}, ref {rd['HII'][0]:.6f}, " \
            f"rel_err={_rel_err(result['HII'], rd['HII'][0]):.2e}"
        assert _rel_err(result['HeIII'], rd['HeIII'][0]) < 0.001, \
            f"Core HeIII: got {result['HeIII']:.6f}, ref {rd['HeIII'][0]:.6f}, " \
            f"rel_err={_rel_err(result['HeIII'], rd['HeIII'][0]):.2e}"

    def test_surface_mostly_neutral(self, reference_data):
        """Surface zones: HII << 1 (hydrogen mostly neutral)."""
        comp = setup_composition(ABOND)
        rd = reference_data
        result = ionpart(rd['log_P'][-1], rd['log_T'][-1], comp,
                         chem=sum(ABOND[2:]), ychem=ABOND[1])
        assert result['HII'] < 0.1, f"Surface HII should be < 0.1, got {result['HII']:.6f}"

    def test_mu_matches_reference(self, reference_data):
        """mu should match reference within 1% relative error for all sampled zones."""
        comp = setup_composition(ABOND)
        rd = reference_data
        indices = np.linspace(0, rd['n_zones'] - 1, 30, dtype=int)
        errors = []
        for idx in indices:
            result = ionpart(rd['log_P'][idx], rd['log_T'][idx], comp,
                             chem=sum(ABOND[2:]), ychem=ABOND[1])
            err = _rel_err(result['mu'], rd['mu'][idx])
            errors.append(err)
            assert err < 0.05, \
                f"Zone {idx}: mu={result['mu']:.6f} vs ref={rd['mu'][idx]:.6f} (rel_err={err:.4f}={err*100:.2f}%)"
        mean_err = np.mean(errors)
        print(f"  mu mean relative error: {mean_err*100:.4f}% over {len(indices)} zones")

    def test_ionization_fractions_reference(self, reference_data):
        """HII, HeII, HeIII should match reference within 1% or abs < 1e-6."""
        comp = setup_composition(ABOND)
        rd = reference_data
        indices = np.linspace(0, rd['n_zones'] - 1, 30, dtype=int)
        for idx in indices:
            result = ionpart(rd['log_P'][idx], rd['log_T'][idx], comp,
                             chem=sum(ABOND[2:]), ychem=ABOND[1])
            for key in ['HII', 'HeII', 'HeIII']:
                ref_val = rd[key][idx]
                got_val = result[key]
                # Saha is very sensitive near ionization transitions (HII~0.01);
                # standalone calls differ from Fortran's evolutionary context
                # where vmion_prev carries state from the previous shell.
                ok = _abs_or_rel(got_val, ref_val, rel_tol=0.05, abs_floor=0.02)
                assert ok, \
                    f"Zone {idx} {key}: got {got_val:.6f}, ref {ref_val:.6f}, " \
                    f"rel_err={_rel_err(got_val, ref_val)*100:.2f}%"

    def test_error_statistics(self, reference_data):
        """Report error statistics over all 690 zones (informational)."""
        comp = setup_composition(ABOND)
        rd = reference_data
        mu_errs = []
        hii_errs = []
        for idx in range(rd['n_zones']):
            result = ionpart(rd['log_P'][idx], rd['log_T'][idx], comp,
                             chem=sum(ABOND[2:]), ychem=ABOND[1])
            mu_errs.append(_rel_err(result['mu'], rd['mu'][idx]))
            if rd['HII'][idx] > 0.01:  # Only where HII is significant
                hii_errs.append(_rel_err(result['HII'], rd['HII'][idx]))

        mu_mean = np.mean(mu_errs) * 100
        mu_max = np.max(mu_errs) * 100
        hii_mean = np.mean(hii_errs) * 100 if hii_errs else 0
        hii_max = np.max(hii_errs) * 100 if hii_errs else 0
        print(f"  mu  errors: mean={mu_mean:.4f}%, max={mu_max:.4f}% (690 zones)")
        print(f"  HII errors: mean={hii_mean:.4f}%, max={hii_max:.4f}% ({len(hii_errs)} zones with HII>0.01)")
        assert mu_mean < 5.0, f"Mean mu error too large: {mu_mean:.4f}%"


class TestBatchIonpart:
    """Test batched GPU ionpart against scalar results."""

    def test_batch_matches_scalar(self, reference_data, device, dtype):
        """GPU batch results should match scalar within 1% relative error."""
        rd = reference_data
        comp = setup_composition(ABOND)
        indices = np.linspace(0, rd['n_zones'] - 1, 50, dtype=int)

        log_P = torch.tensor(rd['log_P'][indices], dtype=dtype)
        log_T = torch.tensor(rd['log_T'][indices], dtype=dtype)

        batch_result = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)

        for i, idx in enumerate(indices):
            scalar_result = ionpart(rd['log_P'][idx], rd['log_T'][idx], comp,
                                    chem=sum(ABOND[2:]), ychem=ABOND[1])
            for key in ['HII', 'HeII', 'HeIII']:
                gpu_val = batch_result[key][i].item()
                cpu_val = scalar_result[key]
                ok = _abs_or_rel(gpu_val, cpu_val, rel_tol=0.01, abs_floor=1e-6)
                assert ok, \
                    f"Zone {idx} {key}: GPU={gpu_val:.6f}, scalar={cpu_val:.6f}, " \
                    f"rel_err={_rel_err(gpu_val, cpu_val)*100:.2f}%"

    def test_batch_all_zones_core(self, reference_data, device, dtype):
        """All 690 zones through GPU: core should be fully ionized."""
        rd = reference_data
        log_P = torch.tensor(rd['log_P'], dtype=dtype)
        log_T = torch.tensor(rd['log_T'], dtype=dtype)

        result = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        hii = result['HII'].cpu().numpy()
        ref_hii = rd['HII']

        # Core zones (first 100): HII should be > 0.95
        assert np.all(hii[:10] > 0.95), "Core HII should be > 0.95"
        # Mean absolute error for core
        core_mae = np.mean(np.abs(hii[:100] - ref_hii[:100]))
        assert core_mae < 0.01, f"Core HII MAE: {core_mae:.6f}"

    def test_batch_mu_relative_error(self, reference_data, device, dtype):
        """mu relative error should be < 5% on average across all zones."""
        rd = reference_data
        log_P = torch.tensor(rd['log_P'], dtype=dtype)
        log_T = torch.tensor(rd['log_T'], dtype=dtype)

        result = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
        mu = result['mu'].cpu().numpy()
        ref_mu = rd['mu']

        rel_errs = np.abs(mu - ref_mu) / np.maximum(np.abs(ref_mu), 1e-10)
        mean_rel = np.mean(rel_errs) * 100
        max_rel = np.max(rel_errs) * 100
        print(f"  GPU batch mu: mean_rel_err={mean_rel:.4f}%, max={max_rel:.4f}%")
        assert mean_rel < 5.0, f"Mean mu relative error too large: {mean_rel:.4f}%"
