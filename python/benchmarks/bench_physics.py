"""Benchmark the full physics pipeline: ionpart + EOS + energy for N shells.

Measures Fortran CPU, Python CPU, Python GPU as complete per-timestep cost.
"""

import time
import subprocess
import os
import numpy as np
import torch
from genec.ionisation import setup_composition, ionpart, batch_ionpart
from genec.eos import eos_ideal, batch_eos_ideal
from genec.energy import total_energy_rate, batch_total_energy

ABOND = [0.72, 0.266, 2.56e-3, 6.42e-3, 1.65e-3, 5.13e-4]
COMP = setup_composition(ABOND)

FORTRAN_BIN = os.path.join(os.path.dirname(__file__), '..', '..', 'bin', 'bench_fortran')


def _make_inputs(n, seed=42):
    rng = np.random.default_rng(seed)
    log_T = rng.uniform(3.7, 7.15, n)
    log_P = rng.uniform(0.5, 17.3, n)
    return log_T, log_P


def run_fortran(func_name, n, timeout=600):
    if not os.path.isfile(FORTRAN_BIN):
        return None
    try:
        result = subprocess.run(
            [FORTRAN_BIN, func_name, str(n)],
            capture_output=True, text=True, timeout=timeout
        )
        if result.returncode != 0:
            return None
        return float(result.stdout.strip())
    except (subprocess.TimeoutExpired, ValueError):
        return None


# ---- Full pipeline: ionpart + EOS + energy per shell ----

def bench_pipeline_fortran(n):
    """Fortran: ionpart + eos + energy sequentially."""
    t_ion = run_fortran('ionpart', n)
    t_eos = run_fortran('eos', n)
    t_nrg = run_fortran('energy', n)
    if t_ion is None or t_eos is None or t_nrg is None:
        return None
    return t_ion + t_eos + t_nrg


def bench_pipeline_cpu(n):
    """Python CPU: ionpart + EOS + energy in scalar loop."""
    log_T, log_P = _make_inputs(n)
    start = time.perf_counter()
    for i in range(n):
        result = ionpart(log_P[i], log_T[i], COMP, chem=sum(ABOND[2:]), ychem=ABOND[1])
        mu = result['mu']
        eos_result = eos_ideal(log_P[i], log_T[i], mu)
        T = 10.0**log_T[i]
        rho = 10.0**eos_result['log_rho']
        total_energy_rate(T, rho, 0.72, 0.266)
    return time.perf_counter() - start


def bench_pipeline_gpu(n, device):
    """Python GPU: ionpart + EOS + energy batched."""
    log_T_np, log_P_np = _make_inputs(n)
    dtype = torch.float64
    log_T = torch.tensor(log_T_np, dtype=dtype, device=device)
    log_P = torch.tensor(log_P_np, dtype=dtype, device=device)

    # Warmup
    if n > 100:
        _warmup_t = log_T[:100]
        _warmup_p = log_P[:100]
        batch_ionpart(_warmup_p, _warmup_t, ABOND, device=device)
    if device.type == 'cuda':
        torch.cuda.synchronize()

    start = time.perf_counter()
    # Ionization
    ion = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
    # EOS
    eos = batch_eos_ideal(log_P, log_T, ion['mu'], device=device, dtype=dtype)
    # Energy
    T = torch.pow(10.0, log_T)
    rho = torch.pow(10.0, eos['log_rho'])
    batch_total_energy(T, rho, 0.72, 0.266, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


# ---- Individual function benchmarks (for per-function breakdown) ----

def bench_ionpart_scalar(n):
    log_T, log_P = _make_inputs(n)
    start = time.perf_counter()
    for i in range(n):
        ionpart(log_P[i], log_T[i], COMP, chem=sum(ABOND[2:]), ychem=ABOND[1])
    return time.perf_counter() - start


def bench_ionpart_gpu(n, device):
    log_T_np, log_P_np = _make_inputs(n)
    dtype = torch.float64
    log_T = torch.tensor(log_T_np, dtype=dtype, device=device)
    log_P = torch.tensor(log_P_np, dtype=dtype, device=device)
    batch_ionpart(log_P[:min(100, n)], log_T[:min(100, n)], ABOND, device=device)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    start = time.perf_counter()
    batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


def bench_eos_scalar(n):
    log_T, log_P = _make_inputs(n)
    rng = np.random.default_rng(42)
    mu = rng.uniform(0.6, 1.28, n)
    start = time.perf_counter()
    for i in range(n):
        eos_ideal(log_P[i], log_T[i], mu[i])
    return time.perf_counter() - start


def bench_eos_gpu(n, device):
    log_T_np, log_P_np = _make_inputs(n)
    rng = np.random.default_rng(42)
    dtype = torch.float64
    log_T = torch.tensor(log_T_np, dtype=dtype, device=device)
    log_P = torch.tensor(log_P_np, dtype=dtype, device=device)
    mu = torch.tensor(rng.uniform(0.6, 1.28, n), dtype=dtype, device=device)
    batch_eos_ideal(log_P[:min(100, n)], log_T[:min(100, n)], mu[:min(100, n)], device=device)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    start = time.perf_counter()
    batch_eos_ideal(log_P, log_T, mu, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


def bench_energy_scalar(n):
    rng = np.random.default_rng(42)
    T = 10.0**rng.uniform(6.5, 7.15, n)
    rho = 10.0**rng.uniform(0, 2, n)
    start = time.perf_counter()
    for i in range(n):
        total_energy_rate(T[i], rho[i], 0.72, 0.266)
    return time.perf_counter() - start


def bench_energy_gpu(n, device):
    rng = np.random.default_rng(42)
    dtype = torch.float64
    T = torch.tensor(10.0**rng.uniform(6.5, 7.15, n), dtype=dtype, device=device)
    rho = torch.tensor(10.0**rng.uniform(0, 2, n), dtype=dtype, device=device)
    batch_total_energy(T[:min(100, n)], rho[:min(100, n)], 0.72, 0.266, device=device)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    start = time.perf_counter()
    batch_total_energy(T, rho, 0.72, 0.266, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


# Registries
BENCHMARKS = {
    'ionpart': (bench_ionpart_scalar, bench_ionpart_gpu),
    'eos': (bench_eos_scalar, bench_eos_gpu),
    'energy': (bench_energy_scalar, bench_energy_gpu),
}

PIPELINE_BENCHMARKS = {
    'full_pipeline': (bench_pipeline_cpu, bench_pipeline_gpu, bench_pipeline_fortran),
}
