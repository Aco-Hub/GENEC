"""Benchmark the full physics pipeline with all optimization backends.

Measures: Fortran CPU, Python CPU, GPU eager, GPU compile, GPU Triton, GPU FP32.
"""

import time
import subprocess
import os
import numpy as np
import torch
from genec.ionisation import setup_composition, ionpart, batch_ionpart, compiled_batch_ionpart
from genec.eos import eos_ideal, batch_eos_ideal, compiled_batch_eos_ideal
from genec.energy import total_energy_rate, batch_total_energy, compiled_batch_total_energy

try:
    from genec.triton_saha import triton_batch_ionpart, HAS_TRITON
except ImportError:
    HAS_TRITON = False

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


# ---- Full pipeline benchmarks per backend ----

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


def _bench_pipeline_gpu_impl(n, device, ionpart_fn, eos_fn, energy_fn, dtype=torch.float64):
    """Generic GPU pipeline benchmark."""
    log_T_np, log_P_np = _make_inputs(n)
    log_T = torch.tensor(log_T_np, dtype=dtype, device=device)
    log_P = torch.tensor(log_P_np, dtype=dtype, device=device)

    # Warmup
    if n > 100:
        ionpart_fn(log_P[:100], log_T[:100], ABOND, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()

    start = time.perf_counter()
    ion = ionpart_fn(log_P, log_T, ABOND, device=device, dtype=dtype)
    eos = eos_fn(log_P, log_T, ion['mu'], device=device, dtype=dtype)
    T = torch.pow(10.0, log_T)
    rho = torch.pow(10.0, eos['log_rho'])
    energy_fn(T, rho, 0.72, 0.266, device=device, dtype=dtype)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


def bench_pipeline_gpu_eager(n, device):
    """GPU eager: standard PyTorch batch ops."""
    return _bench_pipeline_gpu_impl(n, device, batch_ionpart, batch_eos_ideal, batch_total_energy)


def bench_pipeline_gpu_compile(n, device):
    """GPU compile: torch.compile fused kernels."""
    return _bench_pipeline_gpu_impl(n, device, compiled_batch_ionpart,
                                    compiled_batch_eos_ideal, compiled_batch_total_energy)


def bench_pipeline_gpu_triton(n, device):
    """GPU Triton: fused Saha kernel + eager EOS/energy."""
    if not HAS_TRITON:
        return None
    return _bench_pipeline_gpu_impl(n, device, triton_batch_ionpart,
                                    batch_eos_ideal, batch_total_energy)


def bench_pipeline_gpu_fp32(n, device):
    """GPU FP32: ionpart in FP64, EOS+energy in FP32."""
    log_T_np, log_P_np = _make_inputs(n)
    log_T = torch.tensor(log_T_np, dtype=torch.float64, device=device)
    log_P = torch.tensor(log_P_np, dtype=torch.float64, device=device)

    if n > 100:
        batch_ionpart(log_P[:100], log_T[:100], ABOND, device=device, dtype=torch.float64)
    if device.type == 'cuda':
        torch.cuda.synchronize()

    start = time.perf_counter()
    # Ionpart always in FP64
    ion = batch_ionpart(log_P, log_T, ABOND, device=device, dtype=torch.float64)
    # EOS in FP32
    eos = batch_eos_ideal(log_P.float(), log_T.float(), ion['mu'].float(),
                          device=device, dtype=torch.float32)
    # Energy in FP32
    T = torch.pow(10.0, log_T.float())
    rho = torch.pow(10.0, eos['log_rho'])
    batch_total_energy(T, rho, 0.72, 0.266, device=device, dtype=torch.float32)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    return time.perf_counter() - start


# ---- Individual function benchmarks ----

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


def bench_ionpart_triton(n, device):
    if not HAS_TRITON:
        return None
    log_T_np, log_P_np = _make_inputs(n)
    dtype = torch.float64
    log_T = torch.tensor(log_T_np, dtype=dtype, device=device)
    log_P = torch.tensor(log_P_np, dtype=dtype, device=device)
    triton_batch_ionpart(log_P[:min(100, n)], log_T[:min(100, n)], ABOND, device=device)
    if device.type == 'cuda':
        torch.cuda.synchronize()
    start = time.perf_counter()
    triton_batch_ionpart(log_P, log_T, ABOND, device=device, dtype=dtype)
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
    'full_pipeline': (bench_pipeline_cpu, bench_pipeline_gpu_eager, bench_pipeline_fortran),
}

# Extended pipeline benchmarks (all backends)
PIPELINE_BACKENDS = {
    'gpu_eager': bench_pipeline_gpu_eager,
    'gpu_compile': bench_pipeline_gpu_compile,
    'gpu_triton': bench_pipeline_gpu_triton,
    'gpu_fp32': bench_pipeline_gpu_fp32,
}
