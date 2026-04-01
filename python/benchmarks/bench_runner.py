"""Benchmark orchestrator: Fortran CPU vs Python CPU vs Python GPU.

Benchmarks both individual physics functions and the full per-timestep pipeline.
"""

import json
import os
import sys
import statistics
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from benchmarks.hardware_info import get_hardware_info, format_hardware_info
from benchmarks.bench_physics import (
    BENCHMARKS, PIPELINE_BENCHMARKS, run_fortran
)

BATCH_SIZES = [1_000, 5_000, 10_000, 50_000, 100_000, 500_000, 1_000_000, 5_000_000, 10_000_000]
REPEATS = 3

# Fortran function names mapping
FORTRAN_FUNCS = {'ionpart': 'ionpart', 'eos': 'eos', 'energy': 'energy'}


def main():
    results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    os.makedirs(results_dir, exist_ok=True)

    hw = get_hardware_info()
    print("Hardware Info:")
    print(format_hardware_info(hw))
    print()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    print(f"Shells: {BATCH_SIZES}")
    print()

    results = {
        'hardware': hw,
        'batch_sizes': BATCH_SIZES,
        'functions': {}
    }

    # ---- Individual function benchmarks ----
    for name, (scalar_fn, gpu_fn) in sorted(BENCHMARKS.items()):
        func_results = {'batch_sizes': [], 'fortran_cpu': [], 'python_cpu': [], 'python_gpu': []}
        fortran_name = FORTRAN_FUNCS.get(name)

        print(f"--- {name} ---")
        for n in BATCH_SIZES:
            reps = 1 if n >= 5_000_000 else REPEATS
            sys.stdout.write(f"  N={n:>10,} ... ")
            sys.stdout.flush()

            # Fortran
            f_times = []
            for _ in range(reps):
                ft = run_fortran(fortran_name, n) if fortran_name else None
                if ft is not None:
                    f_times.append(ft)
            f_med = statistics.median(f_times) if f_times else None

            # Python CPU
            c_times = [scalar_fn(n) for _ in range(reps)]
            c_med = statistics.median(c_times)

            # GPU
            g_times = [gpu_fn(n, device) for _ in range(reps)]
            g_med = statistics.median(g_times)

            func_results['batch_sizes'].append(n)
            func_results['fortran_cpu'].append(f_med)
            func_results['python_cpu'].append(c_med)
            func_results['python_gpu'].append(g_med)

            f_str = f"{f_med:.4f}s" if f_med else "N/A"
            sp = f_med / g_med if f_med and g_med > 0 else 0
            print(f"F={f_str}  Py={c_med:.4f}s  GPU={g_med:.4f}s  GPU/F={sp:.1f}x")

        results['functions'][name] = func_results

    # ---- Full pipeline benchmark ----
    for name, (cpu_fn, gpu_fn, fortran_fn) in PIPELINE_BENCHMARKS.items():
        func_results = {'batch_sizes': [], 'fortran_cpu': [], 'python_cpu': [], 'python_gpu': []}

        print(f"\n=== {name} (ionpart + EOS + energy) ===")
        for n in BATCH_SIZES:
            reps = 1 if n >= 5_000_000 else REPEATS
            sys.stdout.write(f"  N={n:>10,} ... ")
            sys.stdout.flush()

            # Fortran pipeline
            f_times = []
            for _ in range(reps):
                ft = fortran_fn(n)
                if ft is not None:
                    f_times.append(ft)
            f_med = statistics.median(f_times) if f_times else None

            # Python CPU pipeline
            c_times = [cpu_fn(n) for _ in range(reps)]
            c_med = statistics.median(c_times)

            # GPU pipeline
            g_times = [gpu_fn(n, device) for _ in range(reps)]
            g_med = statistics.median(g_times)

            func_results['batch_sizes'].append(n)
            func_results['fortran_cpu'].append(f_med)
            func_results['python_cpu'].append(c_med)
            func_results['python_gpu'].append(g_med)

            f_str = f"{f_med:.4f}s" if f_med else "N/A"
            sp = f_med / g_med if f_med and g_med > 0 else 0
            print(f"F={f_str}  Py={c_med:.4f}s  GPU={g_med:.4f}s  GPU/F={sp:.1f}x")

        results['functions'][name] = func_results

    output_path = os.path.join(results_dir, 'benchmark_results.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)

    print(f"\nResults saved to {output_path}")


if __name__ == '__main__':
    main()
