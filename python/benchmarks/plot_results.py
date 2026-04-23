"""Generate benchmark plots: Fortran CPU vs Python CPU vs GPU backends."""

import json
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from benchmarks.hardware_info import format_hardware_info

FC = '#d62728'   # Fortran red
PC = '#1f77b4'   # Python CPU blue
GC = '#2ca02c'   # GPU eager green
CC = '#ff7f0e'   # GPU compile orange
TC = '#9467bd'   # GPU Triton purple
FP = '#17becf'   # GPU FP32 cyan

TITLES = {
    'full_pipeline': 'Full Physics Pipeline per Timestep\n(ionpart + EOS + nuclear energy)',
    'ionpart': 'Saha Ionization Solver',
    'eos': 'Equation of State',
    'energy': 'Nuclear Energy Rates (PP + CNO)',
}


def _fmt(n):
    if n >= 1e6: return f'{n/1e6:.0f}M'
    if n >= 1e3: return f'{n/1e3:.0f}k'
    return str(n)


def load_results(path):
    with open(path) as f:
        return json.load(f)


def _safe_array(data, key):
    """Convert list with possible None to numpy array with NaN."""
    if key not in data:
        return None
    return np.array([f if f and f > 0 else np.nan for f in data[key]])


def plot_main(results, hw_text, output_dir):
    """Full pipeline: all backends, absolute time + speedup."""
    data = results['functions'].get('full_pipeline')
    if not data:
        print("  No full_pipeline data found, skipping main plot")
        return None

    sizes = np.array(data['batch_sizes'])
    fortran = _safe_array(data, 'fortran_cpu')
    cpu = np.array(data['python_cpu'])
    gpu = np.array(data['python_gpu'])
    gpu_compile = _safe_array(data, 'gpu_compile')
    gpu_triton = _safe_array(data, 'gpu_triton')
    gpu_fp32 = _safe_array(data, 'gpu_fp32')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7), gridspec_kw={'width_ratios': [1.2, 1]})

    # --- Left: Absolute time ---
    mask_f = ~np.isnan(fortran) if fortran is not None else np.zeros(len(sizes), dtype=bool)
    if mask_f.any():
        ax1.loglog(sizes[mask_f], fortran[mask_f], 'D-', color=FC, lw=2.5, ms=8, label='Fortran CPU')
    ax1.loglog(sizes, cpu, 's-', color=PC, lw=2.5, ms=8, label='Python CPU')
    ax1.loglog(sizes, gpu, '^-', color=GC, lw=2.5, ms=8, label='GPU eager')

    if gpu_compile is not None:
        mask_c = ~np.isnan(gpu_compile)
        if mask_c.any():
            ax1.loglog(sizes[mask_c], gpu_compile[mask_c], 'o-', color=CC, lw=2, ms=7, label='GPU compile')

    if gpu_triton is not None:
        mask_t = ~np.isnan(gpu_triton)
        if mask_t.any():
            ax1.loglog(sizes[mask_t], gpu_triton[mask_t], 'v-', color=TC, lw=2, ms=7, label='GPU Triton')

    if gpu_fp32 is not None:
        mask_32 = ~np.isnan(gpu_fp32)
        if mask_32.any():
            ax1.loglog(sizes[mask_32], gpu_fp32[mask_32], 'P-', color=FP, lw=2, ms=7, label='GPU FP32')

    ax1.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax1.set_ylabel('Execution Time (seconds)', fontsize=13)
    ax1.set_title('Full Physics Pipeline per Timestep', fontsize=15, fontweight='bold')
    ax1.legend(fontsize=10, loc='upper left')
    ax1.grid(True, which='major', alpha=0.3)
    ax1.grid(True, which='minor', alpha=0.1, linestyle=':')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _fmt(int(x))))

    # --- Right: Speedup vs Fortran ---
    if mask_f.any():
        sp_eager = fortran[mask_f] / gpu[mask_f]
        ax2.semilogx(sizes[mask_f], sp_eager, '^-', color=GC, lw=2.5, ms=8, label='Eager vs Fortran')
        ax2.fill_between(sizes[mask_f], 1, sp_eager, alpha=0.1, color=GC)

        if gpu_compile is not None:
            sp_comp = fortran[mask_f] / gpu_compile[mask_f]
            valid = ~np.isnan(sp_comp)
            if valid.any():
                ax2.semilogx(sizes[mask_f][valid], sp_comp[valid], 'o-', color=CC, lw=2, ms=7, label='Compile vs Fortran')

        if gpu_triton is not None:
            sp_tri = fortran[mask_f] / gpu_triton[mask_f]
            valid = ~np.isnan(sp_tri)
            if valid.any():
                ax2.semilogx(sizes[mask_f][valid], sp_tri[valid], 'v-', color=TC, lw=2, ms=7, label='Triton vs Fortran')

        if gpu_fp32 is not None:
            sp_fp32 = fortran[mask_f] / gpu_fp32[mask_f]
            valid = ~np.isnan(sp_fp32)
            if valid.any():
                ax2.semilogx(sizes[mask_f][valid], sp_fp32[valid], 'P-', color=FP, lw=2, ms=7, label='FP32 vs Fortran')

        ax2.axhline(y=1.0, color='gray', ls='--', alpha=0.5, lw=1)

    ax2.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax2.set_ylabel('Speedup Factor (x faster)', fontsize=13)
    ax2.set_title('GPU Speedup vs Fortran', fontsize=15, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, which='both', alpha=0.25)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _fmt(int(x))))

    fig.text(0.5, -0.04, hw_text, ha='center', fontsize=8, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fafad2', alpha=0.9, edgecolor='#ccc'))

    plt.tight_layout()
    path = os.path.join(output_dir, 'bench_full_pipeline.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_summary_bars(results, hw_text, output_dir):
    """Grouped bar chart: all backends at max N."""
    functions = results['functions']
    max_n = results['batch_sizes'][-1]

    data = functions.get('full_pipeline')
    if not data:
        return None

    backends = ['fortran_cpu', 'python_cpu', 'python_gpu', 'gpu_compile', 'gpu_triton', 'gpu_fp32']
    labels = ['Fortran CPU', 'Python CPU', 'GPU eager', 'GPU compile', 'GPU Triton', 'GPU FP32']
    colors = [FC, PC, GC, CC, TC, FP]

    times = []
    valid_labels = []
    valid_colors = []
    for bk, lb, cl in zip(backends, labels, colors):
        if bk in data and data[bk] and data[bk][-1] and data[bk][-1] > 0:
            times.append(data[bk][-1])
            valid_labels.append(lb)
            valid_colors.append(cl)

    x = np.arange(len(valid_labels))
    fig, ax = plt.subplots(figsize=(14, 6))

    bars = ax.bar(x, times, 0.6, color=valid_colors, alpha=0.85, edgecolor='#333', lw=0.5)

    for bar in bars:
        h = bar.get_height()
        if h > 0:
            txt = f'{h:.1f}s' if h >= 1 else (f'{h:.3f}s' if h >= 0.001 else f'{h*1000:.1f}ms')
            ax.text(bar.get_x() + bar.get_width()/2., h * 1.08, txt,
                    ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax.set_yscale('log')
    ax.set_ylabel('Execution Time (seconds, log scale)', fontsize=12)
    ax.set_title(f'Full Pipeline at N={max_n:,} Shells -- All Backends', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(valid_labels, fontsize=11, rotation=20, ha='right')
    ax.grid(True, axis='y', alpha=0.2)

    fig.text(0.5, -0.08, hw_text, ha='center', fontsize=8, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fafad2', alpha=0.9, edgecolor='#ccc'))

    plt.tight_layout()
    path = os.path.join(output_dir, 'bench_summary.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_all_scaling(results, hw_text, output_dir):
    """All functions overlaid: absolute time, all backends."""
    functions = results['functions']
    sizes = np.array(results['batch_sizes'])

    fig, ax = plt.subplots(figsize=(14, 7))

    func_styles = {
        'full_pipeline': ('-', 2.5, 8),
        'ionpart':       ('--', 1.8, 6),
        'eos':           (':', 1.5, 5),
        'energy':        ('-.', 1.5, 5),
    }
    func_colors = {
        'full_pipeline': '#333333',
        'ionpart': '#9467bd',
        'eos': '#ff7f0e',
        'energy': '#17becf',
    }

    for name in ['full_pipeline', 'ionpart', 'eos', 'energy']:
        if name not in functions:
            continue
        d = functions[name]
        ls, lw, ms = func_styles.get(name, ('-', 1.5, 5))
        c = func_colors.get(name, '#333')
        short = TITLES.get(name, name).split('\n')[0].split('(')[0].strip()

        f_arr = _safe_array(d, 'fortran_cpu')
        cpu = np.array(d['python_cpu'])
        gpu = np.array(d['python_gpu'])

        if f_arr is not None:
            mask = ~np.isnan(f_arr)
            if mask.any():
                ax.loglog(sizes[mask], f_arr[mask], 'D', color=c, ms=ms-2, alpha=0.4, linestyle='none')
        ax.loglog(sizes, cpu, 's', color=c, ms=ms-2, alpha=0.4, linestyle='none')
        ax.loglog(sizes, gpu, ls, color=c, lw=lw, ms=ms, marker='^', label=f'{short} (GPU)')

        if name == 'ionpart' and 'gpu_triton' in d:
            tri = _safe_array(d, 'gpu_triton')
            if tri is not None:
                mask_t = ~np.isnan(tri)
                if mask_t.any():
                    ax.loglog(sizes[mask_t], tri[mask_t], '--', color=TC, lw=2, ms=6,
                              marker='v', label=f'{short} (Triton)')

    if 'full_pipeline' in functions:
        d = functions['full_pipeline']
        f_arr = _safe_array(d, 'fortran_cpu')
        if f_arr is not None:
            mask = ~np.isnan(f_arr)
            if mask.any():
                ax.loglog(sizes[mask], f_arr[mask], 'D--', color=FC, lw=2, ms=6, alpha=0.7, label='Pipeline Fortran')
        ax.loglog(sizes, d['python_cpu'], 's--', color=PC, lw=2, ms=6, alpha=0.7, label='Pipeline Python CPU')

    ax.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax.set_ylabel('Execution Time (seconds)', fontsize=13)
    ax.set_title('All Functions Scaling: All Backends', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9, ncol=2)
    ax.grid(True, which='both', alpha=0.2)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _fmt(int(x))))

    fig.text(0.5, -0.03, hw_text, ha='center', fontsize=8, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fafad2', alpha=0.9, edgecolor='#ccc'))

    plt.tight_layout()
    path = os.path.join(output_dir, 'bench_scaling.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path


def main():
    results_dir = os.path.join(os.path.dirname(__file__), '..', 'results')
    results_path = os.path.join(results_dir, 'benchmark_results.json')

    if not os.path.isfile(results_path):
        print(f"No results at {results_path}. Run: python -m benchmarks.bench_runner")
        sys.exit(1)

    results = load_results(results_path)
    hw_text = format_hardware_info(results['hardware'])

    for f in os.listdir(results_dir):
        if f.endswith('.png'):
            os.remove(os.path.join(results_dir, f))

    print("Generating plots...")
    plots = []

    p = plot_main(results, hw_text, results_dir)
    if p:
        plots.append(p)
        print(f"  {os.path.basename(p)}")

    p = plot_summary_bars(results, hw_text, results_dir)
    if p:
        plots.append(p)
        print(f"  {os.path.basename(p)}")

    p = plot_all_scaling(results, hw_text, results_dir)
    if p:
        plots.append(p)
        print(f"  {os.path.basename(p)}")

    print(f"\n{len(plots)} plots saved to {results_dir}/")


if __name__ == '__main__':
    main()
