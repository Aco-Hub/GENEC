"""Generate benchmark plots: Fortran CPU vs Python CPU vs Python GPU."""

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
GC = '#2ca02c'   # GPU green

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


def plot_main(results, hw_text, output_dir):
    """The ONE graph that matters: full pipeline, 3 lines, absolute time."""
    data = results['functions'].get('full_pipeline')
    if not data:
        print("  No full_pipeline data found, skipping main plot")
        return None

    sizes = np.array(data['batch_sizes'])
    fortran = np.array([f if f and f > 0 else np.nan for f in data['fortran_cpu']])
    cpu = np.array(data['python_cpu'])
    gpu = np.array(data['python_gpu'])

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={'width_ratios': [1.2, 1]})

    # --- Left: Absolute time, 3 lines ---
    mask_f = ~np.isnan(fortran)
    ax1.loglog(sizes[mask_f], fortran[mask_f], 'D-', color=FC, lw=2.5, ms=8, label='Fortran CPU (gfortran -O3)')
    ax1.loglog(sizes, cpu, 's-', color=PC, lw=2.5, ms=8, label='Python CPU (scalar loop)')
    ax1.loglog(sizes, gpu, '^-', color=GC, lw=2.5, ms=8, label='Python GPU (PyTorch CUDA)')

    ax1.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax1.set_ylabel('Execution Time (seconds)', fontsize=13)
    ax1.set_title('Full Physics Pipeline per Timestep', fontsize=15, fontweight='bold')
    ax1.legend(fontsize=11, loc='upper left')
    ax1.grid(True, which='major', alpha=0.3)
    ax1.grid(True, which='minor', alpha=0.1, linestyle=':')
    ax1.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _fmt(int(x))))

    # Annotate key points
    for i in [0, 4, -1]:  # 1k, 100k, 10M
        if mask_f[i]:
            ax1.annotate(f'{fortran[i]:.2f}s', (sizes[i], fortran[i]),
                         textcoords='offset points', xytext=(-10, 8), fontsize=8, color=FC)
        ax1.annotate(f'{gpu[i]:.3f}s', (sizes[i], gpu[i]),
                     textcoords='offset points', xytext=(5, -12), fontsize=8, color=GC)

    # --- Right: Speedup vs Fortran ---
    sp_gpu = fortran[mask_f] / gpu[mask_f]
    sp_cpu = fortran[mask_f] / cpu[mask_f]

    ax2.semilogx(sizes[mask_f], sp_gpu, '^-', color=GC, lw=2.5, ms=8, label='GPU vs Fortran')
    ax2.semilogx(sizes[mask_f], sp_cpu, 's-', color=PC, lw=2.5, ms=8, label='Python CPU vs Fortran')
    ax2.axhline(y=1.0, color='gray', ls='--', alpha=0.5, lw=1)
    ax2.fill_between(sizes[mask_f], 1, sp_gpu, alpha=0.1, color=GC)

    ax2.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax2.set_ylabel('Speedup Factor (x faster)', fontsize=13)
    ax2.set_title('GPU Speedup vs Fortran', fontsize=15, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, which='both', alpha=0.25)
    ax2.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: _fmt(int(x))))

    # Peak annotation
    peak = np.argmax(sp_gpu)
    ax2.annotate(f'{sp_gpu[peak]:.0f}x', (sizes[mask_f][peak], sp_gpu[peak]),
                 textcoords='offset points', xytext=(0, 12), fontsize=14,
                 fontweight='bold', color=GC, ha='center')

    fig.text(0.5, -0.04, hw_text, ha='center', fontsize=8, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fafad2', alpha=0.9, edgecolor='#ccc'))

    plt.tight_layout()
    path = os.path.join(output_dir, 'bench_full_pipeline.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_summary_bars(results, hw_text, output_dir):
    """Grouped bar chart: absolute times at 10M shells for all functions."""
    functions = results['functions']
    max_n = results['batch_sizes'][-1]

    names = []
    f_times = []
    c_times = []
    g_times = []

    for name in ['full_pipeline', 'ionpart', 'eos', 'energy']:
        if name not in functions:
            continue
        d = functions[name]
        f = d['fortran_cpu'][-1]
        names.append(TITLES.get(name, name).split('\n')[0])
        f_times.append(f if f and f > 0 else 0)
        c_times.append(d['python_cpu'][-1])
        g_times.append(d['python_gpu'][-1])

    x = np.arange(len(names))
    w = 0.25

    fig, ax = plt.subplots(figsize=(14, 6))
    bars_f = ax.bar(x - w, f_times, w, color=FC, alpha=0.85, edgecolor='#333', lw=0.5, label='Fortran CPU')
    bars_c = ax.bar(x, c_times, w, color=PC, alpha=0.85, edgecolor='#333', lw=0.5, label='Python CPU')
    bars_g = ax.bar(x + w, g_times, w, color=GC, alpha=0.85, edgecolor='#333', lw=0.5, label='Python GPU')

    for bars in [bars_f, bars_c, bars_g]:
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                txt = f'{h:.1f}s' if h >= 1 else (f'{h:.3f}s' if h >= 0.001 else f'{h*1000:.1f}ms')
                ax.text(bar.get_x() + bar.get_width()/2., h * 1.08, txt,
                        ha='center', va='bottom', fontsize=8, fontweight='bold')

    ax.set_yscale('log')
    ax.set_ylabel('Execution Time (seconds, log scale)', fontsize=12)
    ax.set_title(f'Execution Time at N={max_n:,} Stellar Shells', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=10)
    ax.legend(fontsize=11)
    ax.grid(True, axis='y', alpha=0.2)

    fig.text(0.5, -0.05, hw_text, ha='center', fontsize=8, family='monospace',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#fafad2', alpha=0.9, edgecolor='#ccc'))

    plt.tight_layout()
    path = os.path.join(output_dir, 'bench_summary.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    return path


def plot_all_scaling(results, hw_text, output_dir):
    """All functions overlaid: absolute time, 3 backends each."""
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

        f_arr = np.array([f if f and f > 0 else np.nan for f in d['fortran_cpu']])
        cpu = np.array(d['python_cpu'])
        gpu = np.array(d['python_gpu'])

        mask = ~np.isnan(f_arr)
        if mask.any():
            ax.loglog(sizes[mask], f_arr[mask], 'D', color=c, ms=ms-2, alpha=0.4, linestyle='none')
        ax.loglog(sizes, cpu, 's', color=c, ms=ms-2, alpha=0.4, linestyle='none')
        ax.loglog(sizes, gpu, ls, color=c, lw=lw, ms=ms, marker='^', label=f'{short} (GPU)')

    # Add reference lines for Fortran/CPU of full_pipeline
    if 'full_pipeline' in functions:
        d = functions['full_pipeline']
        f_arr = np.array([f if f and f > 0 else np.nan for f in d['fortran_cpu']])
        mask = ~np.isnan(f_arr)
        ax.loglog(sizes[mask], f_arr[mask], 'D--', color=FC, lw=2, ms=6, alpha=0.7, label='Pipeline Fortran')
        ax.loglog(sizes, d['python_cpu'], 's--', color=PC, lw=2, ms=6, alpha=0.7, label='Pipeline Python CPU')

    ax.set_xlabel('Number of Stellar Shells', fontsize=13)
    ax.set_ylabel('Execution Time (seconds)', fontsize=13)
    ax.set_title('All Functions Scaling: Fortran vs Python CPU vs GPU', fontsize=14, fontweight='bold')
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

    # Remove old plots
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
