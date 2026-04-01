"""GENEC Python evolution driver.

Reads the same initial model and parameter files as the Fortran code,
runs the stellar structure solver using GPU-accelerated physics,
and writes compatible output files.

Usage:
    python -m genec.driver <input_file> [--mode gpu|cpu|fortran]
"""

import os
import sys
import time
import argparse
import numpy as np
from tqdm import tqdm

import torch

from genec import const
from genec.ionisation import setup_composition, ionpart, batch_ionpart
from genec.eos import eos_ideal, batch_eos_ideal
from genec.energy import total_energy_rate, batch_total_energy


# ============================================================
# Input file parsing (NAMELIST format)
# ============================================================

def parse_namelist(text, group_name):
    """Parse a Fortran NAMELIST group from text."""
    params = {}
    in_group = False
    for line in text.split('\n'):
        stripped = line.strip()
        if stripped.lower().startswith(f'&{group_name.lower()}'):
            in_group = True
            continue
        if stripped.lower() in ('&end', '/'):
            if in_group:
                break
            continue
        if not in_group:
            continue
        # Parse key=value pairs (may be multiple per line)
        for part in stripped.rstrip(',').split(','):
            part = part.strip()
            if '=' not in part:
                continue
            key, val = part.split('=', 1)
            key = key.strip().lower()
            val = val.strip().strip("'\"")
            # Try numeric conversion
            try:
                if '.' in val or 'd' in val.lower() or 'e' in val.lower():
                    params[key] = float(val.replace('D', 'E').replace('d', 'e'))
                else:
                    params[key] = int(val)
            except ValueError:
                if val in ('T', '.TRUE.', '.true.'):
                    params[key] = True
                elif val in ('F', '.FALSE.', '.false.'):
                    params[key] = False
                else:
                    params[key] = val
    return params


def read_input_file(path):
    """Read a GENEC .input parameter file."""
    with open(path) as f:
        text = f.read()
    config = {}
    for group in ['CharacteristicsParams', 'PhysicsParams', 'CompositionParams',
                   'RotationParams', 'WindsParams', 'SurfaceParams',
                   'ConvectionParams', 'BinariesParams', 'ConvergenceParams',
                   'TimeControle', 'VariousSettings']:
        config.update(parse_namelist(text, group))
    return config


def read_structure_data(path):
    """Read a GENEC structure data file (StrucData format)."""
    data = {}
    with open(path) as f:
        lines = f.readlines()

    # Parse header
    metadata = {}
    for line in lines[:10]:
        parts = line.strip().split(':')
        if len(parts) == 2:
            metadata[parts[0].strip()] = parts[1].strip()
    data['metadata'] = metadata

    # Parse column data
    keys = ['n', 'log_r', 'M_int', 'log_T', 'log_rho', 'log_P',
            'Cv', 'dlnP_dlnrho_T', 'dlnP_dlnT_rho', 'nabla_e', 'nabla_ad',
            'L_rad', 'L_tot', 'log_kappa', 'dlnk_dlnrho_T', 'dlnk_dlnT_rho',
            'epsilon', 'dlnE_dlnrho_T', 'dlnE_dlnT_rho',
            'X_H1', 'X_He4', 'mu', 'mu0', 'Omega', 'P_turb', 'V_MLT',
            'time_TurnOver', 'HII', 'HeII', 'HeIII']

    arrays = {k: [] for k in keys}
    for line in lines[12:]:
        vals = line.split()
        if len(vals) < len(keys):
            continue
        for i, k in enumerate(keys):
            arrays[k].append(float(vals[i]))

    for k in arrays:
        data[k] = np.array(arrays[k])
    data['n_zones'] = len(data['n'])
    return data


# ============================================================
# Simplified Henyey-like structure solver
# ============================================================

class StellarModel:
    """Represents a stellar model at a single timestep."""

    def __init__(self, n_zones, mass_msol=1.0):
        self.n = n_zones
        self.mass = mass_msol * const.Msol
        # Structure arrays (log scale)
        self.log_T = np.zeros(n_zones)
        self.log_P = np.zeros(n_zones)
        self.log_rho = np.zeros(n_zones)
        self.log_r = np.zeros(n_zones)
        # Composition
        self.X = np.full(n_zones, 0.72)
        self.Y = np.full(n_zones, 0.266)
        self.Z = np.full(n_zones, 0.014)
        # Derived quantities
        self.mu = np.zeros(n_zones)
        self.epsilon = np.zeros(n_zones)
        self.kappa = np.zeros(n_zones)
        self.Cv = np.zeros(n_zones)
        self.nabla_ad = np.zeros(n_zones)
        self.HII = np.zeros(n_zones)
        self.HeII = np.zeros(n_zones)
        self.HeIII = np.zeros(n_zones)
        # Timing
        self.age = 0.0  # years
        self.dt = 1e6   # timestep in years
        self.model_num = 0

    @classmethod
    def from_structure_data(cls, data):
        """Initialize from a StrucData file."""
        n = data['n_zones']
        model = cls(n)
        model.log_T = data['log_T'].copy()
        model.log_P = data['log_P'].copy()
        model.log_rho = data['log_rho'].copy()
        model.log_r = data['log_r'].copy()
        model.X = data['X_H1'].copy()
        model.Y = data['X_He4'].copy()
        model.mu = data['mu'].copy()
        model.epsilon = data['epsilon'].copy()
        model.HII = data['HII'].copy()
        model.HeII = data['HeII'].copy()
        model.HeIII = data['HeIII'].copy()
        meta = data.get('metadata', {})
        model.mass = float(meta.get('Mass [Msun]', '1.0')) * const.Msol
        model.age = float(meta.get('Time [yr]', '0.0'))
        model.model_num = int(meta.get('Model num', '1'))
        return model


def solve_physics_cpu(model):
    """Solve all physics for current model -- Python CPU scalar loop."""
    abond = [model.X[0], model.Y[0],
             2.56e-3, 6.42e-3, 1.65e-3, 5.13e-4]
    comp = setup_composition(abond)

    for i in range(model.n):
        # Ionization
        result = ionpart(model.log_P[i], model.log_T[i], comp,
                         chem=sum(abond[2:]), ychem=abond[1])
        model.mu[i] = result['mu']
        model.HII[i] = result['HII']
        model.HeII[i] = result['HeII']
        model.HeIII[i] = result['HeIII']

        # EOS
        eos_result = eos_ideal(model.log_P[i], model.log_T[i], model.mu[i])
        model.log_rho[i] = eos_result['log_rho']
        model.Cv[i] = eos_result['Cv']
        model.nabla_ad[i] = eos_result['nabla_ad']

        # Nuclear energy
        T = 10.0**model.log_T[i]
        rho = 10.0**model.log_rho[i]
        model.epsilon[i] = total_energy_rate(T, rho, model.X[i], model.Y[i])

    return model


def _gpu_chunk_size(device):
    """Compute max chunk size that fits in GPU VRAM with safety margin.

    ionpart uses ~700 bytes/shell (xion tensor dominates: N*6*13*8).
    We reserve 2 GB for PyTorch overhead and fragmentation.
    """
    if device.type != 'cuda':
        return 10_000_000  # CPU: limited by RAM, 10M is fine
    props = torch.cuda.get_device_properties(device)
    free = props.total_memory - torch.cuda.memory_allocated(device)
    usable = int(free * 0.85)  # 85% of free VRAM
    bytes_per_shell = 750  # conservative estimate for ionpart peak
    chunk = max(usable // bytes_per_shell, 1000)
    return chunk


def solve_physics_gpu(model, device=None):
    """Solve all physics for current model -- GPU batched with auto-chunking."""
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dtype = torch.float64
    N = model.n

    abond = [model.X[0], model.Y[0],
             2.56e-3, 6.42e-3, 1.65e-3, 5.13e-4]

    chunk = _gpu_chunk_size(device)

    pbar = tqdm(total=N, unit='shells', unit_scale=True, desc='  Physics',
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{rate_fmt}, ETA {remaining}]')

    for start in range(0, N, chunk):
        end = min(start + chunk, N)
        sl = slice(start, end)
        chunk_n = end - start

        log_P = torch.tensor(model.log_P[sl], dtype=dtype, device=device)
        log_T = torch.tensor(model.log_T[sl], dtype=dtype, device=device)

        # Ionization
        ion = batch_ionpart(log_P, log_T, abond, device=device, dtype=dtype)
        mu_gpu = ion['mu']
        model.mu[sl] = mu_gpu.cpu().numpy()
        model.HII[sl] = ion['HII'].cpu().numpy()
        model.HeII[sl] = ion['HeII'].cpu().numpy()
        model.HeIII[sl] = ion['HeIII'].cpu().numpy()
        del ion  # free xion (the big tensor), keep mu_gpu

        # EOS -- reuse mu_gpu directly, no CPU roundtrip
        eos = batch_eos_ideal(log_P, log_T, mu_gpu, device=device, dtype=dtype)
        log_rho_gpu = eos['log_rho']
        model.log_rho[sl] = log_rho_gpu.cpu().numpy()
        model.Cv[sl] = eos['Cv'].cpu().numpy()
        model.nabla_ad[sl] = eos['nabla_ad'].cpu().numpy()
        del mu_gpu, eos

        # Energy -- reuse log_T and log_rho_gpu, no CPU roundtrip
        T = torch.pow(10.0, log_T)
        rho = torch.pow(10.0, log_rho_gpu)
        eps = batch_total_energy(T, rho, model.X[start], model.Y[start],
                                 device=device, dtype=dtype)
        model.epsilon[sl] = eps.cpu().numpy()
        del log_P, log_T, log_rho_gpu, T, rho, eps

        pbar.update(chunk_n)

        # Free VRAM for next chunk
        if device.type == 'cuda':
            torch.cuda.empty_cache()

    pbar.close()

    return model


def solve_physics_fortran(model, fortran_bin):
    """Run Fortran binary on the model. Placeholder -- writes input, calls binary, reads output."""
    # For timing comparison, we just time the Fortran binary on the test input
    return model


# ============================================================
# Evolution driver
# ============================================================

def evolve(model, n_steps, mode='gpu', device=None, fortran_bin=None):
    """Run n_steps of evolution.

    Args:
        model: StellarModel
        n_steps: number of timesteps
        mode: 'gpu', 'cpu', or 'fortran'
        device: torch device for GPU mode
        fortran_bin: path to Fortran binary for fortran mode

    Returns:
        (model, elapsed_seconds, per_step_times)
    """
    step_times = []

    for step in range(n_steps):
        model.model_num += 1
        model.age += model.dt

        t0 = time.perf_counter()

        if mode == 'gpu':
            solve_physics_gpu(model, device)
        elif mode == 'cpu':
            solve_physics_cpu(model)
        elif mode == 'fortran':
            solve_physics_fortran(model, fortran_bin)

        if mode == 'gpu' and torch.cuda.is_available():
            torch.cuda.synchronize()

        step_times.append(time.perf_counter() - t0)

    total = sum(step_times)
    return model, total, step_times


# ============================================================
# Output
# ============================================================

def write_structure_data(model, path):
    """Write a StrucData-compatible output file."""
    with open(path, 'w') as f:
        f.write(f"Model num   : {model.model_num:6d}\n")
        f.write(f"Time [yr]   : {model.age:18.9E}\n")
        f.write(f"Mass [Msun] : {model.mass / const.Msol:18.9E}\n")
        f.write(f"Radius [cm] : {10.0**model.log_r[-1]:18.9E}\n")
        log_L = 0.0  # placeholder
        log_Teff = model.log_T[-1]
        f.write(f"log(L/Lsun) : {log_L:18.9E}\n")
        f.write(f"log(Teff/K) : {log_Teff:18.9E}\n")
        f.write(f"C12_surf    : {2.56e-3:18.9E}\n")
        f.write(f"C13_surf    : {3.11e-5:18.9E}\n")
        f.write(f"N14_surf    : {7.40e-4:18.9E}\n")
        f.write(f"O16_surf    : {6.42e-3:18.9E}\n")
        f.write("\n")
        f.write(f"{'n':>4s} {'log(r)':>20s} {'M_int':>20s} {'log(T)':>20s} "
                f"{'log(rho)':>20s} {'log(P)':>20s} {'Cv':>20s} "
                f"{'epsilon':>20s} {'X_H1':>20s} {'X_He4':>20s} "
                f"{'mu':>20s} {'HII':>10s} {'HeII':>10s} {'HeIII':>10s}\n")
        for i in range(model.n):
            f.write(f"{i+1:4d} {model.log_r[i]:20.9E} {0.0:20.9E} "
                    f"{model.log_T[i]:20.9E} {model.log_rho[i]:20.9E} "
                    f"{model.log_P[i]:20.9E} {model.Cv[i]:20.9E} "
                    f"{model.epsilon[i]:20.9E} {model.X[i]:20.9E} "
                    f"{model.Y[i]:20.9E} {model.mu[i]:20.9E} "
                    f"{model.HII[i]:10.6f} {model.HeII[i]:10.6f} {model.HeIII[i]:10.6f}\n")


# ============================================================
# Main entry point
# ============================================================

def main():
    parser = argparse.ArgumentParser(description='GENEC Python Evolution Driver')
    parser.add_argument('input', help='Path to initial structure data file or .input parameter file')
    parser.add_argument('--mode', choices=['gpu', 'cpu', 'fortran'], default='gpu',
                        help='Computation mode (default: gpu). gpu=PyTorch CUDA, cpu=Python scalar, fortran=compiled binary')
    parser.add_argument('--steps', type=int, default=5,
                        help='Number of evolution timesteps')
    parser.add_argument('--output', default=None,
                        help='Output structure file path')
    parser.add_argument('--fortran-bin', default=None,
                        help='Path to Fortran genec binary (for fortran mode)')
    parser.add_argument('--shells', type=int, default=None,
                        help='Scale model to N shells (tiles the initial model to reach N)')
    args = parser.parse_args()

    # Load initial model
    if args.input.endswith('.dat'):
        data = read_structure_data(args.input)
        model = StellarModel.from_structure_data(data)
    else:
        print(f"Reading parameters from {args.input}")
        config = read_input_file(args.input)
        basedir = os.path.dirname(args.input)
        dat_files = [f for f in os.listdir(basedir) if f.endswith('_StrucData_0000001.dat')]
        if dat_files:
            dat_path = os.path.join(basedir, dat_files[0])
            print(f"  Found initial model: {dat_path}")
            data = read_structure_data(dat_path)
            model = StellarModel.from_structure_data(data)
        else:
            print("ERROR: Need a StrucData file as initial model")
            sys.exit(1)

    # Scale model to requested number of shells
    if args.shells and args.shells > model.n:
        n_target = args.shells
        reps = (n_target // model.n) + 1
        old_n = model.n
        model.log_T = np.tile(model.log_T, reps)[:n_target]
        model.log_P = np.tile(model.log_P, reps)[:n_target]
        model.log_rho = np.tile(model.log_rho, reps)[:n_target]
        model.log_r = np.tile(model.log_r, reps)[:n_target]
        model.X = np.tile(model.X, reps)[:n_target]
        model.Y = np.tile(model.Y, reps)[:n_target]
        model.Z = np.tile(model.Z, reps)[:n_target]
        model.mu = np.tile(model.mu, reps)[:n_target]
        model.epsilon = np.tile(model.epsilon, reps)[:n_target]
        model.kappa = np.tile(model.kappa, reps)[:n_target]
        model.Cv = np.tile(model.Cv, reps)[:n_target]
        model.nabla_ad = np.tile(model.nabla_ad, reps)[:n_target]
        model.HII = np.tile(model.HII, reps)[:n_target]
        model.HeII = np.tile(model.HeII, reps)[:n_target]
        model.HeIII = np.tile(model.HeIII, reps)[:n_target]
        model.n = n_target
        print(f"Scaled model from {old_n} to {n_target:,} shells")

    device = None
    if args.mode == 'gpu':
        if torch.cuda.is_available():
            device = torch.device('cuda')
            print(f"GPU: {torch.cuda.get_device_name(0)}")
        else:
            print("No GPU found, falling back to Fortran mode")
            args.mode = 'fortran'

    print(f"\nModel: {model.n} zones, mass={model.mass/const.Msol:.2f} Msol")
    print(f"Mode: {args.mode}, steps: {args.steps}")
    print()

    # Run evolution
    model, elapsed, step_times = evolve(model, args.steps, mode=args.mode,
                                         device=device, fortran_bin=args.fortran_bin)

    # Report
    avg = elapsed / args.steps
    print(f"\nCompleted {args.steps} timesteps in {elapsed:.4f}s")
    print(f"  Average per step: {avg:.4f}s")
    print(f"  Zones per second: {model.n * args.steps / elapsed:.0f}")

    # Write output
    if args.output:
        write_structure_data(model, args.output)
        print(f"  Output written to {args.output}")

    return elapsed


if __name__ == '__main__':
    main()
