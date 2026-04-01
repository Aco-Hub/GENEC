# GENEC -- Geneva Stellar Evolution Code

GENEC is a one-dimensional stellar evolution code. It solves the equations of
stellar structure and evolution using the Henyey implicit integration method,
tracking a star from the zero-age main sequence through advanced burning phases.

## One-time setup (copy and paste this whole block)

```bash
# Install Fortran (needed for genec_makeini which creates initial models)
sudo apt-get install -y gfortran

# Build the Fortran tools
cd src
make genec
make genec_makeini
cd ..

# Install Python/GPU (auto-detects your CUDA version)
cd python
bash setup.sh
cd ..
```

Done. You now have:

- `bin/genec_makeini` -- creates initial stellar models
- `bin/genec` -- Fortran evolution binary (used as fallback)
- `python/genec/driver.py` -- GPU-accelerated evolution (**up to 160x faster**)

Everything is launched through one script: `tools/GENEC_launch.py`.

## Run it

```bash
# GPU mode (default -- auto-detects your NVIDIA GPU)
python tools/GENEC_launch.py star_name

# CPU mode (no GPU required)
python tools/GENEC_launch.py star_name --mode cpu

# Fortran mode (original compiled binary)
python tools/GENEC_launch.py star_name --mode fortran
```

No GPU? It falls back to the compiled Fortran binary automatically.

All existing `GENEC_launch.py` options still work (`-p`, `-m`, `-i`, `-z`, `-l`, etc.).

## Tests

```bash
cd python && make test
```

20 tests verify the Python physics against a Fortran 690-zone 1Msol reference model: ionization fractions (HII, HeII, HeIII), mean molecular weight, EOS derivatives, and nuclear energy rates.

## Benchmarks

```bash
cd python && make benchmark   # times GPU vs CPU from 1k to 10M shells (~30 min)
cd python && make plots       # generates comparison graphs in python/results/
```

## What the GPU accelerates

The per-shell physics that runs inside the Henyey solver loop at every timestep:

1. **ionpart / saha** -- Partial ionization equilibrium (iterative Saha equation for H, He, C, O, Ne, Mg)
2. **eos_ideal** -- Ideal gas + radiation equation of state (density, pressure derivatives, adiabatic gradient)
3. **energy** -- Nuclear energy generation rates (PP chain + CNO cycle)

All shells are computed in parallel on the GPU as a single batched tensor operation.

## Physics

- Stellar structure via Henyey integration (up to 5000 shells)
- Nuclear burning with configurable reaction networks (NACRE I/II, CNE, CNEO,
  GENET23, GENET48)
- Equation of state based on the Timmes--Helmholtz tables
- OPAL opacity tables (GN93, AGSS09, and other solar mixtures)
- Convection modelling
- Element diffusion and advection
- Rotation and angular momentum transport
- Magnetorotational instability (MRI)
- Stellar wind mass loss (OB, RSG, Wolf--Rayet prescriptions)
- Binary tidal effects
- Partial ionisation (Saha equation for H, He, C, O, Ne, Mg)

## Repository layout

```
GENEC/
  python/           Python/GPU evolution driver (recommended)
  src/              Fortran source code and Makefile
  src/inputs/       Opacity tables, nuclear network configs, reaction data
  src/taux/         Nuclear reaction rate files (NACRE I and II)
  src/Timmes_EOS/   Helmholtz equation of state tables
  tests/            Unit tests and integration test data
  docs/             User manual (PDF and LaTeX), physics notes
  tools/            Launcher scripts, input generators, analysis utilities
  bin/              genec, genec_makeini (built from src/)
```

## Requirements

- Python 3.9+
- PyTorch 2.0+ (GPU mode needs CUDA, CPU mode works without)
- numpy, pytest, matplotlib, psutil, tabulate
- gfortran (for building genec_makeini)

## Documentation

- `docs/GENEC.pdf` -- main user manual
- `docs/ftfp.pdf` -- nuclear burning and reaction rate format
- `docs/DiffusionCoefficients.pdf` -- diffusion physics
- `docs/Origin2016_Rotation.pdf` -- rotation implementation
- `docs/Origin2016_tidalEffects.pdf` -- binary tidal effects

## License

GENEC is released under the [GNU General Public License v3](LICENSE).
