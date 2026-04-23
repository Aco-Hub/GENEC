"""pytest fixtures for GENEC Python tests."""

import os
import sys

import numpy as np
import pytest
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

REFERENCE_FILE = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "tests",
    "integration",
    "test_1Msol_StrucData_0000001.dat",
)


def _parse_reference_data():
    """Parse the 1Msol reference model into structured arrays."""
    data = {
        "n": [],
        "log_r": [],
        "M_int": [],
        "log_T": [],
        "log_rho": [],
        "log_P": [],
        "Cv": [],
        "dlnP_dlnrho_T": [],
        "dlnP_dlnT_rho": [],
        "nabla_e": [],
        "nabla_ad": [],
        "L_rad": [],
        "L_tot": [],
        "log_kappa": [],
        "dlnk_dlnrho_T": [],
        "dlnk_dlnT_rho": [],
        "epsilon": [],
        "dlnE_dlnrho_T": [],
        "dlnE_dlnT_rho": [],
        "X_H1": [],
        "X_He4": [],
        "mu": [],
        "mu0": [],
        "Omega": [],
        "P_turb": [],
        "V_MLT": [],
        "time_TurnOver": [],
        "HII": [],
        "HeII": [],
        "HeIII": [],
    }
    # Model metadata
    metadata = {}

    with open(REFERENCE_FILE) as f:
        lines = f.readlines()

    # Parse header (lines 0-9)
    for line in lines[:10]:
        parts = line.strip().split(":")
        if len(parts) == 2:
            key = parts[0].strip()
            val = parts[1].strip()
            metadata[key] = val

    # Skip blank line and column header (lines 10-11)
    # Parse data (lines 12+)
    keys = list(data.keys())
    for line in lines[12:]:
        vals = line.split()
        if len(vals) < len(keys):
            continue
        for i, key in enumerate(keys):
            data[key].append(float(vals[i]))

    # Convert to numpy
    for key in data:
        data[key] = np.array(data[key])

    data["metadata"] = metadata
    data["n_zones"] = len(data["n"])
    return data


@pytest.fixture(scope="session")
def reference_data():
    """Load the 690-zone reference model."""
    return _parse_reference_data()


@pytest.fixture
def device():
    """Return CUDA device if available, else CPU."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


@pytest.fixture
def dtype():
    """Return float64 to match Fortran real(8)."""
    return torch.float64


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Save test results to CSV."""
    import csv
    from pathlib import Path

    results_file = Path(config.rootdir) / "test_results.csv"

    stats = terminalreporter.stats
    passed = stats.get("passed", [])
    failed = stats.get("failed", [])

    with open(results_file, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Test", "Status", "Duration_s"])

        for rep in passed + failed:
            writer.writerow([rep.nodeid, rep.outcome, f"{rep.duration:.4f}"])

    print(f"\nTest results saved to {results_file}")
