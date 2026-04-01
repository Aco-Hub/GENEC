#!/bin/bash
# GENEC Python/GPU setup -- auto-detects CUDA and installs the right PyTorch
set -e

echo "=== GENEC Python/GPU Setup ==="

# Install non-torch dependencies first
pip install numpy pytest matplotlib psutil tabulate

# Detect CUDA version
CUDA_VERSION=""
if command -v nvcc &> /dev/null; then
    CUDA_VERSION=$(nvcc --version | grep "release" | sed 's/.*release //' | sed 's/,.*//')
    echo "Found CUDA $CUDA_VERSION via nvcc"
elif command -v nvidia-smi &> /dev/null; then
    CUDA_VERSION=$(nvidia-smi | grep "CUDA Version" | sed 's/.*CUDA Version: //' | sed 's/ .*//')
    echo "Found CUDA $CUDA_VERSION via nvidia-smi"
fi

if [ -z "$CUDA_VERSION" ]; then
    echo "No CUDA found. Installing PyTorch CPU-only."
    pip install torch --index-url https://download.pytorch.org/whl/cpu
else
    MAJOR=$(echo "$CUDA_VERSION" | cut -d. -f1)
    MINOR=$(echo "$CUDA_VERSION" | cut -d. -f2)

    if [ "$MAJOR" -ge 12 ]; then
        echo "Installing PyTorch for CUDA 12.x"
        pip install torch --index-url https://download.pytorch.org/whl/cu124
    elif [ "$MAJOR" -eq 11 ] && [ "$MINOR" -ge 8 ]; then
        echo "Installing PyTorch for CUDA 11.8"
        pip install torch --index-url https://download.pytorch.org/whl/cu118
    else
        echo "CUDA $CUDA_VERSION is old. Installing PyTorch for CUDA 11.8 (best effort)."
        pip install torch --index-url https://download.pytorch.org/whl/cu118
    fi
fi

# Verify
echo ""
echo "=== Verification ==="
python -c "
import torch
if torch.cuda.is_available():
    print(f'GPU ready: {torch.cuda.get_device_name(0)}')
    print(f'CUDA: {torch.version.cuda}')
else:
    print('No GPU detected. Running in CPU mode (still works, just slower).')
print(f'PyTorch: {torch.__version__}')
print()
print('Setup complete. Run:')
print('  python tools/GENEC_launch.py star_name')
"
