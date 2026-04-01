"""Collect hardware information for benchmark reports."""

import platform
import subprocess
import torch
import psutil


def get_hardware_info():
    """Return a dict of hardware information."""
    info = {}

    # CPU
    info['cpu_model'] = platform.processor() or 'Unknown'
    try:
        with open('/proc/cpuinfo') as f:
            for line in f:
                if line.startswith('model name'):
                    info['cpu_model'] = line.split(':')[1].strip()
                    break
    except FileNotFoundError:
        pass
    info['cpu_cores'] = psutil.cpu_count(logical=False) or 0
    info['cpu_threads'] = psutil.cpu_count(logical=True) or 0

    # RAM
    mem = psutil.virtual_memory()
    info['ram_gb'] = round(mem.total / (1024**3), 1)

    # GPU
    info['cuda_available'] = torch.cuda.is_available()
    if info['cuda_available']:
        info['gpu_model'] = torch.cuda.get_device_name(0)
        props = torch.cuda.get_device_properties(0)
        info['gpu_vram_gb'] = round(props.total_memory / (1024**3), 1)
        info['gpu_sm_count'] = props.multi_processor_count
        info['cuda_version'] = torch.version.cuda or 'N/A'
    else:
        info['gpu_model'] = 'None'
        info['gpu_vram_gb'] = 0
        info['gpu_sm_count'] = 0
        info['cuda_version'] = 'N/A'

    # PyTorch
    info['pytorch_version'] = torch.__version__

    # gfortran
    try:
        result = subprocess.run(['gfortran', '--version'], capture_output=True, text=True)
        info['gfortran_version'] = result.stdout.split('\n')[0].strip()
    except FileNotFoundError:
        info['gfortran_version'] = 'Not found'

    return info


def format_hardware_info(info):
    """Format hardware info as a multi-line string for graph annotations."""
    lines = [
        f"CPU: {info['cpu_model']}",
        f"RAM: {info['ram_gb']} GB",
    ]
    if info['cuda_available']:
        lines.append(f"GPU: {info['gpu_model']} ({info['gpu_vram_gb']} GB VRAM)")
        lines.append(f"CUDA: {info['cuda_version']}  PyTorch: {info['pytorch_version']}")
    else:
        lines.append("GPU: None (CPU only)")
        lines.append(f"PyTorch: {info['pytorch_version']}")
    lines.append(f"Fortran: {info['gfortran_version']}")
    return '\n'.join(lines)


if __name__ == '__main__':
    info = get_hardware_info()
    print(format_hardware_info(info))
