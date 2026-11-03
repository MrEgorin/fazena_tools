"""
FAZENA Tools — System Configuration & Environment Discovery
Handles cross-platform paths, GPU acceleration, and engine discovery.
Author: Yehor Franchuk (@MrEgorin), FAZENA
"""

import os
import sys
import platform
import shutil
from pathlib import Path
from typing import Dict, Any

# Base project directory
BASE_DIR = Path(__file__).resolve().parent.parent

# Cache and weights directory
WEIGHTS_DIR = BASE_DIR / "weights"
WEIGHTS_DIR.mkdir(exist_ok=True)

OUTPUT_DIR = BASE_DIR / "output"
OUTPUT_DIR.mkdir(exist_ok=True)

# ═════════════════════════════════════════════════════════════════════
# HARDWARE ACCELERATION DETECTION
# ═════════════════════════════════════════════════════════════════════

def get_hardware_info() -> Dict[str, Any]:
    """Detects CPU, GPU, OS platform, and CUDA capabilities."""
    info = {
        'os': platform.system(),
        'platform_release': platform.release(),
        'python_version': sys.version.split()[0],
        'has_cuda': False,
        'cuda_device_name': "None",
        'cuda_device_count': 0,
        'gpu_vram_gb': 0.0,
        'recommended_device': 'cpu'
    }

    try:
        import torch
        if torch.cuda.is_available():
            info['has_cuda'] = True
            info['cuda_device_count'] = torch.cuda.device_count()
            info['cuda_device_name'] = torch.cuda.get_device_name(0)
            mem_bytes = torch.cuda.get_device_properties(0).total_memory
            info['gpu_vram_gb'] = round(mem_bytes / (1024**3), 2)
            try:
                cap = torch.cuda.get_device_capability(0)
                info['compute_capability'] = f"sm_{cap[0]}{cap[1]}"
            except Exception:
                info['compute_capability'] = None
            # Verify if compiled kernels support this GPU architecture (e.g., sm_120)
            try:
                test_t = torch.zeros(1, device='cuda')
                _ = test_t + 1
                info['recommended_device'] = 'cuda'
                del test_t
            except Exception:
                # Kernel architecture mismatch (e.g. Blackwell sm_120 on cu121)
                info['recommended_device'] = 'cpu'
    except Exception:
        pass

    return info

def get_compute_device():
    """Returns torch.device object for hardware acceleration."""
    try:
        import torch
        if torch.cuda.is_available() and HARDWARE.get('recommended_device') == 'cuda':
            return torch.device('cuda')
    except Exception:
        pass
    return torch.device('cpu')

# ═════════════════════════════════════════════════════════════════════
# DOCKING ENGINE DISCOVERY
# ═════════════════════════════════════════════════════════════════════

def get_docking_engine_info() -> Dict[str, Any]:
    """Detects AutoDock Vina Python API or CLI binaries."""
    has_vina_py = False
    try:
        from vina import Vina
        has_vina_py = True
    except ImportError:
        pass

    has_meeko = False
    try:
        import meeko
        has_meeko = True
    except ImportError:
        pass

    # Look for CLI executables
    vina_cli = shutil.which("vina") or shutil.which("vina_gpu") or shutil.which("AutoDock-Vina")
    
    # Check common linux/windows paths if not in PATH
    if not vina_cli:
        candidates = [
            BASE_DIR / "bin" / "vina",
            BASE_DIR / "bin" / "vina.exe",
            Path("/usr/local/bin/vina"),
            Path("/usr/bin/vina"),
            Path("C:/AutoDock-Vina-develop/build/windows/vina.exe")
        ]
        for c in candidates:
            if c.exists() and os.access(c, os.X_OK):
                vina_cli = str(c)
                break

    return {
        'has_vina_python': has_vina_py,
        'has_meeko': has_meeko,
        'vina_cli_path': vina_cli,
        'docking_ready': has_vina_py or (vina_cli is not None)
    }

HARDWARE = get_hardware_info()
DOCKING_CONFIG = get_docking_engine_info()
