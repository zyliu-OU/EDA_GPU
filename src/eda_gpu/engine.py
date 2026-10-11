"""Equivalent restricted DFT settings, explicit GPU failure, geomeTRIC control."""
from dataclasses import dataclass, asdict
from importlib.metadata import version, PackageNotFoundError
import platform
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from pyscf import gto, dft, lib


@dataclass(frozen=True)
class Settings:
    backend: str = "cpu"
    xc: str = "b3lyp"
    basis: str = "def2-svp"
    auxbasis: str = "def2-universal-jkfit"
    grid_level: int = 3
    conv_tol: float = 1e-10
    max_cycle: int = 100
    threads: int = 4
    density_fitting: bool = True
    dispersion: str = "none"


def hardware():
    try:
        r = subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version,memory.total", "--format=csv,noheader"], capture_output=True, text=True, timeout=10)
        gpu = r.stdout.strip() if r.returncode == 0 else r.stderr.strip()
    except (OSError, subprocess.TimeoutExpired) as e:
        gpu = str(e)
    versions = {}
    for pkg in ("EDA_GPU", "pyscf", "geometric", "numpy", "gpu4pyscf-cuda12x", "cupy-cuda12x", "cutensor-cu12"):
        try:
            versions[pkg] = version(pkg)
        except PackageNotFoundError:
            versions[pkg] = None
    cpu_model = platform.processor()
    mem_kib = None
    try:
        cpu_model = next(line.split(":", 1)[1].strip() for line in Path("/proc/cpuinfo").read_text().splitlines() if line.startswith("model name"))
        mem_kib = int(next(line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines() if line.startswith("MemTotal:")))
    except (OSError, StopIteration, ValueError):
        pass
    return {"cpu_model": cpu_model, "logical_cpus": os.cpu_count(), "memory_total_kib": mem_kib, "platform": platform.platform(), "python": platform.python_version(), "processor": platform.processor(), "gpu": gpu, "versions": versions}


def make_mf(symbols, coords, metadata, settings):
    if metadata["spin"] != 0:
        raise ValueError("initial implementation supports restricted closed-shell spin=0 only")
    if settings.backend not in ("cpu", "gpu"):
        raise ValueError("backend must be cpu or gpu")
    if settings.threads < 1 or settings.max_cycle < 1 or (not np.isfinite(settings.conv_tol) or settings.conv_tol <= 0) or settings.grid_level < 0:
        raise ValueError("invalid numerical settings")
    lib.num_threads(settings.threads)
    mol = gto.M(atom=list(zip(symbols, coords)), unit="Angstrom", charge=metadata["charge"], spin=metadata["spin"], basis=settings.basis, verbose=0)
    if settings.backend == "gpu":
        try:
            import cupy as cp
            if cp.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("no CUDA devices")
            cp.zeros(1).sum().get()
            from gpu4pyscf.dft import rks
            mf = rks.RKS(mol)
        except Exception as e:
            raise RuntimeError(f"requested GPU backend unavailable: {e}; CPU fallback is disabled") from e
    else:
        mf = dft.RKS(mol)
    mf.xc = settings.xc
    # Initial supported functional scope: global hybrids/GGA/LDA, no range separation or NLC.
    if dft.libxc.xc_type(settings.xc) == "HF":
        raise ValueError("use a supported DFT functional; HF is outside this CLI scope")
    if dft.libxc.is_nlc(settings.xc) or abs(dft.libxc.rsh_coeff(settings.xc)[0]) > 1e-14:
        raise ValueError("range-separated and nonlocal-correlation functionals are outside initial scope")
    if settings.density_fitting:
        mf = mf.density_fit(auxbasis=settings.auxbasis)
    mf.grids.level = settings.grid_level
    mf.grids.prune = None
    mf.conv_tol = settings.conv_tol
    mf.max_cycle = settings.max_cycle
    return mf


def host_array(a):
    return np.asarray(a.get() if hasattr(a, "get") else a)


def singlepoint(symbols, coords, metadata, settings):
    start = time.perf_counter()
    mf = make_mf(symbols, coords, metadata, settings)
    energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError("SCF did not converge; no successful result produced")
    grad = host_array(mf.nuc_grad_method().kernel())
    if not np.isfinite(energy) or not np.isfinite(grad).all():
        raise RuntimeError("nonfinite energy or gradient")
    return {"energy_hartree": energy, "gradient_hartree_per_bohr": grad.tolist(), "scf_converged": True, "backend_used": settings.backend, "settings": asdict(settings), "seconds": time.perf_counter() - start}


OPT_THRESHOLDS = dict(convergence_energy=1e-6, convergence_grms=3e-4,
                      convergence_gmax=4.5e-4, convergence_drms=1.2e-3,
                      convergence_dmax=1.8e-3)


def optimize(symbols, coords, metadata, settings, out, maxsteps=100, constraints=None):
    from pathlib import Path
    from tempfile import TemporaryDirectory
    from pyscf.geomopt import geometric_solver, as_pyscf_method
    from .model import xyz_text
    if maxsteps < 1:
        raise ValueError("maxsteps must be positive")
    start = time.perf_counter()
    mf = make_mf(symbols, coords, metadata, settings)
    scanner = mf.nuc_grad_method().as_scanner()
    frames = []
    def energy_gradient(mol):
        energy, gradient = scanner(mol)
        if not scanner.base.converged:
            raise RuntimeError("SCF failed during optimization")
        arr = host_array(gradient)
        if not np.isfinite(energy) or not np.isfinite(arr).all():
            raise RuntimeError("nonfinite optimization evaluation")
        frames.append({"energy_hartree": float(energy), "coordinates_angstrom": mol.atom_coords(unit="Angstrom").tolist()})
        return float(energy), arr
    wrapped = as_pyscf_method(mf.mol, energy_gradient)
    frozen = metadata.get("frozen_atoms", [])
    if constraints and frozen:
        raise ValueError("use frozen_atoms metadata or an external constraints file, not both")
    constraint_text = Path(constraints).read_text() if constraints else ("$freeze\n" + "".join(f"xyz {i+1}\n" for i in frozen) + "$end\n" if frozen else None)
    if constraint_text and "$scan" in constraint_text.lower():
        raise ValueError("constraint scans are outside this single-geometry optimizer scope")
    with TemporaryDirectory(prefix="eda-gpu-") as tmp:
        params = dict(OPT_THRESHOLDS)
        params.update(maxsteps=maxsteps, assert_convergence=True)
        if constraint_text:
            f = Path(tmp) / "constraints.txt"
            f.write_text(constraint_text)
            params.update(constraints=str(f), enforce=1e-6)
        try:
            converged, mol = geometric_solver.kernel(wrapped, **params)
        finally:
            Path(out, "trajectory.xyz").write_text("".join(xyz_text(symbols, f["coordinates_angstrom"], f"evaluation {i}; E={f['energy_hartree']:.12f} Eh") for i, f in enumerate(frames)))
    final_coords = mol.atom_coords(unit="Angstrom")
    final = singlepoint(symbols, final_coords, metadata, settings)
    frozen_error = float(np.max(np.abs(final_coords[frozen] - coords[frozen]))) if frozen else 0.0
    Path(out, "optimized.xyz").write_text(xyz_text(symbols, final_coords))
    final.update(optimization_converged=bool(converged), optimization_thresholds=OPT_THRESHOLDS,
                 constraints=constraint_text, frozen_max_displacement_angstrom=frozen_error,
                 frozen_tolerance_angstrom=1e-5, evaluations=frames,
                 stationary_point_class="constrained stationary point; minimum not established" if constraint_text else "stationary point; minimum not established",
                 seconds=time.perf_counter()-start)
    if frozen_error > 1e-5:
        final["optimization_converged"] = False
    return final
