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
    conv_tol_grad: float | None = None
    max_cycle: int = 100
    threads: int = 4
    density_fitting: bool = True
    dispersion: str = "none"
    cart: bool = False
    solvent: str = "none"
    solvent_eps: float = 78.3553
    solvent_lebedev_order: int = 29
    solvent_vdw_scale: float = 1.2
    solvent_monomer_cavity: str = "complex"
    cpu_core_limit: int = 0
    compute_duty_cycle: float = 1.0
    gpu_temperature_limit: float = 0.0
    gpu_resume_temperature: float = 65.0
    cpu_temperature_limit: float = 0.0
    cpu_resume_temperature: float = 75.0
    max_cooldown_seconds: float = 600.0


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


def make_mf(symbols, coords, metadata, settings,progress=None,label="SCF",cavity_mol=None):
    if settings.backend not in ("cpu", "gpu"):
        raise ValueError("backend must be cpu or gpu")
    if settings.threads < 1 or settings.max_cycle < 1 or (not np.isfinite(settings.conv_tol) or settings.conv_tol <= 0) or settings.grid_level < 0:
        raise ValueError("invalid numerical settings")
    if type(settings.cpu_core_limit) is not int or settings.cpu_core_limit<0:
        raise ValueError("CPU core limit must be a nonnegative integer")
    if settings.cpu_core_limit and hasattr(os,"sched_getaffinity"):
        physical={}
        for cpu in sorted(os.sched_getaffinity(0)):
            topology=Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
            try:
                key=(topology.joinpath("physical_package_id").read_text().strip(),topology.joinpath("core_id").read_text().strip())
            except OSError:
                key=("unknown",cpu)
            physical.setdefault(key,cpu)
        cores=set(list(physical.values())[:settings.cpu_core_limit])
        for task in Path("/proc/self/task").iterdir():
            try:os.sched_setaffinity(int(task.name),cores)
            except ProcessLookupError:pass
    lib.num_threads(min(settings.threads,settings.cpu_core_limit) if settings.cpu_core_limit else settings.threads)
    from threadpoolctl import threadpool_limits
    thread_limits = threadpool_limits(limits=min(settings.threads,settings.cpu_core_limit) if settings.cpu_core_limit else settings.threads)
    mol = gto.M(atom=list(zip(symbols, coords)), unit="Angstrom", charge=metadata["charge"], spin=metadata["spin"], basis=settings.basis, cart=settings.cart, verbose=0)
    if settings.backend == "gpu":
        try:
            import cupy as cp
            if cp.cuda.runtime.getDeviceCount() < 1:
                raise RuntimeError("no CUDA devices")
            cp.zeros(1).sum().get()
            from gpu4pyscf.dft import rks, uks
            mf = uks.UKS(mol) if metadata["spin"] else rks.RKS(mol)
        except Exception as e:
            raise RuntimeError(f"requested GPU backend unavailable: {e}; CPU fallback is disabled") from e
    else:
        mf = dft.UKS(mol) if metadata["spin"] else dft.RKS(mol)
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
    if settings.conv_tol_grad is not None:
        if not np.isfinite(settings.conv_tol_grad) or settings.conv_tol_grad<=0:raise ValueError("SCF gradient tolerance must be positive")
        mf.conv_tol_grad = settings.conv_tol_grad
    mf.max_cycle = settings.max_cycle
    mf.chkfile = None  # no checkpoint files: unused here, and HDF5 locks on them failed during long EDA runs
    from .solvent import attach
    mf = attach(mf, settings, cavity_mol)
    from .runtime import WorkloadPacer
    mf._eda_thread_limits = thread_limits
    device = cp.cuda.Device().pci_bus_id if settings.backend == "gpu" else None
    mf._eda_pacer = WorkloadPacer(settings,device=device)
    mf._eda_pacer.progress=progress
    def callback(env):
        if progress:progress.scf(env,label)
        mf._eda_pacer.checkpoint(env)
    mf.callback = callback
    return mf


def host_array(a):
    return np.asarray(a.get() if hasattr(a, "get") else a)


def singlepoint(symbols, coords, metadata, settings,out=None):
    start = time.perf_counter()
    from .progress import ProgressWriter
    progress=ProgressWriter(out) if out else None
    if progress:
        from .citations import information
        progress.current.update(information(metadata))
    if progress:progress.update("single-point SCF")
    mf = make_mf(symbols, coords, metadata, settings,progress=progress)
    energy = float(mf._eda_pacer.run(mf.kernel))
    if not mf.converged:
        raise RuntimeError("SCF did not converge; no successful result produced")
    if progress:progress.update("analytic gradient")
    grad = host_array(mf._eda_pacer.run(mf.nuc_grad_method().kernel))
    if not np.isfinite(energy) or not np.isfinite(grad).all():
        raise RuntimeError("nonfinite energy or gradient")
    return {"energy_hartree": energy, "gradient_hartree_per_bohr": grad.tolist(), "scf_converged": True, "backend_used": settings.backend, "settings": asdict(settings), "compute_limits": mf._eda_pacer.report(), "seconds": time.perf_counter() - start}


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
    from .progress import ProgressWriter
    progress=ProgressWriter(out)
    from .citations import information
    progress.current.update(information(metadata))
    progress.update("optimization starting")
    mf = make_mf(symbols, coords, metadata, settings,progress=progress,label="optimization SCF")
    scanner = mf.nuc_grad_method().as_scanner()
    frames = []
    def energy_gradient(mol):
        energy, gradient = mf._eda_pacer.run(scanner,mol)
        if not scanner.base.converged:
            raise RuntimeError("SCF failed during optimization")
        arr = host_array(gradient)
        if not np.isfinite(energy) or not np.isfinite(arr).all():
            raise RuntimeError("nonfinite optimization evaluation")
        progress.update("optimization evaluation",optimization_evaluation=len(frames)+1,energy_hartree=float(energy),
                        gradient_rms_hartree_per_bohr=float(np.sqrt(np.mean(arr**2))),
                        gradient_max_hartree_per_bohr=float(np.max(np.abs(arr))))
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
    progress.update("final optimized-geometry verification")
    final = singlepoint(symbols, final_coords, metadata, settings)
    frozen_error = float(np.max(np.abs(final_coords[frozen] - coords[frozen]))) if frozen else 0.0
    Path(out, "optimized.xyz").write_text(xyz_text(symbols, final_coords))
    final.update(optimization_compute_limits=mf._eda_pacer.report(),optimization_converged=bool(converged), optimization_thresholds=OPT_THRESHOLDS,
                 constraints=constraint_text, frozen_max_displacement_angstrom=frozen_error,
                 frozen_tolerance_angstrom=1e-5, evaluations=frames,
                 stationary_point_class="constrained stationary point; minimum not established" if constraint_text else "stationary point; minimum not established",
                 seconds=time.perf_counter()-start)
    if frozen_error > 1e-5:
        final["optimization_converged"] = False
    return final
