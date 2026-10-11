# EDA_GPU

Working scope: restricted closed-shell gas-phase DFT structure validation, energies,
analytic gradients, and geomeTRIC optimization with explicit PySCF CPU or
GPU4PySCF GPU selection. This is not a GAMESS port. Original GKS-EDA is unavailable
until its intermediate-state definitions and matching component references are secured.

## Install

Use Python 3.12 in an isolated environment; do not upgrade system Python/drivers.
CPU: `python -m pip install -r requirements-cpu-lock.txt -e '.[test]'`.
The actual host GPU installation and exact pins are in `requirements-gpu-lock.txt`
after verified local execution. Install with `python -m pip install -r requirements-gpu-lock.txt -e '.[test]'`. CUDA wheel names depend on driver/runtime compatibility.
CUDA 12 runtimes can coexist with the host CUDA 13 toolkit. Source
`scripts/gpu-env.sh` before GPU execution; it exposes the environment's NVIDIA
wheel libraries without modifying the system loader. The script targets the verified
Python 3.12 `.env` layout. Exact host installation:

```bash
/home/zyliu/anaconda3/bin/conda create -p "$PWD/.env" python=3.12 pip -y
.env/bin/python -m pip install -r requirements-gpu-lock.txt -e '.[test]'
source scripts/gpu-env.sh
```

All nine local scientific validation checks and twelve tests passed. GPU acceleration
was actually executed; GKS-EDA components remain blocked. See the recorded evidence.

```bash
source scripts/gpu-env.sh # GPU wheel library paths
.env/bin/eda-gpu inspect
.env/bin/eda-gpu singlepoint examples/water.xyz --metadata examples/water.json --backend cpu --output runs/water_sp
.env/bin/eda-gpu optimize examples/water_dimer.xyz --metadata examples/water_dimer.json --backend gpu --output runs/dimer_gpu
.env/bin/eda-gpu optimize examples/water.xyz --metadata examples/water_frozen.json --backend cpu --output runs/frozen
.env/bin/eda-gpu validate --gpu --output runs/validation
.env/bin/python -m pytest -q
```

Output directories must be new, preserving existing work. Failed calculations write
`result.json` with an error and return nonzero. GPU failures never fall back to CPU.
`validate --gpu` requires GPU execution; without it GPU validation is explicitly absent.

JSON metadata requires identity, protonation, integer charge and spin (Nalpha-Nbeta,
not multiplicity), and a complete fragment partition with charges/spins. Atom indices
are zero based; geomeTRIC constraint files use one based indices. Bonds,
stereochemistry, and other provided chemical metadata are preserved in results.
XYZ alone cannot establish protonation or stereochemistry. No conformer generation
or HPLC retention prediction is implemented.

Results record energies in hartree, gradients in hartree/bohr, angstrom coordinates,
SCF/optimization status, method settings, runtime versions/hardware, constraints,
initial/final XYZ and evaluation trajectories. Convergence does not establish a minimum.

See [methods](docs/METHODS.md), [validation](docs/VALIDATION.md),
[GKS-EDA source gate](docs/GKS_EDA_SPEC.md), and the dated report under `reports/`.
No license is assigned pending the owner's release decision.
