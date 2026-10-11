# EDA_GPU

**Status: under active development — an early research prototype.**

EDA_GPU is a GPU-assisted computational chemistry project for reproducible molecular
structure preparation, electronic-energy calculations, analytic gradients, and geometry
optimization. Its long-term goal is a faithful implementation of the original generalized
Kohn–Sham energy decomposition analysis (GKS-EDA), beginning with a validated CPU
reference and then accelerating the expensive operations on NVIDIA GPUs.

The working foundation uses **PySCF** for CPU calculations, **GPU4PySCF** for GPU
calculations, and **geomeTRIC** for optimization. Current validation covers small,
restricted closed-shell and small unrestricted gas-phase systems. APIs, supported methods, and installation
requirements may change as development proceeds.

## What works today

| Capability | Current status |
| --- | --- |
| Explicit CPU/GPU backend selection | Implemented and tested; unavailable GPU requests fail clearly |
| XYZ plus identity, protonation, charge, spin, and fragment metadata | Implemented with input consistency checks |
| DFT single points and analytic nuclear gradients | Tested with density-fitted B3LYP/def2-SVP |
| Water and hydrogen-bonded water-dimer optimization | Convergence demonstrated on CPU and GPU |
| Frozen-atom optimization | Tested on CPU and GPU; external freeze/set constraints also accepted |
| Machine-readable calculation results | Settings, units, convergence, versions, hardware, timing, and chemical metadata |
| Optimized structures and trajectories | XYZ outputs and evaluated-geometry trajectories |
| Scientific validation and CPU CI | CPU/GPU agreement, finite differences, constraints, and explicit failure handling |
| CPU EDA preparation | Counterpoise fragment references, frozen electrostatics, generalized correlation; independent algebra checks passed, published component validation unresolved |
| Complete experimental GKS-EDA | Six gas-phase terms on CPU/GPU; independent mathematical tests and component-wise CPU/GPU agreement passed; external matched-settings validation pending |

The CLI provides `inspect`, `singlepoint`, `optimize`, `validate`, and experimental
`eda-prepare` and `eda` commands.
This project implements an independently useful optimization foundation; it does not
claim a GAMESS port or a validated EDA decomposition.

## Validation evidence

On the development host (RTX 4080 SUPER, Ryzen 9 5950X, Ubuntu 26.04), **28 local
tests and nine optimization-foundation scientific checks passed**, including actual GPU execution. Water and
water-dimer CPU/GPU single-point energies differed by at most **9.21e-12 hartree**;
maximum gradient differences were **7.35e-7 hartree/bohr**. Water finite-difference
gradients agreed within **3.59e-6 hartree/bohr**, and frozen oxygen displacement was
below **1.79e-11 angstrom**. Both molecules converged on both backends.

These results apply to the recorded settings and small benchmark cases. Optimization
convergence establishes a stationary-point result; a vibrational analysis has not
established a minimum. GitHub CPU CI also passed; its GPU test was explicitly skipped
because that runner has no GPU. Warmed, repeated timing measurements are provided
without extrapolating them to larger systems.

See [validation details](docs/VALIDATION.md),
[full scientific results](reports/validation_2026-10-11/result.json), and
[benchmark measurements](reports/benchmark_2026-10-11.json).

## Development roadmap

1. Extend the structure-preparation and optimization foundation with additional
   carefully validated molecules, installation paths, and supported settings.
2. Secure the original GKS-EDA derivation, supporting material, and trustworthy
   individual component reference calculations under identical settings.
3. Specify every intermediate state and component convention, including orbital
   normalization, exchange/correlation, fragment reference geometry, and BSSE policy
   (experimental closed-shell specification now implemented).
4. Validate the implemented experimental CPU GKS-EDA against external matched-settings
   outputs, component by component.
5. Extend the implemented GPU component path after external validation; improve
   performance beyond the current small-system direct-integral measurements.
6. Explore nucleotide/RNA–triethylammonium–stationary-phase models after the scientific
   foundation is established, with explicit chemical states, solvent/sampling assumptions,
   and experimental validation before interpreting HPLC retention trends.

Complete restricted/unrestricted EDA is implemented experimentally; an externally validated EDA
release and HPLC remain future work. Total-energy closure alone will not count
as component validation, and ALMO-EDA or other decompositions will not be substituted
for the original GKS-EDA definitions. The exact current blockers are documented in
[the GKS-EDA source gate](docs/GKS_EDA_SPEC.md).

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

The optimization foundation passed nine scientific checks with actual GPU execution.
Experimental CPU EDA preparation adds independent electrostatic and correlation checks;
a validated full GKS-EDA release remains unavailable. See the recorded evidence.

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
## Scientific references and project terms

- [PySCF geometry optimization](https://pyscf.org/user/geomopt.html)
- [GPU4PySCF](https://github.com/pyscf/gpu4pyscf)
- [geomeTRIC constraint documentation](https://geometric.readthedocs.io/en/latest/constraints.html)
- [Original GKS-EDA derivation, Su et al. (2014)](https://doi.org/10.1021/jp500405s)
- [GKS-EDA review, Su, Tang, and Wu (2020)](https://doi.org/10.1002/wcms.1460)
- [Official XEDA manual and reference examples](https://xmvb.xmu.edu.cn/docs/xeda/content.html)

This repository is public for development visibility. No project license has been
assigned yet. Published literature is linked rather than redistributed as full-text papers.

## Experimental CPU EDA preparation

The supplied 2020 review was read completely. Its unambiguous definitions now support
frozen electrostatics and generalized correlation, using direct integrals and ghost-basis
counterpoise references at the complex geometry. Exactly two closed-shell, all-electron
fragments and B3LYP/PBE0 are supported by this older CPU-only command.
Use `eda` below for the complete experimental CPU/GPU component calculation.

```bash
.env/bin/eda-gpu eda-prepare examples/water_dimer.xyz --metadata examples/water_dimer.json --backend cpu --output runs/eda_dimer
.env/bin/eda-gpu eda-prepare examples/co2_dimer_xeda_ii.xyz --metadata examples/co2_dimer_xeda_ii.json --backend cpu --basis '6-31g*' --output runs/eda_co2
```

A successful preparation reports `full_eda_available: false`. Exchange, repulsion,
polarization, and total closure remain null until their intermediate-state conventions
and component references are established. Saved densities support further development.
The official XEDA II CO₂ example differs from our partial calculation beyond printed
rounding; its settings are not fully specified, so it is an exploratory comparison. See
[comparison evidence](reports/eda_preparation_2026-10-11/co2_reference_comparison.json),
[mathematical specification](docs/GKS_EDA_SPEC.md), and
[source audit](docs/EDA_SOURCE_AUDIT_2026-10-11.md). Papers remain local.

## Complete experimental CPU/GPU GKS-EDA

`eda` now computes electrostatic, exchange, repulsion, polarization (including charge
transfer), generalized correlation, and explicit dispersion (zero for the supported
non-dispersion-corrected methods). Exactly two all-electron fragments with additive spins,
B3LYP/B3LYP5/PBE0, counterpoise references and direct integrals are supported. It
uses the explicit Su–Li orbital intermediates on KS orbitals, as specified in the
supplied review; no ALMO-EDA substitution is made.

```bash
source scripts/gpu-env.sh
.env/bin/eda-gpu eda examples/water_dimer.xyz --metadata examples/water_dimer.json --backend gpu --output runs/gks_dimer_gpu
.env/bin/eda-gpu eda examples/water_dimer.xyz --metadata examples/water_dimer.json --backend cpu --output runs/gks_dimer_cpu
.env/bin/eda-gpu eda examples/co2_dimer_xeda_ii.xyz --metadata examples/co2_dimer_xeda_ii.json --backend gpu --basis '6-31g*' --cart --output runs/gks_co2_gpu
.env/bin/python scripts/validate_gks.py --gpu --repeats 3 --output runs/gks_checks
```

All terms, intermediate energies, condition checks and reference densities are saved.
Results state `all_components_computed: true`, `scientifically_validated: false`, and
`full_eda_available: false`: the last flag is the external-validation release gate.
A successful process means the experimental calculation converged and its internal
checks passed. It does not imply validated agreement with XEDA.

**28 local tests passed**, including actual GPU EDA, independent integral/PBE0
checks and a many-electron determinant energy check. Repeated CPU/GPU comparisons
on water and CO₂ dimers have a maximum individual-component difference of
**1.82e-8 hartree**. Cartesian functions explain most of the earlier CO₂ reference
discrepancy, but residual component differences remain beyond printed rounding.

GPU4PySCF direct J/K has CPU/GPU hybrid internals; this is not an all-GPU workflow.
For these tiny cases, measured complete EDA GPU medians were slower than CPU
(24.43 versus 3.15 seconds for water dimer; 29.47 versus 9.19 for CO₂). These
measurements do not establish larger-system throughput or speedup. See
[GPU EDA validation and timings](reports/gks_eda_gpu_2026-10-11/result.json),
[methods](docs/GKS_EDA_SPEC.md), and [validation details](docs/VALIDATION.md).

The next scientific gate is individual component output from XEDA with explicitly
matched functional, basis representation, quadrature, integral and counterpoise
settings. Open-shell, solvent, double-hybrid and RNA/HPLC applications remain deferred.


### Tryptamine–water paper benchmark (under development)

The benchmark follows Liu et al., PCCP 2017, DOI
[10.1039/C7CP03491A](https://doi.org/10.1039/C7CP03491A): neutral M06-2X and
radical-cation U-M06-2X with aug-cc-pVDZ and counterpoise EDA. The cation is
created by removing one electron, **not** by protonation. Restricted and unrestricted
CPU/GPU optimization and complete experimental EDA are supported. Open-shell
fragment spins must add to the complex spin; broken-symmetry/antiferromagnetic
coupling is outside the implemented scope. S² is recorded for unrestricted SCFs.

Published Table 3 data are attributed in `references/tryptamine_water_2017.json`.
The supplied four-page SI contains molecular figures and charges, but no full
Cartesian coordinates. `examples/tryptamine_water/constructed_bridge.xyz` is an
explicitly labelled approximate seed constructed from the SI bonding pattern,
not the paper's original geometry. Optional seed generation requires
`rdkit==2026.9.1`; ordinary calculations do not require RDKit.

```bash
source scripts/gpu-env.sh
.env/bin/python scripts/compare_tryptamine.py \
  --xyz examples/tryptamine_water/constructed_bridge.xyz \
  --metadata examples/tryptamine_water/neutral.json \
  --case A_H2O --output runs/tryptamine_comparison
```

This runs CPU/GPU optimization sequentially from identical input, then direct
counterpoise EDA on the same GPU-optimized geometry. Outputs include coordinates,
trajectories, settings, spin diagnostics, aligned geometry RMSD, all EDA terms,
backend differences and differences from the paper. Pauli is exchange plus
repulsion, matching Table 3. Paper agreement remains exploratory until original
coordinates and numerical settings are matched. No ZPE/deformation-corrected D0,
minimum/Hessian certification, intra-EDA or individual-H-bond decomposition is
claimed. A separate H⁺–H₂O starting conformer is required to study that isomer.
