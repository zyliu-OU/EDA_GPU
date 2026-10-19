# EDA_GPU build handoff — 2026-10-11

Repository: https://github.com/zyliu-OU/EDA_GPU (personal account; public visibility authorized by the owner after the initial private build).
Workspace: /home/zyliu/Documents/EDA_GPU. Original supplied handoff is preserved
locally and deliberately not included in published code (contains unrelated project
context and attachment identifiers). No papers, credentials or employer data published.
Initial tested stage committed/pushed as 7dfcad8; this report accompanies the final
validation stage. Use `git log -1` and `git ls-remote origin refs/heads/main` for final SHA.

## Completed and tested

CLI: inspect, singlepoint, optimize, validate. Explicit CPU/GPU selection, restricted
closed-shell gas-phase DFT, density fitting, analytic gradients, geomeTRIC convergence
reporting, XYZ plus explicit identity/protonation/charge/spin/fragment metadata,
frozen atoms and external freeze/set constraints, reproducible settings and JSON results.
Initial/final XYZ and evaluation trajectories are retained. SCF failure, optimizer
iteration limits, missing GPU and existing output directories fail explicitly.

12 local pytest tests passed, including one actual marked GPU test (no GPU skip).
Nine scientific checks passed in reports/validation_2026-10-11/result.json:
- Water and hydrogen-bonded water dimer CPU/GPU energies and gradients agreed.
- Water and dimer optimized to converged stationary points on both backends.
- Water all-nine-component finite differences passed.
- Frozen oxygen examples passed on CPU and GPU.

Largest single-point energy discrepancy: 9.208634e-12 Eh.
Largest single-point gradient discrepancy: 7.342964e-07 Eh/bohr.
Water finite-difference gradient error: 3.582325e-06 Eh/bohr.
Largest frozen atom displacement: 1.781698e-11 angstrom.

Dependency consistency check passed; wheel built as eda_gpu-0.1.0-py3-none-any.whl.
Repeated warmed benchmarks (3 repetitions, 4 CPU threads) and memory snapshots are
in reports/benchmark_2026-10-11.json. These do not establish production speedups.

## Reproducible environment

Ubuntu 26.04 LTS; Ryzen 9 5950X (32 logical CPUs), approximately 30 GiB RAM;
RTX 4080 SUPER 16 GiB, driver 595.91.07. Python 3.12.13 isolated under `.env`.
PySCF 2.14.0, GPU4PySCF 1.8.1, geomeTRIC 1.1.1, CuPy 14.2.0, cuTensor 2.8.1.
Exact CPU/GPU pins: requirements-cpu-lock.txt and requirements-gpu-lock.txt.
CUDA 12 pip wheel runtime libraries coexist with system CUDA 13 toolkit; source
scripts/gpu-env.sh before GPU work. No system Python/driver upgrade occurred.

Run `.env/bin/eda-gpu validate --gpu --output runs/another_validation` to reproduce.
Output directory must be new. GPU status must be confirmed by successful GPU execution.
SciPy 1.18 deprecation warnings from upstream GPU4PySCF do not invalidate current runs.

## Scientific blockers and scope

Faithful original GKS-EDA remains explicitly unavailable. Accessible ACS pages expose
abstract/paywall, not full intermediate-state derivation. Full original paper/supporting
material and reproducible matching global-hybrid individual component reference data
are needed, ideally authorized XEDA execution/source. Official XEDA manual grouped
examples were found but combine exchange/repulsion and correlation/dispersion and use
a range-separated ammonia–water example outside the initial scope. No invented EDA
values, substitute ALMO/DM decomposition or arithmetic-closure validation is present.
See docs/GKS_EDA_SPEC.md for exact missing definitions and source links.

No Hessians/minimum classification, conformer generation, solvent/open-shell/
range-separated/double-hybrid EDA, RNA/TEAH+/C18 models, sampling or HPLC retention
predictions are implemented. Chemical metadata is retained, not inferred from XYZ.

## GitHub verification

Private repository created; initial and validation stages pushed successfully.
Validated code commit: f769e90af9ad7c0b20e6d33afb8af8aacfcf8af7.
Pinned CPU GitHub Actions completed SUCCESS for that commit:
https://github.com/zyliu-OU/EDA_GPU/actions/runs/38170042142.
The job installed pinned CPU dependencies, ran tests and completed CPU scientific validation.
Local hardware validation is independent of remote CI (which lacks a GPU).
This final documentation-only delivery commit uses [skip ci] to avoid repeating unchanged
scientific code tests. Remote/local main equality is checked after pushing it.

## Public development release

The owner explicitly requested public visibility and a detailed development description.
README now documents the current prototype status, tested capabilities, numerical
validation evidence, scientific blockers, roadmap, and references. Repository metadata
identifies the project as under development. The original supplied handoff remains local.
No license choice or new scientific validation claim was introduced by this release.

## EDA preparation stage

The supplied 22-page review was read completely, including visual inspection of the
methodology equations. The official XEDA II manual supplies attributed CO₂ geometry
and individual reference components. Full-text sources stay local. Implemented CPU
`eda-prepare`: common-basis counterpoise complex/fragment SCFs, reference densities,
frozen electrostatics, generalized correlation, and explicit absent dispersion. Full
EDA is explicitly false; exchange/repulsion/polarization and closure remain null.
Eight new algebra/scope tests pass. Partial published comparisons differ beyond rounding
and are exploratory because reference settings are incompletely specified. Original
antisymmetrization/normalization conventions still need authenticated-source verification.
Next gates: secure the original derivation and matched XEDA output, implement and
validate remaining CPU terms, then port validated operations to GPU.

Final stage local verification: 20 tests passed, actual GPU test executed, no skips.

CPU CI passed for EDA code commit `7e2917676cca1c7ef0d5fe22f17d9abb0378c97c`: https://github.com/zyliu-OU/EDA_GPU/actions/runs/38171856321. The clean installation, test suite, CPU scientific validation and EDA CLI smoke check all succeeded; GPU execution is not claimed on this CPU-only runner.

## Complete experimental CPU/GPU GKS-EDA stage

Implemented `eda` with all six gas-phase terms for two closed-shell all-electron
fragments, direct integrals and counterpoise, on CPU and GPU. Su–Li 2009 explicit
orbital energy expressions plus the supplied review resolve the AS/ASN construction
for this experimental implementation; details and source caveats are documented.
GPU SCF/J/K use GPU4PySCF, occupied-space algebra uses CuPy; no CPU fallback or
all-GPU claim. The old CPU preparation API remains intact.

26 tests pass, actual GPU exercised. Independent individual-integral, PBE0
quadrature, QR and many-electron determinant checks pass. Repeated component-wise
CPU/GPU agreement passes on water and CO₂ dimers; maximum difference 1.80e-9 Eh.
Three-repetition timings after one warmup show GPU slower on these tiny cases:
water 19.18 versus CPU 4.27 seconds, CO₂ 28.91 versus CPU 9.51 seconds. Full
settings, component values, differences and memory snapshots are archived.

Cartesian basis functions explain most of the earlier CO₂ discrepancy, but external
component differences remain 1–5e-6 Eh despite grid convergence. XEDA settings
are not fully matched. Results honestly state all_components_computed=true,
scientifically_validated=false and full_eda_available=false (validated-release gate).
Next: obtain explicitly specified XEDA component outputs and confirm original 2014
source wording, validate several cases, then optimize/extend the GPU path. No solvent,
open-shell, double-hybrid or HPLC claims. Papers remain local.

## Unrestricted and tryptamine–water benchmark extension

The latest local full suite passes 28 tests, including actual GPU execution, spin-resolved
four-index energy/exchange, independent QR repulsion/polarization and a meta-GGA gradient
finite difference. Molecular commands now select UKS for nonzero spin; complete EDA
supports additive fragment spins and M06-2X. The older preparation API remains restricted.

The supplied tryptamine paper and SI were read. SI has figures and charges but no Cartesian
coordinates or PDF attachments. The neutral 27-atom constructed seed is explicitly labelled;
CPU/GPU direct M06-2X/aug-cc-pVDZ optimizations are running from the same seed. The reusable
comparison script archives both optimizers and evaluates common-geometry counterpoise EDA.
No completed tryptamine result is claimed by this implementation commit. Original geometries
are needed for matched paper replication; Table 3 is the attributed exploratory target.

## Independent Psi4 energy/gradient cross-check (2026-10-18)

Added `scripts/compare_psi4.py` and `scripts/psi4_reference.py`. Psi4 1.11 (separate
conda env `.psi4-env`, lock `requirements-psi4-conda-lock.txt`, not part of `.env`)
was compared with EDA_GPU CPU and actual GPU single points under matched exact-integral
settings, identical Libxc functional IDs and bohr coordinates. 13/13 comparisons passed
preset tolerances: HF 3e-13 Eh; B3LYP/PBE0 on water, water dimer and CO₂ dimer within
1.6e-8 Eh and 3.9e-7 Eh/bohr; M06-2X/aug-cc-pVDZ tryptamine–water 7.9e-7 Eh and
8.2e-6 Eh/bohr. Archive: `reports/psi4_crosscheck_2026-10-18/`; details in
docs/VALIDATION.md. This is the first independent-program validation of the energy and
gradient layer; it does not validate EDA components (Psi4 has no GKS-EDA). The CO₂
result localizes the XEDA component differences to EDA settings/definitions rather than
the SCF energy. Remaining EDA gate unchanged: matched XEDA component outputs.
