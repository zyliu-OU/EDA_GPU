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
