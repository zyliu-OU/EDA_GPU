# Validation

Scientific validation command: `.env/bin/eda-gpu validate --gpu --output runs/validation`.
It performs CPU/GPU single points and analytic gradients on water and water dimer,
CPU/GPU geometry convergence checks for both, all nine central finite differences
on water (step 1e-3 bohr), and CPU/GPU frozen oxygen optimizations. Numerical tolerances:
energy difference <1e-7 Eh; gradient max difference <1e-6 Eh/bohr;
finite difference max error <2e-5 Eh/bohr; frozen displacement <1e-5 angstrom.
These absolute tolerances accommodate double-precision integral/quadrature ordering
and finite-difference truncation without hiding chemically significant discrepancies.
They are acceptance thresholds, not claimed measured errors; see the dated result JSON.

Input tests check charge/spin errors, duplicate/missing fragment membership, invalid
frozen indices, and an actual SCF nonconvergence failure. GPU unavailability must produce
an explicit failure with no backend substitution. GPU tests on unavailable hosts should
be skipped explicitly; the hardware validation command itself never treats a skip as pass.

No EDA components have been calculated or validated. No Hessian/minimum or HPLC
validation is claimed. Small-system timings are recorded for reproducibility; they are
not reliable throughput benchmarks or measured production speedup.

## Measured on 2026-10-11

All nine checks passed on RTX 4080 SUPER (16 GiB), Ryzen 9 5950X, Ubuntu 26.04,
NVIDIA driver 595.91.07, isolated Python 3.12.13. Exact settings/versions and full
results: `reports/validation_2026-10-11/result.json`. Small final XYZ/evaluation
trajectories and per-optimization JSON are preserved beside it.

| System | CPU/GPU energy difference (Eh) | Maximum gradient difference (Eh/bohr) | Optimization |
| --- | ---: | ---: | --- |
| water | 7.105e-14 | 1.164e-07 | passed CPU/GPU |
| water_dimer | 9.209e-12 | 7.343e-07 | passed CPU/GPU |

Water finite-difference maximum error: 3.582325e-06 Eh/bohr,
all nine Cartesian components checked with central step 1e-3 bohr.
Frozen oxygen maximum displacement: CPU 1.781643e-11 angstrom,
GPU 1.781698e-11 angstrom.
Twelve local pytest tests passed, including actual GPU execution (no skip on this host).
Upstream GPU4PySCF emits SciPy 1.18 deprecation warnings; they did not change validation status.

CPU CI intentionally lacks GPU support; its marked hardware test must be reported as skipped.
Pinned remote CPU CI passed for code commit
`f769e90af9ad7c0b20e6d33afb8af8aacfcf8af7`:
https://github.com/zyliu-OU/EDA_GPU/actions/runs/38170042142.
Remote tests and CPU scientific validation succeeded; its GPU test is skipped honestly.

## Repeated energy/gradient timing

Reproduce after sourcing the CUDA environment:
`.env/bin/python scripts/benchmark.py --output runs/benchmark --repeats 3`.
One energy/gradient warmup per system/backend is excluded, followed by three fresh
SCF+analytic-gradient calls. GPU synchronization brackets timing. CPU threads=4.
Full repetitions, settings, hardware, process cumulative RSS, CuPy pool snapshots and
post-run device memory usage are in `reports/benchmark_2026-10-11.json`. Memory
snapshots are not peak GPU memory, and process RSS includes prior CPU/GPU calls.

| System | CPU median seconds | GPU median seconds |
| --- | ---: | ---: |
| water | 0.6186 | 0.3814 |
| water_dimer | 2.0931 | 0.4843 |

These tiny-system measurements include backend setup and transfers; no larger-system
speedup or GPU EDA benchmark is established. Optimizer iterations/timings can differ
because small gradient differences affect the search path. geomeTRIC, molecule
construction and output remain CPU operations.
