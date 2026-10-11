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

Initial CPU electrostatic and generalized-correlation terms have been calculated and
checked against independent algebra; matched-settings published component validation
has not been established. No Hessian/minimum or HPLC
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

## CPU EDA preparation, 2026-10-11

Eight additional tests cover explicit four-index integral contractions for frozen
electrostatics and diagnostic exchange, fragment permutation invariance, independently
integrated PBE0 exchange/correlation with its exact-exchange fraction, counterpoise
reference cycles, density normalization, and unsupported-scope rejection. Dense ERIs
are used only in tiny test systems; production uses direct J/K evaluation.

Recorded B3LYP/def2-SVP water-dimer and B3LYP/6-31G* CO₂-dimer preparation
results are in `reports/eda_preparation_2026-10-11/`. The CO₂ electrostatic
difference from the XEDA II manual is about 6.28e-5 Eh and generalized-correlation
difference about 7.46e-6 Eh, beyond its 1e-6-Eh printed resolution. Functional variant,
quadrature, integral approximation, thresholds and BSSE policy are not fully matched.
These comparisons do not pass the published component-validation gate. No complete
EDA decomposition or GPU EDA validation is claimed.

Full local suite after this stage: **20 passed**, including actual GPU execution
(7.53 seconds). CPU CI additionally exercises the experimental CLI with PBE0/STO-3G.

CPU CI passed for EDA code commit `7e2917676cca1c7ef0d5fe22f17d9abb0378c97c`: https://github.com/zyliu-OU/EDA_GPU/actions/runs/38171856321. The clean installation, test suite, CPU scientific validation and EDA CLI smoke check all succeeded; GPU execution is not claimed on this CPU-only runner.

## Complete experimental GPU EDA, 2026-10-11

The full local suite now passes **26 tests**, with actual GPU execution and no skips.
Six new tests check every EDA term against independent four-index integrals and QR
orthonormalization; a four-electron He₂ Slater expansion is evaluated using the full
many-electron Hamiltonian to independently test the normalized determinant energy.
Other checks cover invertible occupied-span transformations, fragment swapping,
singular occupied spans, unavailable-GPU failure, experimental CLI flags and density
archives. GPU tests first use fixed CPU orbitals/densities (tolerance 1e-8 Eh), then
independent GPU SCFs (1e-7 Eh) to separate numerical backend and SCF differences.

The complete repeated validation script uses B3LYP, direct integrals, grid 3 without
pruning, SCF tolerance 1e-12 Eh and four CPU threads. Water dimer uses spherical
def2-SVP; CO₂ dimer uses Cartesian 6-31G*. One complete warmup per backend/case
is excluded; three fresh complete reference-SCF decompositions are measured. Every
recorded GPU repetition is compared component by component to the CPU baseline,
with a 1e-7-Eh threshold; closure must be within 1e-8 Eh. The stricter fixed-density
checks are independent of this end-to-end threshold.

| System | Largest CPU/GPU component difference (Eh) | CPU median (s) | GPU median (s) |
| --- | ---: | ---: | ---: |
| Water dimer | 1.059e-9 | 4.267 | 19.182 |
| CO₂ dimer, Cartesian | 1.797e-9 | 9.510 | 28.913 |

Evidence: `reports/gks_eda_gpu_2026-10-11/result.json` includes full settings,
hardware, repetitions, component differences and memory snapshots. Snapshots are
not peak GPU memory; RSS is process cumulative. The GPU is slower on these tiny
direct-integral examples. No speedup or larger-system scaling is established.

External validation is still **not established**. Cartesian CO₂ components differ
from the rounded XEDA II values by approximately +1.10e-6 Eh (electrostatic),
+1.91e-6 (exchange), -5.40e-6 (repulsion), +4.14e-6 (polarization), and -2.43e-6
(generalized correlation). Total differs by only +3.17e-7 Eh, illustrating why total
agreement does not validate individual components. Increasing the grid to levels
5 and 7 does not remove those differences. Named B3LYP and B3LYP5 variants and
spherical/Cartesian bases were checked without changing any component equation.
Full diagnostic differences are in `reports/gks_eda_gpu_2026-10-11/reference_audit.json`.
The reference leaves basis representation, functional variant, quadrature, CP policy
and integral approximation incompletely specified. A matched XEDA output remains
required before marking the implementation scientifically validated.

## Unrestricted M06-2X extension

Two additional tests passed on the actual GPU host: water radical-cation dimer independently contracted four-index total energy and spin-resolved cross exchange, analytic meta-GGA gradient finite difference, CPU/GPU energies, gradients and all EDA components. Tolerances: energy 1e-7 Eh, gradient difference 1e-6 Eh/bohr, EDA 2e-7 Eh, finite-difference derivative 2e-4 Eh/bohr. These small-system tests are internal implementation checks; they do not establish agreement with the tryptamine paper.
