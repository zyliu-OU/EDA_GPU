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

An empty-beta-channel regression passed separately after the 28-test full suite; fully spin-polarized occupied spaces now archive a zero density without an invalid eigenvalue access.

## Independent Psi4 cross-check of energies and gradients, 2026-10-18

Command (after creating the separate Psi4 environment below and sourcing `scripts/gpu-env.sh`):
`.env/bin/python scripts/compare_psi4.py --gpu --output runs/psi4_crosscheck`.
Psi4 1.11 (conda-forge, Libxc 7.1.2, libint 2.13.1) is an independently developed
program with its own integral, quadrature and SCF code. It is installed apart from
`.env` with `conda create -p .psi4-env -c conda-forge --override-channels python=3.12 psi4`;
the exact package set is `requirements-psi4-conda-lock.txt`. EDA_GPU used PySCF 2.14.0
(Libxc 7.0.0) on CPU and GPU4PySCF 1.8.1 on the RTX 4080 SUPER through the normal
`singlepoint` path; HF uses PySCF directly because the CLI is DFT-only.

Matched settings: exact integrals on both sides (no density fitting; Psi4 PK, or direct
for tryptamine–water), identical Libxc functional IDs (B3LYP 402 = VWN-RPA variant,
PBE0 406, M06-2X 450+236), the same spherical/Cartesian basis representation,
unpruned grids (PySCF level 8, Psi4 150 radial × 1202 Lebedev points), no grid-weight
gradients, `no_com`/`no_reorient`, and coordinates passed to Psi4 in bohr using
PySCF's Bohr constant (0.52917721092 Å). The two programs' different Bohr constants
would otherwise shift nuclear repulsion by roughly 1e-6 Eh for the largest case.
Tolerances were fixed before running: HF 1e-8 Eh / 1e-7 Eh/bohr; GGA hybrids
2e-6 Eh / 2e-5 Eh/bohr; M06-2X 2e-5 Eh / 5e-5 Eh/bohr (maximum absolute gradient element),
reflecting the different radial/atomic-partition schemes at finite grid size.
PySCF SCF thresholds: 1e-12 Eh with 1e-7 orbital-gradient norm, except 1e-10/1e-6 for
M06-2X/aug-cc-pVDZ, where the first run showed 1e-12/1e-7 is below the GPU meta-GGA
noise floor (explicit nonconvergence failure; no result was used).

All 13 comparisons passed (`reports/psi4_crosscheck_2026-10-18/result.json`;
EDA_GPU minus Psi4):

| System | Method / basis | Functions | CPU ΔE (Eh) | GPU ΔE (Eh) | CPU / GPU max \|ΔG\| (Eh/bohr) |
| --- | --- | ---: | ---: | ---: | ---: |
| Water | PBE0 / def2-SVP | 24 | −7.3e-10 | −7.5e-10 | 6.3e-9 / 4.3e-9 |
| Water | B3LYP / def2-SVP | 24 | −8.0e-10 | −8.4e-10 | 5.4e-10 / 2.2e-9 |
| Water dimer | HF / def2-SVP | 48 | +3.4e-13 | — | 1.5e-8 / — |
| Water dimer | PBE0 / def2-SVP | 48 | −6.1e-9 | −6.0e-9 | 7.3e-8 / 7.1e-8 |
| Water dimer | B3LYP / def2-SVP | 48 | −2.1e-9 | −2.1e-9 | 1.8e-8 / 1.6e-8 |
| CO₂ dimer (XEDA II geometry) | B3LYP / 6-31G* Cartesian | 90 | −1.6e-8 | −1.6e-8 | 3.9e-7 / 3.9e-7 |
| Tryptamine–water, neutral constructed geometry | M06-2X / aug-cc-pVDZ | 425 | +7.9e-7 | +7.9e-7 | 8.2e-6 / 8.0e-6 |

Interpretation. The HF agreement (3e-13 Eh) confirms basis definitions, integrals,
nuclear repulsion and geometry handling; the DFT differences are at the expected
quadrature level, and CPU and GPU deviate from Psi4 almost identically. The CO₂ dimer
total energy agrees with Psi4 to 1.6e-8 Eh, so the 1–5e-6 Eh component differences
from the rounded XEDA II values (above) do not originate in the SCF energy layer.

Scope. This validates restricted closed-shell energies and analytic gradients only.
Psi4 has no GKS-EDA, so no EDA component, intermediate state or decomposition is
validated here; open-shell/UKS gradients and optimized geometries were not compared.
Timings in the JSON are single cold runs with different grids and programs, not
benchmarks; they establish no speedup.

## Matched-settings GAMESS-XEDA kit, 2026-10-18 (not run; deferred by project decision 2026-10-19)

`scripts/xeda_matched.py prepare --gpu --output runs/xeda_matched` writes GAMESS-XEDA
inputs and computes EDA_GPU components for the same five cases; `compare` parses the
`ES/EX/REP/POL/CORR/E` table from GAMESS logs. Inputs use only keywords documented in
the official [XEDA GAMESS manual](https://xmvb.xmu.edu.cn/docs/xeda/_downloads/222e6e0bd8ccfe3ce8fe39d17fc3eef3/XEDA_gms.pdf)
and fix every setting the earlier manual comparison left open: functional variant,
`ISPHER`, bohr coordinates (PySCF constant), `NRAD=150 NLEB=1202`, `SUPBAS=.T.`
counterpoise, tight integral/SCF thresholds. PBE0 cases are primary because they have no
VWN-variant ambiguity; B3LYP cases assume GAMESS `B3LYP` = VWN5 (`b3lyp5`), which each
log must confirm. Preset per-component tolerances: 2e-6 Eh (PBE0/B3LYP), 1e-5 Eh
(M06-2X); XEDA prints 1e-6 Eh. Kit, instructions and EDA_GPU values:
`reports/xeda_matched_2026-10-18/` (see its README).

EDA_GPU side (level-8 unpruned grid, exact integrals): CPU and GPU agree within
2.0e-8 Eh per component on all four small cases; closure residuals ≤7e-12 Eh. No
GAMESS/XEDA installation exists on the development host, so **no external component
validation is claimed from this kit until the GAMESS logs are compared.**

Diagnostic: CO₂ dimer B3LYP5 (VWN5) components differ from the rounded XEDA II manual
values by 1–3e-5 Eh, versus 1–5e-6 Eh for the VWN-RPA `b3lyp` variant reported above,
indicating the manual example used the VWN-RPA form. Tryptamine–water M06-2X components
on the level-8 grid differ from the archived level-3 GPU EDA by ≤2.4e-6 Eh (≤0.0015
kcal/mol), so grid choice cannot explain the 0.01–0.15 kcal/mol differences from the paper.

## Original tryptamine–water data search and Table 1 geometry, 2026-10-18

No original Cartesian coordinates or calculation outputs from Liu et al. (2017) were found
on the development host: the accepted manuscript (30 pages) and four-page SI contain no
coordinate tables. The manuscript states the EDA used a local GAMESS version with
counterpoise correction and M06-2X/U-M06-2X/aug-cc-pVDZ. Table 1 hydrogen-bond geometry
is the available geometric reference (`scripts/tryptamine_hbond_geometry.py`;
`reports/tryptamine_neutral_gpu_2026-10-11/hbond_geometry_comparison.json`):

| Neutral A–H₂O H-bond | EDA_GPU GPU-optimized (Å / °) | Paper Table 1 (Å / °) | Difference |
| --- | --- | --- | --- |
| 1: water O–H···N (amine) | 1.908 / 166.3 | 1.905 / 165.5 | +0.003 Å / +0.8° |
| 2: water O···H–C (indole C2–H) | 2.367 / 128.9 | 2.390 / 128.3 | −0.023 Å / +0.6° |
| 3: water O···H–C (side-chain CH₂) | 2.659 / 122.6 | 2.612 / 123.4 | +0.047 Å / −0.8° |

The constructed seed optimizes to the paper's binding motif within 0.05 Å and 1°. This is
geometric consistency, not identity: the comparison remains exploratory until original
coordinates are available. The paper's text gives H-bond-3 as 2.621 Å (Table 1: 2.612 Å)
and labels H-bonds 2/3 inconsistently between text and caption; bonds are mapped by chemistry.

## Optimized-geometry cross-check against Psi4, 2026-10-18

`scripts/compare_psi4_geometry.py` optimizes water and water dimer with EDA_GPU (geomeTRIC,
default production thresholds) and Psi4 1.11 (optking, `gau_tight`) from the same seeded
symmetry-broken start (0.05 Å uniform perturbation, seed 20261018), PBE0/def2-SVP, exact
integrals, unpruned large grids. Every final structure is classified by an analytic PySCF
Hessian at the same settings. Preset tolerances: aligned RMSD 2e-3 Å (water) / 2e-2 Å
(dimer), O–H 1e-3 Å, O···O 5e-3 Å, energy 1e-6 Eh, and both structures true minima.

| System | Aligned RMSD (Å) | ΔE (Eh) | Key distance, EDA_GPU / Psi4 (Å) | Lowest frequency, EDA_GPU / Psi4 (cm⁻¹) | Result |
| --- | ---: | ---: | --- | --- | --- |
| Water | 1.9e-5 | +3.4e-9 | O–H 0.96267 / 0.96268 | 1643.3 / 1643.4 | PASS |
| Water dimer | 1.1e-3 | +7.9e-9 | O···O 2.75039 / 2.75054 | 100.7 / 100.8 | PASS |

Archive: `reports/psi4_geometry_2026-10-18/`. Finding recorded during development: from the
exactly mirror-symmetric example dimer start, geomeTRIC converged (as reported, "minimum not
established") to a trans-linear first-order saddle (one 53.7i cm⁻¹ mode, 3.2e-4 Eh above the
cyclic Ci minimum at PBE0/def2-SVP without counterpoise), because symmetric gradients preserve
the mirror plane (`symmetric_start_saddle_diagnostic.json`). The archived 2026-10-11 validation
dimer (B3LYP/def2-SVP, density-fitted, grid 3) was checked separately and is a true minimum on
CPU and GPU (lowest frequency 105.1 cm⁻¹, no imaginary modes). Optimizer convergence still does
not establish a minimum; Hessian classification is required where minima are claimed.

## Neutral tryptamine–water CPU/GPU EDA at the common geometry, 2026-10-18

CPU direct GKS-EDA at exactly the archived GPU-optimized geometry, with the archived GPU EDA
settings (M06-2X/aug-cc-pVDZ, 425 functions, grid 3, conv 1e-9 / 1e-6; only backend and workload
pacing differ). CPU minus GPU: electrostatic 4.7e-10, exchange −1.1e-8, repulsion 1.1e-8,
polarization 7.1e-9, generalized correlation −7.3e-9, counterpoise total 2.4e-11 Eh; closure
9.7e-12 Eh. Maximum component difference 1.1e-8 Eh (7e-6 kcal/mol). Archive:
`reports/tryptamine_neutral_gpu_2026-10-11/cpu_gpu_common_geometry.json` and
`cpu_common_geometry_eda/`. Single CPU and GPU runs only; timings are not benchmarks.

## Radical-cation tryptamine–water (A⁺–H₂O), 2026-10-18

Following the paper's procedure, A⁺–H₂O was optimized from the neutral GPU-optimized geometry
with one electron removed (doublet, spin on tryptamine): GPU U-M06-2X/aug-cc-pVDZ, density fitting
with aug-cc-pvdz-jkfit, grid 3, 20 evaluations, converged; the A⁺–H₂O binding motif is retained
(no collapse to H⁺–H₂O). GPU direct counterpoise GKS-EDA (grid 3, conv 1e-9/1e-6), closure
1.7e-10 Eh. Archive: `reports/tryptamine_cation_gpu_2026-10-18/`.

| A⁺–H₂O H-bond | EDA_GPU (Å / °) | Paper Table 1 (Å / °) | Difference |
| --- | --- | --- | --- |
| 1: water O–H···N | 1.942 / 156.5 | 1.950 / 155.3 | −0.008 Å / +1.2° |
| 2: water O···H–C (indole) | 2.078 / 130.8 | 2.057 / 134.0 | +0.021 Å / −3.2° |
| 3: water O···H–C (side chain) | 2.494 / 122.0 | 2.500 / 122.6 | −0.006 Å / −0.6° |

| A⁺–H₂O, kcal/mol | EDA_GPU GPU | Paper Table 3 | Difference |
| --- | ---: | ---: | ---: |
| Electrostatic | −22.33 | −22.24 | −0.09 |
| Pauli (exchange + repulsion) | 23.68 | 23.48 | +0.20 |
| Polarization | −9.62 | −9.57 | −0.05 |
| Generalized correlation | −5.30 | −5.22 | −0.08 |
| Total | −13.56 | −13.56 | −0.00 |

Neutral and cation totals both agree with the paper within 0.012 kcal/mol, and individual
components within 0.20 kcal/mol. The component offsets have opposite signs in the two charge
states (neutral: Pauli −0.15, attractive terms 0.02–0.06 less negative; cation: Pauli +0.20,
attractive terms 0.05–0.09 more negative), which points to geometry differences between the
reconstructed and original structures rather than a systematic convention difference; grid
(≤0.0015 kcal/mol) and backend (≤7e-6 kcal/mol) are excluded. (Corrected 2026-10-19: an earlier
version stated the same sign pattern in both charge states.) Accepted by the project owner on 2026-10-19 as the external component benchmark for gas-phase
GKS-EDA (`scientifically_validated: true` in gas-phase results from that date; earlier archived
results keep `false`). Agreement tighter than ~0.2 kcal/mol per component is not demonstrated.

## Psi4 cross-check extended to density fitting and UKS, 2026-10-18

`.env/bin/python scripts/compare_psi4.py --gpu --set extended --output runs/psi4_crosscheck_extended`
uses the same matched settings and preset tolerances as the core set. Density-fitted cases use
the same auxiliary basis in both programs (Psi4 `scf_type df`, `df_basis_scf`). All 8 comparisons
passed (`reports/psi4_crosscheck_extended_2026-10-18/result.json`; EDA_GPU minus Psi4):

| System | Method / basis (auxiliary) | Reference | CPU ΔE (Eh) | GPU ΔE (Eh) | CPU / GPU max \|ΔG\| (Eh/bohr) |
| --- | --- | --- | ---: | ---: | ---: |
| Water dimer | B3LYP / def2-SVP (def2-universal-jkfit) | RKS | −2.1e-9 | −2.1e-9 | 1.9e-8 / 1.4e-8 |
| Tryptamine–water | M06-2X / aug-cc-pVDZ (aug-cc-pvdz-jkfit) | RKS | +7.5e-7 | +7.5e-7 | 3.3e-5 / 2.7e-5 |
| Water radical cation | PBE0 / def2-SVP | UKS | −6.0e-10 | −6.3e-10 | 4.0e-9 / 3.4e-9 |
| Water radical cation | M06-2X / def2-SVP | UKS | −2.7e-11 | −5.8e-11 | 1.5e-7 / 2.4e-7 |

The first two rows are the density-fitted paths actually used by `validate` and by the tryptamine
optimizations. The tryptamine density-fitted gradient difference (3.3e-5 Eh/bohr) is four times
the exact-integral value (8.2e-6) but within the preset 5e-5 tolerance; it is attributed to
differences in density-fitted gradient implementations, not investigated further. UKS energies
and gradients, including meta-GGA, agree at the quadrature level. Still not compared: EDA
components (no Psi4 GKS-EDA) and open-shell optimized geometries.

## CPU/GPU scaling benchmark (gate 6), 2026-10-18

`scripts/benchmark_scaling.py` times the two production paths with one excluded warm-up and
three repetitions per case and backend; every repetition reruns all SCFs. 16 CPU threads, no
workload pacing, grid 3; energy+gradient uses density fitting (optimizer evaluation), EDA uses
direct integrals (all reference and intermediate SCFs). Peak memory is whole-device nvidia-smi
sampling (0.25 s; idle baseline ~360 MiB). Archive: `reports/benchmark_scaling_2026-10-18.json`.

| System (basis functions) | Task | CPU median (s) | GPU median (s) | CPU/GPU | CPU−GPU value (Eh) | GPU peak (MiB) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Water dimer, B3LYP/def2-SVP (48) | energy+gradient | 1.2 | 0.5 | 2.4 | 9.4e-12 | 826 |
| Water dimer, B3LYP/def2-SVP (48) | GKS-EDA | 2.4 | 14.1 | 0.17 | 5.2e-12 | 1487 |
| Tryptamine–water, M06-2X/def2-SVP (252) | energy+gradient | 56.9 | 8.4 | 6.8 | 8.4e-12 | 3055 |
| Tryptamine–water, M06-2X/def2-SVP (252) | GKS-EDA | 169.5 | 53.9 | 3.1 | 2.1e-11 | 1837 |
| Tryptamine–water, M06-2X/aug-cc-pVDZ (425) | energy+gradient | 137.5 | 22.0 | 6.3 | 7.7e-11 | 8887 |
| Tryptamine–water, M06-2X/aug-cc-pVDZ (425) | GKS-EDA | 548.6 | 121.1 | 4.5 | 1.1e-11 | 1591 |

Repetition spread is small (GPU EDA 121.1 s in all three repetitions). For tiny systems the GPU
EDA is slower because fixed per-SCF overheads dominate; from ~250 basis functions the GPU is
3–7× faster on this host. Remaining CPU work: molecule construction, one-electron integrals for
the common-basis check, PCM/output bookkeeping and parts of GPU4PySCF's hybrid direct J/K.
The aug-cc-pVDZ EDA row was rerun separately (same settings, exclusive host) after an HDF5
checkpoint-file lock failure; checkpoint files are now disabled (`mf.chkfile = None`) because
nothing uses them. These numbers describe this host only and do not establish general speedups.

## Experimental GKS-EDA(sol), 2026-10-18

Specification and interpretations: `docs/GKS_EDA_SOL_SPEC.md`. Validation:
`.env/bin/python scripts/validate_gks_sol.py --gpu --output runs/gks_sol_validation`
(`reports/gks_sol_validation_2026-10-18/`). Water dimer, PBE0/def2-SVP, exact integrals.

| Check (tolerance fixed before running) | Result |
| --- | --- |
| Vacuum limit ε = 1.000001 vs gas-phase GKS-EDA (1e-6 Eh) | max component difference 1.6e-9 Eh; desolvation 0 — PASS |
| CPU/GPU, C-PCM (1e-7 Eh) | max component difference 9.6e-9 Eh — PASS |
| CPU/GPU, IEF-PCM (1e-7 Eh) | max component difference 9.3e-9 Eh — PASS |
| PCM analytic gradient vs central differences, water (2e-5 Eh/bohr) | 3.8e-6 (CPU), 3.7e-6 (GPU) — PASS |
| Separation screening: total → 0, −ΔG_desol/ΔE_ele → 1 − 1/ε = 0.987 (±0.02) | totals −4.0e-6, −5.4e-7, −2.4e-8 Eh at +4/+8/+12 Å; ratios 0.988, 0.992 — PASS |
| Preregistered "desolvation → 0 at separation" | FAIL, recorded and non-counting: criterion physically wrong |

The preregistered separation criterion was an error in the test design: the gas-phase-form
electrostatic term is the unscreened vacuum interaction, so at long range desolvation must
cancel the fraction (1 − 1/ε) of it while the total decays. The observed ratios confirm the
dielectric screening behaviour; the corrected criterion was fixed before the archived rerun.

Exploratory XEDA II BH₃–NH₃ C-PCM example (PBE0/6-31G, water; XEDA radii, discretization and
RI settings unknown): differences are large (desolvation −17.8 mEh, components up to 12 mEh).
Diagnostics (`diagnostics/`): the monomer-own-cavity variant is further from XEDA; enlarging
the cavity (modified Bondi × 1.0→1.8) moves every component monotonically toward XEDA (largest
difference 25 → 3.7 mEh) without converging, and our gas-phase EDA alone already differs from
XEDA's solvated total by 3.6 mEh. The comparison is unresolved; no settings were tuned to it.
The scan shows desolvation for this polar complex changes by ~13 kcal/mol across ordinary radius
choices, so application results must report cavity-radius sensitivity.

## HPLC model stages 7.1 and 7.2 (application, not validation), 2026-10-19

Models and limits: `docs/HPLC_MODELS.md`. M06-2X/def2-SVP, DF (def2-universal-jkfit), grid 3,
C-PCM water; GPU optimization then direct-integral counterpoise GKS-EDA(sol) at each optimized
geometry. Archive: `reports/hplc_ionpair_2026-10-19/`, `reports/hplc_surface_2026-10-19/`
(summary, `analysis.json`, per-case optimized XYZ and optimization/EDA results). Table from
`scripts/analyze_hplc.py`: lowest-PCM-free-energy placement per group, kcal/mol.

| Group | n | elec | exch | rep | pol | gcorr | desolv | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| residue_A···TEAH⁺ | 2/2 | -108.31 | -57.63 | 102.38 | -25.15 | -10.87 | 80.83 | -18.77 |
| residue_C···TEAH⁺ | 2/2 | -113.42 | -61.76 | 110.23 | -27.44 | -11.09 | 83.74 | -19.74 |
| residue_G···TEAH⁺ | 1/2 | -112.01 | -62.18 | 109.30 | -25.90 | -12.99 | 84.15 | -19.64 |
| residue_U···TEAH⁺ | 2/2 | -109.62 | -62.00 | 110.28 | -26.47 | -12.71 | 80.54 | -19.98 |
| C8 proxy···TEAH⁺ | 2/2 | -25.84 | -32.65 | 56.12 | -6.27 | -12.09 | 9.40 | -11.33 |
| residue_A···C8 | 1/1 | -27.95 | -57.14 | 97.05 | -14.35 | -15.51 | 2.99 | -14.92 |
| residue_C···C8 | 1/1 | -39.65 | -67.40 | 115.80 | -18.85 | -16.28 | 7.52 | -18.87 |
| residue_G···C8 | 1/1 | -35.09 | -65.79 | 110.64 | -13.93 | -19.32 | 8.06 | -15.43 |
| residue_U···C8 | 1/1 | -36.49 | -60.95 | 104.44 | -16.46 | -15.69 | 10.29 | -14.86 |

Component closure ≤ 3e-10 Eh; bonding topology preserved in every converged case. residue_G···TEAH⁺
placement 13 hit the 200-step cap on a flat surface (rms gradient 6e-5 Eh/bohr; 2.9 kcal/mol above
placement 14) and is excluded. Single local placements at def2-SVP without sampling: base differences
of ~1–4 kcal/mol are within placement/basis/cavity uncertainty and are not established.

## HPLC stage 7.3 ternary optimization (application), 2026-10-20

Residue···TEAH⁺···C8-proxy complexes (104–108 atoms; TEAH⁺ on the lower-energy phosphate O of each
optimized ion pair, C8 silanol O–H on the other non-bridging O), GFN2-xTB/ALPB pre-refinement, then GPU
M06-2X/def2-SVP, DF, C-PCM water, frozen silane anchor. Archive: `reports/hplc_ternary_2026-10-20/`.

| Complex | DFT evaluations | Wall time (h) | N–H···O–P (Å) | Si–O–H···O–P (Å) |
| --- | ---: | ---: | ---: | ---: |
| residue_A···TEAH⁺···C8 | 38 | 2.31 | 1.58 | 1.61 |
| residue_C···TEAH⁺···C8 | 35 | 1.84 | 1.58 | 1.61 |
| residue_G···TEAH⁺···C8 | 63 | 3.97 | 1.60 | 1.60 |
| residue_U···TEAH⁺···C8 | 53 | 2.83 | 1.54 | 1.62 |

All four converged with bonding topology preserved; frozen-anchor displacement ≤ 6.2e-8 Å; maximum GPU
temperature 61 °C. Constrained stationary points (no Hessian at this size; minima not established).
Single seeds (selected before the 2026-10-19 placement-ranking fix). Three-body GKS-EDA(sol) follows.

### Three-body GKS-EDA(sol) at the ternary geometries

`scripts/three_body_hplc.py` (M06-2X/def2-SVP, direct integrals, C-PCM water; R = residue, T = TEAH⁺,
S = C8 proxy). Pairs in their own dimer basis, (RT)|S and R|(TS) in the trimer basis;
ΔE₃ = ΔE((RT)|S) − ΔE(T|S) − ΔE(R|S). Closure ≤ 2.1e-10 Eh in all 20 decompositions. kcal/mol.

| Base | R···T | T···S | R···S | (RT)···S | R···(TS) | ΔE₃ elec | ΔE₃ pol | ΔE₃ desolv | ΔE₃ total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| A | -18.28 | -2.02 | -13.78 | -14.80 | -31.31 | +1.35 | +2.41 | -2.13 | +1.00 |
| C | -17.70 | -1.15 | -15.14 | -15.16 | -31.77 | +1.14 | +4.41 | -3.35 | +1.13 |
| G | -18.95 | -2.07 | -16.74 | -17.84 | -34.97 | +0.93 | +2.46 | -1.64 | +0.97 |
| U | -18.60 | -0.33 | -15.73 | -14.87 | -33.29 | +1.43 | +3.02 | -3.09 | +1.19 |

The three-body term is small and anti-cooperative (+0.97 to +1.19 kcal/mol) for every base:
polarization and electrostatics are less favourable in the trimer, partly offset by desolvation.
In this arrangement the residue binds the surface directly through the silanol O–H···O–P bond
(R···S −13.8 to −16.7) while TEAH⁺ barely touches the proxy (T···S −0.3 to −2.1), so the larger
R···(TS) than (RT)···S reflects the geometry built here (silanol on the second phosphate O), not a
general preference between the two adsorption pictures. One placement per base; the base spread
(≈3 kcal/mol in (RT)···S, G most favourable) is not established.

## Ion-pair reagent series (application), 2026-10-21

Same protocol as stages 7.1/7.2 (M06-2X/def2-SVP, C-PCM water, GPU optimization then GKS-EDA(sol)); seeds
from `scripts/build_hplc_reagents.py` (placement ranking after the 2026-10-19 fix; C8 pairs include
chain-contact seeds). Archive: `reports/hplc_reagents_2026-10-21/` (TEAH⁺ residue_A value from stage 7.1;
TEAH⁺···C8 from the two chain-contact seeds). All 25 cases converged with topology preserved and EDA
closure ≤ 3e-10 Eh. c8_hexylammonium_13_0: the xTB refinement was rejected (anchor distortion 0.020 Å
against the 0.02 Å gate) and DFT started from the RDKit seed. kcal/mol; "selected" = lowest PCM free
energy; range over all placements. Buried nonpolar SASA from `scripts/buried_sasa.py` (Å², C and C/Si-bound H;
FreeSASA cross-check on one complex: −0.6 %).

| Reagent | residue_A··· total (selected) | range | ···C8 total (selected) | range | ···C8 Pauli | ···C8 gcorr | buried nonpolar SASA, selected | all C8 placements |
| --- | ---: | --- | ---: | --- | ---: | ---: | ---: | --- |
| TEAH⁺ (TEAA) | -18.77 | -18.8 to -18.2 | -10.55 | -10.7 to -10.6 | 23.6 | -11.4 | 125 | 125–139 |
| hexylammonium (HAA) | -17.74 | -17.7 to -17.4 | -12.52 | -12.5 to -8.8 | 24.1 | -12.3 | 128 | 49–138 |
| DIPEA-H⁺ | -19.08 | -19.1 to -18.8 | -12.08 | -12.1 to -10.7 | 30.3 | -16.7 | 170 | 134–170 |
| dibutylammonium (DBAA) | -23.42 | -23.4 to -19.1 | -10.82 | -14.2 to -8.9 | 22.7 | -13.1 | 181 | 88–183 |
| tributylammonium (TBuAA) | -23.52 | -24.3 to -23.5 | -12.67 | -12.7 to -8.0 | 32.8 | -18.7 | 223 | 62–223 |

Reagent···C8 interaction totals are nearly flat across reagents (−10.6 to −12.7) and the placement range
within one reagent (up to 5.3 kcal/mol) exceeds the spread between reagents. This is expected: C-PCM has
no cavitation term, so the hydrophobic effect that dominates alkylammonium adsorption on C18 in water is
absent from GKS-EDA(sol). The buried nonpolar area tracks reagent hydrophobicity (tributylammonium largest)
and is reported as a separate empirical descriptor, without a surface-tension coefficient (literature
values span ~0.005–0.03 kcal/mol/Å²). A single C8 ligand cannot wrap a bulky reagent as a C18 brush can,
so buried areas of large reagents are likely underestimated. One residue (A); single local placements.

## Bridged ternaries (application), 2026-10-22

residue_A···reagent⁺···C8-proxy with the C8 ligand docked on the reagent alkyl chains of each optimized residue_A
ion pair (`scripts/build_hplc_bridged.py`, most proxy-C···reagent-C contacts within 4.5 Å) and no direct
residue–silane contact: phosphate on the reagent N–H, alkyl chains on C8. GFN2-xTB/ALPB pre-refinement, then
GPU M06-2X/def2-SVP, DF, C-PCM water, frozen silane anchor; three-body GKS-EDA(sol) as in stage 7.3.
Archive: `reports/hplc_bridged_2026-10-22/`. All five converged with topology preserved (frozen-anchor
displacement ≤ 2.7e-8 Å); closure ≤ 4.0e-10 Eh in all 25 decompositions. Constrained stationary points, one
seed per reagent. kcal/mol.

| Reagent | DFT evals | Wall (h) | N–H···O–P (Å) | R···T | T···S | R···S | (RT)···S | R···(TS) | ΔE₃ elec | ΔE₃ pol | ΔE₃ desolv | ΔE₃ total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TEAH⁺ | 126 | 6.04 | 1.54 | -19.22 | -2.15 | +0.02 | -2.11 | -19.30 | -0.12 | +0.32 | -0.33 | +0.02 |
| hexylammonium | 28 | 1.52 | 1.55 | -18.27 | -3.75 | -0.09 | -3.83 | -18.40 | +0.01 | +0.27 | -0.29 | +0.01 |
| DIPEA-H⁺ | 59 | 3.71 | 1.59 | -19.99 | -2.59 | -1.04 | -3.92 | -21.47 | -0.05 | +0.53 | -0.38 | -0.30 |
| dibutylammonium | 33 | 2.01 | 1.57 | -23.44 | -2.37 | +0.02 | -2.35 | -23.47 | -0.04 | +0.25 | -0.24 | +0.01 |
| tributylammonium | 67 | 5.21 | 1.64 | -24.70 | -3.17 | +0.04 | -3.11 | -24.74 | +0.27 | +0.36 | -0.76 | +0.03 |

In the bridged arrangement the three-body term is negligible (|ΔE₃| ≤ 0.3, against +1.0 to +1.2 for the
silanol-bound TEAH⁺ ternaries), so the complex is pairwise additive: the residue binds the reagent
(R···T −18.3 to −24.7, ordered as in the ion-pair series) and the reagent binds C8 only weakly through its
chains (T···S −2.1 to −3.7, mostly exchange + generalized correlation), with no direct residue–C8 term.

Reagent···C8 is 7–10 kcal/mol weaker here than in the isolated reagent···C8 pairs of the reagent series
(−10.6 to −12.7). Re-checking those 17 optimized pairs shows that every one, including the chain-contact
seeds, ended with the reagent N–H hydrogen-bonded to a silanol O (N–H···O 1.75–2.04 Å); in the bridged
complexes this N–H is used by the phosphate (N–H···O(silanol) 7.2–10.3 Å) while C···C contacts within 4.5 Å are
as many or more (10–20 vs 0–16). The isolated-pair values therefore measure mainly the N–H···silanol
hydrogen bond, not hydrophobic retention; the chain-only contribution is ≈ −2 to −4 kcal/mol in C-PCM, which
still lacks a cavitation term (see the buried-SASA discussion above). For the bridged model, the
reagent-dependent part of retention at this level is R···T (DBA/TBA ≈ 4–6 kcal/mol stronger than TEAH⁺/HA)
plus the hydrophobic term C-PCM omits. DIPEA-H⁺ is the only case with a small direct residue–C8 contact
(R···S −1.0, closest heavy atoms 3.2 Å). One residue, one seed per reagent; differences of ~1–2 kcal/mol are
not established.

## SMD CDS term (application), 2026-10-26

C-PCM GKS-EDA(sol) has electrostatic solvation only. To test whether the missing non-electrostatic term
explains the reagent and base orders, the SMD cavitation–dispersion–solvent-structure term (Marenich, Cramer &
Truhlar 2009, via PySCF/MNSOL; `eda_gpu.cds`, `scripts/nonpolar_hplc.py`) was evaluated on every converged
complex: dG_CDS = G_CDS(AB) − ΣG_CDS(fragments) at the complex geometry. It needs no SCF, so no new DFT was run.
SMD CDS was fitted together with SMD electrostatics and SMD radii, so EDA + dG_CDS is a mixed-model sensitivity
check, not a consistent free energy, and it is not part of the EDA closure. Archive:
`reports/hplc_nonpolar_2026-10-26/`. Lowest-energy placement per group; kcal/mol; np = buried nonpolar SASA (Å²).

| Complex | EDA total | dG_CDS water | EDA + CDS | dG_CDS acetonitrile | np | CDS per np (cal/mol/Å²) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| residue_A···TEAH⁺ | -18.77 | -2.36 | -21.13 | +0.20 | 109 | -22 |
| residue_C···TEAH⁺ | -19.74 | -1.85 | -21.59 | +0.61 | 113 | -16 |
| residue_G···TEAH⁺ | -19.64 | -2.44 | -22.07 | +0.30 | 121 | -20 |
| residue_U···TEAH⁺ | -19.98 | -2.51 | -22.49 | +0.46 | 153 | -16 |
| residue_A···C8 | -14.92 | -2.17 | -17.09 | +1.20 | 129 | -17 |
| residue_C···C8 | -18.87 | -2.63 | -21.49 | +0.83 | 89 | -30 |
| residue_G···C8 | -15.43 | -2.74 | -18.17 | +1.77 | 149 | -18 |
| residue_U···C8 | -14.86 | -2.43 | -17.29 | +0.84 | 106 | -23 |
| bridged T···S, TEAH⁺ | -2.15 | -0.72 | -2.87 | +0.81 | 171 | -4 |
| bridged T···S, DIPEA-H⁺ | -2.59 | -0.92 | -3.50 | +0.76 | 168 | -5 |
| bridged T···S, DBAA | -2.37 | -0.85 | -3.22 | +0.70 | 144 | -6 |
| bridged T···S, HAA | -3.75 | -1.03 | -4.77 | +1.32 | 198 | -5 |
| bridged T···S, TBuAA | -3.17 | -0.98 | -4.15 | +1.09 | 180 | -5 |
| bridged R···T, TEAH⁺ / DIPEA / DBAA / HAA / TBuAA | -19.2 / -20.0 / -23.4 / -18.3 / -24.7 | -2.3 / -2.2 / -2.7 / -1.9 / -3.0 | | +0.6 / +0.5 / +1.9 / +0.3 / +2.2 | | |

Findings:

- dG_CDS(water) is small (0 to −3.0 kcal/mol per pair; −3.5 to −5.2 for R···(TS) in the silanol-bound
  ternaries) and nearly additive (three-body CDS ≤ 0.1).
- It is weak exactly where the hydrophobic effect should act. Alkyl-chain contact (bridged T···S) gains only
  4–6 cal/mol per Å² of buried nonpolar area, against ≈ 25 cal/mol/Å² usually derived from alkane transfer free energies;
  polar contacts (ion pairs, residue···silanol) gain 16–30. SMD CDS in water therefore does not supply the
  missing hydrophobic driving force for reagent adsorption; its reagent spread on C8 is 0.3 kcal/mol.
- Reagent order on C8 (bridged T···S): TEAH⁺ < DBAA < DIPEA-H⁺ < TBuAA < HAA for EDA, CDS and their sum alike,
  Spearman ρ = 0.8 (p = 0.10, n = 5) against the published TEA < DIPEA < DBA < HA < TBA. CDS does not change the
  rank; the chain-contact EDA term alone already gives it. Buried nonpolar SASA alone gives ρ = 0.5. A single C8
  ligand cannot wrap the reagent like a C18 brush, so areas and contact energies of the bulky reagents are
  likely underestimated.
- Base order is not changed: ion pair + CDS U < G < C < A (−22.5 to −21.1, spread 1.4 kcal/mol), and
  residue···C8 + CDS still puts C strongest; measured T > A > G > C is not reproduced.
- In acetonitrile dG_CDS is positive for every pair (+0.1 to +3.0): the organic modifier opposes association
  through the non-electrostatic term, qualitatively as in gradient elution.

Conclusion: adding SMD CDS is not enough; the next checks are a silanol-free surface proxy, a larger surface
(several chains) that can bury the reagent, and placement sampling.

## Silanol-free C8 patch (application), 2026-10-27

Reagent⁺···patch pairs on a silanol-free three-chain C8 patch (`examples/hplc/c8x3_patch`: three tethered all-anti
n-octanes on a 7.6 Å triangle, C1/C2 of each chain frozen, no polar site; `scripts/build_hplc_patch.py`), two
placements per reagent. GFN2-xTB/ALPB pre-opt, then GPU M06-2X/def2-SVP, DF, C-PCM water; GKS-EDA(sol) and SMD
CDS/buried SASA as above. Archive: `reports/hplc_patch_2026-10-27/`. All ten converged with topology preserved
(frozen-anchor displacement ≤ 2.7e-8 Å, closure ≤ 7.3e-10 Eh). Constrained stationary points. kcal/mol;
ΔG_rel = relative optimized C-PCM energy of the two placements; N h / r = reagent N height above the C1 anchor
plane / lateral distance from the triangle centre (Å; chain tips at h ≈ 8.9, triangle inradius 2.2);
CC = reagent-C···patch-C contacts within 4.5 Å; Pauli = exchange + repulsion; np = buried nonpolar SASA (Å²).

| Reagent, placement | ΔE_rel | DFT evals | Wall (h) | N h / r | CC | elec | Pauli | pol | gen. corr | desolv | EDA total | dG_CDS water | EDA + CDS | dG_CDS MeCN | np |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TEAH⁺ 0 | 1.13 | 19 | 0.44 | 5.1 / 5.8 | 19 | -3.81 | 12.43 | -1.64 | -10.55 | +0.45 | -3.12 | -1.53 | -4.66 | +1.14 | 255 |
| TEAH⁺ 1 | 0.00 | 29 | 0.64 | 7.9 / 4.9 | 20 | -4.71 | 12.45 | -2.65 | -11.56 | +1.53 | -4.94 | -2.01 | -6.95 | +0.87 | 272 |
| DIPEA-H⁺ 0 | 0.00 | 21 | 0.61 | 5.5 / 5.0 | 22 | -4.25 | 12.66 | -1.96 | -13.12 | +1.18 | -5.49 | -1.57 | -7.06 | +1.42 | 323 |
| DIPEA-H⁺ 1 | 0.94 | 31 | 0.85 | 4.3 / 5.9 | 23 | -5.50 | 16.58 | -1.71 | -15.09 | +0.99 | -4.72 | -1.82 | -6.54 | +1.45 | 310 |
| DBAA 0 | 0.00 | 85 | 2.49 | 7.0 / 1.6 | 37 | -6.38 | 21.46 | -3.02 | -21.77 | +0.86 | -8.85 | -3.09 | -11.94 | +2.42 | 404 |
| DBAA 1 | 1.50 | 53 | 1.52 | 7.5 / 2.2 | 30 | -5.26 | 16.88 | -2.47 | -17.45 | +1.06 | -7.25 | -2.22 | -9.47 | +2.02 | 374 |
| HAA 0 | 0.00 | 74 | 1.58 | 10.0 / 3.7 | 27 | -5.83 | 16.51 | -1.23 | -16.24 | +0.82 | -5.97 | -2.03 | -8.00 | +1.69 | 305 |
| HAA 1 | 0.53 | 37 | 0.83 | 8.9 / 6.1 | 26 | -4.94 | 15.66 | -0.96 | -15.17 | +0.13 | -5.27 | -1.91 | -7.18 | +1.79 | 284 |
| TBuAA 0 | 0.85 | 69 | 2.17 | 9.3 / 6.2 | 19 | -4.57 | 13.83 | -1.18 | -12.95 | +0.77 | -4.11 | -1.85 | -5.96 | +1.27 | 309 |
| TBuAA 1 | 0.00 | 20 | 0.73 | 8.8 / 4.8 | 24 | -5.35 | 16.33 | -1.26 | -15.58 | +0.65 | -5.21 | -2.00 | -7.22 | +1.67 | 324 |

Findings:

- No reagent stayed polar-bound: with no silanol, every reagent binds through chain contact, and the energy is
  dispersion-like (generalized correlation −10.6 to −21.8, electrostatics −3.8 to −6.4, polarization ≤ 3.0).
  The patch roughly doubles the single-C8 chain contact (bridged T···S −2.2 to −3.8) to −4.9 to −8.9 (lower-energy
  placement), still below the silanol-H-bonded reagent···C8 pairs (−10.6 to −12.7).
- Reagent order (lower-energy placement; best-EDA placement gives the same): EDA TEAH⁺ < TBuAA < DIPEA-H⁺ < HAA <
  DBAA, Spearman ρ = 0.3 (p = 0.62, n = 5) against the published TEA < DIPEA < DBA < HA < TBA; EDA + CDS(water)
  TEAH⁺ < DIPEA-H⁺ < TBuAA < HAA < DBAA, ρ = 0.6 (p = 0.28); buried nonpolar SASA ρ = 0.5. TEAH⁺ is weakest in
  every measure, as measured; TBuAA is the outlier (most retained experimentally, 4th here).
- DBAA is strongest because both placements inserted its linear chains into the triangle (N at r = 1.6–2.2 Å,
  h = 7.0–7.5 Å, ≈ 1.5–2 Å below the chain tips; 30–37 C···C contacts). The branched/tertiary reagents (TBuAA,
  DIPEA-H⁺, TEAH⁺) stayed outside (r = 4.8–6.2 Å). With three chains the patch can be entered by a linear
  dialkylammonium but not wrapped around a bulky tertiary one, which biases against TBuAA; whether a denser or
  larger brush reverses this is not tested.
- Placement matters as much as reagent identity: the two placements differ by up to 1.8 kcal/mol in EDA (TEAH⁺,
  DBAA), comparable to the spread between reagents (excluding DBAA, −4.9 to −6.0). Two placements are not a
  converged sampling.
- SMD CDS(water) adds −1.5 to −3.1 (5–8 cal/mol/Å² of buried nonpolar area, again well below the ≈ 25 from alkane
  transfer); in acetonitrile it is positive (+0.9 to +2.4) for every pair.

Conclusion: removing the silanol gives the qualitatively right binding mode (chain contact, TEA weakest), but the
reagent order is not reproduced (ρ = 0.3–0.6, n = 5); the three-chain patch is too small for the bulky reagents
and placement sampling is too sparse to resolve 1 kcal/mol differences.

## Silanol-free C8 patch — placement sampling, 2026-10-28

Three extra placements per reagent on the same patch, chosen by a GFN2-xTB/ALPB screen (lowest dE_xTB of eight new
seeds, ranks 2–9 of `scripts/build_hplc_patch.py --keep 10`), same protocol as above; combined with placements 0/1
(5 per reagent). Archive: `reports/hplc_patch_sampling_2026-10-28/`. All 15 new cases converged with topology
preserved (no broken or new bonds; frozen-anchor displacement ≤ 1.6e-8 Å, closure ≤ 1.2e-9 Eh, 17–96 DFT
evaluations, 0.6–2.3 h each, 15.7 h total); no reagent heavy atom went below 2.7 Å above the anchor plane.
Constrained stationary points. Columns as in the previous section; ΔE_rel is now relative to the lowest of all five
placements; * = new placement.

| Reagent, placement | ΔE_rel | N h / r | CC | elec | Pauli | pol | gen. corr | desolv | EDA total | dG_CDS water | EDA + CDS | np |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| TEAH⁺ 0 | 2.07 | 5.1 / 5.8 | 19 | -3.81 | 12.43 | -1.64 | -10.55 | +0.45 | -3.12 | -1.53 | -4.66 | 255 |
| TEAH⁺ 1 | 0.94 | 7.9 / 4.9 | 20 | -4.71 | 12.45 | -2.65 | -11.56 | +1.53 | -4.94 | -2.01 | -6.95 | 272 |
| TEAH⁺ 2* | 0.00 | 4.7 / 5.4 | 25 | -3.73 | 14.34 | -1.84 | -14.16 | +0.07 | -5.33 | -2.01 | -7.33 | 265 |
| TEAH⁺ 3* | 1.25 | 4.7 / 5.8 | 17 | -4.45 | 12.90 | -1.85 | -11.53 | +1.08 | -3.86 | -1.49 | -5.35 | 274 |
| TEAH⁺ 4* | 1.97 | 10.2 / 2.2 | 15 | -2.91 | 8.83 | -1.11 | -8.09 | +0.52 | -2.76 | -1.91 | -4.67 | 260 |
| DIPEA-H⁺ 0 | 0.00 | 5.5 / 5.0 | 22 | -4.25 | 12.66 | -1.96 | -13.12 | +1.18 | -5.49 | -1.57 | -7.06 | 323 |
| DIPEA-H⁺ 1 | 0.94 | 4.3 / 5.9 | 23 | -5.50 | 16.58 | -1.71 | -15.09 | +0.99 | -4.72 | -1.82 | -6.54 | 310 |
| DIPEA-H⁺ 5* | 0.92 | 4.6 / 5.9 | 21 | -4.95 | 15.60 | -1.58 | -14.36 | +0.69 | -4.60 | -1.77 | -6.38 | 317 |
| DIPEA-H⁺ 6* | 2.99 | 11.1 / 1.2 | 14 | -3.62 | 10.10 | -1.10 | -9.45 | +0.82 | -3.25 | -1.81 | -5.06 | 303 |
| DIPEA-H⁺ 8* | 2.10 | 10.6 / 1.5 | 16 | -3.70 | 10.41 | -1.17 | -9.90 | +0.87 | -3.49 | -1.67 | -5.16 | 284 |
| DBAA 0 | 0.00 | 7.0 / 1.6 | 37 | -6.38 | 21.46 | -3.02 | -21.77 | +0.86 | -8.85 | -3.09 | -11.94 | 404 |
| DBAA 1 | 1.50 | 7.5 / 2.2 | 30 | -5.26 | 16.88 | -2.47 | -17.45 | +1.06 | -7.25 | -2.22 | -9.47 | 374 |
| DBAA 2* | 2.87 | 7.6 / 5.3 | 21 | -3.46 | 12.55 | -1.70 | -12.81 | +0.09 | -5.34 | -1.92 | -7.26 | 275 |
| DBAA 3* | 3.53 | 9.6 / 3.8 | 19 | -4.07 | 11.49 | -1.41 | -11.57 | +0.84 | -4.73 | -1.85 | -6.58 | 266 |
| DBAA 8* | 4.24 | 8.6 / 4.7 | 14 | -3.67 | 10.54 | -2.29 | -10.43 | +1.15 | -4.71 | -1.99 | -6.69 | 247 |
| HAA 0 | 0.00 | 10.0 / 3.7 | 27 | -5.83 | 16.51 | -1.23 | -16.24 | +0.82 | -5.97 | -2.03 | -8.00 | 305 |
| HAA 1 | 0.53 | 8.9 / 6.1 | 26 | -4.94 | 15.66 | -0.96 | -15.17 | +0.13 | -5.27 | -1.91 | -7.18 | 284 |
| HAA 3* | 1.59 | 12.2 / 4.6 | 20 | -4.97 | 14.03 | -0.58 | -13.26 | +0.49 | -4.29 | -1.52 | -5.82 | 268 |
| HAA 4* | 1.86 | 11.2 / 5.9 | 24 | -5.08 | 14.92 | -0.71 | -13.62 | +0.42 | -4.07 | -1.52 | -5.59 | 260 |
| HAA 5* | 1.68 | 12.1 / 3.2 | 17 | -4.74 | 12.75 | -0.74 | -12.43 | +0.91 | -4.26 | -1.59 | -5.84 | 307 |
| TBuAA 0 | 3.52 | 9.3 / 6.2 | 19 | -4.57 | 13.83 | -1.18 | -12.95 | +0.77 | -4.11 | -1.85 | -5.96 | 309 |
| TBuAA 1 | 2.67 | 8.8 / 4.8 | 24 | -5.35 | 16.33 | -1.26 | -15.58 | +0.65 | -5.21 | -2.00 | -7.22 | 324 |
| TBuAA 3* | 0.00 | 7.7 / 3.0 | 31 | -6.49 | 20.37 | -1.92 | -19.99 | +0.90 | -7.14 | -2.55 | -9.68 | 434 |
| TBuAA 5* | 1.66 | 8.4 / 2.5 | 32 | -5.99 | 19.01 | -2.25 | -18.47 | +0.98 | -6.72 | -2.80 | -9.51 | 428 |
| TBuAA 9* | 3.67 | 10.2 / 4.0 | 21 | -4.40 | 12.12 | -1.60 | -11.86 | +1.05 | -4.69 | -2.11 | -6.80 | 312 |

Reagent order (least → most retained published: TEA < DIPEA < DBA < HA < TBA), n = 5:

| Measure | Order (weakest → strongest) | Spearman ρ (p) | 10-27 (2 placements) |
| --- | --- | ---: | ---: |
| EDA, lowest-energy = best placement = Boltzmann | TEAH⁺ < DIPEA-H⁺ < HAA < TBuAA < DBAA | 0.7 (0.19) | 0.3 |
| EDA + CDS(water), same three | DIPEA-H⁺ < TEAH⁺ < HAA < TBuAA < DBAA | 0.6 (0.28) | 0.6 |
| EDA + CDS(water), mean of 5 | TEAH⁺ < DIPEA-H⁺ < HAA < TBuAA < DBAA | 0.7 (0.19) | — |
| buried nonpolar SASA, lowest-energy | TEAH⁺ < HAA < DIPEA-H⁺ < DBAA < TBuAA | 0.7 (0.19) | 0.5 |

Findings:

- In every reagent the lowest-energy optimized placement is also the most strongly bound one in EDA and EDA + CDS,
  and it carries 63–93 % of the 298 K Boltzmann weight (electronic C-PCM energy only), so lowest-energy, best and
  Boltzmann-weighted give the same order. Lowest-energy placement changed for TEAH⁺ (→ 2, −5.33) and TBuAA
  (→ 3, −7.14); DIPEA-H⁺, DBAA and HAA kept placement 0.
- TBuAA moves from 4th to 2nd: placements 3 and 5 insert it over the triangle (N r = 2.5–3.0 Å, h = 7.7–8.4 Å;
  31–32 C···C contacts, the largest buried area of any pair, 428–434 Å²), the same binding mode as DBAA 0/1. The
  10-27 statement that the patch cannot bind a bulky tertiary reagent this way was a sampling artefact.
- Insertion over the triangle (r ≤ 3 Å, h ≤ 8.5 Å, ≥ 30 contacts) gives the four strongest pairs (−6.7 to −8.9);
  every other placement, including those sitting on the chain tips (r ≤ 2.2 Å but h ≥ 10 Å: TEAH⁺ 4, DIPEA-H⁺ 6/8),
  is −2.8 to −6.0. TEAH⁺, DIPEA-H⁺ and HAA found no inserted placement among five.
- The remaining disagreement is DBAA > TBuAA (−8.85 vs −7.14) and HAA below both. With TBuAA and DBAA swapped,
  EDA would give ρ = 0.9; the gap (1.7 kcal/mol) is about one placement spread, so it is not resolved.
- xTB screening is a weak filter: pooled over the 15 new placements dE_xTB vs EDA has ρ = 0.78, mostly from reagent
  identity (TBuAA strongest in both); within reagents the xTB best pick was also the DFT best in 2 of 5 (DIPEA-H⁺,
  HAA). xTB overbinds by 1.1–3.5 kcal/mol.

Conclusion: five placements per reagent make the ordering stable within each reagent (lowest-energy = best), raise
the EDA rank correlation to ρ = 0.7 (n = 5, not significant) and remove the TBuAA artefact, but DBAA is still
strongest. Whether DBA vs TBA and HA follow the published order depends on the brush: a larger/denser patch
(5–7 chains) that both linear and tertiary reagents can insert into is the next test.
