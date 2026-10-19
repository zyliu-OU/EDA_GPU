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

## Matched-settings GAMESS-XEDA kit, 2026-10-18 (external run pending)

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
components within 0.20 kcal/mol, with the same sign pattern (Pauli high by 0.15–0.20 kcal/mol,
attractive terms slightly more negative). A systematic pattern across charge states suggests a
geometry or convention difference rather than numerical noise; grid (≤0.0015 kcal/mol) and
backend (≤7e-6 kcal/mol) are excluded. These remain exploratory comparisons
(`scientifically_validated=false`) until original coordinates or matched XEDA output exist.

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
