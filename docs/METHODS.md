# Methods and conventions

The implemented model is restricted closed-shell Kohn–Sham DFT in gas phase.
Default: B3LYP/def2-SVP, density fitting with def2-universal-jkfit, unpruned level-3
quadrature, SCF energy tolerance 1e-10 Eh, at most 100 SCF cycles, four CPU threads,
and no explicit dispersion. CPU and GPU use identical specified settings. Open-shell,
range-separated, nonlocal-correlation, double-hybrid, solvent and EDA methods are
outside the initial scope. Alternate global hybrid/GGA/LDA functionals are exposed
but only the reported B3LYP case is validated. No guessed interaction components.

PySCF constructs the molecule and CPU RKS. GPU requests check CuPy device execution
and construct GPU4PySCF RKS directly. Integrals/J/K and DFT evaluations are handled by
the selected backend; geomeTRIC and structure bookkeeping run on CPU. Final gradients
are transferred to NumPy for optimizer control/reporting. This is not an all-GPU workflow.

Optimization wraps the backend's analytic-gradient scanner with the documented PySCF
custom energy-gradient interface so geomeTRIC receives host arrays and checks SCF
convergence at each evaluation. Optimization status uses `geometric_solver.kernel`'s
actual boolean, rather than assuming successful return means convergence. Thresholds:
energy change 1e-6 Eh, RMS/max gradient 3e-4/4.5e-4 Eh/bohr,
RMS/max displacement 1.2e-3/1.8e-3 angstrom. Final energy and gradient are recomputed.
The trajectory contains evaluated geometries including trial/rejected steps.

Frozen atoms use `$freeze` Cartesian `xyz` constraints, translated from metadata's
zero based indices to geomeTRIC's one based indices. Enforce threshold is 1e-6 in
geomeTRIC internal units; final displacement must be below 1e-5 angstrom. External
constraint files are also accepted but cannot be mixed with frozen metadata.
A constrained stationary point can carry nonzero gradients along frozen directions.
No vibrational Hessian has established any reported stationary point as a minimum.

Input validation checks finite coordinates, electron parity, fragment electron/spin
consistency, exact partition and charge sums. Fragment geometry/order remains inherited
from the complex; no fragment optimizations or BSSE correction are performed by this
optimizer. Bond/stereochemistry metadata are retained, not inferred or chemically repaired.

Official API sources checked 2026-10-11:
- https://pyscf.org/user/gpu.html
- https://pyscf.org/user/geomopt.html
- https://github.com/pyscf/gpu4pyscf
- https://geometric.readthedocs.io/en/latest/constraints.html
- Installed PySCF/geomeTRIC/GPU4PySCF source and package versions.

## Experimental EDA preparation

`eda-prepare` uses CPU direct integrals, B3LYP/PBE0, two restricted closed-shell
all-electron fragments, and Boys–Bernardi ghost-basis counterpoise at the unchanged
complex geometry. Every density uses the same AO basis with spin-summed occupancy 2.
Frozen electrostatics follows review Eq. 2; generalized correlation follows Eq. 9,
evaluated as KS energy minus the full HF-like functional on the same KS density.
Explicit dispersion is zero for these settings. The diagnostic bare cross HF exchange
is not assigned to the original EDA exchange component. Exchange, repulsion, polarization
and closure are unavailable. This preparation is experimental and does not satisfy the
full EDA reference gate. See GKS_EDA_SPEC.md for equations and normalization blockers.

## Complete experimental GKS-EDA

The new `eda` command computes all six gas-phase components on CPU or GPU, with
the same direct-integral/closed-shell/counterpoise policy. It adds the unnormalized
exchange energy expression and normalized occupied-span determinant from Su–Li
Eqs. 6–12, evaluated on KS orbitals following the supplied review. Repulsion is the
normalized-minus-unnormalized HF-like energy difference; polarization is the
relaxed-KS-orbital HF-like energy minus the normalized frozen determinant energy.
Generalized correlation contains the remaining hybrid exchange/correlation difference.
The smallest occupied-metric eigenvalue must exceed 1e-8; electron count, density
idempotency/symmetry, electrostatic/exchange consistency and closure are checked.

CPU uses PySCF; GPU uses GPU4PySCF direct RKS/J/K and CuPy for density contractions
and occupied-space orthonormalization. GPU4PySCF direct J/K has CPU/GPU hybrid
internals. CPU molecule/basis bookkeeping and output remain. No fallback labels CPU
as GPU. `--cart` explicitly selects Cartesian basis functions, while the default is
spherical. The named B3LYP5 variant is also available experimentally. No density fitting
or explicit dispersion correction is added. All terms are computed, but external
matched-settings component validation and a validated full-EDA release remain pending.
The mathematical specification in GKS_EDA_SPEC.md supersedes earlier EDA blockers.

The complete `eda` CLI defaults to SCF energy tolerance 1e-12 Eh, matching the
recorded component benchmarks; other calculation commands retain 1e-10 Eh.

The complete `eda` and molecular calculation commands additionally support UKS for nonzero spin and the global hybrid meta-GGA M06-2X. See GKS_EDA_SPEC.md for spin-resolved formulas and the additive-spin restriction. This extension is experimental.
