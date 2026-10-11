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
