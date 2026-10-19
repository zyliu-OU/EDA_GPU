# EDA_GPU session handoff — 2026-10-18

Repository: https://github.com/zyliu-OU/EDA_GPU (public). Workspace: /home/zyliu/Documents/EDA_GPU.
Prior cumulative history: `reports/EDA_GPU_HANDOFF_2026-10-11.md` (build, validation, EDA stages).
The original design brief `EDA_GPU_CODEX_HANDOFF_2026-10-11.md` stays local and untracked by design.

## State at stop

- `main` = `origin/main`; code/results commit `a6255b3` (Psi4 cross-check). This handoff is a
  separate documentation commit; check `git log -1` for its SHA.
- CPU CI passed for `a6255b3`: https://github.com/zyliu-OU/EDA_GPU/actions/runs/39181092250.
- Local suite: 37 tests passed (`.env/bin/python -m pytest -q`), GPU exercised.
- No calculations are running. Nothing uncommitted except the local-only original brief.

## Done this session

1. Reviewed the original brief against the repository; the brief predates the build and is stale.
2. Answered "is it verified by another data source": before today only internal checks
   (CPU/GPU, finite differences, independent algebra, closure) plus two unmatched published
   comparisons (XEDA II CO₂ manual values; Liu et al. PCCP 2017 tryptamine–water Table 3).
3. Added the first independent-program validation of the energy/gradient layer:
   - `scripts/compare_psi4.py` (run with `.env`) and `scripts/psi4_reference.py` (run by
     `.psi4-env/bin/python -I`). Psi4 1.11, Libxc 7.1.2, libint 2.13.1; env lock
     `requirements-psi4-conda-lock.txt`; recreate with
     `conda create -p .psi4-env -c conda-forge --override-channels python=3.12 psi4`.
   - Matched exact integrals, identical Libxc IDs (B3LYP 402, PBE0 406, M06-2X 450+236),
     unpruned grids (PySCF level 8; Psi4 150×1202), bohr coordinates via PySCF's constant.
   - 13/13 pass: HF 3.4e-13 Eh; GGA hybrids ≤1.6e-8 Eh, ≤3.9e-7 Eh/bohr;
     M06-2X/aug-cc-pVDZ tryptamine–water +7.9e-7 Eh, 8.2e-6 Eh/bohr (CPU and GPU).
   - Archive: `reports/psi4_crosscheck_2026-10-18/`; write-up in `docs/VALIDATION.md`.
   - First attempt (`runs/psi4_crosscheck_full`, not archived): GPU M06-2X SCF did not reach
     1e-12/1e-7 and failed explicitly; large-case SCF thresholds now 1e-10/1e-6.
4. README: Psi4 status row and citation (references renumbered); rebased over remote commit
   `c5541b2` (identical README bibliography) without force-push.

## Interpretation and limits

Energy/gradient layer agrees with an independent program on CPU and GPU. EDA components are
NOT validated (Psi4 has no GKS-EDA). CO₂ total energy agrees with Psi4 to 1.6e-8 Eh, so the
1–5e-6 Eh component gaps to XEDA II are from EDA settings/definitions, not the SCF energy.
Not compared: UKS/open-shell gradients, optimized geometries, density-fitted paths. Timings in
the archive are single cold runs, not benchmarks.

## Next steps (priority)

1. Matched XEDA runs (GAMESS+XEDA) for water dimer and CO₂ dimer with recorded basis
   representation, B3LYP variant, grid, CP policy; compare each component to ~1e-6 Eh.
   Offer to prepare inputs once the user confirms XEDA access.
2. Ask whether original tryptamine–water coordinates/outputs from the 2017 work exist.
3. Optional: extend the Psi4 check to UKS (water radical-cation dimer) and DF paths.
4. Pending from earlier: CPU common-geometry tryptamine comparison and cation benchmarks.

## Session 2 update (same day) — XEDA kit, Table 1 geometry, Psi4 geometry; 3a–3c running

Committed and pushed `626622f`:
- Item 1 (matched XEDA): `scripts/xeda_matched.py` + `reports/xeda_matched_2026-10-18/`
  (5 GAMESS-XEDA inputs from official manual keywords, EDA_GPU values, README with run/check
  steps). **Blocked on user:** no GAMESS/XEDA on this host. User must run the five inputs and
  supply logs; then `scripts/xeda_matched.py compare ...`. Check B3LYP VWN form in GAMESS logs.
- Item 2: no original 2017 coordinates/outputs on host (only manuscript + SI PDFs in
  ~/Downloads). Table 1 H-bond geometry matches within 0.05 Å / 1°
  (`scripts/tryptamine_hbond_geometry.py`).
- Item 3d: Psi4 geometry cross-check passes (water, water dimer; Hessian-verified minima).
  Symmetric example dimer start → geomeTRIC Cs saddle at PBE0/def2-SVP (documented);
  archived B3LYP validation dimer verified to be a true minimum.

Running at time of writing (outputs under `runs/`, not yet archived):
- 3a `runs/tryptamine_neutral_cpu_common_eda` — CPU EDA at the GPU common geometry with the
  archived GPU EDA settings; compare components with
  `reports/tryptamine_neutral_gpu_2026-10-11/result.json`.
- 3b `runs/tryptamine_cation_gpu_df_opt` — GPU UKS M06-2X/aug-cc-pVDZ (DF aug-cc-pvdz-jkfit,
  grid 3) optimization from the neutral optimized geometry with one electron removed
  (paper procedure). Next: GPU direct EDA with the neutral EDA settings, then compare with
  Table 3 `A_cation_H2O` and Table 1 cation H-bonds (`--case A_cation_H2O`).
- 3c not started: `.env/bin/python scripts/compare_psi4.py --gpu --set extended --output
  runs/psi4_crosscheck_extended` (DF and UKS cases), after 3a/3b free the CPU.
