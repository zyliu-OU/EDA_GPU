# EDA_GPU handoff — 2026-11-16

Previous: `reports/EDA_GPU_HANDOFF_2026-11-11.md`.

## Running: step 4 — bases A/C/G/U on the C8 patch (systemd user unit `hplc-patch-bases`, started 2026-11-16)
- Seeds: `scripts/build_hplc_patch.py --reagents residue_A residue_C residue_G residue_U --keep 10 --output
  examples/hplc/patch_bases` (120/118/121/117 atoms; `c8x3_patch.*` unchanged). Residues are the capped RNA
  monoanions (charge -1).
- `runs/queue_hplc_patch_bases.sh`: same protocol as the reagents — placements 0/1 → xTB pre-opt + DFT opt +
  GKS-EDA(sol); xTB screen of 2–9 (`--known` = 0/1), best 3 → DFT; then nonpolar (CDS/SASA) and analysis.
  Thermal guards at `run_hplc_stage.py` defaults (GPU 80 °C, CPU 85 °C). Resumable (rerun the same command).
- Output `runs/hplc_patch_bases_cpcm_def2svp/`, log `runs/hplc_patch_bases.log`, ends with "queue done".
  20 DFT cases; reagent cases of this size took 0.6–2.5 h each → roughly 1–2 days.
- Linger is off: logging out of the desktop session stops it. Restart:
  `systemd-run --user --unit=hplc-patch-bases --working-directory=$PWD /bin/bash -c 'exec /bin/bash runs/queue_hplc_patch_bases.sh >> runs/hplc_patch_bases.log 2>&1'`

## Done today — chain-splay check: `reports/hplc_patch_splay_2026-11-16/`
Geometry-only analysis of the 40 optimized reagent placements. The chains do **not** splay around inserted
reagents: C5–C8 chain–chain distance 6.61 Å inserted vs 6.81 Å other vs 7.23 Å in the constructed patch; free
radius on the axis 2.28 / 2.50 / 2.78 Å. The central site exists by construction (7.6 Å spacing), so option 1a
(restrain/freeze outer chain carbons) would not close it and is dropped.

## Discussed with the user today (no decision yet)
- "Hardware limit" = GPU memory (16 GB) at the DFT gradient for the 7-chain brush (223 atoms), with host RAM
  (32 GB) also tight. A 4-chain patch (~145 atoms) was never tested; estimated gradient buffers ~3–4 GB, so it
  probably fits. A memory test (single point + gradient) would be needed first.
- SSD as overflow: works for host RAM (swap; 8 GB swap file now, 1.2 TB free, needs sudo to enlarge), not directly
  for GPU memory. Untested routes for the 7-chain brush: cupy managed (unified) memory, or batching the gradient
  buffers in gpu4pyscf.
- Note: a 4th chain on a rhombus adds a second triangular hollow rather than filling the first one.

## Next when the queue finishes
- Check convergence, closures, broken/new bonds (adapt `check_patch_cases.py`; r/h there are defined by the first
  guest N — for residues use the base ring centroid, and report insertion in the hollow per placement).
- Boltzmann EDA per base, components, CDS; compare A/C/G/U order with the user's base-ratio retention data and
  with reagent values. Archive to `reports/hplc_patch_bases_<date>/`, VALIDATION.md section, README.
- Then: off-centre re-sampling for DEA/BA/DBA (option 2, ~10–15 h) if the reagent ranking is still wanted, or the
  4-chain memory test if the user wants a larger surface.

## Uncommitted
11-11 work (retention-fit archive, `scripts/fit_patch_retention.py`, Donegan CSV, VALIDATION/README edits, progress
report, 11-11 handoff), today's splay archive, `examples/hplc/patch_bases/`, this handoff. 7-chain files
(`c8x7_patch.*`, `reports/c8x7_memtest_2026-11-09/`) uncommitted on purpose. User has not said whether to commit/push.
