# Silanol-free three-chain C8 patch — reagent···patch pairs (2026-10-27)

Five ion-pair reagents × 2 placements on `examples/hplc/c8x3_patch` (three tethered all-anti n-octanes, 7.6 Å
triangle, C1/C2 of each chain frozen, no polar site). GFN2-xTB/ALPB pre-opt → M06-2X/def2-SVP DF C-PCM water
optimization → GKS-EDA(sol); SMD CDS and buried SASA at the optimized geometry. Table and discussion:
`docs/VALIDATION.md`, section "Silanol-free C8 patch (application), 2026-10-27".

- `summary.json`, `analysis.json`, `nonpolar.json`: copies from `runs/hplc_patch_cpcm_def2svp/`.
- `<case>/optimized.xyz`, `<case>/optimization_result.json`, `<case>/eda_result.json`: final geometry, optimizer
  result, full GKS-EDA(sol) result (fragment 0 = patch, fragment 1 = reagent). Reference densities, xTB files and
  trajectories are not archived (in `runs/`).
- `queue.log`: per-case progress and EDA lines from `runs/hplc_patch.log` (warnings removed).

Command: `runs/queue_hplc_patch.sh` (run as `systemd-run --user --unit=hplc-patch`, 2026-10-26 18:12 → 2026-10-27 06:58).
