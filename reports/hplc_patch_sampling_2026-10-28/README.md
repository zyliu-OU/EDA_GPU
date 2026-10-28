# Silanol-free C8 patch — placement sampling (2026-10-28)

Three extra placements per reagent on `examples/hplc/c8x3_patch`, added to the two of
`reports/hplc_patch_2026-10-27/` (5 per reagent, 25 in all). Seeds: `scripts/build_hplc_patch.py --keep 10 --output
examples/hplc/patch_sampling` (ranks 2–9); GFN2-xTB/ALPB screen of the 40 new seeds (`scripts/screen_hplc_patch.py`),
lowest three dE_xTB per reagent kept. Then the same protocol as 10-27: GFN2-xTB pre-opt → M06-2X/def2-SVP DF C-PCM
water optimization → GKS-EDA(sol); SMD CDS and buried SASA at the optimized geometry. Table and discussion:
`docs/VALIDATION.md`, section "Silanol-free C8 patch — placement sampling, 2026-10-28".

- `summary.json`, `analysis.json`, `nonpolar.json`: copies from `runs/hplc_patch_sampling_cpcm_def2svp/` (the 15
  new cases only).
- `screen/screen.json`, `screen/selection.json`: xTB screen of placements 2–9 and the picks.
- `<case>/optimized.xyz`, `<case>/optimization_result.json`, `<case>/eda_result.json`: as in the 10-27 archive
  (fragment 0 = patch, fragment 1 = reagent).
- `cases.json`: per-case checks and values for all 25 placements (convergence, broken/new bonds, frozen-anchor
  displacement, EDA closure, N height `h` / lateral distance `r` from the C1 triangle, lowest reagent heavy atom
  above the C1 plane, C···C contacts within 4.5 Å, EDA components, CDS, buried SASA, ΔE relative within reagent).
  Produced by `check_patch_cases.py` (run from the repo root on both `runs/` directories).
- `ranking.json`: per-reagent lowest-energy / best / Boltzmann-weighted (298 K, on the optimized C-PCM electronic
  energy only; no ZPE/entropy) / mean values, Spearman ρ against the published order, xTB-vs-DFT comparison.
  Produced by `rank_patch_reagents.py cases.json ranking.json`.
- `queue.log`: progress and EDA lines from `runs/hplc_patch_sampling.log` (warnings removed).

Command: `runs/queue_hplc_patch_sampling.sh` (as `systemd-run --user --unit=hplc-patch-sampling`,
2026-10-27 09:21 → 2026-10-28 02:23).
