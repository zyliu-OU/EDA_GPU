# Silanol-free C8 patch — solvent-model sensitivity (2026-11-10)

GKS-EDA(sol) single points on the 25 patch placements of `reports/hplc_patch_2026-10-27/` and
`reports/hplc_patch_sampling_2026-10-28/`. The geometries are the M06-2X/def2-SVP C-PCM(water, Bondi×1.2) optimized
structures and are not re-optimized. Only one field of the recorded `eda_settings` changes in each variant:

- `vdw14`: `solvent_vdw_scale` 1.2 → 1.4 (larger cavity), water ε = 78.3553
- `acn`: `solvent_eps` 78.3553 → 35.688 (acetonitrile dielectric), Bondi × 1.2
- `ref`: unchanged, on patch_teah_2 only (reproduction check: −5.33 kcal/mol, same as the original)

Files:
- `summary.json`: per case and variant, the components (kcal/mol), complex energy, closure residual, time and
  solvent record. Copy of `runs/hplc_patch_sensitivity/summary.json`.
- `ranking.json`: per-variant lowest-energy / best / Boltzmann (298 K, from the variant's own complex energies) /
  mean values per reagent, Spearman ρ against the published order, and per-placement shifts relative to the base run.
  Produced by `python rank_sensitivity.py runs/hplc_patch_sensitivity/summary.json
  reports/hplc_patch_sampling_2026-10-28/cases.json ranking.json`. EDA + CDS adds the SMD CDS term of the matching
  solvent (water for base/vdw14, acetonitrile for acn), computed at the same geometry in the 10-26/10-28 runs.
- `queue.log`: the queue log with warnings removed. `queue_hplc_sensitivity.sh`: the queue itself, run as
  `systemd-run --user --unit=hplc-sensitivity` from 2026-11-09 22:09 to 2026-11-10 02:28.

Discussion: `docs/VALIDATION.md`, section "Silanol-free C8 patch — solvent-model sensitivity, 2026-11-10".
