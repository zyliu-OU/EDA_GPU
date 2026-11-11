# EDA_GPU handoff — 2026-11-09

Previous: `reports/EDA_GPU_HANDOFF_2026-10-28.md` (its 11-09 sections: 5/7-chain dropped, 3-chain patch only).

## Plan agreed with the user (3-chain patch only)

1. Solvent sensitivity on the 25 existing patch placements (single points, no re-optimization) — RUNNING.
2. aug-cc-pVDZ single points on the best placement per reagent; estimate memory from the basis size first
   (no CPU PCM/large-matrix builds).
3. Add BA, DEA, TPA to the reagent series (≈100–110 atoms with the patch; full xTB → DFT workflow), n 5 → 8.
4. Bases on the patch: residue···patch and residue···reagent pairs at DFT; full ternary (140–160 atoms) only as
   xTB geometry + DFT EDA single point.
5. Retention fit (`scripts/hplc_retention.py`) after 3, then after 4.
Housekeeping: commit the engine GPU-cache release fix; 7-chain builder / `c8x7_patch.*` not needed.

## Running: step 1

- New `scripts/sensitivity_hplc.py`: repeats GKS-EDA(sol) at `<run>/<case>/optimization/optimized.xyz` with the
  run's recorded `eda_settings`, changing only the variant fields: `vdw14` (solvent_vdw_scale 1.4), `acn`
  (solvent_eps 35.688), `ref` (unchanged, reproduction check). Records components, complex energy (for Boltzmann
  weights at the variant), closure. Resumable.
- `runs/queue_hplc_sensitivity.sh` as systemd user unit `hplc-sensitivity` (started 2026-11-09 22:09):
  `ref` on patch_teah_2 (must reproduce tot = −5.33), then all 25 × `vdw14`, then all 25 × `acn`.
  ~6 min per EDA → ~5–6 h. Output `runs/hplc_patch_sensitivity/` (summary.json, `<variant>/<case>/`), log
  `runs/hplc_sensitivity.log`, ends with "queue done".
- Survives closing the terminal window. Linger is off, so logging out of the desktop session would stop it.

## Next when it finishes

- Check ref reproduction, closures, failures. Re-rank with `reports/hplc_patch_sampling_2026-10-28/rank_patch_reagents.py`
  logic per variant (tot; Boltzmann weights from the variant complex energies); compare ρ and whether DBAA > TBuAA
  and HAA-low persist. Archive to `reports/hplc_patch_sensitivity_<date>/`, add a `docs/VALIDATION.md` section,
  update README, commit. Then step 2.
