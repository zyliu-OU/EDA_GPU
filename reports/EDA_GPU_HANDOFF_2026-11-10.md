# EDA_GPU handoff — 2026-11-10

Previous: `reports/EDA_GPU_HANDOFF_2026-11-09.md` (plan: 3-chain patch only, steps 1–5).

## Done: step 1 (solvent sensitivity), committed b3a498d

- 50/50 EDAs OK (closure ≤ 1.2e-9); ref reproduces −5.33. Archive `reports/hplc_patch_sensitivity_2026-11-10/`,
  `docs/VALIDATION.md` section, README line. Bondi ×1.4: all pairs −0.57 ± 0.32 more bound, ρ 0.5–0.7. ACN ε: EDA
  shifts < 0.03 kcal/mol (C-PCM (ε−1)/ε); ACN CDS halves EDA+CDS, no rank change. DBAA > TBuAA persists.
- Engine GPU-cache release fix committed (46/46 tests).

## Running: step 2 — aug-cc-pVDZ (systemd unit `hplc-augdz`, started 2026-11-10)

- Memory: the EDA uses **direct integrals** (density_fitting False), so no DF tensor (a DF one would be 28–46 GB, not
  feasible: matters only if aug-cc-pVDZ *optimizations* were ever wanted). Observed ~3.7 GB GPU in the first SCF.
- `scripts/sensitivity_hplc.py` variant `augdz`; `runs/queue_hplc_augdz.sh`: lowest-energy placement per reagent
  (teah_2, hexyl_0, dipeah_0, DBAA_0, TBuAA_3) at the def2-SVP geometries. Output `runs/hplc_patch_augdz/`
  (summary.json, mem.log every 10 s), log `runs/hplc_augdz.log`, ends "queue done". Time per case unknown (expect ~1 h+).
- Next: compare with def2-SVP per reagent (tot and components; ρ), archive `reports/hplc_patch_augdz_<date>/`, docs.

## Queued: step 3 — BA⁺, DEA⁺, TPA⁺ (systemd unit `hplc-patch-extended`, waits for `hplc-augdz` to finish)

- Monomers: `scripts/build_hplc_reagents.py` (+3 entries; also wrote residue_A ion-pair and c8 surface seeds under
  `examples/hplc/reagents/`). `scripts/hplc_retention.py` mapping extended (deaa/baa/tpaa).
- Seeds: `build_hplc_patch.py --keep 10 --output examples/hplc/patch_extended` (95/95/110 atoms; c8x3_patch unchanged).
- `runs/queue_hplc_patch_extended.sh`: placements 0/1 → xTB + DFT opt + EDA; xTB screen of 2–9 (`--known` = 0/1),
  best 3 → DFT; then nonpolar (CDS/SASA) and analysis. Output `runs/hplc_patch_extended_cpcm_def2svp/`, log
  `runs/hplc_patch_extended.log`. 15 DFT cases, ~0.6–2.3 h each → roughly 15–25 h.
- Published order (Donegan Table 4, 100 mM acetate, %ACN 20T): DEA 10.2 < TEA 11.9 < DIPEA 15.0 < BA 16.3 <
  TPA 27.2 < DBA 32.9 < HA 35.3 (< TBA). n = 8 once done; rank with the 10-28 `rank_patch_reagents.py` logic.
- Not yet committed: the reagent/retention script edits, new seeds, queue scripts (runs/ is ignored; archive copies).
- Still uncommitted on purpose: 7-chain builder diff, `c8x7_patch.*`, `reports/c8x7_memtest_2026-11-09/`.

## After 3: step 5 retention fit; then step 4 (bases on the patch).
Linger is off: closing the terminal is fine, logging out of the desktop session would stop both units.

## Status check 22:58
- augdz: 4/5 done (43–52 min each, closure ≤ 5e-9, peak ~7.8 GB GPU / ~4.6 GB RAM, GPU ≤ 61 °C). TBuAA_3 running since 21:49.
- Δ(aDZ − SVP) total: TEA −2.27, HA −2.33, DIPEA −1.86, DBAA −3.34; nearly all from generalized correlation (−3.4 to −6.2),
  polarization shrinks by 0.8–1.6 (HA pol turns +0.17). Order DBAA > HA > TEA ≈ DIPEA (TEA/DIPEA swap, 0.25 apart; 0.16 at SVP).
