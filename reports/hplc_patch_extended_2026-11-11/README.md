# Silanol-free C8 patch — BA⁺, DEA⁺, TPA⁺ added (2026-11-11)

Same protocol as 10-27/10-28 on `examples/hplc/c8x3_patch` (seeds `examples/hplc/patch_extended`): placements 0/1
plus the best 3 of an xTB screen of 2–9 → GFN2-xTB pre-opt → M06-2X/def2-SVP C-PCM water opt → GKS-EDA(sol); SMD CDS
and buried SASA. 15 cases, all converged, no broken/new bonds, frozen anchors fixed, closure ≤ 1e-9.

- `summary.json`, `analysis.json`, `nonpolar.json`, `screen/`: from `runs/hplc_patch_extended_cpcm_def2svp/`.
- `cases.json`: all 40 placements (8 reagents × 5) from the three run dirs (`check_patch_cases.py`, run from repo root).
- `ranking.json` / `ranking8.json`: `rank_patch_reagents.py cases.json ranking8.json` — n = 8 vs Donegan Table 4 order
  (100 mM acetate, 20T %ACN): DEA < TEA < DIPEA < BA < TPA < DBA < HA, TBA placed last.

Spearman ρ (n = 8): EDA total 0.45 (best/Boltzmann, p 0.26), EDA+CDS 0.43, buried npSASA 0.38–0.52;
generalized correlation 0.71–0.81 (p 0.02–0.05), electrostatic 0.71–0.83 (p 0.01–0.05). 20 metric/aggregate
combinations tested, so component p-values are not multiplicity-corrected. Without TBA (n = 7) total ρ = 0.25.
Main outlier: DEA⁺ (least retained) binds −6.3 (Boltzmann), 3rd strongest; 4/5 placements sit at the patch centre
(lateral r 1.2–1.8 Å), likely helped by the open 3-chain pocket.

Command: `queue_hplc_patch_extended.sh` (systemd unit `hplc-patch-extended`, ended 2026-11-11 16:22).
