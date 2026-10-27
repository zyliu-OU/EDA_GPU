# EDA_GPU handoff — 2026-10-27

Previous: `reports/EDA_GPU_HANDOFF_2026-10-26.md` (SMD CDS term, silanol-free patch built, queue launched).

## Done this session

- Patch queue (`hplc-patch.service`) finished 2026-10-27 06:58; unit gone, nothing running.
- Checked all 10 reagent···patch cases: converged, topology preserved, frozen anchors ≤ 2.7e-8 Å, EDA closure
  ≤ 7.3e-10 Eh. No reagent went below the anchor plane; DBAA (both placements) inserted between the three chains
  (N 1.6–2.2 Å from triangle centre, ~1.5–2 Å below the chain tips); all others sit outside or on top.
- Archive `reports/hplc_patch_2026-10-27/`; `docs/VALIDATION.md` section "Silanol-free C8 patch (application),
  2026-10-27"; README "gets wrong" + next steps updated.
- Result: chain contact −4.9 to −8.9 kcal/mol (≈ 2× single-C8 T···S), TEAH⁺ weakest; order TEAH⁺ < TBuAA <
  DIPEA-H⁺ < HAA < DBAA, Spearman ρ = 0.3 (EDA) / 0.6 (EDA + CDS water) / 0.5 (buried SASA), n = 5, vs published
  TEA < DIPEA < DBA < HA < TBA. Placement spread up to 1.8 kcal/mol ≈ reagent spread.

## Next (in order)

1. Placement sampling on the patch (≥ 5 seeds per reagent, e.g. xTB screening then DFT on the best 2–3), and/or a
   larger/denser brush (5–7 chains) that can wrap TBuAA/DIPEA. Not started.
2. Sensitivities (handoff 10-19): aug-cc-pVDZ, vdw_scale 1.2 vs 1.4, acetonitrile ε.
3. Bases on the patch only if the reagent order is reproduced; retention data via `scripts/hplc_retention.py`.
