# EDA_GPU handoff — 2026-10-28

Previous: `reports/EDA_GPU_HANDOFF_2026-10-27.md` (patch results, placement sampling launched).

## Done this session

- Placement-sampling queue (`hplc-patch-sampling.service`) finished 2026-10-28 02:23 (17 h wall); unit gone,
  nothing running.
- Checked all 15 new cases (and re-derived the 10 from 10-27, which match the earlier table): converged, no
  broken/new bonds, frozen anchors ≤ 1.6e-8 Å, closure ≤ 1.2e-9 Eh, no reagent below the anchor plane.
- Archive `reports/hplc_patch_sampling_2026-10-28/` (incl. `check_patch_cases.py`, `rank_patch_reagents.py`,
  `cases.json` for all 25 placements, `ranking.json`); `docs/VALIDATION.md` section "Silanol-free C8 patch —
  placement sampling, 2026-10-28"; README "gets wrong" + next steps updated.
- Result: lowest-energy placement = strongest-bound in every reagent; order TEAH⁺ < DIPEA-H⁺ < HAA < TBuAA <
  DBAA, ρ = 0.7 (EDA; was 0.3), 0.6 (EDA + CDS), 0.7 (buried SASA), n = 5. TBuAA 3/5 insert over the triangle
  like DBAA 0/1 (the 10-27 "TBuAA can't be wrapped" was sampling). Remaining misses: DBAA > TBuAA by 1.7, HAA low.
  xTB screen: within-reagent best pick right 2/5; overbinds 1.1–3.5 kcal/mol.

## Next (in order)

1. Larger/denser brush (5–7 chains, e.g. hexagon + centre; pick the spacing from a bonded-phase density, the
   current triangle side is 7.6 Å): extend `scripts/build_hplc_patch.py`; seeds for the 5 reagents with
   insertion-favouring placement (N over the brush centre); xTB screen; DFT the best 2–3 per reagent. Check
   cost first (≥ 150 atoms with TBuAA).
2. Sensitivities (handoff 10-19): aug-cc-pVDZ, vdw_scale 1.2 vs 1.4, acetonitrile ε.
3. Bases on the patch only if the reagent order is reproduced; retention data via `scripts/hplc_retention.py`.

## Later 2026-10-28 — 7-chain brush timing test (not committed, nothing running)

- `scripts/build_hplc_patch.py --chains 7` added (centre + hexagon at 7.6 Å; `--chains 3` output byte-identical);
  wrote `examples/hplc/c8x7_patch.{xyz,json}` (182 atoms). Uncommitted.
- Timing test (scratchpad, TBuAA seed, 223 atoms, M06-2X/def2-SVP DF C-PCM, as production): setup 4.5 min, then
  54 s per SCF cycle (vs ≈ 4 s at 101 atoms), then **CUDA out of memory at the first gradient** (16 GB card).
  Production run as configured is not possible. Waiting on user to choose: xTB geometry + DFT EDA single point,
  5-chain brush timing test, or a memory fix in the engine.
