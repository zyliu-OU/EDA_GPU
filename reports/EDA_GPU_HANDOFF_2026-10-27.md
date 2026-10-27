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

## RUNNING (started 2026-10-27 09:21) — patch placement sampling

- Seeds: `scripts/build_hplc_patch.py --keep 10 --output examples/hplc/patch_sampling` (same RNG; ranks 0–1 were
  identical to the finished runs and were deleted) → 8 extra placements per reagent, `examples/hplc/patch_sampling/`.
- xTB screen (done, 40/40 ok, no duplicates of each other or of the DFT'd placements 0/1):
  `scripts/screen_hplc_patch.py` → `runs/hplc_patch_sampling_cpcm_def2svp/screen/{screen,selection}.json`.
  Picked (lowest dE_xTB): DBAA 8,3,2; DIPEA 5,8,6; HAA 3,5,4; TEAH 3,2,4; TBuAA 5,3,9. TBuAA_5/_3 (−10.3/−9.3
  kcal/mol xTB) are stronger than any other reagent's placements.
- Unit `hplc-patch-sampling.service` (`systemd-run --user`) runs `runs/queue_hplc_patch_sampling.sh`; log
  `runs/hplc_patch_sampling.log`; output `runs/hplc_patch_sampling_cpcm_def2svp/`. 15 cases, best pick of every
  reagent first; then nonpolar.json, analysis.json, "queue done". Thermal guards on (defaults).
- Estimate ≈ 1.3 h per case (0.5–2.6) → ≈ 20 h (12–30 h), done ≈ 2026-10-28 morning.
- Check: `systemctl --user status hplc-patch-sampling`, `grep -v -i warn runs/hplc_patch_sampling.log | grep '^patch_'`.
- Linger=no: survives closing the terminal, stops on logout/reboot. Resume (skips finished cases):
  `systemd-run --user --unit=hplc-patch-sampling --working-directory=$PWD /bin/bash -c 'exec /bin/bash runs/queue_hplc_patch_sampling.sh >> runs/hplc_patch_sampling.log 2>&1'`

## Next (in order)

1. When done: same checks as the 10-27 patch set (convergence, topology, N height/lateral position); combine with
   placements 0/1 (5 per reagent), tabulate best/Boltzmann-weighted EDA and EDA + CDS, Spearman vs published order;
   check xTB-vs-DFT ranking; archive `reports/hplc_patch_sampling_<date>/`, VALIDATION section, README, commit/push.
   If TBuAA still fails to rank high, build a larger/denser brush (5–7 chains).
2. Sensitivities (handoff 10-19): aug-cc-pVDZ, vdw_scale 1.2 vs 1.4, acetonitrile ε.
3. Bases on the patch only if the reagent order is reproduced; retention data via `scripts/hplc_retention.py`.
