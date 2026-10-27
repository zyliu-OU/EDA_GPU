# EDA_GPU handoff — 2026-10-26

Previous: `reports/EDA_GPU_HANDOFF_2026-10-19.md` (stages 7.1–7.3, reagent series, bridged ternaries, literature comparison).

## Done this session (committed + pushed: 23b717c, 1903050)

1. **SMD CDS non-electrostatic term.** `src/eda_gpu/cds.py` (PySCF/MNSOL, geometry only, CPU), `tests/test_cds.py`,
   `scripts/nonpolar_hplc.py` (dG_CDS water/acetonitrile + buried SASA + EDA + CDS per pair and three-body split).
   Archive `reports/hplc_nonpolar_2026-10-26/`; `docs/VALIDATION.md` section "SMD CDS term"; README updated.
   Result: does NOT supply the hydrophobic term. Chain contact gains only 4–6 cal/mol/Å² of buried nonpolar area
   (vs ≈25 from alkane transfer); no reagent or base rank changes. Bridged T···S order TEAH < DBA < DIPEA < TBA < HA,
   ρ = 0.8 (p 0.10, n 5) vs published TEA < DIPEA < DBA < HA < TBA, with or without CDS.
2. **Silanol-free three-chain patch.** `scripts/build_hplc_patch.py` → `examples/hplc/c8x3_patch.{xyz,json}` (3 tethered
   all-anti n-octanes on a 7.6 Å triangle, C1/C2 of each frozen, 78 atoms, no polar site) and 10 seeds
   `examples/hplc/patch/patch_<reagent>_{0,1}` (101–119 atoms; 2 placements per reagent ≥3 Å apart). Rationale and
   model: `docs/HPLC_MODELS.md` (last section). No residue: bridged ternaries were additive with R···S ≈ 0.

## RUNNING (started 2026-10-26 ~18:00; restarted 18:12 as a systemd user service)

- Unit `hplc-patch.service` (transient, `systemd-run --user`) runs `runs/queue_hplc_patch.sh`; log `runs/hplc_patch.log`;
  output `runs/hplc_patch_cpcm_def2svp/` (`summary.json`, per-case `xtb/`, `optimization/`, `eda/`).
  Pipeline per seed: GFN2-xTB/ALPB pre-opt → M06-2X/def2-SVP DF C-PCM water opt (maxsteps 300) → GKS-EDA(sol).
  Then `nonpolar.json` and `analysis.json` in the same directory, then "queue done" in the log.
  Order: all `_0` seeds (DBA, DIPEA, HA, TEAH, TBA), then all `_1`. Thermal guards on (GPU 80/65 °C, CPU 85/75 °C).
- Speed ~2.1 min per DFT step at 107 atoms. Estimate 2–5 h per case → `_0` set ~10–25 h, all ~20–50 h (≈ 10-27 to 10-28).
- `patch_dibutylammonium_0/optimization_failed_0/` is the ~5-step run killed when moving the job into systemd (not a real failure).
- Check: `systemctl --user status hplc-patch`, `grep -v -i warn runs/hplc_patch.log | grep '^patch_'`.
- Survives closing the terminal / Claude. Linger=no, so it stops on logout or reboot. Resume (skips finished cases):
  `cd /home/zyliu/Documents/EDA_GPU && systemd-run --user --unit=hplc-patch --working-directory=$PWD /bin/bash -c 'exec /bin/bash runs/queue_hplc_patch.sh >> runs/hplc_patch.log 2>&1'`
  (if an optimization was cut off, it restarts that case from the xTB geometry).

## Next (in order)

1. When the queue finishes: check every case converged with topology preserved, and that the reagent stayed on the
   chains (reagent-C···patch-C contacts, N height above C1 plane; no burial of N–H inside the chains or below the anchor plane).
2. Tabulate reagent···patch EDA (components), EDA + dG_CDS, buried nonpolar SASA per reagent (best and both placements);
   Spearman vs published TEA < DIPEA < DBA < HA < TBA; compare with bridged T···S (single C8) and with the silanol-bound pairs.
3. Archive `reports/hplc_patch_<date>/`, add a `docs/VALIDATION.md` section, update README "gets wrong"/next steps, commit/push.
4. Then (handoff 10-19 list): aug-cc-pVDZ / vdw_scale 1.2 vs 1.4 / acetonitrile ε sensitivity; placement sampling;
   bases on the patch if the reagent order is reproduced; retention data via `scripts/hplc_retention.py` when provided.
