# EDA_GPU handoff — 2026-10-19

Previous: `reports/EDA_GPU_HANDOFF_2026-10-18.md` (validation status, Psi4/XEDA, GKS-EDA(sol)).

## Stage 7.1 + 7.2 results (complete; M06-2X/def2-SVP, DF+C-PCM water optimization, direct GKS-EDA(sol))

Lowest-PCM-free-energy placement per group (`scripts/analyze_hplc.py`; kcal/mol; closure ≤ 3e-10 Eh, topology preserved in all converged cases):

| Group | n | elec | exch | rep | pol | gcorr | desolv | total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| residue_A···TEAH+ | 2/2 | -108.31 | -57.63 | 102.38 | -25.15 | -10.87 | 80.83 | -18.77 |
| residue_C···TEAH+ | 2/2 | -113.42 | -61.76 | 110.23 | -27.44 | -11.09 | 83.74 | -19.74 |
| residue_G···TEAH+ | 1/2 | -112.01 | -62.18 | 109.30 | -25.90 | -12.99 | 84.15 | -19.64 |
| residue_U···TEAH+ | 2/2 | -109.62 | -62.00 | 110.28 | -26.47 | -12.71 | 80.54 | -19.98 |
| C8 proxy···TEAH+ | 2/2 | -25.84 | -32.65 | 56.12 | -6.27 | -12.09 | 9.40 | -11.33 |
| residue_A···C8 | 1/1 | -27.95 | -57.14 | 97.05 | -14.35 | -15.51 | 2.99 | -14.92 |
| residue_C···C8 | 1/1 | -39.65 | -67.40 | 115.80 | -18.85 | -16.28 | 7.52 | -18.87 |
| residue_G···C8 | 1/1 | -35.09 | -65.79 | 110.64 | -13.93 | -19.32 | 8.06 | -15.43 |
| residue_U···C8 | 1/1 | -36.49 | -60.95 | 104.44 | -16.46 | -15.69 | 10.29 | -14.86 |

Single local placements, def2-SVP, no sampling: base differences (~1–4 kcal/mol) are within plausible
placement/basis/cavity uncertainty and are NOT established. G_13 ion pair: 200-step cap reached (flat creep).
Per-seed values: `runs/hplc_{ionpair,surface}_cpcm_def2svp/summary.json`.

## Session 2 (2026-10-19, after reboot)

- G_13 continuation skipped: G_14 is 2.9 kcal/mol lower and G_13 is flat, so it cannot change the ternary host.
- Stage 7.1/7.2 archived in `reports/hplc_{ionpair,surface}_2026-10-19/` and documented in `docs/VALIDATION.md`;
  `docs/HPLC_MODELS.md` status updated. `scripts/analyze_hplc.py` now creates its output directory.
- Committed and pushed: the archive plus `--continue-unconverged`, `--host-xyz`, `build_hplc_ternary.py`, and the
  four ternary seeds in `examples/hplc/ternary/` (A_18, C_17, G_13, U_16 = free phosphate O; 104–108 atoms;
  charges −1/+1/0; min heavy-atom clearance 2.9–3.7 Å).
- Stage 7.3 optimization RUNNING in the background: `runs/queue_hplc_ternary.sh` → `runs/hplc_ternary.log`,
  output `runs/hplc_ternary_cpcm_def2svp/`. Re-running the same script resumes (converged cases are skipped).

## Next (in order)

1. When stage 7.3 finishes, check each case for convergence and preserved topology in `runs/hplc_ternary_cpcm_def2svp/summary.json`, then
   per base: `.env/bin/python scripts/three_body_hplc.py runs/hplc_ternary_cpcm_def2svp/<case>/optimization/optimized.xyz --metadata examples/hplc/ternary/<case>.json --names residue_X teah c8_silane_proxy --output runs/hplc_three_body/<case>`;
   archive and document like 7.1/7.2.
2. Sensitivity: aug-cc-pVDZ EDA(sol) on the lowest pair per base; vdw_scale 1.2 vs 1.4; acetonitrile ε 35.688.
3. (Dropped 2026-10-19 by user decision: external XEDA run.)

## Earlier stopped state (session 1, superseded)

- Stopped for a computer restart. `runs/queue_hplc_followup.sh` was killed before stage 7.3; the G_13
  continuation was killed at step ~114 (gradient rms 6e-5, max 3.5e-4 Eh/bohr; energy change since the
  capped run ~2e-6 Eh; geomeTRIC trust radius at its 1e-4 floor, i.e. a flat surface).
  Its directory is kept as `runs/hplc_ionpair_cpcm_def2svp/residue_G_teah_13_0/optimization_interrupted_reboot/`
  (no coordinates saved mid-run); the 200-step capped run was moved back to `optimization/` so
  `--continue-unconverged` resumes from its last geometry. summary.json still records G_13 as not converged.
- Resume after reboot:
  `cd /home/zyliu/Documents/EDA_GPU && nohup /bin/bash runs/queue_hplc_followup.sh > runs/hplc_followup.log 2>&1 &`
  (G_13 continuation, then ternary build and optimization). Or skip G_13 (G_14 already provides the G ion pair)
  by running only the last two lines of that script.
- Stage 7.3 was NOT started: no ternary seeds in `examples/hplc/ternary/`, no `runs/hplc_ternary_cpcm_def2svp/`.
- Nothing is running. Nothing new was committed or pushed in this session after `3e90de9`.

## Session-1 uncommitted work (now committed in session 2)

- `scripts/run_hplc_stage.py`: `--continue-unconverged` (restart a capped optimization from its last geometry) — used for G_13.
- `scripts/build_hplc_complexes.py`: `--host-xyz` (optimized host geometry) — dry-run tested in scratch.
- `scripts/build_hplc_ternary.py` (new): ternary seeds on the lowest ion pair per base — dry-run built 4 seeds
  (107 atoms, charges −1/+1/0, anchors shifted correctly), not committed to `examples/`.
- `runs/queue_hplc_followup.sh` (runs/ is local).
- This handoff.


## Correction (session 2)

`docs/VALIDATION.md` (and the 2026-10-18 handoff, line 93) said the tryptamine–water Pauli offset vs the
PCCP 2017 Table 3 was +0.15–0.20 kcal/mol in both charge states. Actually neutral −0.15, cation +0.20
(opposite signs). VALIDATION.md corrected. Asked the user (the paper's author) for original 2017
coordinates/GAMESS outputs, which would allow an identical-geometry comparison without installing XEDA.

## Decision (2026-10-19)

The user (project owner) accepts the PCCP 2017 tryptamine–water agreement (totals ≤0.012, components ≤0.20 kcal/mol)
as sufficient external validation. The matched XEDA run and the original-coordinate comparison are dropped from the
next steps. README/VALIDATION updated; `scientifically_validated` left false in code (strict definition = matched XEDA output).
- Follow-up: per the user, gas-phase GKS-EDA results now set `scientifically_validated: true` (basis: the published benchmark,
  stated in `reference_component_validation`); solvated results stay `false`; `full_eda_available` stays false
  (reserved for matched XEDA). Tests updated (`tests/test_gks.py`, `tests/test_solvent.py`); 13 CPU tests pass.
  The running stage 7.3 process loaded the old code but uses `--skip-eda`, so it is unaffected.

## Session 2, later: reagent series, retention comparison, thermal guard

- `scripts/build_hplc_reagents.py`: hexylammonium, dibutylammonium, tributylammonium, DIPEA-H+ monomers
  (`examples/hplc/<name>.xyz/json`) and 23 seeds in `examples/hplc/reagents/` (8 residue_A ion pairs; 15 C8 surface
  pairs incl. chain-contact seeds; one duplicate contact seed removed). Not yet optimized.
- Fixed `build_hplc_complexes.py` placement ranking: clearance included the acceptor/donor-X atoms and was identical
  for every placement (choice fell to sampling order). TEAH+ 7.1–7.3 seeds predate the fix (valid, not meaningfully
  ranked). Added `--select contact`. `analyze_hplc.py` folds `_contact` seeds into their pair group.
- `scripts/hplc_retention.py` + `examples/hplc/retention_template.csv`: correlate measured retention with
  composition-weighted EDA descriptors per condition group (LOO RMSE vs group mean). Tested on synthetic data only.
- `run_hplc_stage.py`: `--cpu-temperature-limit` (default 85 C, resume 75) and `--gpu-temperature-limit` (80, resume 65).
  The running ternary job predates this (GPU guard 80 C only, CPU guard off).
- User asked for CPU/GPU power caps (no overheating). Hardware caps need sudo (user must run; see chat).
- `runs/queue_hplc_reagents.sh` prepared (waits for the ternary job, then reagent ion pairs + surface pairs with EDA);
  REPLACED (PID 23198, log `runs/hplc_reagents.log`): waits for the ternary optimization, then three-body EDA
  (`runs/hplc_three_body/<case>/three_body.json`) on each converged ternary, then the reagent ion pairs and surface
  pairs with EDA. No hardware-cap gate (would stall unattended); software guards GPU 80/65 C, CPU 85/75 C on all runs;
  cap status is logged. User authorized continuing the build without asking. After a reboot: re-apply caps if wanted,
  restart `runs/queue_hplc_ternary.sh` and then `runs/queue_hplc_reagents.sh` (both resume).

## 2026-10-20: stage 7.3 complete

- Ternary optimizations: all four converged, topology preserved (archive `reports/hplc_ternary_2026-10-20/`).
- Three-body GKS-EDA(sol): ΔE₃ = +0.97 to +1.19 kcal/mol (small, anti-cooperative) for all bases; closure ≤ 2.1e-10 Eh.
  Table in `docs/VALIDATION.md`.
- Model caveat: this ternary arrangement has the residue H-bonded directly to the silanol (residual-silanol
  picture) and TEAH⁺ barely touching the C8 proxy (T···S −0.3 to −2.1). The usual IP-RP picture is a *bridged*
  arrangement: reagent alkyl chains on the C8 ligand, phosphate paired with the reagent N–H, no direct R–S contact.
- Reagent series started automatically at ~05:50 (queue `runs/queue_hplc_reagents.sh`, log `runs/hplc_reagents.log`).

## Next (updated 2026-10-20)

1. Reagent series finishes → `analyze_hplc.py` on `runs/hplc_reagents_{ionpair,surface}_cpcm_def2svp`, archive, document.
2. Bridged ternaries: residue phosphate on the N–H of each optimized reagent···C8 *contact* pair (TEAH⁺ included:
   needs a TEAH⁺ contact seed on C8), optimize, three-body EDA. Compare with the silanol-bound arrangement above.
3. Sensitivity: aug-cc-pVDZ EDA(sol) on the lowest pair per base; vdw_scale 1.2 vs 1.4; acetonitrile ε 35.688.
4. When the user provides retention data: `scripts/hplc_retention.py`.

## 2026-10-21: reagent series complete; bridged ternaries running

- Reagent series: 25/25 converged + EDA (c8_hexylammonium_13_0 re-run with `--preopt none` after the xTB anchor gate
  rejected it at 0.020 Å). Archive `reports/hplc_reagents_2026-10-21/`; table in `docs/VALIDATION.md`.
- Finding: reagent···C8 EDA totals are flat (−10.6 to −12.7) and within placement spread; no cavitation term in
  C-PCM → hydrophobic effect missing. Added `scripts/buried_sasa.py` (FreeSASA cross-check −0.6 %); buried nonpolar
  area tracks hydrophobicity (TBA 223, DBA 183, DIPEA 170, HA 138, TEAH 139 Å² max).
- Bridged ternaries: `scripts/build_hplc_bridged.py` → `examples/hplc/bridged/` (5 reagents, residue_A).
  RUNNING: `runs/queue_hplc_bridged.sh` → `runs/hplc_bridged.log` (optimize, then three-body into
  `runs/hplc_three_body_bridged/`).

## Next (updated 2026-10-21)

1. Bridged ternaries finish → archive, compare with the silanol-bound ternary (TEAH⁺) and across reagents.
2. Sensitivity: aug-cc-pVDZ EDA(sol) on selected pairs; vdw_scale 1.2 vs 1.4; acetonitrile ε 35.688.
3. Placement sampling for the base comparison (more seeds per base, xTB-ranked).
4. When the user provides retention data: `scripts/hplc_retention.py` (add buried-SASA descriptors).

## 2026-10-21 08:40 status — user: stop after the current queue, save handoff

- Running: `runs/queue_hplc_bridged.sh` (PID 393702). Bridged optimizations: dibutylammonium converged (33 evals,
  1.9 h); DIPEA-H⁺ at step ~39 (2.6 h); hexylammonium, TEAH⁺, tributylammonium (107–125 atoms) still to run; then
  three-body EDA on each converged case (≈40–60 min each). Estimated finish ≈ late evening 2026-10-21 to early 10-22.
- Plan agreed with the user: when this queue finishes, archive bridged results to `reports/hplc_bridged_<date>/`,
  document in `docs/VALIDATION.md`, update this handoff, commit/push, and STOP — do not start the sensitivity stage.

## 2026-10-22: bridged ternaries complete — STOPPED (per user plan)

- `runs/queue_hplc_bridged.sh` finished 2026-10-22 01:31. All 5 converged, topology preserved; three-body EDA on all
  5 (closure ≤ 4.0e-10 Eh). Archive `reports/hplc_bridged_2026-10-22/`; table in `docs/VALIDATION.md`.
- Findings: ΔE₃ ≈ 0 (|ΔE₃| ≤ 0.3) — bridged complex is pairwise additive. R···T −18.3 to −24.7 (DBA/TBA strongest),
  T···S only −2.2 to −3.8, R···S ≈ 0 (DIPEA −1.0).
- Re-interpretation: all 17 isolated reagent···C8 pairs of the reagent series optimized to an N–H···O(silanol)
  H-bond (1.75–2.04 Å), so their −10.6 to −12.7 is mostly that H-bond, not chain contact. Chain-only ≈ −2 to −4.
- Nothing running. Sensitivity stage NOT started (user instruction).

## Next (updated 2026-10-22)

1. Sensitivity: aug-cc-pVDZ EDA(sol) on selected pairs; vdw_scale 1.2 vs 1.4; acetonitrile ε 35.688.
2. Decide the surface model for retention: silanol-bound vs bridged (or both, weighted); consider endcapped proxy
   (no free silanol) to test the N–H···silanol contribution directly.
3. Placement sampling (more seeds per base/reagent, xTB-ranked).
4. When the user provides retention data: `scripts/hplc_retention.py` (add buried-SASA descriptors).
