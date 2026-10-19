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
3. External XEDA run still waits on the user (`reports/xeda_matched_2026-10-18/README.md`).

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
