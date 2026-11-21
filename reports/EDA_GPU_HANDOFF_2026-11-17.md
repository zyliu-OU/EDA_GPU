# EDA_GPU handoff — 2026-11-17

Previous: `reports/EDA_GPU_HANDOFF_2026-11-16.md` (queue design, splay check, open options).

## User instructions today
- Stop the bases queue after U_1 (end of the placement 0/1 stage). The user is rebooting the computer on 2026-11-18.
- **Do not start any new run** after it ends until the user says so (memory `no-new-runs-after-bases`). Analysis of
  finished results is fine.

## State at 2026-11-17 ~19:00 PST
- systemd user unit `hplc-patch-bases` (bash PID 23169, python PID 23176) is running stage 1 of
  `runs/queue_hplc_patch_bases.sh`: placements 0/1 of A/C/G/U, xTB pre-opt → DFT opt → GKS-EDA(sol).
- Done, all converged, closure ≤ 4e-10 Eh, topology preserved (total EDA, kcal/mol):
  A_0 −9.72, A_1 −11.58, C_0 −7.04, C_1 −7.50, G_0 −10.65, G_1 −8.78.
- Running: U_0 DFT opt (step ~49, close to converged), expected done ~20:00–22:00 tonight. Then U_1 (~3–4 h),
  expected done ~01:00–03:00 on 11-18.
- **Queue truncation:** `runs/queue_hplc_patch_bases.sh` was edited in place (same inode, first 651 bytes, which is
  bash's read offset, left unchanged). After the stage-1 python exits, the queue prints
  "stopped after placements 0/1 ..." and exits, so the xTB screen and the 12 extra DFT cases do **not** run.
  The full original queue is kept as `runs/queue_hplc_patch_bases.full.sh`.

## After the reboot — check first
1. `tail -3 runs/hplc_patch_bases.log`. It should end with "stopped after placements 0/1".
2. `python3 -c "import json;print({k:v.get('eda') for k,v in json.load(open('runs/hplc_patch_bases_cpcm_def2svp/summary.json'))['cases'].items()})"`
   should show 8 cases with `eda: success`.
3. If the reboot happened before U_1 finished: the case is incomplete. Completed stages are skipped on restart,
   so rerunning resumes it. **Only do this if the user asks** (no-new-runs rule).

## Resume commands (only when the user approves new runs)
- Finish the interrupted stage 1 only:
  `systemd-run --user --unit=hplc-patch-bases --working-directory=$PWD /bin/bash -c 'exec /bin/bash runs/queue_hplc_patch_bases.sh >> runs/hplc_patch_bases.log 2>&1'`
- Run the full protocol (stage 1 skips done cases, then xTB screen of placements 2–9, best 3 per base → 12 DFT cases
  of ~3.5 h each ≈ 2 days, then nonpolar + analysis):
  `systemd-run --user --unit=hplc-patch-bases --working-directory=$PWD /bin/bash -c 'exec /bin/bash runs/queue_hplc_patch_bases.full.sh >> runs/hplc_patch_bases.log 2>&1'`
- Linger is off. Logging out stops the unit.
- Thermal guards: `run_hplc_stage.py` defaults (GPU 80 °C, CPU 85 °C). Last seen: GPU 58 °C, CPU 75.5 °C.

## Allowed without new runs
- Analyze the 8 placement-0/1 cases: convergence and closure, base-ring position in the hollow (adapt
  `check_patch_cases.py`, using the base ring centroid instead of the first guest N), Boltzmann EDA per base over 2
  placements, and the components. Compare the A/C/G/U order with the user's base-ratio retention data and with the
  reagents. Mark it as preliminary (2 placements per base). `nonpolar_hplc.py` / `analyze_hplc.py` are CPU-side
  post-processing. Check with the user before running them if in doubt.

## Open options (unchanged from 11-16, waiting on the user)
Full base sampling (above), off-centre re-sampling for DEA/BA/DBA, 4-chain patch memory test.

## Uncommitted
Everything listed in the 11-16 handoff, plus this handoff and `runs/queue_hplc_patch_bases.full.sh` (runs/ is
probably gitignored). The user has not said whether to commit/push.
