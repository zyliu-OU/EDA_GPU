# EDA_GPU handoff — 2026-11-18

Previous: `reports/EDA_GPU_HANDOFF_2026-11-17.md`.

## State
- Stage 1 (placements 0/1) finished before the reboot. All 8 cases have `eda: success`. U_1 total is −8.12 kcal/mol
  (U_0 is in summary.json).
- 2026-11-18 23:24 PST: the user said "start the run". The full protocol was launched as systemd user unit
  `hplc-patch-bases` (bash PID 29554) running `runs/queue_hplc_patch_bases.full.sh` and appending to
  `runs/hplc_patch_bases.log`. Stage 1 skips the done cases. Next is the xTB screen of placements 2–9, then the best
  3 per base → 12 DFT opt + GKS-EDA(sol) cases (~3.5 h each, ≈2 days), then nonpolar + analysis.
- Thermal guards: run_hplc_stage.py defaults (GPU 80 °C, CPU 85 °C). GPU was at 50 °C at launch. Hardware power caps
  (`nvidia-smi -pl`) reset at reboot and need user sudo; the limit is currently 320 W.
- Linger is off. Logging out stops the unit. Completed stages are skipped on restart, so relaunch with the same command.

## After this queue
Do not start another run without asking (memory `no-new-runs-after-bases`). Open options are unchanged: off-centre
re-sampling for DEA/BA/DBA, 4-chain memory test.
