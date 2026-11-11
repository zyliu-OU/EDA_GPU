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

## 2026-11-09 — larger brush dropped

- User: the machine (32 GB RAM, 16 GB GPU) cannot run larger systems. The 5–7-chain brush (next step 1) is dropped;
  the 3-chain C8 patch stays the production surface model. 7-chain builder and `c8x7_patch.*` remain uncommitted.
- Remaining: sensitivities (memory-check aug-cc-pVDZ first), bases on the 3-chain patch, retention fit.

## 2026-11-09 (later) — GPU memory fix for the 7-chain brush (option B)

- Correction: the 10-28 OOM was **GPU** memory (16 GB), not host RAM; the SCF fit. 7-chain brush back on the table.
- Cause (read in gpu4pyscf 1.8.1 `df/grad/rhf.py`): the DF gradient calls `with_df.reset()` then sizes its
  batches from `cudaMemGetInfo` (`get_avail_mem(exclude_memory_pool=True)`), which does not count the blocks cupy
  keeps cached in its pool after the SCF.
- Fix (uncommitted): `src/eda_gpu/engine.py` `release_gpu_cache()` + `gradient_method()`; every GPU gradient
  (single point and optimizer scanner) first resets the DF tensor and frees cupy's device + pinned pools;
  `make_mf` frees the pools before each new SCF (EDA too; DF tensor not reset there). 46/46 tests pass.
- Estimate for 7-chain + TBuAA (nao 1736, naux 7949, nocc 284): fixed K-gradient buffers j3c_oo + dm_oo
  ≈ 10.3 GB peak; rest auto-batched → ~11.5–12 GB expected.
- Memory test running: production-settings single point (SCF + gradient) of the 7-chain TBuAA seed, scratchpad
  `c8x7_memtest/` (`mem.log` = time, GPU MiB, GPU °C, host MiB, swap MiB every 5 s).
- Next if it passes: commit builder + fix; seeds for all 5 reagents (`--chains 7`), xTB screen, DFT best 2–3 each.
- Moved to systemd user unit `c8x7-memtest` (survives closing the window); files in
  `reports/c8x7_memtest_2026-11-09/` (`run.out` ends with rc=, `mem.log`, `run.err` has `time -v` peak RSS).
  First attempt showed host RAM 28.5 GB + 2.5 GB swap during SCF setup: host RAM is tight too.
- **Memory test result (21:18–21:33): failed, GPU OOM again in the gradient.** SCF converged (cycle 11,
  E = −2735.103088 Eh, ~13.7 min). At the gradient, GPU use went 8.5 → 15.3 GB, then OOM in
  `df/grad/rhf.py _jk_energy_per_atom → int2c2e(auxmol) → apply_CT_mat_C`, i.e. after the large K-gradient
  buffers were already allocated. So the cache-release fix works (the SCF pool was freed first) but the fixed
  buffers (~10+ GB) plus the aux 2c2e step don't fit in 16 GB. Host peak RSS 25.7 GB (of 32). rc=1, unit gone.
  Options: xTB geometry + DFT EDA single point only (no DFT gradient), or a 5-chain brush.

## 2026-11-09 (late) — 7- and 5-chain dropped for good; 3-chain only

- Size estimate (def2-SVP, gradient buffers are built in the Cartesian basis, which is why the earlier estimate was low):
  3 / 5 / 7 chains + TBuAA = 119 / 171 / 223 atoms, naux(cart) 4937 / 7081 / 9225, nocc 152 / 218 / 284;
  K-gradient j3c_oo 0.9 / 2.7 / 6.0 GB, and the gradient holds it twice (j3c_oo + dm_oo).
- A CPU PCM-grid estimate (pyscf `PCM.build()` builds dense ngrid² matrices) hit host OOM at 26 GB and was
  killed (22:03), taking the terminal with it. Don't build PCM/large matrices on CPU for 7-chain-sized systems.
- User decision: focus on the 3-chain patch data only. No 5- or 7-chain work. 7-chain builder + `c8x7_patch.*`
  stay uncommitted (not needed); the engine GPU-cache release fix is general and can be kept.
- Next (3-chain): sensitivities (aug-cc-pVDZ memory check first, vdw_scale 1.2 vs 1.4, acetonitrile ε),
  bases on the patch, retention fit (`scripts/hplc_retention.py`).
