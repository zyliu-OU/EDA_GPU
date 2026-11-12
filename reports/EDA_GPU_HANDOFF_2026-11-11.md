# EDA_GPU handoff — 2026-11-11

Previous: `reports/EDA_GPU_HANDOFF_2026-11-10.md`. Nothing running (machine rebooted 21:30 after both queues ended).

## Step 2 done — aug-cc-pVDZ: `reports/hplc_patch_augdz_2026-11-11/`
5/5 OK. Uniform −1.9 to −3.3 shift (generalized correlation), ρ(SVP,aDZ) 0.9; def2-SVP ranking is basis-robust
except TEA/DIPEA (0.25 apart).

## Step 3 done — BA⁺, DEA⁺, TPA⁺: `reports/hplc_patch_extended_2026-11-11/`
15/15 clean. n = 8 Spearman vs Donegan: total 0.45 (p 0.26), +CDS 0.43; GC 0.71–0.81, elec 0.71–0.83 (uncorrected).
Adding the three reagents lowered total ρ from 0.7 (n = 5) — DEA⁺ is the outlier (binds like a long-chain reagent;
sits in the centre pocket of the 3-chain patch). TPA⁺ and BA⁺ land mid-pack, plausible.

## Not done / pending
- Not committed: archives above, reagent/retention script edits, new seeds, queue scripts; docs/VALIDATION.md and
  README sections for steps 2–3 not yet written.
- Still uncommitted on purpose: 7-chain builder diff, `c8x7_patch.*`, `reports/c8x7_memtest_2026-11-09/`.
- Next per plan: step 5 retention fit (consider fitting on components, and a DEA pocket-occupancy check), then step 4
  (bases on the patch).
