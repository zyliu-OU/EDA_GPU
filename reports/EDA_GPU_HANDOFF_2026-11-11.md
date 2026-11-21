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

## Step 5 done — retention fit: `reports/hplc_retention_fit_2026-11-11/`
New `scripts/fit_patch_retention.py`, data `examples/hplc/donegan2022_table4.csv` (Table 4 %ACN 20T/60T, Table 1
logP/logD). logD alone LOO 7.5 %ACN (mean 11.2); best EDA (off-centre-hollow Boltzmann total) 8.1, all-placement 13.1.
No EDA descriptor beats logD (permutation p_select 0.33). logD+pol (LOO 3.7) is HAA-leveraged → spurious.
Docs: VALIDATION.md section 2026-11-11, README lines. Not committed yet.
Open question for the user before step 4: the 3-chain patch's open hollow is the failure mode; the 7-chain brush
OOMed (10-28). Options: freeze/cap a hollow-blocking chain, constrain placements off-centre, or accept and move to bases.

## Plain summary of the last finished runs (as explained to the user, end of 2026-11-11)
- Run 1, aug-cc-pVDZ (5 best placements, original reagents): every reagent 2–3 kcal/mol more bound, nearly all
  generalized correlation (chain–chain dispersion); order unchanged (ρ 0.9; TEA/DIPEA swap, 0.25 apart).
  def2-SVP is adequate for reagent comparison; basis is not the cause of the mismatch with experiment.
- Run 2, DEA/BA/TPA (15 cases, all clean): Boltzmann EDA DEA −6.3, TEA −5.0, DIPEA −5.2, BA −5.8, TPA −5.4,
  DBA −8.7, HA −5.6, TBA −7.1. ρ vs published fell 0.7 (n 5) → 0.45 (n 8). DEA (least retained) is 3rd strongest
  because it sinks into the open central gap of the 3-chain patch; bulky TEA/DIPEA cannot.
- Fit: logD beats every EDA descriptor; removing centre-gap placements brings EDA close to, not past, logD.
- Takeaway: numerics are solid (convergence, basis, solvent tests); the model surface is the weak point.

## Options discussed for next session (user to decide; nothing started)
1. Stop the hollow opening: (a) restrain/freeze outer chain carbons (no extra atoms; rerun the 7 inserted placements
   DEA 4, BA 1, DBA 2, ~10 h); (b) a 4th chain/stub — full octane ≈145 atoms exceeds the hardware size limit.
2. Re-sample off-centre only (r ≥ 2.5 Å) for DEA/BA/DBA, ~10–15 h; post-hoc, at best matches logD.
3. Accept and go to step 4 (bases A/C/G/U on patch, 20 cases, ~1–2 days).
Recommended: first a minutes-long check whether chains splay apart in the inserted placements (from optimized
geometries), then 1a if they do, then step 4.
Still uncommitted: retention-fit archive, fit script, Donegan CSV, VALIDATION/README edits, this handoff
(and, on purpose, the 7-chain files). User has not yet said whether to commit/push.

## Late 2026-11-11 — PDF progress report
`reports/progress_report_2026-11-11/EDA_GPU_progress_report_2026-11-11.pdf` (13 pp., research-paper style: methods,
validation, HPLC results, hollow failure mode, retention fit, progress table, open issues/decisions). Sources: `report.html`
(rendered with headless Chrome), `make_figs.py` (figures from archived JSON/XYZ only), `fig/`. Uncommitted. No runs started.
