# Chain-splay check on the 40 optimized patch placements (2026-11-16)

`python3 reports/hplc_patch_splay_2026-11-16/splay_check.py` (geometry only, no new calculations; output `splay.log`,
`splay.json`). Reads `optimized.xyz` of all 40 reagent-on-patch cases (archives 10-27, 10-28, 11-11) and `r`, `h`,
`tot` from `reports/hplc_patch_extended_2026-11-11/cases.json`. "Inserted" = r < 2.5 A and h < 8 A, as in the
retention fit (DEA 4, BA 1, DBA 2).

Per case: mean chain-chain distance at each carbon level (C1..C8), the three C8-C8 distances, the smallest distance
of any upper patch atom to the patch axis (`free`), reagent heavy atoms within 3 A of the axis and below the chain
tops (`gin`), gauche count.

| Group | n | chain-chain, C5-C8 mean (A) | C8-C8 mean (A) | widest C8-C8 (A) | free radius (A) | closest chain-chain atoms (A) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| constructed seed | - | 7.23 | 7.60 | 7.60 | 2.78 | - |
| inserted | 7 | 6.61 +- 0.36 | 6.58 | 8.01 | 2.28 | 2.51 |
| other | 33 | 6.81 +- 0.43 | 7.01 | 8.08 | 2.50 | 3.10 |

- The chains do not splay apart around an inserted reagent. In the inserted placements they stand closer together
  than in the others and closer than in the constructed patch; they lean in and wrap the guest.
- The triangle does distort (two chains pair up at 4.3-5.1 A, the third sits about 8 A away), but equally in
  placements where the reagent is outside (e.g. butylammonium_1: 5.0 / 8.4 / 10.4 A). It is not caused by insertion.
- All chains stay all-anti (0 gauche in 40 cases); the deformation is bending/leaning from the frozen C1-C2 tether.
- The inserted reagents are truly between the chains: 3-6 heavy atoms within 3 A of the axis below the chain tops
  (all 5 of DEA in two placements).

Conclusion: the central site exists by construction (7.6 A spacing, 4.4 A from each chain axis; the builder
docstring says "room for one inserted alkyl chain"), not because the chains open. Restraining or freezing the outer
chain carbons (option 1a of the 11-11 handoff) would not close it: the frozen seed is more open (free radius 2.78 A)
than the relaxed inserted geometries (2.28 A). 1a is dropped. What the 3-chain patch lacks is the cost of the site:
neighbouring chains and solvent that would compete for it on a real phase.
