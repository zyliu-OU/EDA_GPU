# Reagent retention fit on the C8 patch (2026-11-11)

`scripts/fit_patch_retention.py --cases reports/hplc_patch_extended_2026-11-11/cases.json --data
examples/hplc/donegan2022_table4.csv --output fit.json` (stdout in `fit.log`). Target: % ACN at elution of dT20 and
dT60, 100 mM alkylammonium acetate (Donegan 2022 Table 4), 7 reagents with values; TBAA (not eluted) only as most
retained in rank tests. Baselines: logP, logD (ionized), carbon count (Table 1). Descriptors per reagent from 5
placements: Boltzmann and best (EDA, EDA+CDS, buried npSASA) or lowest-energy value (single terms), over all
placements or `off_centre` (drops r < 2.5 Å and h < 8 Å, i.e. inserted in the patch's central hollow: 4/5 DEA, 1 BA,
2 DBA). LOO RMSE of linear fits; permutation check (2000) that repeats the whole descriptor search.

dT20 (group-mean LOO RMSE 11.2 % ACN):

| Descriptor | r | LOO RMSE | ρ (n 7) | ρ (+TBA, n 8) |
| --- | ---: | ---: | ---: | ---: |
| logD (baseline) | 0.77 | 7.5 | 0.61 | 0.74 |
| EDA total, Boltzmann, off-centre | 0.74 | 8.1 | 0.75 | 0.83 |
| EDA + CDS, Boltzmann, off-centre | 0.75 | 8.1 | 0.82 | 0.88 |
| electrostatic, Boltzmann, all | 0.65 | 9.0 | 0.71 | 0.81 |
| EDA total, Boltzmann, all | 0.43 | 13.1 | 0.25 | 0.45 |

- No EDA descriptor beats logD; the best single fit of the whole search is logD itself (p_select 0.33, i.e. a
  shuffled target finds an equally good fit a third of the time). logD alone: p 0.05.
- Dropping centre-pocket placements turns the EDA total from no fit (r 0.43) into one about as good as logD (r 0.74,
  ρ 0.75–0.83). The DEA outlier is a pocket effect of the open 3-chain patch; but off-centre DEA rests on 1 placement.
- logD + off-centre polarization: LOO 3.7 (p_select 0.025) but driven by HAA alone (least negative pol, largest
  logD residual +11.9); without HAA it improves logD only 5.5 → 4.8, and the sign (weaker pol = more retained) has
  no physical reading. Treated as spurious.
- Largest logD residuals: HAA +11.9, DIPEAA −8.8 (primary/linear amines retained more, branched less, than logD
  says). EDA reproduces neither: off-centre EDA has HAA at −5.6 vs DBAA −7.5 (if any, the shape effect points to
  chain packing that a 3-chain patch cannot represent).
- dT60 gives the same picture (logD LOO 9.0 vs mean 14.7; off-centre EDA+CDS 10.2).
