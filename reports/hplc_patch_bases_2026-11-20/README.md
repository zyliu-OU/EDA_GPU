# Capped RNA residues on the C8 patch vs base retention (2026-11-20)

Campaign: `runs/hplc_patch_bases_cpcm_def2svp/`. Capped residues A/C/G/U (5′-O-methylphosphodiester, 3′-O-methyl,
charge −1) on the silanol-free three-chain C8 patch. 5 placements per base: seeds 0/1 plus the best 3 of an xTB screen
of 2–9. M06-2X/def2-SVP C-PCM optimisation, then GKS-EDA(sol). All 20 succeeded with topology preserved and
|closure| ≤ 5.7e-10 Eh.

Commands:
```
.env/bin/python scripts/analyze_hplc.py runs/hplc_patch_bases_cpcm_def2svp --output runs/hplc_patch_bases_cpcm_def2svp/analysis.json
.env/bin/python scripts/compare_base_retention.py runs/hplc_patch_bases_cpcm_def2svp --data examples/hplc/maurer2025_fig9_ntp.csv \
    --seeds examples/hplc/patch_bases --output reports/hplc_patch_bases_2026-11-20/base_retention.json
```
The output is in `base_retention.json` and `base_retention.log`. The robustness checks in the log were run inline. The
same session fixed the grouping in `analyze_hplc.py`. Placement names with one index (`patch_<guest>_<k>`) are now
grouped per guest. Before the fix, each placement of a patch campaign was its own group. Groupings for the other
campaigns are unchanged (checked on all 12 `runs/hplc_*` summaries).

## Retention data
Maurer/Guillarme 2025 (J Chromatogr A 1740:465574) Fig. 9 shows NTPs on a BEH C18 column with 100 mM BA + 50 mM TPA at
1 % ACN and 65 °C. Peak apex times, read off the figure: CTP 1.50, UTP 1.87, GTP 1.99, ATP 2.82 min, so
**C < U < G < A**. The same paper's text gives the nucleotide hydrophobicity order T > A > G > C, which agrees on
A > G > C (T is not modelled).

## Result (kcal/mol, more negative = more strongly bound)
| Descriptor | C | U | G | A | ρ vs t_R | Order weak → strong |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| EDA total, Boltzmann (298 K) | −8.25 | −9.64 | −10.61 | −11.52 | −1.00 | C U G A |
| EDA total, lowest-energy placement | −7.50 | −10.83 | −10.61 | −11.58 | −0.80 | C G U A |
| EDA + CDS(water), Boltzmann | −12.10 | −14.35 | −13.84 | −15.87 | −0.80 | C G U A |
| generalized correlation, Boltzmann | −29.51 | −32.34 | −34.79 | −32.82 | −0.80 | C U A G |
| electrostatic, Boltzmann | −13.24 | −14.62 | −16.71 | −14.04 | −0.40 | C A U G |

- The Boltzmann EDA total reproduces the measured order exactly (exact permutation p = 1/24 = 0.04, one-sided).
  Unlike the silanol proxy (C strongest, −18.9), C is now the weakest binder and A the strongest. The inversion came
  from the silanol H-bonds.
- **Robustness.** Dropping any one of the 20 placements keeps ρ = −1 in 19/20 cases (worst −0.8). A bootstrap of
  placements within each base gives the exact order in 49 % of resamples. With CDS added, U and G swap, giving
  ρ −0.8 (exact order in 22 % of bootstrap resamples). The U/G gap (1.0 kcal/mol) and the G/A gap (0.9 kcal/mol) are
  smaller than the spread between placements (4–5 kcal/mol), so A > G > U > C is the most likely order, not a firm
  result.
- **Selection.** Of 40 energy descriptors looked at, at least one would reach ρ = −1 for 67 % of shuffled orders. For
  only the four headline descriptors (EDA total and EDA + CDS, each as best placement and Boltzmann) the figure is 8 %.
  The EDA-total Boltzmann was the headline descriptor in the earlier reagent fits, so it was not chosen after seeing
  these data. Even so, n = 4 cannot establish significance.
- **Mechanism.** No single component gives the order. A and G put more base-ring atoms in contact with the chains
  (55–64 contacts vs 34–40 for C/U) and fewer backbone atoms. Purines lie flatter on the alkyl chains and gain
  correlation (dispersion-like) energy, and A pays the least desolvation and Pauli repulsion. This is the expected
  purine > pyrimidine hydrophobicity pattern.

## Caveats
- The model is not the experiment. The NTPs carry a −4 triphosphate and are ion-paired with BA/TPA on C18; the model
  is a −1 capped monophosphate residue on a bare three-chain C8 patch with no reagent. The comparison tests only
  whether base–chain contact alone can rank the bases.
- The retention times were read off a figure (±0.03 min). U and G are 0.12 min apart and not baseline-resolved on
  shorter columns.
- U stands in for RNA; T (the DNA base in "T > A > G > C") has an extra methyl that the model does not include.
