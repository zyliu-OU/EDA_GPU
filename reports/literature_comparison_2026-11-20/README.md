# Published base-composition data vs the per-base EDA (2026-11-20)

Published data only. The papers and SI are local copies in the repo root and are not committed (`*.pdf`, `*.zip`
are ignored).

## Data
- **Sturm et al. 2007** (*Nucleic Acids Res.* 35:4195). Supplementary Table A.1 is transcribed to
  `examples/hplc/sturm2007_tableS1.csv`: 72 DNA oligonucleotides, 15–48 nt. Base counts were checked against every
  sequence. Conditions: 100 mM TEAA pH 7, PS-DVB monolith (no silanols), normalized RT at 30–80 °C.
- **Gilar et al. 2002** (*J. Chromatogr. A* 958:167), Eqs. 2–5. % ACN at elution of homo-oligonucleotides, TEAA,
  XTerra C18, 60 °C: A 2.109 ln N + 1.47; T 2.125 ln N + 3.77; C 1.887 ln N + 0.71; G 1.887 ln N + 0.89. G was not
  measured: it is the C fit plus 0.18 % from a CCG-motif ladder. The text gives the hydrophobicity order C < G < A < T.

## Fit to Sturm (`scripts/fit_base_composition.py` → `sturm_composition.json`, `.log`)
At 80 °C, where secondary structure is minimal, the Gilar-type model fits well (R² 0.988, LOO RMSE 0.18 min). The
additive contribution of one base in a 24-mer, relative to C, is:

| | G | C | A | T |
| --- | ---: | ---: | ---: | ---: |
| ΔRT per base (min), 80 °C | −0.02 | 0 | +0.22 | +0.36 |
| Mean shift of single-base exchanges, 80 °C (n pairs) | C→G 0.00 (43) | — | C→A +0.22 (89), G→A +0.26 (33) | G→T +0.40 (18) |

Order: **T > A > C ≈ G**, in agreement with Gilar (T > A > G ≈ C). Below 60 °C the fits degrade (R² 0.64–0.80) because
hairpins and self-complementary sequences dominate, so per-base values there are not meaningful.

## Against the per-base EDA (RNA residues on the 3-chain C8 patch, Boltzmann, kcal/mol)
A −11.52 < G −10.61 < U −9.64 < C −8.25 (more negative = more strongly bound).

| Comparison | Literature | EDA | Agreement |
| --- | --- | --- | --- |
| A vs C | A more retained (Sturm, Gilar, Maurer NTP) | A 3.3 kcal/mol stronger | yes |
| A vs G | A more retained (Sturm, Gilar, Maurer NTP) | A 0.9 stronger | yes |
| G vs C | Oligos: equal (Sturm 0.00 min, Gilar assumed +0.18 %); NTPs: G more retained | G 2.4 stronger | partly: matches the NTPs, too large for oligos |
| T/U | T most retained (DNA, has 5-methyl) | U (no methyl) between C and G | not comparable; the model has no T |

## Caveats
- DNA oligonucleotides vs the RNA residue model; PS-DVB/C18 vs the C8 patch; TEAA is present in experiment but not in
  the base model.
- The EDA order depends on bases inserted into the patch's central hollow (`reports/hplc_patch_bases_2026-11-20/`).
- A deoxythymidine residue was not computed, so the strongest literature trend (T > A) is untested.
