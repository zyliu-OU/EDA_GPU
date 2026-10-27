# SMD CDS non-electrostatic association terms (2026-10-26)

Post-processing of the converged HPLC complexes (stages 7.1–7.3, reagent series, bridged ternaries); no new
DFT. `scripts/nonpolar_hplc.py` (module `eda_gpu.cds`): dG_CDS = G_CDS(AB) − ΣG_CDS(fragments) at the
complex geometry, SMD water and acetonitrile, beside the GKS-EDA(sol) total and buried SASA. Table and
discussion: `docs/VALIDATION.md`, section "SMD CDS term (application), 2026-10-26".

- `nonpolar.json`: per case and split, dG_cds_water, dG_cds_acetonitrile, buried total/nonpolar SASA,
  EDA total and EDA + dG_cds_water where an EDA exists; `three_body` = RT|S − T|S − R|S.
- `nonpolar.log`: one line per case.

Command:
`.env/bin/python scripts/nonpolar_hplc.py runs/hplc_{ionpair,surface,reagents_ionpair,reagents_surface,ternary,bridged}_cpcm_def2svp --three-body runs/hplc_three_body --three-body runs/hplc_three_body_bridged --output reports/hplc_nonpolar_2026-10-26/nonpolar.json`
