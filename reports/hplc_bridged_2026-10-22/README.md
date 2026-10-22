# Bridged ternaries — residue_A···reagent⁺···C8 (2026-10-22)

Complete set (supersedes `reports/hplc_bridged_interim_2026-10-21/`). M06-2X/def2-SVP, C-PCM water, frozen
silane anchor; three-body GKS-EDA(sol). Table and discussion: `docs/VALIDATION.md`, section
"Bridged ternaries (application), 2026-10-22".

- `summary.json`: optimization summary for all five cases (copy of `runs/hplc_bridged_cpcm_def2svp/summary.json`).
- `<case>/optimized.xyz`, `<case>/optimization_result.json`: final DFT geometry and optimizer result.
- `<case>/three_body/three_body.json`: pair, (RT)|S, R|(TS) components and ΔE₃ (kcal/mol).
- `<case>/three_body/<job>/eda_result.json`: full GKS-EDA(sol) result per decomposition
  (R = residue_A, T = reagent, S = C8 proxy). Reference densities and logs are not archived (in `runs/`).
