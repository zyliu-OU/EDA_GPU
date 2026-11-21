# EDA_GPU handoff — 2026-11-20

Previous: `reports/EDA_GPU_HANDOFF_2026-11-18.md`.

## State
- The bases queue (`runs/queue_hplc_patch_bases.full.sh`) finished at 2026-11-20 08:19 PST ("queue done" in
  `runs/hplc_patch_bases.log`). It took about 33 h, against the ~2 days estimated.
- Results are in `runs/hplc_patch_bases_cpcm_def2svp/`: 20 cases (5 placements × A/C/G/U), all `eda: success`. Bonding
  topology was preserved in all of them and max |closure residual| is 5.7e-10 Eh. `nonpolar.json` (SMD CDS
  water/acetonitrile, buried SASA) and `analysis.json` were written.
- Nothing is running. The GPU is idle (44 °C). The systemd unit has exited.

## Per-base results (lowest optimized PCM energy placement, GKS-EDA(sol) total, kcal/mol)
| Base | Best placement | EDA total (best) | EDA range over 5 placements |
| --- | --- | ---: | --- |
| A | A_1 | −11.58 | −11.58 … −6.23 |
| G | G_2 (G_0 +0.06) | −10.61 | −10.65 … −8.65 |
| U | U_4 (U_1 +0.24) | −10.83 | −10.83 … −8.12 |
| C | C_1 (C_4 +0.18) | −7.50 | −9.34 … −7.02 |

The grouping is now fixed: `scripts/analyze_hplc.py` groups single-index placement names per guest, and
`analysis.json` was regenerated. Near-degenerate minima (G, U, C within 0.25 kcal/mol) make the best-placement
values sensitive, so a Boltzmann or range-based comparison is more robust.

The log has ~9.5k "Basis not found for atom N H" lines from the nonpolar (SMD CDS) step. The earlier extended run
shows the same warnings, and the CDS values are populated.

## Base retention comparison (done 2026-11-20)
`reports/hplc_patch_bases_2026-11-20/` (new script `scripts/compare_base_retention.py`, data
`examples/hplc/maurer2025_fig9_ntp.csv` = Maurer 2025 Fig. 9 NTP times read off the figure: C 1.50 < U 1.87 < G 1.99 <
A 2.82 min). The Boltzmann EDA total gives C −8.25 > U −9.64 > G −10.61 > A −11.52, so ρ = −1 (p = 1/24). Dropping
one placement keeps ρ = −1 in 19/20 cases. The bootstrap gives the exact order in 49 % of resamples. With CDS, U and G
swap (ρ −0.8). The base inversion seen with the silanol proxy is gone on the patch. Purines put more ring atoms on the
chains.

Hollow check (also 2026-11-20, `check_base_hollow.py`): 15/20 base placements have the ring inside the patch's central
hollow (r < 2.5 Å, h < 8 Å), and those carry nearly all the Boltzmann weight. The 5 off-centre placements alone give
ρ +0.4. So the order depends on the hollow, the same artefact as DEA⁺. The VALIDATION.md section and the README were
updated with this caveat.

## Next
Do not start another run without asking (memory `no-new-runs-after-bases`). Open options: off-centre re-sampling of the bases
(tests the hollow dependence) and of DEA/BA/DBA. The 4-chain test is dropped (memory `hardware-size-limit`).
