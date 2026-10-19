# Ion-pair RP-HPLC model study (gate 7): models, reference states, limits

Status: stage 7.1 (residue···TEAH⁺ ion pairs) and 7.2 (pairs with a stationary-phase proxy)
in progress; 7.3 (ternary complexes and three-body cycle) planned. Nothing here predicts
retention. EDA quantities are electronic interaction descriptors of labelled small models.

## Chemical models (`scripts/build_hplc_models.py`, `examples/hplc/`)

| Model | Definition | Charge | Rationale |
| --- | --- | ---: | --- |
| Capped residue A/G/C/U | β-D-ribonucleoside, 5′-O-methyl phosphodiester, 3′-O-methyl | −1 | One phosphodiester per residue as in RNA; identical termini for every base; neutral canonical bases |
| TEAH⁺ | protonated triethylamine | +1 | Ion-pair reagent at typical TEAA buffer pH (pKa ≈ 10.7) |
| C8 silane proxy | octyldimethylsilyl–O–Si(OH)₃, all-anti chain | 0 | Short stand-in for a bonded C18 ligand with residual silanols; anchor Si/O atoms frozen |
| C4 silane proxy | butyldimethylsilyl–O–Si(OH)₃ | 0 | Smaller fallback proxy |

Residue conformers are A-form-like by construction: anti χ, C3′-endo-like δ (55–105°), trans β,
gauche⁺ γ, and no intramolecular phosphate O···H–N/O contact below 2.5 Å (gas-phase force
fields otherwise fold the anion onto base or sugar donors). Seeds are RDKit ETKDG/MMFF94
geometries; every complex is re-optimized quantum mechanically before analysis.

Not modelled: neutral TEA, acetate counterions, multiple TEAH⁺ (no concentration dependence),
residual silanolate (SiO⁻), ligand packing or multiple tethered chains, endcapping, explicit
water/acetonitrile, finite-temperature sampling.

## Complexes and method

Seeds (`scripts/build_hplc_complexes.py`) place a donor X–H on a host acceptor at H···A 1.8 Å,
sampling approach cones and rotations with motif-aware clash rejection:

- 7.1 residue···TEAH⁺: TEAH⁺ N–H → each non-bridging phosphate O (2 seeds per base).
- 7.2 TEAH⁺···proxy: TEAH⁺ N–H → silanol O (2 seeds); residue···proxy: silanol O–H → phosphate O
  (1 seed per base).

`scripts/run_hplc_stage.py --preopt xtb`: each seed is first refined with GFN2-xTB/ALPB(water)
(xtb 6.7.1, separate conda environment `.xtb-env`; anchor internal geometry constrained), and the
refined structure is rejected if any metadata bond breaks or a hydrogen changes partner (e.g.
proton transfer). This only shortens the DFT optimization (the first unrefined ion-pair run
needed ~60 s per DFT step); all reported structures are DFT-optimized. Then GPU M06-2X/def2-SVP
optimization (def2-universal-jkfit density
fitting, grid 3, C-PCM water ε = 78.3553, frozen anchors), then direct-integral counterpoise
GKS-EDA(sol) at the optimized geometry with the same functional, basis and solvent. The same
topology gate is applied to every DFT-optimized structure; EDA is skipped if bonding changed,
since the fragment charges would no longer hold. Planned checks on a subset: aug-cc-pVDZ EDA,
cavity radii (modified Bondi × 1.2 vs 1.4), and dielectric (acetonitrile ε = 35.688, a
high-organic mobile-phase bound).

## Reference states and cycle (stage 7.3)

All interaction energies use complex geometries and counterpoise fragments (no deformation
term); fragment deformation energies will be reported separately. For residue R, reagent T⁺ and
proxy S the planned cycle compares, at a common ternary geometry:

- ΔG(R···T⁺), ΔG(T⁺···S), ΔG(R···S): pair terms;
- ΔG([R T⁺]···S): adsorption of the pre-associated ion pair;
- ΔG(R···[T⁺ S]): adsorption on a TEAH⁺-modified surface;
- three-body nonadditivity ΔG₃ = ΔG(RTS) − ΔG(RT) − ΔG(TS) − ΔG(RS).

The two adsorption pictures are different reference states and are never mixed in one comparison.

## Interpretation limits

Single optimized placements are local stationary points (no Hessian classification for these
sizes), not ensembles. PCM electrostatics omit cavitation and dispersion–repulsion solvation
terms; desolvation depends strongly on cavity radii (see the BH₃–NH₃ radius scan in
`docs/VALIDATION.md`). Base-dependent EDA differences are hypotheses for comparison with
controlled retention data on same-length, same-termini oligonucleotides, not retention factors.
