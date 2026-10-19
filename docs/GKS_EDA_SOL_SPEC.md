# GKS-EDA(sol) specification (experimental)

Status: implemented on CPU and GPU for two closed-shell or additive-spin fragments,
C-PCM or IEF-PCM, the same functionals and all-electron bases as gas-phase GKS-EDA.
**Not externally validated.** Definitions below follow the supplied review; points the
review leaves open are marked *interpretation* and must be confirmed against the original
derivation (Su et al. 2014, 10.1021/jp500405s; EDA-PCM, Su, Liu & Wu,
J. Chem. Phys. 2012, 137, 034111, as cited by the XEDA GAMESS manual) or matched XEDA
output before any validated claim.

## Source equations

Su, Tang & Wu, WIREs Comput. Mol. Sci. 2020, 10, e1460, Section 2.3:

- Eq. 12: ΔG_TOT = G_S − Σ_M G_M = ΔG_ele + ΔG_ex + ΔG_rep + ΔG_pol + ΔG_corr + ΔG_disp + ΔG_desol.
- "the first six terms are the same as the corresponding terms in [the gas-phase
  decomposition] except that the orbitals used in these terms are optimized by SCRF".
- Eq. 13: ΔG_desol = ⟨Φ_S|v_SOL|Φ_S⟩ − Σ_M ⟨Φ_M|v_SOL|Φ_M⟩.
- The supermolecule cavity is constructed according to the interacting distance; at short
  distances "monomers are located in the same cavity", because monomer solvation computed
  in each monomer's own cavity is "too negative".
- Procedure: create cavities; assign dielectric constant(s); obtain v_SOL; perform the
  modified GKS-EDA with the SCRF wave functions.

## Implemented definitions

For every state X ∈ {S, M₁, M₂} an SCRF Kohn–Sham calculation in the common
(counterpoise) complex basis gives the density D_X and the free energy

  G_X = E_gas[D_X] + E_solv[D_X],  E_solv[D] = ½ qᵀ V[D],

where V is the total (electronic + nuclear) solute potential on the cavity surface and q the
apparent surface charges of the chosen PCM (PySCF/GPU4PySCF conventions). E_gas[D_X] is
evaluated independently by a vacuum calculation on D_X and must equal G_X − E_solv within
1e-8 Eh, otherwise the run fails.

1. **Six gas-phase-form terms**: the existing GKS-EDA expressions (Su–Li 2009 intermediate
   states with KS orbitals; review Eqs. 2–5, 8–9, 11) evaluated with the SCRF orbitals and
   E_gas in place of the KS energy. Their sum equals E_gas[D_S] − Σ E_gas[D_M]
   (reported as `gas_phase_functional_interaction_hartree`).
2. **Desolvation** (*interpretation of Eq. 13*): ⟨Φ|v_SOL|Φ⟩ is taken as the PCM
   electrostatic solvation free energy E_solv[D_X] (including the nuclear contribution and the
   ½ factor), so ΔG_desol = E_solv[D_S] − Σ E_solv[D_M]. This is the only reading under which
   Eq. 12 holds with the six terms defined as stated.
3. **Monomer cavity** (*interpretation of Section 2.3*): each monomer SCRF uses the full
   supermolecule cavity (spheres on all complex atoms) with only its own nuclei generating the
   surface nuclear potential, i.e. a cavity analogue of counterpoise. The library ghost-atom
   cavity handling is bypassed because it differs between PySCF (element-0 sphere) and
   GPU4PySCF (ghosts dropped); the reconstructed nuclear potential is checked against the
   library value for the complex charges to 1e-10.
4. **Total**: ΔG_TOT = G_S − Σ G_M (reported as `counterpoise_interaction_energy_hartree` in
   solvent runs). Seven-term closure is an identity and is checked to 1e-8 Eh; it is not
   evidence of correctness.

Not included: nonelectrostatic solvation (cavitation, dispersion–repulsion), heterogeneous
dielectric assignment, non-equilibrium solvation, explicit dispersion, and the monomer-own-cavity
variant.

## Numerical settings recorded with every result

PCM model (`C-PCM` with f(ε) = (ε−1)/ε, or `IEF-PCM`), ε (default 78.3553, water, as in the
XEDA solvent table), Lebedev order (default 29, 302 points per sphere), SWIG discretization,
radii = modified Bondi × 1.2 without probe (PySCF/GPU4PySCF default table). These differ
from GAMESS (GEPOL tesserae) and XEDA II defaults, so cross-program agreement is limited by
cavity construction even when the decomposition is identical.

## Validation design (`scripts/validate_gks_sol.py`)

1. Vacuum limit ε → 1 must reproduce every gas-phase component (1e-6 Eh) with ΔG_desol → 0.
2. Separation limit: ΔG_desol and ΔG_TOT must decay monotonically toward zero.
3. CPU/GPU component agreement for C-PCM and IEF-PCM (1e-7 Eh).
4. PCM analytic gradients against central finite differences (2e-5 Eh/bohr).
5. Exploratory comparison with the XEDA II manual BH₃–NH₃ C-PCM example (settings unknown).

External matched validation still requires XEDA output with explicitly stated cavity radii,
discretization and integral settings, plus confirmation of interpretations 2–3 from the
original derivation.
