# Experimental complete restricted/unrestricted GKS-EDA: CPU and GPU

All six gas-phase terms are now computed by the experimental `eda` command. This
supersedes the earlier two-component preparation stage. Independent algebra and
many-electron determinant checks support the implementation. External validation
against XEDA under explicitly identical settings remains unresolved; therefore
`scientifically_validated` and the validated-release flag `full_eda_available` remain
false, while `all_components_computed` is true. No solvent extension is claimed. The unrestricted extension is described below. The older `eda-prepare` command retains its CPU-only partial result format.

## Sources and provenance

The supplied 22-page review, Su, Tang & Wu (2020),
[10.1002/wcms.1460](https://doi.org/10.1002/wcms.1460), was read completely and
methodology pages 3–5 inspected visually. Local PDF SHA-256:
`9172c6f69ff059ea122b5f309fe287da9b4a740996f63e734b9f949bc80a7983`.
Its Section 2.2 applies the orbital intermediates of Eqs. 2–5 to KS orbitals; Eqs.
8, 9 and 11 specify the component sum, generalized correlation and dispersion.

The underlying Su & Li (2009) derivation,
[10.1063/1.3159673](https://doi.org/10.1063/1.3159673), supplies explicit orbital
energy formulas in Eqs. 6–12. Its
[institutional record](https://digitalcommons.unl.edu/chemfacpub/22/) establishes
provenance, but the linked download returned HTTP 403. The
[author-uploaded article text](https://www.researchgate.net/publication/26655483_Energy_Decomposition_Analysis_of_Covalent_Bonds_and_Intermolecular_Interactions)
was accessible and the relevant formulas read. In particular, Eq. 8 uses monomer
orbitals without the inverse interfragment overlap metric; Eq. 10 introduces that
metric. Those explicit approximate energy formulas resolve the review's AS/ASN
terminology sufficiently to implement its stated construction experimentally.
AS here denotes an energy formula, not a second normalized physical determinant.

The 2014 original [10.1021/jp500405s](https://doi.org/10.1021/jp500405s) remains
unavailable as an authenticated PDF. The university bibliography's original PDF
links returned 404. The
[public transcription](https://doczz.net/doc/8055485/energy-decomposition-scheme-based-on-the-generalized-kohn)
has interleaved equations and different step numbering; it is not used to replace
the explicit Su–Li formulas plus the supplied review. Confirmation against the
original PDF or reference implementation remains a scientific audit item.

The [official XEDA II tutorial](https://xmvb.xmu.edu.cn/docs/xeda/_downloads/415e0f39b016129f5cf2eb24dfe28b4f/XEDA_II_en.pdf)
supplies the CO₂ geometry (PDF page 24) and individual rounded results (page 39).
Local SHA-256:
`269131a8383e117af77fd32a806d22ec0c8cd5a08213477c237ba59a76d89952`.
Attributed small reference data are in `references/xeda_ii_co2_b3lyp.json`.
All full-text papers/manuals stay local and ignored by git.

## Scope

Exactly two all-electron fragments with additive spins; RKS for spin zero and UKS
for nonzero spin; real global-hybrid KS orbitals; B3LYP, the explicitly named
B3LYP5 variant, PBE0, or M06-2X; gas phase; no explicit dispersion; no density fitting;
double precision. Both spherical (default) and Cartesian (`--cart`) basis
representations are recorded. The `eda` CLI defaults to SCF tolerance 1e-12 Eh. Nonadditive spin coupling, ECPs, range separation, double hybrids,
nonlocal correlation, PCM and intramolecular extensions are excluded.

Fragment references retain the complex geometry and entire AO basis using ghost
nuclei. This is Boys–Bernardi counterpoise, with no fragment deformation or
separate geometry optimization. Indices preserve input ordering. Spin is Nalpha
minus Nbeta, and all differences are complex minus fragment sums. Energies are
hartree, input coordinates angstrom, internal integrals atomic units.

## Executable orbital construction

Let $C_A,C_B$ be occupied spatial KS coefficients in the same complex AO basis,
$S$ the AO overlap, and $P_M=2C_MC_M^T$. Each monomer satisfies
$C_M^TSC_M=I$, $\mathrm{Tr}(P_MS)=N_M$, and $P_MSP_M=2P_M$.
The simple sum $P_0=P_A+P_B$ generally is not idempotent.

Define $C=(C_A,C_B)$, $T=C^TSC$. The normalized antisymmetrized-product density is

\[
P_{\mathrm{ASN}}=2CT^{-1}C^T.
\]

Implementation diagonalizes the symmetric occupied metric, builds an orthonormal
occupied representation and forms its density. A minimum metric eigenvalue at or
below $10^{-8}$ is rejected: no pseudoinverse, electron dropping or regularization.
The inverse-metric expression and any orthonormal representation of this same
occupied span yield the same determinant. Density count, idempotency and symmetry
are checked. Condition number and smallest eigenvalue are recorded.

For a spin-summed restricted real density, define the HF-like functional on KS
orbitals, with full exact exchange rather than the SCF hybrid fraction:

\[
H[P]=\mathrm{Tr}(Ph)+\tfrac12\mathrm{Tr}(PJ[P])
-\tfrac14\mathrm{Tr}(PK[P])+E_{NN},
\]

where $h$ is the complex kinetic plus nuclear-attraction operator,
$J[P]_{\mu\nu}=\sum_{\lambda\sigma}(\mu\nu|\lambda\sigma)P_{\lambda\sigma}$,
$K[P]_{\mu\nu}=\sum_{\lambda\sigma}(\mu\lambda|\nu\sigma)P_{\lambda\sigma}$.
Monomer HF-like energies use their own real nuclear potentials and nuclear
repulsions, with the unchanged full AO basis. This is not a separate HF SCF.

The three approximate complex energy expressions are:

\[
E_{\mathrm{DP}}=\mathrm{Tr}(P_0h)+\tfrac12\mathrm{Tr}(P_0J[P_0])
+E_X^{\mathrm{HF}}[P_A]+E_X^{\mathrm{HF}}[P_B]+E_{NN},
\]
\[
E_{\mathrm{AS}}=H[P_0],\qquad E_{\mathrm{ASN}}=H[P_{\mathrm{ASN}}].
\]

These correspond to Su–Li Eqs. 6, 8 and 10. DP omits interfragment exchange; AS
adds it without replacing the monomer metric with $T^{-1}$; ASN uses the
normalized determinant. This avoids mistaking two equivalent normalized
determinants for different physical intermediate states.

## Individual terms

\[
\Delta E_{\mathrm{ele}}=E_{\mathrm{DP}}-H_A[P_A]-H_B[P_B],
\quad \Delta E_{\mathrm{ex}}=E_{\mathrm{AS}}-E_{\mathrm{DP}},
\]
\[
\Delta E_{\mathrm{rep}}=E_{\mathrm{ASN}}-E_{\mathrm{AS}},
\quad \Delta E_{\mathrm{pol}}=H[P_S]-E_{\mathrm{ASN}}.
\]

Electrostatic and exchange terms are independently contracted as

\[
\Delta E_{\mathrm{ele}}=\mathrm{Tr}(P_AV_B)+\mathrm{Tr}(P_BV_A)
+\mathrm{Tr}(P_AJ[P_B])+E_{NN}^{AB},
\qquad \Delta E_{\mathrm{ex}}=-\tfrac12\mathrm{Tr}(P_AK[P_B]).
\]

Their agreement with the intermediate differences and Coulomb reciprocity are
checked. The cross-HF exchange previously labeled diagnostic is assigned to the
experimental exchange term only after reading the explicit Su–Li Eq. 9. The old
preparation command retains its diagnostic label. Repulsion and polarization are
evaluated from their defined intermediates, never assigned by fitting a residual.
Polarization includes charge transfer and all orbital relaxation; no ALMO split.

Generalized correlation is review Eq. 9:

\[
G_{\mathrm{corr}}[P]=E_{\mathrm{KS}}[P]-H[P],\qquad
\Delta E_{\mathrm{corr}}=G_{\mathrm{corr}}[P_S]-G_{\mathrm{corr}}[P_A]-G_{\mathrm{corr}}[P_B].
\]

For a simple global hybrid this is
$(1-\gamma)(\Delta E_X^{\mathrm{DFT}}-\Delta E_X^{\mathrm{HF}})+\Delta E_C$.
Each result records its installed LibXC functional definition. Explicit dispersion
(review Eq. 11) is zero for the selected methods; no D3/D4 extension is implied.

All six terms must sum to the counterpoise KS interaction energy within $10^{-8}$
Eh. Closure is an internal consistency check, never external component validation.

## Backend and validation policy

CPU uses PySCF RKS and J/K; GPU uses GPU4PySCF RKS/direct J/K plus CuPy double
precision for occupied metrics, orthonormalization, densities and contractions.
Molecule construction, basis verification, reporting and archive output run on CPU.
GPU4PySCF itself describes direct J/K as a CPU/GPU hybrid algorithm; there is no
all-GPU or speedup claim. Requested unavailable GPUs fail without CPU substitution.
Host transfers occur for scalar checks and final archives; full intermediate matrices
stay on the selected backend. Application code does not build a full four-index ERI
array; tiny independent tests use one, while backend SCF may manage integral caches.

Tests independently evaluate each component with explicit ERI contractions and QR
orthonormalization, numerically integrate separate PBE exchange/correlation, and
compare ASN energy to the full many-electron Hamiltonian acting on a Slater expansion.
They also check occupied-span invariance, fragment order, singular-span rejection,
missing-GPU behavior and honest experimental output. GPU tests compare fixed CPU
densities first to isolate integral/matrix errors, then independent GPU SCFs.

Published CO₂ comparisons cover spherical/Cartesian bases and explicitly named
B3LYP variants; grid convergence is checked separately. Basis representation explains
much of the earlier discrepancy. Unknown XEDA quadrature, integral approximation
(default LR), CP policy and functional variant still prevent a matched-settings
reference-validation claim. The next release gate is explicitly specified XEDA output
for several individual-component cases. Solvent, nonadditive spin coupling and HPLC follow later.


## Unrestricted and M06-2X extension

For each spin σ, Dσ=Cσ Cσᵀ. With D=Dα+Dβ, the HF-like functional is
Tr(Dh)+½Tr(DJ[D])−½ΣσTr(DσK[Dσ])+Enn. Cross exchange between fragments is
−ΣσTr(DaσK[Dbσ]). Each spin's normalized product density is
Cσ(CσᵀSCσ)⁻¹Cσᵀ; occupied spans are independently checked. A closed-shell
fragment contributes half its spin-summed density to each channel. All other
component definitions and the independently verified closure are unchanged.

M06-2X is a global hybrid meta-GGA (54% SCF exact exchange); EDA still uses full
HF-like exchange at KS orbitals, not 54% exchange. Generalized correlation remains
EKS−EHF-like. Tests check open-shell four-index energy and cross exchange, analytic
meta-GGA gradients against energy differences, and actual CPU/GPU component and
gradient agreement. UKS spin contamination is measured using S², not removed.
Only additive parallel fragment-spin coupling is supported. The old `eda-prepare`
command retains its original closed-shell B3LYP/PBE0 restrictions.
