# Original GKS-EDA: source-backed CPU preparation; full decomposition gated

Status as of 2026-10-11: experimental CPU reference preparation and two partial
components are implemented. Full original GKS-EDA, including exchange, repulsion,
and polarization, is unavailable. No GPU EDA or component-reference validation is claimed.

## Sources and provenance

The owner supplied the complete 22-page review: Su, Tang, Wu (2020),
[doi:10.1002/wcms.1460](https://doi.org/10.1002/wcms.1460).
SHA-256: `9172c6f69ff059ea122b5f309fe287da9b4a740996f63e734b9f949bc80a7983`.
All pages were read; methodology pages 3–5 were also visually inspected to check signs
and equation notation. The PDF stays local and is ignored by git.

The original derivation is Su, Jiang, Chen, Wu (2014),
[doi:10.1021/jp500405s](https://doi.org/10.1021/jp500405s).
The publisher exposes an abstract/paywall. The institution's 2014 publication list
links to `http://pcss.xmu.edu.cn/download/paper2014/JPhysChemA-118-2531.pdf`, which
returned 404 (HTTPS and the current pcoss host also returned 404).
A [public transcription](https://doczz.net/doc/8055485/energy-decomposition-scheme-based-on-the-generalized-kohn)
contains the original title/authors/DOI and methodology; its equation layout is interleaved
and its downloadable original PDF has not been secured. It is a lead for resolving
conventions, not an authenticated replacement PDF or implementation authority alone.

The [official XEDA II tutorial](https://xmvb.xmu.edu.cn/docs/xeda/_downloads/415e0f39b016129f5cf2eb24dfe28b4f/XEDA_II_en.pdf)
was downloaded locally. SHA-256:
`269131a8383e117af77fd32a806d22ec0c8cd5a08213477c237ba59a76d89952`.
Sections 1.3.5/3.1.5 provide a closed-shell CO2 dimer, B3LYP/6-31G*, and individually
reported rounded components. This improves the previous reference-data gate, but several
settings are unspecified. Small coordinates/numbers with source attribution are in
`examples/co2_dimer_xeda_ii.*` and `references/xeda_ii_co2_b3lyp.json`; full manuals
and papers are excluded from git.

## Scope and conventions

Initial preparation: exactly two restricted closed-shell fragments and a restricted
closed-shell complex; real all-electron B3LYP or PBE0 orbitals, gas phase, no explicit
dispersion, direct two-electron integrals, double precision. The complex geometry is
fixed and fragment coordinates inherit original atom ordering. No deformation energies
are included. Counterpoise is explicit: every fragment uses the whole complex AO basis,
with zero-charge ghost nuclei on the other fragment. No own-basis/full-basis reference
cycles are mixed. Integral results are contracted directly; full four-index arrays are
used only in tiny independent STO-3G tests.

All differences use complex minus sum of fragment quantities. Energies are hartree;
coordinates are read in angstrom and integrals use atomic units. Spin is Nalpha−Nbeta
and must be zero for every state. Density matrices are spin summed:

\[
P_M=2C_{M,\mathrm{occ}}C_{M,\mathrm{occ}}^T,
\quad \mathrm{Tr}(P_MS)=N_M,
\quad P_MSP_M=2P_M.
\]

The factor two is occupancy, not an arbitrary spin factor. Each fragment occupies
orbitals in the common basis; fragment orbital spaces may overlap. Their simple summed
density is not assumed to be an idempotent supermolecule density.

## Target definitions from the supplied review

Review Eq.8 defines

\[
\Delta E_{\mathrm{tot}}=E_S-\sum_M E_M
=\Delta E_{\mathrm{ele}}+\Delta E_{\mathrm{ex}}+\Delta E_{\mathrm{rep}}
+\Delta E_{\mathrm{pol}}+\Delta E_{\mathrm{corr}}+\Delta E_{\mathrm{disp}}.
\]

Review Eqs.2–5 describe direct-product (DP), antisymmetrized (AS), and
antisymmetrized/orthonormalized (ASN) intermediates; Section 2.2 replaces HF orbitals
with KS orbitals. Polarization includes charge transfer and orbital relaxation.
The executable specification below covers the unambiguous preparation quantities only.

### Frozen electrostatic term: implemented, not reference-validated

Review Eq.2 specifies the classical Coulomb interaction of frozen monomer densities.
With $V_A,V_B$ the nuclear-attraction matrices of real nuclei in A/B,
$J(P)_{\mu\nu}=\sum_{\lambda\sigma}(\mu\nu|\lambda\sigma)P_{\lambda\sigma}$,
and common-basis spin-summed densities,

\[
\Delta E_{\mathrm{ele}}=
\mathrm{Tr}(P_AV_B)+\mathrm{Tr}(P_BV_A)
+\mathrm{Tr}(P_AJ(P_B))
+\sum_{a\in A,b\in B}\frac{Z_aZ_b}{R_{ab}}.
\]

No cross-exchange appears in this classical quantity. Each nucleus contributes once,
and the electron–electron term has no extra factor 1/2 because it is the cross term.
Independent tests use explicit electron-repulsion integrals and fragment nuclear
potentials rather than repeating the production direct-J evaluation. Reciprocity under
A/B exchange is tested. Numerical correctness of these operations does not establish
agreement with an external original GKS-EDA implementation.

### Generalized correlation term: implemented, not reference-validated

For a spin-summed restricted density, define the HF-like energy at **KS orbitals**:

\[
E_{\mathrm{HF-like}}[P]=\mathrm{Tr}(Ph)+\tfrac12\mathrm{Tr}(PJ(P))
-\tfrac14\mathrm{Tr}(PK(P))+E_{NN},
\quad K(P)_{\mu\nu}=\sum_{\lambda\sigma}(\mu\lambda|\nu\sigma)P_{\lambda\sigma}.
\]

This is an energy evaluation, not a separate optimized HF calculation. Review Eq.9
is evaluated without guessing a composite functional's separately normalized exchange:

\[
G_{\mathrm{corr}}[P]=E_{\mathrm{KS}}[P]-E_{\mathrm{HF-like}}[P],
\quad\Delta E_{\mathrm{corr}}=G_{\mathrm{corr}}[P_S]-\sum_M G_{\mathrm{corr}}[P_M].
\]

For a simple global hybrid this equals
$(1-\gamma)(\Delta E_X^{\mathrm{DFT}}-\Delta E_X^{\mathrm{HF}})+\Delta E_C$.
For PBE0 the code is independently tested against separate PBE exchange and correlation
quadrature with gamma=1/4. B3LYP uses the installed LibXC definition, which is recorded
in each result; the tutorial's B3LYP variant is unknown. No MP2/double-hybrid correction,
range separation, or nonlocal correlation is included.

Review Eq.11 gives the explicit dispersion difference. It is exactly zero for the
selected non-dispersion-corrected methods; that does not mean correlation cannot
contribute to dispersion-like attraction.

### Exchange, repulsion, polarization: unresolved and unavailable

A diagnostic bare interfragment HF exchange can be computed as
$-\tfrac12\mathrm{Tr}(P_AK(P_B))$. It is deliberately **not labeled** the original
GKS-EDA exchange component. The original transcription Eq.22 describes a difference
of exact exchange between steps 1 and 2; review Eq.3 describes AS minus DP before its
ASN step. The orbital/normalization convention must be verified against original
page 2533 and its cited LMO-EDA/Stone construction before choosing an executable split.
A normalized determinant and an orthonormalized representation of the same occupied
span cannot simply be treated as different physical states by fiat.

No repulsion/polarization equations are fabricated by subtracting partial components
from the total. Their JSON values, and full closure residual, remain null.
A completed `eda-prepare` means reference preparation succeeded, not full EDA success.

## Validation gate and next work

Independent integral and PBE0-functional tests establish the implemented algebra.
The XEDA II CO2 numbers are compared explicitly as exploratory data, never as a
passing matched-settings test: the tutorial leaves B3LYP variant, grids/pruning, SCF
thresholds, exact program build and explicit CP setting unspecified; its integral
algorithm defaults to LR rather than this direct-integral calculation.

Next: secure original page 2533 in a reliable PDF/source, resolve intermediate exchange
normalization, obtain an XEDA output with explicit settings/orbitals, then implement
and validate exchange, repulsion and polarization component by component. GPU EDA
porting follows a validated CPU reference. Solvent (review Eqs.12–13), open shells,
double hybrids, intramolecular fragmentation and HPLC interpretation remain deferred.
