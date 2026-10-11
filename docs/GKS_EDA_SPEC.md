# Original GKS-EDA specification: BLOCKED / NOT IMPLEMENTED

Authoritative derivation: Su, Jiang, Chen, Wu (2014), J. Phys. Chem. A 118,
2531–2542, https://doi.org/10.1021/jp500405s. The Su–Tang–Wu (2020) review
is https://doi.org/10.1002/wcms.1460. The handoff's supplied review PDF is not
present in this project; its former scratch attachment is unavailable locally.

The ACS accessible pages expose the abstract and an "Available to Purchase" notice,
not the complete derivation. No complete authorized local derivation or XEDA
executable/source has been found. Therefore this document is an explicit source gate,
not a completed mathematical specification or a validated CPU EDA implementation.

The inherited review describes interaction energy as complex energy minus fragment
energies, decomposed into electrostatic, exchange, repulsion, polarization, generalized
correlation and explicit dispersion terms. Polarization includes orbital relaxation/
charge transfer. The inherited generalized correlation expression mixes the DFT/HF
exchange difference and DFT correlation. These statements need equation-by-equation
verification in the original paper before executable component definitions are added.

Unresolved definitions required before implementation:
1. Fragment direct-product, antisymmetrized-product, and orthonormalized intermediate
   energies, their normalization, spin factors and one-/two-electron density conventions.
2. Exchange/repulsion separation and global hybrid exchange convention, including
   the exact B3LYP correlation variant and integration settings of each reference.
3. Exact reference geometry and ghost-basis/BSSE treatment for intermediate states.
4. Individual component references using identical orbitals, basis, grids, thresholds,
   and dispersion settings. Arithmetic closure alone is insufficient.

The official XEDA manual (https://xmvb.xmu.edu.cn/docs/xeda/content.html) provides
ammonia–water coordinates and rounded grouped components for wB97XD/aug-cc-pVTZ.
It combines exchange/repulsion and correlation/dispersion; this range-separated case
is outside the initial global-hybrid scope and lacks separate component precision.
This is a useful future comparison lead, not a fixture validating our implementation.
The manual also describes direct and RIJ-COSK evaluation, which cannot be presumed
equivalent to the optimizer's density fitting settings.

Next authorized scientific stage requires the full original derivation/supporting
material and reproducible XEDA/reference output for restricted closed-shell gas-phase
global-hybrid water dimer/ammonia–water/benzene dimer. Keep papers local unless
redistribution is authorized. No ALMO-EDA or DM-EDA substitution is made.

Future reference policy: use fragment coordinates from the complex, report deformation
separately; choose and verify a consistent counterpoise policy; exclude solvent,
open shells, range separation and double hybrids until their definitions are secured.
HPLC models and retention inference remain deferred under the handoff's sampling gates.
