# EDA_GPU

**Status: under active development — an early research prototype.**

EDA_GPU is a GPU-assisted computational chemistry project for reproducible molecular
structure preparation, electronic-energy calculations, analytic gradients, and geometry
optimization. Its long-term goal is a faithful implementation of the original generalized
Kohn–Sham energy decomposition analysis (GKS-EDA), beginning with a validated CPU
reference and then accelerating the expensive operations on NVIDIA GPUs.

The working foundation uses **PySCF** for CPU calculations, **GPU4PySCF** for GPU
calculations, and **geomeTRIC** for optimization. Current validation covers small,
restricted closed-shell and small unrestricted gas-phase systems. APIs, supported methods, and installation
requirements may change as development proceeds.

See [scientific references, research papers, and websites](#scientific-references-research-papers-and-websites)
for method attribution, benchmark sources, and software documentation.

## What works today

| Capability | Current status |
| --- | --- |
| Explicit CPU/GPU backend selection | Implemented and tested; unavailable GPU requests fail clearly |
| XYZ plus identity, protonation, charge, spin, and fragment metadata | Implemented with input consistency checks |
| DFT single points and analytic nuclear gradients | Tested with density-fitted B3LYP/def2-SVP |
| Water and hydrogen-bonded water-dimer optimization | Convergence demonstrated on CPU and GPU |
| Frozen-atom optimization | Tested on CPU and GPU; external freeze/set constraints also accepted |
| Machine-readable calculation results | Settings, units, convergence, versions, hardware, timing, and chemical metadata |
| Optimized structures and trajectories | XYZ outputs and evaluated-geometry trajectories |
| Scientific validation and CPU CI | CPU/GPU agreement, finite differences, constraints, and explicit failure handling |
| CPU EDA preparation | Counterpoise fragment references, frozen electrostatics, generalized correlation; independent algebra checks passed, published component validation unresolved |
| Complete experimental GKS-EDA | Six gas-phase terms on CPU/GPU; independent mathematical tests and component-wise CPU/GPU agreement passed; external matched-settings validation pending |

The CLI provides `inspect`, `singlepoint`, `optimize`, `validate`, and experimental
`eda-prepare` and `eda` commands.
This project implements an independently useful optimization foundation; it does not
claim a GAMESS port or a validated EDA decomposition.

## Validation evidence

### Test computer and software environment

The local tests and neutral tryptamine–water benchmark were run on the following
computer. These are the recorded test specifications, not minimum requirements.

| Item | Test environment |
| --- | --- |
| CPU | AMD Ryzen 9 5950X; 16 physical cores, 32 logical threads |
| GPU | NVIDIA GeForce RTX 4080 SUPER; 16 GB VRAM (16,376 MiB reported) |
| System memory | 30.25 GiB usable RAM reported by Linux |
| Operating system | Ubuntu 26.04 LTS, x86-64 |
| NVIDIA driver | 595.91.07 |
| Python | 3.12.13 |
| PySCF / GPU4PySCF | 2.14.0 / 1.8.1 |
| CuPy / cuTensor | 14.2.0 / 2.8.1 |
| geomeTRIC | 1.1.1 |

The archived neutral EDA run used 16 physical cores/compute threads, 50% workload
pacing, and CPU/GPU thermal cooldown thresholds of 85°C/75°C. Earlier small-system
timings used their own recorded settings and are not directly comparable with this
paced run. GPU power limits and temporary CPU frequency limits are shared hardware
settings; no fixed frequency or power setting is claimed for all tests.
Exact per-run hardware, versions, and calculation settings are saved in result JSON.

On the development host (RTX 4080 SUPER, Ryzen 9 5950X, Ubuntu 26.04), **37 local
tests passed in the latest full suite, including actual GPU tests**. Separately,
nine optimization-foundation scientific checks passed with actual GPU execution. Water and
water-dimer CPU/GPU single-point energies differed by at most **9.21e-12 hartree**;
maximum gradient differences were **7.35e-7 hartree/bohr**. Water finite-difference
gradients agreed within **3.59e-6 hartree/bohr**, and frozen oxygen displacement was
below **1.79e-11 angstrom**. Both molecules converged on both backends.

These results apply to the recorded settings and small benchmark cases. Optimization
convergence establishes a stationary-point result; a vibrational analysis has not
established a minimum. GitHub CPU CI also passed; its GPU test was explicitly skipped
because that runner has no GPU. Warmed, repeated timing measurements are provided
without extrapolating them to larger systems.

See [validation details](docs/VALIDATION.md),
[full scientific results](reports/validation_2026-10-11/result.json), and
[benchmark measurements](reports/benchmark_2026-10-11.json).

## Development roadmap

1. Extend the structure-preparation and optimization foundation with additional
   carefully validated molecules, installation paths, and supported settings.
2. Secure the original GKS-EDA derivation, supporting material, and trustworthy
   individual component reference calculations under identical settings.
3. Specify every intermediate state and component convention, including orbital
   normalization, exchange/correlation, fragment reference geometry, and BSSE policy
   (experimental closed-shell specification now implemented).
4. Validate the implemented experimental CPU GKS-EDA against external matched-settings
   outputs, component by component.
5. Extend the implemented GPU component path after external validation; improve
   performance beyond the current small-system direct-integral measurements.
6. Explore nucleotide/RNA–triethylammonium–stationary-phase models after the scientific
   foundation is established, with explicit chemical states, solvent/sampling assumptions,
   and experimental validation before interpreting HPLC retention trends.

Complete restricted/unrestricted EDA is implemented experimentally; an externally validated EDA
release and HPLC remain future work. Total-energy closure alone will not count
as component validation, and ALMO-EDA or other decompositions will not be substituted
for the original GKS-EDA definitions. The exact current blockers are documented in
[the GKS-EDA source gate](docs/GKS_EDA_SPEC.md).

## Install

Use Python 3.12 in an isolated environment; do not upgrade system Python/drivers.
CPU: `python -m pip install -r requirements-cpu-lock.txt -e '.[test]'`.
The actual host GPU installation and exact pins are in `requirements-gpu-lock.txt`
after verified local execution. Install with `python -m pip install -r requirements-gpu-lock.txt -e '.[test]'`. CUDA wheel names depend on driver/runtime compatibility.
CUDA 12 runtimes can coexist with the host CUDA 13 toolkit. Source
`scripts/gpu-env.sh` before GPU execution; it exposes the environment's NVIDIA
wheel libraries without modifying the system loader. The script targets the verified
Python 3.12 `.env` layout. Exact host installation:

```bash
/home/zyliu/anaconda3/bin/conda create -p "$PWD/.env" python=3.12 pip -y
.env/bin/python -m pip install -r requirements-gpu-lock.txt -e '.[test]'
source scripts/gpu-env.sh
```

The optimization foundation passed nine scientific checks with actual GPU execution.
Experimental CPU EDA preparation adds independent electrostatic and correlation checks;
a validated full GKS-EDA release remains unavailable. See the recorded evidence.

```bash
source scripts/gpu-env.sh # GPU wheel library paths
.env/bin/eda-gpu inspect
.env/bin/eda-gpu singlepoint examples/water.xyz --metadata examples/water.json --backend cpu --output runs/water_sp
.env/bin/eda-gpu optimize examples/water_dimer.xyz --metadata examples/water_dimer.json --backend gpu --output runs/dimer_gpu
.env/bin/eda-gpu optimize examples/water.xyz --metadata examples/water_frozen.json --backend cpu --output runs/frozen
.env/bin/eda-gpu validate --gpu --output runs/validation
.env/bin/python -m pytest -q
```

Output directories must be new, preserving existing work. Failed calculations write
`result.json` with an error and return nonzero. GPU failures never fall back to CPU.
`validate --gpu` requires GPU execution; without it GPU validation is explicitly absent.

JSON metadata requires identity, protonation, integer charge and spin (Nalpha-Nbeta,
not multiplicity), and a complete fragment partition with charges/spins. Atom indices
are zero based; geomeTRIC constraint files use one based indices. Bonds,
stereochemistry, and other provided chemical metadata are preserved in results.
XYZ alone cannot establish protonation or stereochemistry. No conformer generation
or HPLC retention prediction is implemented.

Results record energies in hartree, gradients in hartree/bohr, angstrom coordinates,
SCF/optimization status, method settings, runtime versions/hardware, constraints,
initial/final XYZ and evaluation trajectories. Convergence does not establish a minimum.

See [methods](docs/METHODS.md), [validation](docs/VALIDATION.md),
[GKS-EDA source gate](docs/GKS_EDA_SPEC.md), and the dated report under `reports/`.
## Scientific references, research papers, and websites

If EDA_GPU has benefited your research, please cite
[this project](https://github.com/zyliu-OU/EDA_GPU) with the version or commit used,
the GKS-EDA publications, and the software, functional, and basis-set papers relevant
to your calculation. The bibliography below covers the current implementation,
documented benchmarks, and supporting software. Each entry explains its role.

### GKS-EDA theory and reference program

1. **Su, P.; Li, H. (2009).** “Energy decomposition analysis of covalent bonds and
   intermolecular interactions.” *Journal of Chemical Physics* **131**, 014102.
   [DOI: 10.1063/1.3159673](https://doi.org/10.1063/1.3159673).
   Supplies the orbital intermediate energy expressions used by the experimental implementation.
   [University of Nebraska institutional record](https://digitalcommons.unl.edu/chemfacpub/22/).
2. **Su, P.; Jiang, Z.; Chen, Z.; Wu, W. (2014).** “Energy Decomposition Scheme
   Based on the Generalized Kohn–Sham Scheme.” *Journal of Physical Chemistry A*
   **118**, 2531–2542. [DOI: 10.1021/jp500405s](https://doi.org/10.1021/jp500405s).
   Original GKS-EDA method paper and the project's intended decomposition.
3. **Su, P.; Tang, Z.; Wu, W. (2020).** “Generalized Kohn-Sham energy
   decomposition analysis and its applications.” *WIREs Computational Molecular
   Science* **10**, e1460. [DOI: 10.1002/wcms.1460](https://doi.org/10.1002/wcms.1460).
   Review used to map KS orbital intermediates, component definitions, and generalized correlation.
4. **Tang, Z.; Song, Y.; Zhang, S.; Wang, W.; Xu, Y.; Wu, D.; Wu, W.; Su, P.
   (2021).** “XEDA, a fast and multipurpose energy decomposition analysis program.”
   *Journal of Computational Chemistry* **42**, 2341–2351.
   [DOI: 10.1002/jcc.26765](https://doi.org/10.1002/jcc.26765).
   Reference-program publication; the current comparisons use published manual data.

Official reference resources:

- [XEDA documentation and tutorials](https://xmvb.xmu.edu.cn/docs/xeda/),
  including the [web manual, examples, and citation guidance](https://xmvb.xmu.edu.cn/docs/xeda/content.html).
- [XEDA II English manual](https://xmvb.xmu.edu.cn/docs/xeda/_downloads/415e0f39b016129f5cf2eb24dfe28b4f/XEDA_II_en.pdf):
  source of the CO₂ dimer geometry and rounded components in
  [the attributed reference fixture](references/xeda_ii_co2_b3lyp.json).
- [GAMESS XEDA manual](https://xmvb.xmu.edu.cn/docs/xeda/_downloads/222e6e0bd8ccfe3ce8fe39d17fc3eef3/XEDA_gms.pdf):
  reference interface and counterpoise keyword conventions.
- [Xiamen University 2014 publication list](https://pcoss.xmu.edu.cn/2014-paperlist.pdf):
  bibliographic provenance recorded during the source audit.

Source-access history and equation provenance are in the
[source audit](docs/EDA_SOURCE_AUDIT_2026-10-11.md) and
[mathematical specification](docs/GKS_EDA_SPEC.md). The original 2014 full text
has not been authenticated locally; the supplied review and the Su–Li formulas
support the experimental construction. No local XEDA reference execution has
been established. Citations and internal energy closure do not establish external
component validation. Historical transcription and author-upload links are preserved
in those audit documents with their access and reliability limitations.

### Published molecular benchmark

5. **Liu, Z.; Trindle, C. O.; Gu, Q.; Wu, W.; Su, P. (2017).** “Unravelling
   hydrogen bonding interactions of tryptamine–water dimer from neutral to cation.”
   *Physical Chemistry Chemical Physics* **19**, 25260–25269.
   [DOI: 10.1039/C7CP03491A](https://doi.org/10.1039/C7CP03491A).
   Source of the neutral/radical-cation M06-2X/aug-cc-pVDZ comparison and Table 3
   component values in [the attributed benchmark fixture](references/tryptamine_water_2017.json).
   [Publisher article and supporting-information navigation](https://pubs.rsc.org/en/content/articlelanding/2017/cp/c7cp03491a);
   [four-page supporting information](https://www.rsc.org/suppdata/c7/cp/c7cp03491a/c7cp03491a2.pdf).
   The supplied SI has no complete Cartesian geometries; comparisons using the
   constructed seed remain exploratory.

### Electronic structure, optimization, and numerical software

6. **Sun, Q., et al. (2018).** “PySCF: the Python-based simulations of chemistry
   framework.” *WIREs Computational Molecular Science* **8**, e1340.
   [DOI: 10.1002/wcms.1340](https://doi.org/10.1002/wcms.1340).
   CPU electronic-structure foundation.
7. **Sun, Q., et al. (2020).** “Recent developments in the PySCF program package.”
   *Journal of Chemical Physics* **153**, 024109.
   [DOI: 10.1063/5.0006074](https://doi.org/10.1063/5.0006074).
   PySCF methods and program development. Also see **Sun, Q., et al. (2026)**,
   “The Python simulations of chemistry framework: 10 years of an open-source
   quantum chemistry project,” *Journal of Chemical Physics* **165**, 102502,
   listed on the [official citation page](https://pyscf.org/about.html#how-to-cite).
8. **Sun, Q. (2015).** “Libcint: An efficient general integral library for Gaussian
   basis functions.” *Journal of Computational Chemistry* **36**, 1664–1671.
   [DOI: 10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981).
   Gaussian-integral library underlying PySCF.
9. **Li, R.; Sun, Q.; Zhang, X.; Chan, G. K.-L. (2024).** “Introducing
   GPU-acceleration into the Python-based Simulations of Chemistry Framework.”
   [arXiv:2407.09700](https://arxiv.org/abs/2407.09700).
   GPU4PySCF reference listed by its official repository.
10. **Wu, X., et al. (2024).** “Enhancing GPU-acceleration in the Python-based
    Simulations of Chemistry Framework.”
    [arXiv:2404.09452](https://arxiv.org/abs/2404.09452).
    Additional GPU4PySCF reference listed by its official repository.
11. **Wang, L.-P.; Song, C. (2016).** “Geometry optimization made simple with
    translation and rotation coordinates.” *Journal of Chemical Physics* **144**,
    214108. [DOI: 10.1063/1.4952956](https://doi.org/10.1063/1.4952956).
    geomeTRIC optimizer and TRIC coordinates;
    [2022 erratum](https://doi.org/10.1063/5.0102029).
12. **Lehtola, S.; Steigemann, C.; Oliveira, M. J. T.; Marques, M. A. L. (2018).**
    “Recent developments in Libxc — A comprehensive library of functionals for
    density functional theory.” *SoftwareX* **7**, 1–5.
    [DOI: 10.1016/j.softx.2017.11.002](https://doi.org/10.1016/j.softx.2017.11.002).
    Exchange-correlation functional library used through the calculation backends.
13. **Harris, C. R., et al. (2020).** “Array programming with NumPy.” *Nature*
    **585**, 357–362. [DOI: 10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2).
    CPU arrays, density algebra, and result processing.
14. **Virtanen, P., et al. (2020).** “SciPy 1.0: fundamental algorithms for
    scientific computing in Python.” *Nature Methods* **17**, 261–272.
    [DOI: 10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).
    Numerical infrastructure used by the scientific software stack.

Official software websites and implementation documentation:

- **PySCF:** [website](https://pyscf.org/), [source](https://github.com/pyscf/pyscf),
  [geometry optimization](https://pyscf.org/user/geomopt.html),
  [GPU interface](https://pyscf.org/user/gpu.html), and
  [density fitting](https://pyscf.org/user/df.html).
- **GPU4PySCF:** [source and citation guidance](https://github.com/pyscf/gpu4pyscf);
  [v1.8.1 direct J/K implementation](https://github.com/pyscf/gpu4pyscf/blob/v1.8.1/gpu4pyscf/scf/hf.py)
  inspected for the CPU/GPU hybrid execution description.
- **geomeTRIC:** [documentation](https://geometric.readthedocs.io/en/latest/),
  [constraints](https://geometric.readthedocs.io/en/latest/constraints.html), and
  [source and citation guidance](https://github.com/leeping/geomeTRIC).
- **Libxc / libcint:** [Libxc website and citations](https://libxc.gitlab.io/),
  [functional definitions](https://libxc.gitlab.io/functionals/), and
  [libcint source](https://github.com/sunqm/libcint).
- **CuPy:** [website](https://cupy.dev/), [documentation](https://docs.cupy.dev/en/stable/),
  and [source](https://github.com/cupy/cupy); GPU arrays and linear algebra.
- **NVIDIA:** [CUDA documentation](https://docs.nvidia.com/cuda/),
  [cuTENSOR documentation](https://docs.nvidia.com/cuda/cutensor/latest/), and
  [NVIDIA System Management Interface](https://docs.nvidia.com/deploy/nvidia-smi/);
  runtime, tensor contractions, and hardware monitoring.
- **NumPy / SciPy:** [NumPy citation guidance](https://numpy.org/citing-numpy/)
  and [SciPy citation guidance](https://scipy.org/citing-scipy/).
- **RDKit (optional seed generation):** [website](https://www.rdkit.org/),
  [overview and citation guidance](https://www.rdkit.org/docs/Overview.html), and
  [source](https://github.com/rdkit/rdkit). Used to construct the labelled approximate
  tryptamine seed; ordinary electronic-structure calculations do not require it.
- **Supporting tools:** [h5py](https://www.h5py.org/),
  [NetworkX](https://networkx.org/), [threadpoolctl](https://github.com/joblib/threadpoolctl),
  and [pytest](https://docs.pytest.org/en/stable/); data handling, optimizer graph
  infrastructure, thread limits, and tests. Exact installed versions are recorded
  in the lock files and calculation results.

### Density functionals and counterpoise convention

15. **Boys, S. F.; Bernardi, F. (1970).** “The calculation of small molecular
    interactions by the differences of separate total energies. Some procedures
    with reduced errors.” *Molecular Physics* **19**, 553–566.
    [DOI: 10.1080/00268977000101561](https://doi.org/10.1080/00268977000101561).
    Ghost-basis counterpoise convention for fragment references and BSSE treatment.
16. **Becke, A. D. (1993).** “Density-functional thermochemistry. III. The role of
    exact exchange.” *Journal of Chemical Physics* **98**, 5648–5652.
    [DOI: 10.1063/1.464913](https://doi.org/10.1063/1.464913).
    Three-parameter hybrid construction relevant to B3LYP.
17. **Stephens, P. J.; Devlin, F. J.; Chabalowski, C. F.; Frisch, M. J. (1994).**
    “Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra
    Using Density Functional Force Fields.” *Journal of Physical Chemistry* **98**,
    11623–11627. [DOI: 10.1021/j100096a001](https://doi.org/10.1021/j100096a001).
    B3LYP reference returned by the installed Libxc interface.
18. **Becke, A. D. (1988).** “Density-functional exchange-energy approximation
    with correct asymptotic behavior.” *Physical Review A* **38**, 3098–3100.
    [DOI: 10.1103/PhysRevA.38.3098](https://doi.org/10.1103/PhysRevA.38.3098).
    B88 exchange contribution to B3LYP/B3LYP5.
19. **Lee, C.; Yang, W.; Parr, R. G. (1988).** “Development of the Colle-Salvetti
    correlation-energy formula into a functional of the electron density.”
    *Physical Review B* **37**, 785–789.
    [DOI: 10.1103/PhysRevB.37.785](https://doi.org/10.1103/PhysRevB.37.785).
    LYP correlation contribution to B3LYP/B3LYP5.
20. **Vosko, S. H.; Wilk, L.; Nusair, M. (1980).** “Accurate spin-dependent
    electron liquid correlation energies for local spin density calculations:
    a critical analysis.” *Canadian Journal of Physics* **58**, 1200–1211.
    [DOI: 10.1139/p80-159](https://doi.org/10.1139/p80-159).
    Local correlation and the explicitly recorded B3LYP variant convention.
21. **Perdew, J. P.; Burke, K.; Ernzerhof, M. (1996).** “Generalized Gradient
    Approximation Made Simple.” *Physical Review Letters* **77**, 3865–3868.
    [DOI: 10.1103/PhysRevLett.77.3865](https://doi.org/10.1103/PhysRevLett.77.3865).
    PBE exchange/correlation underlying PBE0 and independent PBE component checks.
22. **Adamo, C.; Barone, V. (1999).** “Toward reliable density functional methods
    without adjustable parameters: The PBE0 model.” *Journal of Chemical Physics*
    **110**, 6158–6170. [DOI: 10.1063/1.478522](https://doi.org/10.1063/1.478522).
    PBE0 global hybrid functional.
23. **Ernzerhof, M.; Scuseria, G. E. (1999).** “Assessment of the
    Perdew–Burke–Ernzerhof exchange-correlation functional.” *Journal of Chemical
    Physics* **110**, 5029–5036.
    [DOI: 10.1063/1.478401](https://doi.org/10.1063/1.478401).
    Companion PBE0 reference returned by the installed Libxc interface.
24. **Zhao, Y.; Truhlar, D. G. (2008).** “The M06 suite of density functionals
    for main group thermochemistry, thermochemical kinetics, noncovalent
    interactions, excited states, and transition elements: two new functionals
    and systematic testing of four M06-class functionals and 12 other functionals.”
    *Theoretical Chemistry Accounts* **120**, 215–241.
    [DOI: 10.1007/s00214-007-0310-x](https://doi.org/10.1007/s00214-007-0310-x).
    M06-2X functional used for the tryptamine–water benchmark.

### Orbital and auxiliary basis sets

25. **Weigend, F.; Ahlrichs, R. (2005).** “Balanced basis sets of split valence,
    triple zeta valence and quadruple zeta valence quality for H to Rn: Design and
    assessment of accuracy.” *Physical Chemistry Chemical Physics* **7**, 3297–3305.
    [DOI: 10.1039/B508541A](https://doi.org/10.1039/B508541A).
    def2-SVP basis used in the optimization-foundation validation.
26. **Weigend, F. (2002).** “A fully direct RI-HF algorithm: Implementation,
    optimised auxiliary basis sets, demonstration of accuracy and efficiency.”
    *Physical Chemistry Chemical Physics* **4**, 4285–4291.
    [DOI: 10.1039/B204199P](https://doi.org/10.1039/B204199P).
    Correlation-consistent J/K-fitting auxiliary bases relevant to the
    aug-cc-pVDZ-JKFIT density-fitted optimization path.
27. **Dunning, T. H., Jr. (1989).** “Gaussian basis sets for use in correlated
    molecular calculations. I. The atoms boron through neon and hydrogen.”
    *Journal of Chemical Physics* **90**, 1007–1023.
    [DOI: 10.1063/1.456153](https://doi.org/10.1063/1.456153).
    Correlation-consistent cc-pVDZ foundation.
28. **Kendall, R. A.; Dunning, T. H., Jr.; Harrison, R. J. (1992).** “Electron
    affinities of the first-row atoms revisited. Systematic basis sets and wave
    functions.” *Journal of Chemical Physics* **96**, 6796–6806.
    [DOI: 10.1063/1.462569](https://doi.org/10.1063/1.462569).
    Diffuse augmentation underlying aug-cc-pVDZ in the molecular benchmark.
29. **Hehre, W. J.; Ditchfield, R.; Pople, J. A. (1972).** “Self—Consistent
    Molecular Orbital Methods. XII. Further Extensions of Gaussian—Type Basis
    Sets for Use in Molecular Orbital Studies of Organic Molecules.”
    *Journal of Chemical Physics* **56**, 2257–2261.
    [DOI: 10.1063/1.1677527](https://doi.org/10.1063/1.1677527).
    6-31G basis underlying the XEDA CO₂ example.
30. **Hariharan, P. C.; Pople, J. A. (1973).** “The influence of polarization
    functions on molecular orbital hydrogenation energies.” *Theoretica Chimica
    Acta* **28**, 213–222.
    [DOI: 10.1007/BF00533485](https://doi.org/10.1007/BF00533485).
    Polarization functions relevant to the example's 6-31G* basis.
31. **Weigend, F. (2008).** “Hartree–Fock exchange fitting basis sets for H to Rn.”
    *Journal of Computational Chemistry* **29**, 167–175.
    [DOI: 10.1002/jcc.20702](https://doi.org/10.1002/jcc.20702).
    Universal J/K-fitting auxiliary basis family used as def2-universal-jkfit
    in the optimization-foundation calculations.

Basis definitions and their source references can also be inspected through
[PySCF basis documentation](https://pyscf.org/user/gto.html#basis-set) and the
[Basis Set Exchange](https://www.basissetexchange.org/). Spherical/Cartesian
representation, auxiliary basis, and installed functional definitions must be
reported with each result because these choices affect comparison with published data.

### Project terms

This repository is public for development visibility. No project license has been
assigned yet. Published literature and manuals are linked; full-text papers are
not redistributed in this repository. Software projects retain their own licenses
and citation guidance.

## Experimental CPU EDA preparation

The supplied 2020 review was read completely. Its unambiguous definitions now support
frozen electrostatics and generalized correlation, using direct integrals and ghost-basis
counterpoise references at the complex geometry. Exactly two closed-shell, all-electron
fragments and B3LYP/PBE0 are supported by this older CPU-only command.
Use `eda` below for the complete experimental CPU/GPU component calculation.

```bash
.env/bin/eda-gpu eda-prepare examples/water_dimer.xyz --metadata examples/water_dimer.json --backend cpu --output runs/eda_dimer
.env/bin/eda-gpu eda-prepare examples/co2_dimer_xeda_ii.xyz --metadata examples/co2_dimer_xeda_ii.json --backend cpu --basis '6-31g*' --output runs/eda_co2
```

A successful preparation reports `full_eda_available: false`. Exchange, repulsion,
polarization, and total closure remain null until their intermediate-state conventions
and component references are established. Saved densities support further development.
The official XEDA II CO₂ example differs from our partial calculation beyond printed
rounding; its settings are not fully specified, so it is an exploratory comparison. See
[comparison evidence](reports/eda_preparation_2026-10-11/co2_reference_comparison.json),
[mathematical specification](docs/GKS_EDA_SPEC.md), and
[source audit](docs/EDA_SOURCE_AUDIT_2026-10-11.md). Papers remain local.

## Complete experimental CPU/GPU GKS-EDA

`eda` now computes electrostatic, exchange, repulsion, polarization (including charge
transfer), generalized correlation, and explicit dispersion (zero for the supported
non-dispersion-corrected methods). Exactly two all-electron fragments with additive spins,
B3LYP/B3LYP5/PBE0, counterpoise references and direct integrals are supported. It
uses the explicit Su–Li orbital intermediates on KS orbitals, as specified in the
supplied review; no ALMO-EDA substitution is made.

```bash
source scripts/gpu-env.sh
.env/bin/eda-gpu eda examples/water_dimer.xyz --metadata examples/water_dimer.json --backend gpu --output runs/gks_dimer_gpu
.env/bin/eda-gpu eda examples/water_dimer.xyz --metadata examples/water_dimer.json --backend cpu --output runs/gks_dimer_cpu
.env/bin/eda-gpu eda examples/co2_dimer_xeda_ii.xyz --metadata examples/co2_dimer_xeda_ii.json --backend gpu --basis '6-31g*' --cart --output runs/gks_co2_gpu
.env/bin/python scripts/validate_gks.py --gpu --repeats 3 --output runs/gks_checks
```

All terms, intermediate energies, condition checks and reference densities are saved.
Results state `all_components_computed: true`, `scientifically_validated: false`, and
`full_eda_available: false`: the last flag is the external-validation release gate.
A successful process means the experimental calculation converged and its internal
checks passed. It does not imply validated agreement with XEDA.

**28 local tests passed**, including actual GPU EDA, independent integral/PBE0
checks and a many-electron determinant energy check. Repeated CPU/GPU comparisons
on water and CO₂ dimers have a maximum individual-component difference of
**1.82e-8 hartree**. Cartesian functions explain most of the earlier CO₂ reference
discrepancy, but residual component differences remain beyond printed rounding.

GPU4PySCF direct J/K has CPU/GPU hybrid internals; this is not an all-GPU workflow.
For these tiny cases, measured complete EDA GPU medians were slower than CPU
(24.43 versus 3.15 seconds for water dimer; 29.47 versus 9.19 for CO₂). These
measurements do not establish larger-system throughput or speedup. See
[GPU EDA validation and timings](reports/gks_eda_gpu_2026-10-11/result.json),
[methods](docs/GKS_EDA_SPEC.md), and [validation details](docs/VALIDATION.md).

The next scientific gate is individual component output from XEDA with explicitly
matched functional, basis representation, quadrature, integral and counterpoise
settings. Open-shell, solvent, double-hybrid and RNA/HPLC applications remain deferred.


### Tryptamine–water paper benchmark (under development)

The benchmark follows Liu et al., PCCP 2017, DOI
[10.1039/C7CP03491A](https://doi.org/10.1039/C7CP03491A): neutral M06-2X and
radical-cation U-M06-2X with aug-cc-pVDZ and counterpoise EDA. The cation is
created by removing one electron, **not** by protonation. Restricted and unrestricted
CPU/GPU optimization and complete experimental EDA are supported. Open-shell
fragment spins must add to the complex spin; broken-symmetry/antiferromagnetic
coupling is outside the implemented scope. S² is recorded for unrestricted SCFs.

Published Table 3 data are attributed in `references/tryptamine_water_2017.json`.
The supplied four-page SI contains molecular figures and charges, but no full
Cartesian coordinates. `examples/tryptamine_water/constructed_bridge.xyz` is an
explicitly labelled approximate seed constructed from the SI bonding pattern,
not the paper's original geometry. Optional seed generation requires
`rdkit==2026.9.1`; ordinary calculations do not require RDKit.

```bash
source scripts/gpu-env.sh
.env/bin/python scripts/compare_tryptamine.py \
  --xyz examples/tryptamine_water/constructed_bridge.xyz \
  --metadata examples/tryptamine_water/neutral.json \
  --case A_H2O --opt-density-fitting --output runs/tryptamine_comparison
```

This runs CPU/GPU density-fitted optimization (aug-cc-pVDZ-JKFIT) sequentially from identical input, then direct
counterpoise EDA on the same GPU-optimized geometry. Outputs include coordinates,
trajectories, settings, spin diagnostics, aligned geometry RMSD, all EDA terms,
backend differences and differences from the paper. Pauli is exchange plus
repulsion, matching Table 3. Paper agreement remains exploratory until original
coordinates and numerical settings are matched. No ZPE/deformation-corrected D0,
minimum/Hessian certification, intra-EDA or individual-H-bond decomposition is
claimed. A separate H⁺–H₂O starting conformer is required to study that isomer.

### Live calculation progress and research attribution

```bash
.env/bin/eda-gpu monitor runs/my-calculation
.env/bin/eda-gpu status runs/my-calculation
```

`monitor` follows SCF cycles, optimization evaluations, energies, gradient norms,
temperature cooldowns and final EDA components. Ctrl+C exits the viewer while the
calculation continues. `status` prints one snapshot; `--json` gives machine-readable
output. For older running jobs, add `--log /path/to/calculation.log`. New calculations
record atomic `progress.json` snapshots and an `events.jsonl` history beside results.

Calculation output, monitoring and JSON results carry this acknowledgement:
“If EDA_GPU has benefited your research, please support us by citing this project
and the GKS-EDA publications below.” Method attribution includes Peifeng Su et al.,
[original GKS-EDA, 2014](https://doi.org/10.1021/jp500405s),
[Su, Tang and Wu review, 2020](https://doi.org/10.1002/wcms.1460), and
[Su and Li orbital derivation, 2009](https://doi.org/10.1063/1.3159673).
Tryptamine runs additionally cite [Liu et al., 2017](https://doi.org/10.1039/C7CP03491A).

### Computing limits

Molecular CLI commands default to 16 physical cores/compute threads, 50% workload
pacing, GPU cooldown at 75°C until 65°C and CPU cooldown at 85°C until 75°C.
These settings affect this calculation process. `--threads`, `--cpu-core-limit`,
`--compute-duty-cycle`, `--gpu-temperature-limit` and `--cpu-temperature-limit`
are configurable; temperature limit 0 disables that sensor guard. Library Settings
retain unpaced defaults for explicitly controlled benchmarks. CPU sensor availability
and actual pacing/temperature statistics are reported. Software checks occur at SCF
cycles and workload boundaries; long kernels cannot be interrupted by these controls.
Pacing is an average workload target, not a fixed instantaneous utilization or watt cap.

GPU power/clocks are never changed automatically. `nvidia-smi -pl` changes the shared
GPU and affects other projects. The developer's original RTX 4080 SUPER power limit
was 320 W; restoring it is an explicit administrator action.

For an optional **temporary system-wide CPU frequency cap**:

```bash
sudo /usr/bin/python3 scripts/cpu_frequency_cap.py apply --max-mhz 3000 --minutes 60
sudo /usr/bin/python3 scripts/cpu_frequency_cap.py restore
```

The helper saves the original limits, applies the ceiling across CPU policies and
starts an automatic restore timer. Manual restore can end it early; later manual
policy changes are preserved. This hardware setting affects other CPU applications
while active and is not installed as a boot service. Neither this helper nor process
pacing modifies firmware or Windows configuration.

### Completed exploratory neutral tryptamine–water benchmark

GPU M06-2X/aug-cc-pVDZ direct counterpoise GKS-EDA at the GPU
density-fitted optimized geometry converged with all six components.
Values below are kcal/mol, compared with PCCP Table 3:

| Component | GPU | Paper |
|---|---:|---:|
| Electrostatic | -17.279 | -17.34 |
| Pauli | 19.535 | 19.69 |
| Polarization | -7.048 | -7.07 |
| Generalized correlation | -5.489 | -5.54 |
| Total | -10.282 | -10.27 |

This is encouraging exploratory agreement, not a verified reproduction.
The maximum component difference is 0.155 kcal/mol. Original Cartesian
coordinates are unavailable; the seed was reconstructed. Paper grid settings
are unspecified. CPU comparison and cation calculations remain pending.
SCF thresholds were 1e-9 Eh and orbital-gradient 1e-6; an earlier 1e-12 Eh
attempt did not converge and was excluded. Energy closure residual was
8.45e-11 Eh. See [archived results](reports/tryptamine_neutral_gpu_2026-10-11/paper_comparison.json).
