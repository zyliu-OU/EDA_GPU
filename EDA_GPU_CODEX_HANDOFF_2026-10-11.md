# EDA_GPU — Codex terminal build handoff

Date: 2026-10-11 (America/Los_Angeles)
Project name: `EDA_GPU`
Recommended target workspace: `/home/zyliu/Documents/EDA_GPU`
Status: Design only. No project repository, optimizer, EDA code, GPU tests, or GitHub publication has been created or verified here.

## Paste this prompt into Codex terminal

```text
Build my new project EDA_GPU. Read EDA_GPU_CODEX_HANDOFF_2026-10-11.md completely first.

Use /home/zyliu/Documents/EDA_GPU unless an existing EDA_GPU repository is found. Inspect the current directory, applicable AGENTS.md files, existing files, git status/remotes, hardware, NVIDIA driver, Python, and available environments before editing or installing. Preserve existing work. Keep this project separate from RNA-StableAI.

Start with a working structure-preparation and geometry-optimization module using PySCF as the CPU backend, GPU4PySCF as the GPU backend, and geomeTRIC as the optimization controller. This backend choice is accepted for the initial implementation; exact GAMESS compatibility is not required. Do not describe this as a full GAMESS port.

Use an isolated environment and verify current official documentation and installed APIs. Build a CLI with inspect, singlepoint, optimize, and validate commands. Start with water and a small hydrogen-bonded dimer. Support explicit backend selection, XYZ plus charge/spin/fragment metadata, constraints where supported, reproducible settings, and machine-readable outputs. A requested GPU run must fail clearly if GPU support is unavailable; never silently run CPU and label it GPU. Validate CPU/GPU energy and gradient agreement under equivalent settings, optimization convergence, finite-difference gradients on a small case, and one frozen-atom constraint example.

Then secure the original Su–Tang–Wu GKS-EDA derivation and component reference data, write a precise mathematical specification, implement a CPU reference, validate every component, and port costly operations to GPU. Do not substitute ALMO-EDA, invent component equations, or claim validation from total-energy closure alone. If source access or reference data block faithful EDA, finish the working optimizer, document the exact blocker, and keep EDA explicitly unavailable rather than fabricate a decomposition.

I authorize you to commit and push the project to GitHub. Inspect gh authentication and the authenticated account without printing credentials. Reuse an existing verified EDA_GPU remote. If no matching repository exists, create EDA_GPU as a private repository under the authenticated personal account and push verified work. Do not assume the account from prior conversation, overwrite unrelated repositories, change an existing repository's visibility, or force-push. If the account/remote is ambiguous or authentication is unavailable, complete local work and tell me precisely what is needed. Do not publish employer SOPs, confidential data, credentials, or copyrighted full-text papers.

Maintain README.md, docs/METHODS.md, docs/VALIDATION.md, and reports/EDA_GPU_HANDOFF_<date>.md. Commit coherent tested stages and update GitHub as work progresses. Distinguish implemented, tested, experimental, and blocked features. End with the repository URL, commit SHA, completed capabilities, validation evidence, and remaining blockers. Continue autonomously through feasible stages; this handoff authorizes the build and GitHub updates.
```

## First usable milestone

A package with a documented CLI and examples that can run small CPU calculations and, on verified compatible hardware, GPU single points, analytic gradients, and geometry optimization. The optimizer is independently useful while the EDA reference is being secured.

Suggested layout (adapt after inspecting the workspace):

- `pyproject.toml`: packaging, CLI entry point, optional backend dependencies.
- `src/eda_gpu/`: input validation, CPU/GPU backends, optimization, outputs, CLI.
- `examples/`: openly shareable small benchmark XYZ and metadata files.
- `tests/`: input validation and scientifically meaningful calculation checks.
- `docs/METHODS.md`: methods, units, conventions, constraints, approximation choices.
- `docs/GKS_EDA_SPEC.md`: source-backed equations and intermediate-state definitions before EDA implementation.
- `docs/VALIDATION.md`: component and backend comparisons with reproducible commands.
- `reports/`: environment summary, benchmark tables, build/resume handoffs.

Avoid promising fixed dependency versions or CUDA package names until the actual machine and current compatibility requirements have been checked. Record exact working versions and separate CPU-only installation from the verified GPU installation. Do not upgrade system Python or NVIDIA drivers merely to get started.

## Inputs and outputs

Define charge explicitly; define spin as PySCF's Nalpha minus Nbeta, not multiplicity. Check electron-count/spin consistency. First EDA scope is restricted closed-shell only. Track atomic ordering and fragment membership through optimization, validate that each atom has exactly one fragment assignment for the chosen partition, and record fragment charges/spins. Check fragment-charge sums against the complex.

Inputs should retain chemical identity/protonation and relevant bonding/stereochemistry metadata beyond XYZ. Do not guess nucleotide protonation from an XYZ file or a pH value alone. Use simple provided structures for initial tests; add conformer generation as a separate documented capability later.

Outputs should include initial/optimized XYZ, an optimization trajectory when available, final energy in hartree, gradients in explicit units, final convergence status, SCF status, constraints, method/basis/grid/density-fitting settings, backend actually used, timings, versions, and hardware summary. Handle unsuccessful SCF or optimization explicitly; a completed process is not proof of convergence. EDA outputs later need component units, total closure residual, fragment reference geometry, and BSSE policy.

## Acceptance gates

1. **Environment:** verify the real host and record a reproducible isolated installation; no GPU claim without successful GPU execution.
2. **Single points/gradients:** compare equivalent CPU/GPU calculations, including density fitting, grids, thresholds, and dispersion settings; document numerical tolerances and justify them. Check analytic gradients against finite differences on a small system.
3. **Optimization:** demonstrate convergence for water and a small dimer, and verify frozen atoms remain fixed within a stated tolerance. Report constrained stationary points accurately; optimizer convergence alone does not establish a minimum.
4. **EDA specification:** source-backed definitions for all intermediate states, spin factors, exchange/correlation conventions, normalization, reference geometries, and BSSE treatment.
5. **CPU EDA:** compare individual components with trustworthy reference calculations under identical settings. Arithmetic closure is necessary but insufficient.
6. **GPU EDA:** component-by-component agreement with the CPU reference and reproducible benchmarks on the actual host. Include warm-up policy, CPU thread settings, repeated timings, memory, and remaining CPU operations.
7. **HPLC application:** begin with explicitly labeled small models and consistent reference states; quantify sampling and model limitations before interpreting retention trends.

Pure CPU CI can validate input/schema behavior and small reference calculations. GPU tests must be marked and skipped honestly on hosts without hardware. Never report a skipped GPU test as passed validation. Do not turn guessed EDA numbers into test fixtures.

## GitHub delivery rules

The GitHub build and updates are delegated to Codex on the user's PC; no repository has been created by this handoff. Private is the default for a new repository because visibility was unspecified. Preserve an existing repository's visibility. Public release can be chosen later.

Use the project name `EDA_GPU` for the directory and repository; `eda_gpu` is a suitable Python import name. Inspect staged changes before committing, use a focused .gitignore for environments, caches, credentials, and bulky run artifacts, and include small reproducible validation results. Do not add a license with assumed ownership terms; ask if a license choice becomes necessary for release. Cite the literature and link to papers; keep supplied PDFs local unless redistribution rights are established.

If GitHub authentication requires the user, do not expose tokens or claim the push succeeded. Preserve local commits and provide the next exact action. Each final status should distinguish local commit, remote push, and actual hardware validation.

## Scientific design inherited from the earlier handoff

The original handoff below remains the scientific context. Its earlier open question about optimizer backend is resolved for the first build: use GPU4PySCF + geomeTRIC without exact GAMESS compatibility. Initial benchmark chemistry is fixed; nucleotide/RNA choice and experimental HPLC conditions can be resolved when that application stage begins. References below were checked in the previous discussion and must be verified against current documentation before implementation.

---

# GPU GKS-EDA and HPLC interaction modeling — handoff

Date: 2026-10-11 (America/Los_Angeles)
Status: Research/design only. No optimizer, GKS-EDA implementation, repository, GPU environment, or computational benchmark has been built or verified in this conversation.

## Resume prompt

> Resume my GPU GKS-EDA / RNA ion-pair HPLC project. Read this handoff first. Build the structure-preparation and GPU geometry-optimization module first, then implement a CPU reference and GPU version of the original Su–Tang–Wu GKS-EDA definitions. Use GPU4PySCF and geomeTRIC where appropriate; do not call an ALMO-EDA implementation GKS-EDA or claim to have ported all of GAMESS. Begin with small closed-shell benchmark complexes, validate individual components, and then study nucleotide or short RNA fragment interactions with TEAH+ and a model bonded C18 surface. Preserve existing files and clearly distinguish tested behavior from proposed features. Inspect my actual hardware, environment, and any existing repository before installing or editing. Do not predict RNA retention directly from isolated NTP binding energies without experimental validation.

## User goal

Develop a GPU-assisted molecular workflow to understand how A/G/C/U composition influences ion-pair reversed-phase HPLC retention under triethylamine conditions. The initial question is the interaction of each nucleotide with stationary-phase ligands in the presence of protonated triethylammonium (TEAH+).

Two connected modules are wanted:

1. Prepare and optimize 3D structures before EDA.
2. Calculate a faithful GPU GKS-EDA decomposition for representative configurations.

The user asked about making a GPU version of GAMESS. The recommendation was GPU4PySCF for energies/gradients plus geomeTRIC for optimization, rather than a full GAMESS rewrite. This recommendation has not been implemented; clarify whether exact GAMESS compatibility is required before substituting any GAMESS-specific functionality.

## Target computer and workspace

User-reported target: Ubuntu 26.04 LTS; AMD Ryzen 9 5950X; approximately 30 GiB RAM; NVIDIA RTX 4080 SUPER; approximately 1.6 TB root disk.

These are previously reported specifications, not verified during this project. The current ChatGPT scratch workspace is not the user's Ubuntu PC. No project directory on that PC was chosen. Suggested new directory: `/home/zyliu/Documents/gpu-gks-eda-hplc` (proposal only). Keep separate from the existing RNA-StableAI project.

Check `nvidia-smi`, driver/CUDA compatibility, Python version, available memory, and package versions on the actual target. Use an isolated environment; do not change system Python. No speedup has been measured. Data-center GPU benchmarks are not predictions for the RTX 4080 SUPER.

## Authoritative supplied paper

Su, P.; Tang, Z.; Wu, W. *Generalized Kohn-Sham energy decomposition analysis and its applications*. WIREs Computational Molecular Science, 2020, 10, e1460. DOI: https://doi.org/10.1002/wcms.1460

User supplied the complete PDF. It was read directly and methodology pages were inspected.

- Filename: `Generalized Kohn-Sham energy decomposition analysis and its applications.pdf`
- Current attachment path: `/workspace/scratch/7ef617d70ffb/upload/Generalized Kohn-Sham energy decomposition analysis and its applications.pdf`
- Persistent file ID: `libfile_5152eafbbc808191afafceaee70f5d08`

Do not assume the scratch path survives a new session. Retrieve the exact supplied file if needed. This is a review with useful equations, not a complete implementation specification. Consult the original derivation and supporting material before finalizing conventions:

Su et al., *Energy Decomposition Scheme Based on the Generalized Kohn–Sham Scheme*, J. Phys. Chem. A 2014, 118, 2531–2542. https://doi.org/10.1021/jp500405s

Original full text was not retrieved during this conversation. XEDA is the reference implementation; access to its source/executable has not been established.

## Mathematical requirements

The review's Eq. 8 decomposes complex-minus-fragment interaction energy:

`ΔEtotal = Ecomplex − Σ Efragment = ΔEele + ΔEex + ΔErep + ΔEpol + ΔEcorr + ΔEdisp`.

- Eqs. 2–5 define the electrostatic, exchange, repulsion, and polarization construction in LMO-EDA; Section 2.2 carries those definitions to GKS-EDA using KS orbitals.
- Intermediate constructions are the fragment direct product, antisymmetrized product, and antisymmetrized/orthonormalized product. Implement their energy evaluations exactly; do not approximate them by arbitrary constrained DFT states.
- Exchange and repulsion may be combined as exchange–repulsion in reporting, but keep individual values internally for validation.
- Polarization includes charge transfer/orbital relaxation; it is not the separately constrained polarization term of ALMO-EDA.
- Eq. 9: `ΔEcorr = (1−γ)(ΔEX − ΔEX_HF) + ΔEC`, where γ is the exact-exchange fraction. Extract signs and definitions consistently from the original derivation. The label “correlation” includes an exchange correction; do not substitute only the DFT correlation functional.
- Eq. 10 covers double hybrids. Exclude these from the first implementation.
- Eq. 11 defines explicit dispersion as complex dispersion correction minus the fragment corrections. Zero explicit dispersion does not imply no dispersion-like contribution in the functional.
- Eqs. 12–13 and Section 2.3 define GKS-EDA(sol), including desolvation and solvent-cavity/reaction-field treatment. A generic PCM energy difference is not automatically equivalent.
- The review discusses counterpoise/BSSE handling. Document ghost basis treatment and reference geometry explicitly.

Start with restricted, closed-shell, gas-phase, two-fragment systems and a global hybrid functional. Delay range-separated functionals, open shells, double hybrids, and the solvent extension until their definitions and reference results are secured.

## Proposed software architecture

- Structure preparation: accept XYZ plus separate charge, spin, fragment membership, protonation, bonding/stereochemistry, and frozen-atom metadata. XYZ alone does not encode chemical identity sufficiently.
- Optimization: GPU4PySCF energies/analytic gradients; geomeTRIC optimization control. Verify exact method/gradient/constraint support in installed versions.
- EDA: transparent CPU reference first; GPU integral/J/K evaluation and CuPy matrix operations second.
- Output: optimized XYZ, trajectories, convergence flags, energies/units, atom-to-fragment mapping, complete method settings, package versions, and EDA component tables.
- Preserve GPU arrays where possible, avoid unnecessary CPU transfers, and use double precision initially. Do not store the full four-index integral tensor for larger systems if direct/block or validated density-fitting evaluation is available.

GAMESS already has selected GPU paths through LibCChem; universal GPU DFT/gradient support was not established. Recompiling CPU code with CUDA is not a full GPU port. GPU4PySCF's experimental ALMO-EDA2 is a different decomposition and can only be a separately labeled comparator.

## Structure optimization requirements

Generate multiple starting conformations and ion/surface placements. Optimization gives local minima, not a complete finite-temperature adsorption ensemble. Define protonation/charge before structure generation; an optimizer does not select the correct pH-dependent species automatically.

For a bonded surface patch, freeze attachment/support atoms as appropriate and relax exposed groups. A free octadecane molecule is only a preliminary ligand model, not a complete C18 stationary phase.

Perform complex optimization, then use fragment coordinates from that complex for the primary interaction decomposition. Separately optimized fragment energies describe additional deformation contributions and must be reported separately. Validate minima or constrained stationary points using suitable gradient and, where warranted, vibrational checks.

## HPLC model and interpretation

Distinguish neutral TEA from TEAH+. Specify pH, solvent composition, counterions, temperature, and ion population. One placed TEAH+ does not represent a bulk TEA concentration; concentration-dependent conclusions require sampled occupancy/speciation or an experimentally calibrated model.

Include both solution ion association and TEAH+ adsorption on the stationary phase. The surface can be dynamically modified by adsorbed ion-pair reagent. Consider multiple tethered chains, ligand packing, residual silanols, explicit solvent contacts, and solvent competition as scope grows.

NTPs have triphosphates; internal RNA residues have phosphodiester linkages. Therefore:

- ATP/GTP/CTP/UTP models directly address free NTP behavior.
- Nucleosides or AMP/GMP/CMP/UMP can isolate preliminary chemical differences, but are imperfect RNA-residue proxies.
- Matched short RNA fragments such as CAC/CGC/CCC/CUC, with identical terminal chemistry, provide a more relevant comparison for RNA composition.

Use consistent charge, ion counts, surface model, solvent treatment, and configuration sampling for comparisons. Examine both pair interactions and three-body nonadditivity with a properly specified reference cycle. Compare adsorption of a preassociated nucleotide–TEAH+ ensemble to adsorption on a TEAH+-modified surface; avoid mixing reference states.

EDA electronic interaction energies are descriptors, not retention factors or transfer free energies. Solvent displacement, conformational entropy, adsorption kinetics, surface heterogeneity, and pore transport remain relevant. The solvent extension alone does not provide the full chromatographic free energy.

An exploratory additive descriptor is `S = N Σ fi εi`; it is a hypothesis, not a validated RNA energy or retention law. A fixed-length, fixed-condition experimental model can use `ln k = b0 + bA fA + bG fG + bC fC`, with U as reference because fractions sum to one. Coefficients must be calibrated; electronic energies are not directly those coefficients. Use appropriate gradient models for gradient retention rather than treating gradient retention time as isocratic k.

## Validation and build order

1. Inspect the actual machine and establish an isolated, pinned environment. Run small CPU/GPU energy and gradient comparisons.
2. Implement optimization on water, then a small hydrogen-bonded complex. Confirm constraints, convergence reporting, and deterministic outputs.
3. Obtain the original GKS-EDA derivation and trustworthy reference component data. Specify equations, spin factors, orbital normalization, cavity conventions, and BSSE policy before implementation.
4. Implement a CPU GKS-EDA reference; validate water dimer, ammonia–water, and benzene dimer with identical settings against XEDA/reference data.
5. Port the costly operations to GPU. Require individual component agreement and total-energy closure; cancellation can hide errors. Select numerical tolerances after accounting for reference rounding and integration/integral approximations.
6. Benchmark CPU/GPU using equivalent methods and hardware settings, report wall time, memory, and any remaining CPU operations.
7. Add nucleotide/TEAH+/surface models, then solvent treatment and sampling. Label simplified models explicitly.
8. Validate HPLC trends with same-length RNAs of different composition and same-composition sequence permutations. Keep column, pH, temperature, ion-pair reagent, and solvent program controlled. Evaluate held-out sequences; do not extrapolate short-oligo behavior to 2,000-nt RNA without evidence.

## Reference links checked during discussion

- GPU4PySCF: https://github.com/pyscf/gpu4pyscf
- GPU documentation: https://pyscf.org/user/gpu.html
- Geometry optimization: https://pyscf.org/user/geomopt.html
- XEDA manual: https://xmvb.xmu.edu.cn/docs/xeda/content.html
- XEDA paper: https://doi.org/10.1002/jcc.26765
- GAMESS release notes: https://www.msg.chem.iastate.edu/gamess/versions.html
- Oligonucleotide retention/composition model: https://pubmed.ncbi.nlm.nih.gov/12134814/
- Composition and sequence dependence: https://doi.org/10.1002/jssc.201600701
- Retention modeling using sequence/secondary structure: https://academic.oup.com/nar/article/35/12/4195/1134171
- Ion-pair reagent adsorption study: https://pubmed.ncbi.nlm.nih.gov/36470076/

These are references, not a frozen support matrix. Check current official documentation/source and installed versions before choosing APIs.

## Open decisions

- Is exact GAMESS compatibility required, or is the recommended GPU4PySCF optimizer acceptable?
- Is the first application free NTP chromatography or RNA-residue composition effects?
- Which column/surface chemistry and mobile-phase conditions should be modeled first?
- Which protonation states, counterions, and starting coordinates are available?
- Can XEDA reference calculations and the original derivation be obtained?

The only completed work is literature review, PDF methodology inspection, workflow design, and this handoff. No molecular interaction energies, optimized geometries, predictive coefficients, or GPU speedups have been produced.
