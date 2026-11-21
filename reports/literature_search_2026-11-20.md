# Literature search for published comparison data (2026-11-20)

Constraint (2026-11-20): only published data are used; no proprietary or company retention data.

## Base-composition retention data (HPLC demonstration)

| Source | What it has | Usable? |
| --- | --- | --- |
| **Sturm, Quinten, Huber, Kohlbacher 2007**, *Nucleic Acids Res.* 35:4195, doi:10.1093/nar/gkm338 (open access, PMC1919494) | 72 DNA oligonucleotides (15–48 nt), 100 mM TEAA pH 7, poly(styrene/divinylbenzene) monolith (no silanols), 30/40/50/60/80 °C, 0–16 % ACN. Normalized retention times in Supplementary **Table S1**. Includes two groups with the 24mer composition where one G is exchanged for A or for C | **Best candidate.** Silanol-free hydrophobic phase and TEAA are close to our patch and reagent models. The G→A / G→C swaps at fixed length isolate base effects; the 80 °C data minimise secondary structure. The SI (`nar_gkm338_1.html`) is behind a browser check, so it needs a manual download |
| Gilar et al. 2002, *J. Chromatogr. A* 958:167, doi:10.1016/S0021-9673(02)00306-0 | Retention model built from base counts and length, with two coefficients per base (aᵢ, bᵢ) fitted on homo-oligonucleotides. Valid at 60–80 °C | Coefficients would give a direct per-base contribution, but the paper is paywalled and the values were not found in open sources |
| Maurer et al. 2025, *J. Chromatogr. A* 1740:465574, Fig. 9 | NTP retention times (C < U < G < A), already used | In use (`examples/hplc/maurer2025_fig9_ntp.csv`) |
| Rahal et al. 2025, arXiv:2511.15753; Enmark et al. 2022, *J. Chromatogr. A* 1671 | 876 and ~100 antisense oligonucleotides. Base frequencies are the most important features | Not usable: no public data, and conditions not given |

## Solvated GKS-EDA benchmark (method paper)

| Source | What it has | Usable? |
| --- | --- | --- |
| **Shen, Su, Wu 2018**, *Phys. Chem. Chem. Phys.* 20:26126, doi:10.1039/C8CP05358H | GKS-EDA of neutral halogen bonds X′–X···Y (X′X = BrF, ClF, I₂, Br₂, Cl₂; Y = pyridine, NH₃, H₂S, HCN, H₂O, dimethyl ether) in gas phase and at several dielectric constants up to 78.36. ωB97X-D and M06-2X with aug-cc-pVDZ. The SI has component tables including desolvation | **Best candidate for an external GKS-EDA(sol) check**, and the systems are small. Unverified: whether coordinates and cavity radii are given. The RSC pages refused automated access, so the article and ESI need a manual download. I₂ and Br-containing systems need ECP/all-electron basis choices to match |
| Su, Liu, Wu 2012, *J. Chem. Phys.* 137:034111 (EDA-PCM) | Original solvated EDA (LMO-EDA based) | Possible secondary reference; not checked |

## Next actions
1. Manually download the Sturm 2007 supplementary data (Table S1) and the Shen 2018 article + ESI into the repo root
   (like the two HPLC papers), so they can be read locally.
2. Sturm: extract retention at 80 °C. Compare the G→A / G→C exchange series with the per-base EDA (sign and order of
   the shift). Fit per-base contributions from the composition where the design allows.
3. Shen 2018: if coordinates and radii are present, reproduce 2–3 complexes (e.g. ClF···NH₃, Br₂···pyridine) in gas
   and water as the external GKS-EDA(sol) validation.
