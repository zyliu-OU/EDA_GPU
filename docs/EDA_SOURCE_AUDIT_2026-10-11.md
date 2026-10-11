# EDA source audit — 2026-10-11

- Owner-supplied review: present, 22 pages, fully read; equations visually checked on
  pages 3–5. SHA-256 and equation mapping are in GKS_EDA_SPEC.md.
- Original 2014 paper: publisher full text unavailable. A university bibliography
  identifies the original download; original pcss HTTP/HTTPS and current pcoss HTTPS
  paths return 404. Public transcription found but equations are interleaved, download
  inaccessible, and authenticity/normalization details cannot be verified from a PDF.
- New official XEDA II manual: retrieved locally, 52 pages. The relevant input,
  closed-shell/global-hybrid CO2 exercise and reference-result sections were inspected.
  Unlike the earlier web manual, it supplies separate EX/REP/POL/CORR values.
- GAMESS XEDA manual: retrieved locally; SUPBAS requests Boys–Bernardi CP and defaults
  true in that interface. XEDA II's implicit CP defaults need explicit confirmation.
- XEDA execution: no local executable or authenticated chemistry XEDA cloud session
  established. The unrelated PyPI/RTL package named xeda is not the chemical program.
- No contact messages or account registration were performed. No PDFs or full-text
  transcriptions are staged for public redistribution.

Source URLs:
- https://doi.org/10.1002/wcms.1460
- https://doi.org/10.1021/jp500405s
- https://pcoss.xmu.edu.cn/2014-paperlist.pdf
- https://xmvb.xmu.edu.cn/docs/xeda/
- https://xmvb.xmu.edu.cn/docs/xeda/_downloads/415e0f39b016129f5cf2eb24dfe28b4f/XEDA_II_en.pdf
- https://xmvb.xmu.edu.cn/docs/xeda/_downloads/222e6e0bd8ccfe3ce8fe39d17fc3eef3/XEDA_gms.pdf
- https://doczz.net/doc/8055485/energy-decomposition-scheme-based-on-the-generalized-kohn

## Complete experimental CPU/GPU stage

The Su–Li 2009 institutional record and author-uploaded article were found. The
institutional PDF download returns 403; author-hosted article text supplies explicit
Eqs. 6–12. Combined with the supplied review Section 2.2, these distinguish the
un-normalized AS energy expression from the inverse-metric normalized determinant,
resolving the earlier implementation ambiguity for an experimental construction.
The original 2014 PDF remains unavailable, and its transcribed step labels still
require confirmation. Official download navigation currently leads to a 404 page;
no local chemical XEDA reference executable or authenticated account was established.

Additional sources:
- https://doi.org/10.1063/1.3159673
- https://digitalcommons.unl.edu/chemfacpub/22/
- https://www.researchgate.net/publication/26655483_Energy_Decomposition_Analysis_of_Covalent_Bonds_and_Intermolecular_Interactions
- https://raw.githubusercontent.com/pyscf/gpu4pyscf/v1.8.1/gpu4pyscf/scf/hf.py

## Tryptamine–water source audit

Read the user-supplied 30-page accepted manuscript (29 printed pages plus cover),
DOI 10.1039/C7CP03491A, and all four pages of `Unravelling hydrogen bonding_SI.pdf`.
SI page 2 was also inspected visually; it contains conformer/AIM figures, not Cartesian
coordinates. PDF attachment inspection reports zero embedded files. Manuscript SHA-256:
`ea00ca2f48c7860440818ac54cd1628f7f2fba16f465b5035912ca1676585e3a`;
SI SHA-256: `7cb9c4d7d8772abd1e698038fa500e45da27ae9346bd60a66732ea93c29f1ddb`.
The alternative downloaded SI filename has the same hash. Table 3 reference numbers
and Table 1 distances are attributed separately; no full paper is published.

Methods are M06-2X/U-M06-2X, aug-cc-pVDZ and counterpoise. The supplied manuscript
contains inconsistencies (e.g. Table 1 vs introductory third H-bond distance and
individual-H-bond Table 4 vs prose); the benchmark uses Tables 1 and 3 explicitly.
The four-page SI supplies charges and qualitative figures but no complete input
geometries. RDKit/SI-pattern starting coordinates therefore carry constructed provenance,
and comparison with Table 3 is exploratory. Original geometries and program numerical
settings are still needed to establish a matched reference calculation.
