# Patch basis sensitivity — aug-cc-pVDZ single points (2026-11-11)

GKS-EDA(sol), M06-2X/aug-cc-pVDZ, C-PCM water, direct integrals, at the def2-SVP geometries of the lowest-energy
placement per reagent. `scripts/sensitivity_hplc.py` variant `augdz`. 43–78 min each, closure ≤ 5e-9,
peak ~7.8 GB GPU / ~4.6 GB RAM (`mem.log`), GPU ≤ 65 °C.

| case | SVP | aDZ | Δ | ΔPol | ΔGC |
| --- | ---: | ---: | ---: | ---: | ---: |
| teah_2 | −5.33 | −7.60 | −2.27 | +1.11 | −3.91 |
| dipeah_0 | −5.49 | −7.35 | −1.86 | +0.83 | −3.41 |
| dibutylammonium_0 | −8.85 | −12.19 | −3.34 | +1.58 | −6.15 |
| hexylammonium_0 | −5.97 | −8.30 | −2.33 | +1.40 | −4.19 |
| tributylammonium_3 | −7.14 | −10.06 | −2.92 | +1.59 | −5.17 |

All more bound by 1.9–3.3 kcal/mol, almost all generalized correlation; ρ(SVP, aDZ) = 0.9; ρ vs published 0.7 → 0.6
(TEA/DIPEA swap, 0.25 apart). DBAA > TBuAA > HA persists.
