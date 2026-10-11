"""Source-backed, experimental building blocks for original GKS-EDA.

Review doi:10.1002/wcms.1460, Eqs. 2,8,9,11. No repulsion/polarization
intermediate is guessed, and a preparation result is never a full decomposition.
"""
from dataclasses import asdict
import time
import numpy as np
from pyscf import dft, scf
from .engine import make_mf
from .model import electron_check

REVIEW_DOI = "10.1002/wcms.1460"
BLOCKERS = [
    "Verify original antisymmetrized/intermediate-state exchange and normalization conventions against the original PDF/source.",
    "Implement repulsion and polarization only after those conventions are secured.",
    "Obtain reference output with matched functional variant, quadrature, SCF thresholds, integral approximation, and BSSE policy.",
]


def trace_product(a, b):
    return float(np.einsum("ij,ji->", a, b))


def validate_scope(symbols, coords, metadata, settings):
    if settings.backend != "cpu":
        raise ValueError("EDA preparation currently supports CPU only; GPU EDA is unavailable")
    if settings.density_fitting:
        raise ValueError("EDA reference preparation requires direct integrals; density fitting is not supported yet")
    if settings.dispersion != "none":
        raise ValueError("explicit dispersion is outside initial EDA preparation scope")
    if settings.xc.lower() not in ("b3lyp", "pbe0"):
        raise ValueError("EDA preparation supports only the explicit global hybrids b3lyp and pbe0")
    if len(metadata["fragments"]) != 2:
        raise ValueError("initial EDA preparation requires exactly two fragments")
    if not np.isfinite(coords).all() or np.shape(coords) != (len(symbols), 3):
        raise ValueError("invalid EDA coordinates")
    electron_check(symbols, metadata["charge"], metadata["spin"])
    if metadata["spin"] != 0 or any(f["spin"] != 0 for f in metadata["fragments"]):
        raise ValueError("EDA preparation requires a closed-shell complex and closed-shell fragments")
    atoms = [i for f in metadata["fragments"] for i in f["atoms"]]
    if sorted(atoms) != list(range(len(symbols))):
        raise ValueError("EDA fragment partition must cover every atom exactly once")
    if sum(f["charge"] for f in metadata["fragments"]) != metadata["charge"]:
        raise ValueError("fragment charges must sum to the complex charge")
    for f in metadata["fragments"]:
        electron_check([symbols[i] for i in f["atoms"]], f["charge"], f["spin"])


def state_energy(mol, dm):
    """HF-like functional evaluated at supplied KS density, not an HF SCF."""
    j, k = scf.hf.get_jk(mol, dm, hermi=1)
    h = mol.intor_symmetric("int1e_kin") + mol.intor_symmetric("int1e_nuc")
    one = trace_product(dm, h)
    coulomb = .5 * trace_product(dm, j)
    exchange = -.25 * trace_product(dm, k) # spin-summed restricted DM = 2 Cocc Cocc^T
    return {"one_electron_hartree": one, "coulomb_hartree": coulomb,
            "exact_exchange_hartree": exchange, "nuclear_hartree": float(mol.energy_nuc()),
            "hf_like_total_hartree": one + coulomb + exchange + mol.energy_nuc()}


def calculate_state(symbols, coords, metadata, settings):
    mf = make_mf(symbols, coords, metadata, settings)
    if mf.mol.has_ecp() or mf.mol.pseudo:
        raise ValueError("EDA reference preparation supports all-electron bases only")
    energy = float(mf._eda_pacer.run(mf.kernel))
    if not mf.converged:
        raise RuntimeError("EDA reference SCF did not converge")
    dm = np.asarray(mf.make_rdm1())
    if np.iscomplexobj(dm) or not np.isfinite(dm).all() or not np.isfinite(energy):
        raise RuntimeError("EDA reference requires finite real restricted densities")
    overlap = mf.get_ovlp()
    count = trace_product(dm, overlap)
    idempotency = float(np.max(np.abs(dm @ overlap @ dm - 2*dm)))
    if abs(count - mf.mol.nelectron) > 1e-7 or idempotency > 1e-7:
        raise RuntimeError("invalid closed-shell reference density normalization")
    parts = state_energy(mf.mol, dm)
    # Review Eq.9: generalized correlation = mixed-XC contribution - full HF exchange.
    # Taking the difference of independently evaluated totals avoids guessing the
    # normalization of a composite hybrid's separately weighted exchange functional.
    gcorr = energy - parts["hf_like_total_hartree"]
    summary = dict(parts, energy_hartree=energy, scf_converged=True,
                   electron_count=count, occupied_orbitals=mf.mol.nelectron//2,
                   restricted_idempotency_max_error=idempotency,
                   generalized_correlation_hartree=gcorr,
                   libxc_definition=repr(dft.libxc.parse_xc(settings.xc)))
    return mf, dm, summary


def frozen_electrostatics(mol, da, db, atom_indices_a, atom_indices_b):
    """Classical frozen-density Coulomb interaction, review Eq.2.

    Nuclear potentials include only real nuclei; both DMs are in the common
    complex AO basis. Direct J evaluation avoids storing four-index integrals.
    """
    coords = mol.atom_coords()
    z = mol.atom_charges()
    va, vb = np.zeros_like(da), np.zeros_like(db)
    for indices, potential in ((atom_indices_a, va), (atom_indices_b, vb)):
        for atom in indices:
            with mol.with_rinv_origin(coords[atom]):
                potential -= z[atom] * mol.intor_symmetric("int1e_rinv")
    ja, ka = scf.hf.get_jk(mol, da, hermi=1)
    jb, kb = scf.hf.get_jk(mol, db, hermi=1)
    ab = trace_product(da, jb)
    ba = trace_product(db, ja)
    en_ab = trace_product(da, vb)
    en_ba = trace_product(db, va)
    nn = 0.0
    for a in atom_indices_a:
        for b in atom_indices_b:
            distance = np.linalg.norm(coords[a]-coords[b])
            if distance <= 1e-12:
                raise ValueError("coincident nuclei across fragments")
            nn += z[a]*z[b]/distance
    bare_exchange = -.5*trace_product(da, kb)
    if abs(ab-ba) > 1e-8:
        raise RuntimeError("cross-Coulomb reciprocity failed")
    return {"electrostatic_hartree": ab + en_ab + en_ba + nn,
            "electron_electron_hartree": ab, "electron_a_nuclear_b_hartree": en_ab,
            "electron_b_nuclear_a_hartree": en_ba, "nuclear_nuclear_hartree": float(nn),
            "coulomb_reciprocity_error_hartree": abs(ab-ba),
            "bare_cross_hf_exchange_hartree": bare_exchange,
            "bare_exchange_note": "Diagnostic between nonorthogonal fragment orbitals; not reported as the original GKS-EDA exchange component."}


def prepare(symbols, coords, metadata, settings, out):
    from pathlib import Path
    from .model import xyz_text
    start = time.perf_counter()
    validate_scope(symbols, coords, metadata, settings)
    complex_mf, complex_dm, complex_summary = calculate_state(symbols, coords, metadata, settings)
    fragment_states, summaries = [], []
    for i, f in enumerate(metadata["fragments"]):
        actual = set(f["atoms"])
        # Retain all atoms/basis centers in original order; ghost atoms carry no nuclear charge.
        ghost_symbols = [s if n in actual else "ghost-"+s for n,s in enumerate(symbols)]
        fm, dm, summary = calculate_state(ghost_symbols, coords, f, settings)
        if fm.mol.ao_labels() != complex_mf.mol.ao_labels():
            # PySCF labels ghost symbols differently, so check shell structure and overlap instead.
            same_shells = fm.mol.nbas == complex_mf.mol.nbas and np.array_equal(fm.mol._bas[:, :5], complex_mf.mol._bas[:, :5])
            same_primitives = same_shells and all(
                np.array_equal(fm.mol.bas_exp(n), complex_mf.mol.bas_exp(n)) and
                np.array_equal(fm.mol.bas_ctr_coeff(n), complex_mf.mol.bas_ctr_coeff(n))
                for n in range(fm.mol.nbas))
            if not same_primitives or not np.allclose(fm.get_ovlp(), complex_mf.get_ovlp(), atol=1e-12, rtol=0):
                raise RuntimeError("fragment AO basis does not match complex AO basis")
        fragment_states.append((fm,dm))
        summary.update(name=f.get("name",f"fragment_{i}"), atoms=f["atoms"], charge=f["charge"], spin=f["spin"],
                       ghost_atoms=[n for n in range(len(symbols)) if n not in actual],
                       reference_geometry="unrelaxed fragment coordinates inherited from complex")
        summaries.append(summary)
        Path(out, f"fragment_{i}.xyz").write_text(xyz_text([symbols[n] for n in f["atoms"]], [coords[n] for n in f["atoms"]]))
    frozen = frozen_electrostatics(complex_mf.mol, fragment_states[0][1], fragment_states[1][1],
                                  metadata["fragments"][0]["atoms"], metadata["fragments"][1]["atoms"])
    gcorr = complex_summary["generalized_correlation_hartree"] - sum(f["generalized_correlation_hartree"] for f in summaries)
    total = complex_summary["energy_hartree"] - sum(f["energy_hartree"] for f in summaries)
    # Save reference densities for reproducible intermediate-state development.
    np.savez_compressed(Path(out,"reference_densities.npz"), overlap=complex_mf.get_ovlp(), complex_density=complex_dm,
                        fragment_0_density=fragment_states[0][1], fragment_1_density=fragment_states[1][1])
    return {"stage": "experimental_source_backed_preparation", "full_eda_available": False,
            "reference_component_validation": "not established; published-reference comparisons are exploratory",
            "backend_used": "cpu", "scf_converged": True, "settings": asdict(settings),
            "bsse_policy": "Boys-Bernardi counterpoise; complex basis on every fragment via ghost nuclei",
            "fragment_reference_geometry": "complex coordinates; no separate fragment optimization or deformation term",
            "reference_states": {"complex": complex_summary, "fragments": summaries},
            "counterpoise_interaction_energy_hartree": total,
            "partial_components_hartree": {"electrostatic": frozen["electrostatic_hartree"], "exchange": None,
                                          "repulsion": None, "polarization": None, "generalized_correlation": gcorr,
                                          "explicit_dispersion": 0.0},
            "frozen_electrostatic_details": frozen, "total_closure_residual_hartree": None,
            "component_sources": {"electrostatic": REVIEW_DOI+" Eq.2; KS substitution in Section 2.2",
                                  "generalized_correlation": REVIEW_DOI+" Eq.9; mixed-XC minus full exact exchange",
                                  "explicit_dispersion": REVIEW_DOI+" Eq.11; absent for selected B3LYP/PBE0"},
            "blockers": BLOCKERS, "seconds": time.perf_counter()-start}
