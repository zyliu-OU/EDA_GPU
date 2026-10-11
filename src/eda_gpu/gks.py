"""Experimental complete restricted/unrestricted GKS-EDA on CPU or NVIDIA GPU.

Su & Li 2009 (10.1063/1.3159673), Eqs. 6-12, evaluated on KS
orbitals as specified by Su, Tang & Wu 2020 (10.1002/wcms.1460),
Eqs. 2-5,8-9,11. External matched-settings validation remains pending.
"""
from dataclasses import asdict
from pathlib import Path
import time
import numpy as np
from pyscf import dft
from .engine import make_mf, host_array
from .model import xyz_text, electron_check


class Algebra:
    """Backend arrays; GPU requests never select CPU J/K or SCF."""
    def __init__(self, backend):
        self.backend = backend
        if backend == 'gpu':
            try:
                import cupy
                if cupy.cuda.runtime.getDeviceCount() < 1:
                    raise RuntimeError('no CUDA device')
                cupy.zeros(1).sum().get()
                self.xp = cupy
            except Exception as exc:
                raise RuntimeError(f'requested GPU EDA backend unavailable: {exc}; CPU fallback is disabled') from exc
        elif backend == 'cpu':
            self.xp = np
        else:
            raise ValueError('backend must be cpu or gpu')

    def scalar(self, value):
        return float(value.item() if hasattr(value, 'item') else value)

    def trace(self, a, b):
        return self.scalar(self.xp.einsum('ij,ji->', a, b))

    def synchronize(self):
        if self.backend == 'gpu':
            self.xp.cuda.get_current_stream().synchronize()


def hf_like(mf, density, algebra):
    """Full exact exchange, even when SCF uses a hybrid fraction."""
    xp = algebra.xp
    density = xp.asarray(density)
    j, k = mf.get_jk(mf.mol, density, hermi=1)
    h = xp.asarray(mf.get_hcore())
    total = density.sum(axis=0) if density.ndim == 3 else density
    one = algebra.trace(total, h)
    coulomb = .5 * algebra.trace(total, j.sum(axis=0) if density.ndim == 3 else j)
    exchange = (-.5 * sum(algebra.trace(density[i], k[i]) for i in range(2))
                if density.ndim == 3 else -.25 * algebra.trace(density, k))
    nuclear = float(mf.mol.energy_nuc())
    return dict(one_electron_hartree=one, coulomb_hartree=coulomb,
                exact_exchange_hartree=exchange, nuclear_hartree=nuclear,
                hf_like_total_hartree=one+coulomb+exchange+nuclear)


def density_checks(density, overlap, electrons, algebra):
    xp = algebra.xp
    if xp.iscomplexobj(density) or not bool(xp.isfinite(density).all()):
        raise RuntimeError('density must be finite and real')
    count = algebra.trace(density, overlap)
    error = algebra.scalar(xp.max(xp.abs(density @ overlap @ density - 2*density)))
    symmetry = algebra.scalar(xp.max(xp.abs(density-density.T)))
    if abs(count-electrons) > 1e-7 or error > 1e-7 or symmetry > 1e-10:
        raise RuntimeError('invalid restricted density normalization, idempotency or symmetry')
    return dict(electron_count=count, restricted_idempotency_max_error=error,
                density_symmetry_max_error=symmetry)


def orthonormalized_density(occupied, overlap, algebra, minimum_eigenvalue=1e-8):
    """Normalized antisymmetrized product: 2 C (C^T S C)^-1 C^T.

    Fail on a dependent occupied span; never drop electrons or regularize it.
    Symmetric orthonormalization leaves this same determinant density invariant.
    """
    xp = algebra.xp
    if not np.isfinite(minimum_eigenvalue) or minimum_eigenvalue <= 0:
        raise ValueError('occupied-overlap threshold must be finite and positive')
    occupied, overlap = xp.asarray(occupied), xp.asarray(overlap)
    if xp.iscomplexobj(occupied) or not bool(xp.isfinite(occupied).all()):
        raise ValueError('occupied orbitals must be finite and real')
    metric = occupied.T @ overlap @ occupied
    metric = .5*(metric+metric.T)
    eigenvalues, vectors = xp.linalg.eigh(metric)
    smallest = algebra.scalar(eigenvalues[0])
    if smallest <= minimum_eigenvalue:
        raise ValueError(f'fragment occupied spaces are linearly dependent or ill-conditioned: minimum overlap eigenvalue {smallest:.3e}')
    orthogonal = occupied @ (vectors / xp.sqrt(eigenvalues)[None,:])
    density = 2 * orthogonal @ orthogonal.T
    checks = density_checks(density, overlap, 2*occupied.shape[1], algebra)
    checks.update(occupied_overlap_min_eigenvalue=smallest,
                  occupied_overlap_condition_number=algebra.scalar(eigenvalues[-1])/smallest,
                  occupied_overlap_threshold=minimum_eigenvalue)
    return density, orthogonal, checks


def reference_state(symbols, coords, metadata, settings, algebra):
    mf = make_mf(symbols, coords, metadata, settings)
    if mf.mol.has_ecp() or mf.mol.pseudo:
        raise ValueError('EDA supports all-electron bases only')
    energy = float(mf.kernel())
    if not mf.converged or not np.isfinite(energy):
        raise RuntimeError('EDA reference SCF did not converge to a finite energy')
    xp = algebra.xp
    density = xp.asarray(mf.make_rdm1())
    overlap = xp.asarray(mf.get_ovlp())
    if density.ndim == 3:
        occupied = tuple(xp.asarray(mf.mo_coeff)[i][:,xp.asarray(mf.mo_occ)[i]>0] for i in range(2))
        checks = spin_density_checks(density, overlap, mf.mol.nelec, algebra)
        ss, multiplicity = mf.spin_square()
        checks.update(spin_square=float(ss), spin_multiplicity=float(multiplicity))
    else:
        occupied = xp.asarray(mf.mo_coeff)[:,xp.asarray(mf.mo_occ)>0]
        checks = density_checks(density, overlap, mf.mol.nelectron, algebra)
    parts = hf_like(mf, density, algebra)
    summary = dict(parts, **checks, energy_hartree=energy, scf_converged=True,
                   occupied_orbitals=[c.shape[1] for c in occupied] if isinstance(occupied,tuple) else occupied.shape[1],
                   generalized_correlation_hartree=energy-parts['hf_like_total_hartree'],
                   libxc_definition=repr(dft.libxc.parse_xc(settings.xc)))
    return mf, density, occupied, summary


def assert_common_basis(complex_mf, fragment_mf):
    a,b = complex_mf.mol, fragment_mf.mol
    same = (a.cart == b.cart and a.nao == b.nao and a.nbas == b.nbas
            and np.array_equal(a._bas[:,:5], b._bas[:,:5]))
    same = same and all(np.array_equal(a.bas_exp(n), b.bas_exp(n))
                        and np.array_equal(a.bas_ctr_coeff(n), b.bas_ctr_coeff(n))
                        for n in range(a.nbas))
    if not same or not np.allclose(a.intor_symmetric('int1e_ovlp'), b.intor_symmetric('int1e_ovlp'), atol=1e-12, rtol=0):
        raise RuntimeError('fragment AO basis does not match the complex basis')


def decompose_states(complex_state, fragment_states, metadata, algebra):
    """Evaluate all intermediate energies, not a residual component fit."""
    xp = algebra.xp
    mf, density, _, summary = complex_state
    (_,da,ca,sa),(_,db,cb,sb) = fragment_states
    unrestricted = density.ndim == 3
    if unrestricted:
        da,ca = spin_form(da,ca,algebra)
        db,cb = spin_form(db,cb,algebra)
    # Direct-product and AS formula densities are intentionally NOT idempotent.
    summed = da+db
    as_energy = hf_like(mf, summed, algebra)
    overlap = xp.asarray(mf.get_ovlp())
    if unrestricted:
        spin_results = [orthonormalized_density(xp.concatenate((ca[i],cb[i]),axis=1), overlap, algebra) for i in range(2)]
        asn_density = xp.stack([r[0]/2 for r in spin_results])
        asn_orbitals = tuple(r[1] for r in spin_results)
        checks = dict(spin_channels=[r[2] for r in spin_results], **spin_density_checks(asn_density,overlap,mf.mol.nelec,algebra))
    else:
        asn_density, asn_orbitals, checks = orthonormalized_density(xp.concatenate((ca,cb), axis=1), overlap, algebra)
    asn_energy = hf_like(mf, asn_density, algebra)
    # DP expression keeps intrafragment exchange but omits cross exchange.
    dp = (as_energy['one_electron_hartree'] + as_energy['coulomb_hartree']
          + sa['exact_exchange_hartree'] + sb['exact_exchange_hartree']
          + as_energy['nuclear_hartree'])
    monomer_sum = sa['hf_like_total_hartree']+sb['hf_like_total_hartree']
    # Independently contract classical interaction to check the DP difference.
    h = xp.asarray(mf.get_hcore())
    ha = xp.asarray(fragment_states[0][0].get_hcore())
    hb = xp.asarray(fragment_states[1][0].get_hcore())
    jb,kb = mf.get_jk(mf.mol, db, hermi=1)
    ja,_ = mf.get_jk(mf.mol, da, hermi=1)
    ta,tb = (da.sum(axis=0),db.sum(axis=0)) if unrestricted else (da,db)
    en_ab, en_ba = algebra.trace(ta,h-ha), algebra.trace(tb,h-hb)
    ee = algebra.trace(ta,jb.sum(axis=0) if unrestricted else jb)
    reciprocity = abs(ee-algebra.trace(tb,ja.sum(axis=0) if unrestricted else ja))
    nn = float(mf.mol.energy_nuc()-fragment_states[0][0].mol.energy_nuc()-fragment_states[1][0].mol.energy_nuc())
    electrostatic = ee+en_ab+en_ba+nn
    exchange = -sum(algebra.trace(da[i],kb[i]) for i in range(2)) if unrestricted else -.5*algebra.trace(da,kb)
    components = dict(electrostatic=electrostatic, exchange=exchange,
                      repulsion=asn_energy['hf_like_total_hartree']-as_energy['hf_like_total_hartree'],
                      polarization=summary['hf_like_total_hartree']-asn_energy['hf_like_total_hartree'],
                      generalized_correlation=summary['generalized_correlation_hartree']-sa['generalized_correlation_hartree']-sb['generalized_correlation_hartree'],
                      explicit_dispersion=0.0)
    total = summary['energy_hartree']-sa['energy_hartree']-sb['energy_hartree']
    residual = total-sum(components.values())
    dp_error = abs(electrostatic-(dp-monomer_sum))
    exchange_error = abs(exchange-(as_energy['hf_like_total_hartree']-dp))
    if max(dp_error, exchange_error, reciprocity, abs(residual)) > 1e-8:
        raise RuntimeError('EDA intermediate consistency or closure check failed')
    return dict(components_hartree=components, counterpoise_interaction_energy_hartree=total,
                total_closure_residual_hartree=residual,
                intermediate_states={'fragment_hf_like_sum_hartree':monomer_sum,
                    'direct_product_energy_hartree':dp,
                    'antisymmetrized_formula':as_energy,
                    'normalized_antisymmetrized_product':dict(asn_energy,**checks),
                    'direct_product_electrostatic_consistency_error_hartree':dp_error,
                    'exchange_consistency_error_hartree':exchange_error,
                    'coulomb_reciprocity_error_hartree':reciprocity},
                frozen_electrostatic_details=dict(electron_electron_hartree=ee,
                    electron_a_nuclear_b_hartree=en_ab,electron_b_nuclear_a_hartree=en_ba,
                    nuclear_nuclear_hartree=nn)), asn_density, asn_orbitals


def decompose(symbols, coords, metadata, settings, out):
    # CPU scope checker validates chemistry/numerics before requesting hardware.
    validate_gks_scope(symbols,coords,metadata,settings)
    algebra = Algebra(settings.backend)
    algebra.synchronize()
    start = time.perf_counter()
    complex_state = reference_state(symbols,coords,metadata,settings,algebra)
    fragments=[]
    summaries=[]
    for i,f in enumerate(metadata['fragments']):
        actual = set(f['atoms'])
        labels = [s if n in actual else 'ghost-'+s for n,s in enumerate(symbols)]
        state = reference_state(labels,coords,f,settings,algebra)
        assert_common_basis(complex_state[0],state[0])
        state[3].update(name=f.get('name',f'fragment_{i}'),atoms=f['atoms'],charge=f['charge'],spin=f['spin'],
                        ghost_atoms=[n for n in range(len(symbols)) if n not in actual])
        summaries.append(state[3]); fragments.append(state)
        Path(out,f'fragment_{i}.xyz').write_text(xyz_text([symbols[n] for n in f['atoms']],[coords[n] for n in f['atoms']]))
    result,asn_density,asn_orbitals = decompose_states(complex_state,fragments,metadata,algebra)
    algebra.synchronize()
    calculation_seconds=time.perf_counter()-start
    # Host transfer only at archive/report boundary; SCF and matrix operations stay on selected backend.
    np.savez_compressed(Path(out,'reference_densities.npz'),
        overlap=host_array(complex_state[0].get_ovlp()),complex_density=host_array(complex_state[1]),
        fragment_0_density=host_array(fragments[0][1]),fragment_1_density=host_array(fragments[1][1]),
        **archive_occupied(fragments),
        asn_density=host_array(asn_density),
        **({'asn_occupied_alpha':host_array(asn_orbitals[0]),'asn_occupied_beta':host_array(asn_orbitals[1])} if isinstance(asn_orbitals,tuple) else {'asn_occupied':host_array(asn_orbitals)}))
    result.update(stage=('experimental_complete_unrestricted_gks_eda' if metadata['spin'] else 'experimental_complete_closed_shell_gks_eda'),
        all_components_computed=True,full_eda_available=False,scientifically_validated=False,
        reference_component_validation='pending external matched-settings validation; complete experimental component calculation only',
        backend_used=settings.backend,scf_converged=True,settings=asdict(settings),
        bsse_policy='Boys-Bernardi counterpoise; whole complex basis on every fragment',
        fragment_reference_geometry='complex coordinates; no separate optimization or deformation term',
        reference_states=dict(complex=complex_state[3],fragments=summaries),
        component_sources={'orbital_intermediates':'10.1063/1.3159673 Eqs.6-12, evaluated with KS orbitals per 10.1002/wcms.1460 Eqs.2-5 and Section 2.2',
            'generalized_correlation':'10.1002/wcms.1460 Eq.9', 'explicit_dispersion':'10.1002/wcms.1460 Eq.11; absent for selected methods'},
        execution={'scf':settings.backend,'direct_jk':settings.backend,'orbital_matrix_algebra':settings.backend,
            'one_electron_integrals':'backend implementation; common-basis verification on CPU',
            'molecule_construction_and_output':'cpu',
            'gpu_jk_note':'GPU4PySCF direct J/K uses CPU/GPU hybrid internals; no all-GPU claim' if settings.backend=='gpu' else None},
        validation_blockers=['External XEDA component output with explicitly matched basis representation, functional variant, grid, integral settings and CP policy remains needed.',
                             'Authenticated original 2014 PDF and source-level confirmation of its step-number/exchange wording remain desirable.'],
        calculation_seconds=calculation_seconds,seconds=time.perf_counter()-start)
    return result


def spin_form(density, occupied, algebra):
    if density.ndim == 3:
        return density, occupied
    return algebra.xp.stack((density/2,density/2)), (occupied,occupied)


def spin_density_checks(density, overlap, electrons, algebra):
    checks=[]
    for i in range(2):
        c=density_checks(2*density[i],overlap,2*electrons[i],algebra)
        checks.append(dict(electron_count=c['electron_count']/2,
                           idempotency_max_error=c['restricted_idempotency_max_error']/2))
    return dict(electron_count=sum(c['electron_count'] for c in checks),spin_density_channels=checks)


def archive_occupied(fragments):
    result={}
    for i,state in enumerate(fragments):
        occupied=state[2]
        if isinstance(occupied,tuple):
            for label,c in zip(('alpha','beta'),occupied):
                result[f'fragment_{i}_occupied_{label}']=host_array(c)
        else:
            result[f'fragment_{i}_occupied']=host_array(occupied)
    return result


def validate_gks_scope(symbols,coords,metadata,settings):
    if settings.density_fitting or settings.dispersion != 'none':
        raise ValueError('complete EDA requires direct integrals and no explicit dispersion')
    if settings.xc.lower() not in ('b3lyp','b3lyp5','pbe0','m06-2x'):
        raise ValueError('complete EDA supports b3lyp, b3lyp5, pbe0 and m06-2x')
    if len(metadata['fragments']) != 2:
        raise ValueError('complete EDA requires exactly two fragments')
    if np.shape(coords) != (len(symbols),3) or not np.isfinite(coords).all():
        raise ValueError('invalid EDA coordinates')
    electron_check(symbols,metadata['charge'],metadata['spin'])
    atoms=[i for f in metadata['fragments'] for i in f['atoms']]
    if sorted(atoms) != list(range(len(symbols))):
        raise ValueError('fragment partition must cover every atom exactly once')
    if sum(f['charge'] for f in metadata['fragments']) != metadata['charge']:
        raise ValueError('fragment charges must sum to complex charge')
    if sum(f['spin'] for f in metadata['fragments']) != metadata['spin']:
        raise ValueError('only parallel spin coupling with additive fragment spins is supported')
    for f in metadata['fragments']:
        electron_check([symbols[i] for i in f['atoms']],f['charge'],f['spin'])
