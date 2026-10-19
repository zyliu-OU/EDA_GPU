"""PCM reaction field for experimental GKS-EDA(sol).

Su, Tang & Wu 2020 (10.1002/wcms.1460) Eqs. 12-13, Section 2.3. Monomers are
placed in the supermolecule cavity: the cavity and AO basis come from the
complex, while the surface nuclear potential uses only the monomer's real nuclei.
This avoids backend-specific ghost-atom cavity handling (PySCF assigns ghost
atoms an element-0 sphere; GPU4PySCF drops them) so CPU and GPU match.
"""
import numpy as np
from pyscf import gto

MODELS = {'cpcm': 'C-PCM', 'iefpcm': 'IEF-PCM'}


def surface_nuclear_potential(cavity_mol, surface, charges):
    """Nuclear potential on Gaussian surface charges, as in PySCF/GPU4PySCF PCM.build."""
    host = lambda a: np.asarray(a.get() if hasattr(a, 'get') else a)
    fakemol = gto.fakemol_for_charges(host(surface['grid_coords']), expnt=host(surface['charge_exp'])**2)
    fakemol_nuc = gto.fakemol_for_charges(cavity_mol.atom_coords())
    return np.asarray(charges, dtype=float) @ gto.mole.intor_cross(cavity_mol._add_suffix('int2c2e'), fakemol_nuc, fakemol)


def attach(mf, settings, cavity_mol=None):
    """Return mf with a PCM reaction field; cavity_mol puts a ghost-basis monomer in the complex cavity."""
    if settings.solvent == 'none':
        return mf
    if settings.solvent not in MODELS:
        raise ValueError(f'solvent must be none, {", ".join(MODELS)}')
    if not np.isfinite(settings.solvent_eps) or settings.solvent_eps <= 1:
        raise ValueError('solvent dielectric constant must exceed 1')
    mf = mf.PCM()
    pcm = mf.with_solvent
    if cavity_mol is not None:
        if cavity_mol.natm != mf.mol.natm or not np.allclose(cavity_mol.atom_coords(), mf.mol.atom_coords(), rtol=0, atol=1e-12):
            raise ValueError('cavity molecule must share atom order and coordinates with the monomer')
        if any(e.upper().startswith('GHOST') for e in cavity_mol.elements):
            raise ValueError('cavity molecule must contain real atoms only')
        pcm.reset(cavity_mol)
    pcm.method = MODELS[settings.solvent]
    pcm.eps = settings.solvent_eps
    pcm.lebedev_order = settings.solvent_lebedev_order
    pcm.vdw_scale = settings.solvent_vdw_scale
    pcm.build()
    if cavity_mol is not None:
        reproduced = surface_nuclear_potential(cavity_mol, pcm.surface, cavity_mol.atom_charges())
        library = np.asarray(pcm.v_grids_n.get() if hasattr(pcm.v_grids_n, 'get') else pcm.v_grids_n)
        if np.max(np.abs(reproduced - library)) > 1e-10:
            raise RuntimeError('surface nuclear potential does not reproduce the PCM library')
        monomer = surface_nuclear_potential(cavity_mol, pcm.surface, mf.mol.atom_charges())
        if hasattr(pcm.v_grids_n, 'get'):
            import cupy
            monomer = cupy.asarray(monomer)
        pcm.v_grids_n = monomer
    mf._eda_cavity = 'complex' if cavity_mol is not None else 'own'
    return mf


def describe(settings):
    if settings.solvent == 'none':
        return None
    return {'model': MODELS[settings.solvent], 'eps': settings.solvent_eps, 'lebedev_order': settings.solvent_lebedev_order,
            'radii': f'modified Bondi x {settings.solvent_vdw_scale} (PySCF/GPU4PySCF table), no probe',
            'surface_discretization': 'SWIG', 'energy': 'electrostatic solvation 0.5 q.V only; no cavitation/dispersion terms',
            'monomer_cavity': ('supermolecule cavity (complex atoms), monomer nuclei only' if settings.solvent_monomer_cavity == 'complex'
                               else 'diagnostic: monomer own cavity (GPU4PySCF drops ghost atoms)')}
