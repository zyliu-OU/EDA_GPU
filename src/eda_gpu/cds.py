"""SMD non-electrostatic (cavitation-dispersion-solvent structure, CDS) solvation term.

Marenich, Cramer & Truhlar 2009 (10.1021/jp810292n) Eqs. 21-27, via the PySCF/MNSOL routine. G_CDS depends
only on element types and coordinates (atomic surface tensions on SASA with SMD radii + 0.4 A solvent
radius), so the association term needs no SCF: dG_CDS = G_CDS(AB) - sum G_CDS(fragment), all at the
complex geometry. It complements C-PCM GKS-EDA(sol), whose desolvation is electrostatic only. It is not
part of the seven-term EDA closure: SMD CDS was parametrized with SMD electrostatics and SMD radii, not
C-PCM with modified Bondi radii, so the sum is a mixed-model estimate.
"""
import numpy as np
from pyscf import gto
from pyscf.solvent import smd

HARTREE_KCAL = 627.509474
SOLVENTS = ('water', 'acetonitrile')


def cds_kcal(symbols, coords, solvent='water'):
    """G_CDS (kcal/mol) of atoms at coords (A)."""
    if solvent not in smd.solvent_db:
        raise ValueError(f'unknown SMD solvent {solvent}')
    mol = gto.M(atom=[(s, tuple(map(float, r))) for s, r in zip(symbols, coords)], basis={}, spin=None, verbose=0)
    return smd.get_cds_legacy(smd.SMD(mol, solvent=solvent))[0] * HARTREE_KCAL


def association_cds(symbols, coords, groups, solvent='water'):
    """dG_CDS (kcal/mol) for binding atom-index groups, each taken at the complex geometry."""
    coords = np.asarray(coords, float)
    if sorted(i for g in groups for i in g) != sorted(set(i for g in groups for i in g)):
        raise ValueError('groups must not overlap')
    atoms = [i for g in groups for i in g]
    whole = cds_kcal([symbols[i] for i in atoms], coords[atoms], solvent)
    return whole - sum(cds_kcal([symbols[i] for i in g], coords[g], solvent) for g in groups)
