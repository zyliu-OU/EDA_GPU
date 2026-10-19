#!/usr/bin/env python
"""Build labelled monomer seeds for the ion-pair RP-HPLC model study (optional rdkit==2026.9.1).

Capped RNA residues: beta-D-ribonucleosides (PubChem isomeric SMILES) with a 5'-O methyl
phosphodiester (-1) and a 3'-O methyl cap, i.e. one phosphodiester per residue and identical
termini for A/G/C/U; neutral canonical bases. Conformer: lowest MMFF94 energy among ETKDG
conformers with an anti glycosidic angle (A-form chi ~ -160 deg), C3'-endo-like delta (55-105 deg;
A-form ~82), trans beta (140-220) and gauche+ gamma (30-90; A-form ~55), which keeps the
5'-phosphate off the base as in A-form RNA, plus no intramolecular phosphate O...H-N/O contact
below 2.5 A (gas-phase MMFF otherwise folds the anion onto base or sugar donors). TEAH+: protonated triethylamine. Stationary-phase proxy:
octyl- (and butyl-) dimethylsilyl ether of Si(OH)3, a short alkylsilane stand-in for a bonded
C18 ligand with residual silanols; anchor Si(OH)3 heavy atoms are frozen in optimizations.
Seeds are force-field geometries, not quantum-optimized structures.
"""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
from eda_gpu.model import xyz_text

OUT = Path(__file__).resolve().parents[1] / 'examples/hplc'
NUCLEOSIDES = {  # PubChem isomeric SMILES
    'A': ('adenosine', 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N', 'purine'),
    'G': ('guanosine', 'C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=C(NC2=O)N', 'purine'),
    'C': ('cytidine', 'C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O', 'pyrimidine'),
    'U': ('uridine', 'C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O', 'pyrimidine'),
}
OTHERS = {
    'teah': ('triethylammonium (TEAH+)', 'CC[NH+](CC)CC', 'protonated tertiary amine; N-H donor', None),
    'c8_silane_proxy': ('octyldimethylsilyl-O-Si(OH)3 (C18 stationary-phase proxy)', 'CCCCCCCC[Si](C)(C)O[Si](O)(O)O',
                        'neutral; three residual silanol OH', '[Si](O)(O)(O)O[Si]'),
    'c4_silane_proxy': ('butyldimethylsilyl-O-Si(OH)3 (C18 stationary-phase proxy)', 'CCCC[Si](C)(C)O[Si](O)(O)O',
                        'neutral; three residual silanol OH', '[Si](O)(O)(O)O[Si]'),
}
RESIDUE_PATTERN = Chem.MolFromSmarts('[CH3]O[P](=O)([O-])O[CH2][CH]1O[CH]([n,N])[CH](O)[CH]1O[CH3]')


def capped(smiles):
    for ring in '23':
        old = f'(O{ring})CO)O)'
        if old in smiles:
            return smiles.replace(old, f'(O{ring})COP(=O)([O-])OC)OC)')
    raise ValueError('unexpected nucleoside SMILES')


def embed(mol, n, seed):
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n, randomSeed=seed))
    if not ids: raise RuntimeError('embedding failed')
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None: raise RuntimeError('MMFF94 parameters unavailable')
    energies = [e for _, e in AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=2000)]
    return ids, energies


def residue_torsions(mol, conf, kind):
    m = mol.GetSubstructMatch(RESIDUE_PATTERN)
    if not m: raise RuntimeError('capped residue pattern not found')
    c5, c4, o4, c1, nglyc, c3, o3 = m[6], m[7], m[8], m[9], m[10], m[13], m[14]
    ring = [a.GetIdx() for a in mol.GetAtomWithIdx(nglyc).GetNeighbors() if a.GetIdx() != c1]
    # chi: O4'-C1'-N9-C4 (purine) / O4'-C1'-N1-C2 (pyrimidine): the ring neighbour bonded to N of the other ring/C=O
    if kind == 'purine':
        ref = next(i for i in ring if mol.GetAtomWithIdx(i).IsInRingSize(5) and mol.GetAtomWithIdx(i).IsInRingSize(6))
    else:
        ref = next(i for i in ring if any(b.GetBondType() == Chem.BondType.DOUBLE and b.GetOtherAtom(mol.GetAtomWithIdx(i)).GetSymbol() == 'O'
                                          for b in mol.GetAtomWithIdx(i).GetBonds()))
    chi = T.GetDihedralDeg(conf, o4, c1, nglyc, ref)
    delta = T.GetDihedralDeg(conf, c5, c4, c3, o3)
    beta = T.GetDihedralDeg(conf, m[2], m[5], c5, c4) % 360
    gamma = T.GetDihedralDeg(conf, m[5], c5, c4, c3)
    terminal = [m[3], m[4]]
    polar_h = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H' and a.GetNeighbors()[0].GetSymbol() in ('N', 'O')]
    xyz = conf.GetPositions()
    contact = min(float(np.linalg.norm(xyz[o] - xyz[h])) for o in terminal for h in polar_h)
    return chi, delta, beta, gamma, contact, {'C5p': c5, 'C4p': c4, 'O4p': o4, 'C1p': c1, 'N_glycosidic': nglyc, 'C3p': c3, 'O3p': o3}


def bonds(mol):
    return [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()] for b in mol.GetBonds()]


def write(name, mol, conf_id, meta):
    coords = mol.GetConformer(conf_id).GetPositions()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    meta.update(charge=Chem.GetFormalCharge(mol), spin=0, bonds=bonds(mol),
                fragments=[{'name': name, 'atoms': list(range(len(symbols))), 'charge': Chem.GetFormalCharge(mol), 'spin': 0}],
                smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)),
                geometry_provenance='RDKit ETKDG + MMFF94 seed (rdkit 2026.9.1); not quantum-optimized')
    (OUT / f'{name}.xyz').write_text(xyz_text(symbols, coords, f'{name}; RDKit/MMFF94 seed; angstrom'))
    (OUT / f'{name}.json').write_text(json.dumps(meta, indent=2) + '\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for base, (parent, smiles, kind) in NUCLEOSIDES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(capped(smiles)))
        ids, energies = embed(mol, 1000, 2026)
        ok = []
        for cid, e in zip(ids, energies):
            chi, delta, beta, gamma, contact, atoms = residue_torsions(mol, mol.GetConformer(cid), kind)
            if (chi <= -120 or chi >= 150) and 55 <= delta <= 105 and 140 <= beta <= 220 and 30 <= gamma <= 90 and contact >= 2.5:
                ok.append((e, cid, chi, delta, beta, gamma, contact, atoms))
        if not ok: raise RuntimeError(f'{base}: no anti/C3prime-endo-like conformer')
        e, cid, chi, delta, beta, gamma, contact, atoms = min(ok)
        stereo = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
        write(f'residue_{base}', mol, cid, {
            'identity': f'capped RNA residue {base}: 5prime-O-methylphosphodiester, 3prime-O-methyl {parent}',
            'protonation': 'neutral canonical base; phosphodiester monoanion (-1); no added counterion',
            'stereochemistry': f'beta-D-ribofuranosyl; CIP {stereo}',
            'conformer_selection': {'criterion': 'lowest MMFF94 with anti chi, delta 55-105, beta 140-220, gamma 30-90 deg, no intramolecular phosphate O...H-N/O < 2.5 A',
                                    'chi_deg': chi, 'delta_deg': delta, 'beta_deg': beta, 'gamma_deg': gamma,
                                    'min_intramolecular_phosphate_O_polar_H_angstrom': contact,
                                    'mmff94_kcal_mol': e, 'candidates': len(ok), 'embedded': len(ids)},
            'key_atoms': atoms})
        print(f'residue_{base}: chi {chi:.1f} delta {delta:.1f} beta {beta:.1f} gamma {gamma:.1f} from {len(ok)}/{len(ids)} conformers')
    for name, (identity, smiles, protonation, anchor) in OTHERS.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        ids, energies = embed(mol, 200, 2026)
        meta = {'identity': identity, 'protonation': protonation, 'stereochemistry': 'no stereocenters'}
        cid = ids[int(np.argmin(energies))]
        if anchor:
            # Bonded-phase proxy: lowest-energy conformer whose C...C-Si backbone is all anti (|dihedral| >= 150 deg).
            carbons = smiles.index('[Si]')  # proxy SMILES start with an unbranched run of single-letter C atoms
            chain = list(range(carbons + 1))
            if mol.GetAtomWithIdx(chain[-1]).GetSymbol() != 'Si': raise RuntimeError('unexpected proxy atom order')
            # Random embedding rarely yields a fully extended chain: set each backbone dihedral to 180 deg on the
            # lowest-energy conformer, re-minimize with MMFF94 (all-anti is a local minimum) and verify.
            conf = mol.GetConformer(cid)
            for k in range(len(chain) - 3):
                T.SetDihedralDeg(conf, *chain[k:k + 4], 180.0)
            field = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid)
            field.Minimize(maxIts=5000)
            dihedrals = [T.GetDihedralDeg(conf, *chain[k:k + 4]) for k in range(len(chain) - 3)]
            if not all(abs(d) >= 170 for d in dihedrals): raise RuntimeError(f'{name}: extended chain did not stay all-anti')
            meta['conformer_selection'] = {'criterion': 'lowest MMFF94 ETKDG conformer, backbone set all-anti and re-minimized',
                                           'backbone_dihedrals_deg': dihedrals, 'mmff94_kcal_mol': field.CalcEnergy(), 'embedded': len(ids)}
            meta['frozen_atoms'] = sorted(mol.GetSubstructMatch(Chem.MolFromSmarts(anchor)))
            meta['frozen_atoms_note'] = 'silanol Si, its three hydroxyl O, the bridging O and the ligand Si: tethered-surface anchor'
        write(name, mol, cid, meta)
        print(f'{name}: {mol.GetNumAtoms()} atoms')


if __name__ == '__main__':
    main()
