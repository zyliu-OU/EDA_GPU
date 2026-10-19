#!/usr/bin/env python
"""Ion-pair reagent series for the HPLC model study (optional rdkit==2026.9.1).

Alkylammonium cations of common IP-RP reagents (acetate/HFIP counter-anions are not modelled), built
like TEAH+ in build_hplc_models.py: lowest-MMFF94 ETKDG conformer, then every C-C bond between two
non-terminal heavy atoms is set anti (180 deg) and re-minimized, so alkyl chains are extended
seeds rather than gas-phase folds. Each reagent N-H is then placed on the residue phosphate O
(ion pair) and on the C8 proxy silanol O (surface pair) with build_hplc_complexes.py, matching
the TEAH+ stage 7.1/7.2 seeds; surface pairs get a second, chain-contact seed (--select contact)
because alkyl-chain contact with the C8 ligand is the interaction of interest there. Seeds only; geometries are refined by run_hplc_stage.py.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms as T
from build_hplc_models import OUT, embed, write

REAGENTS = {  # name: (identity, SMILES of the cation, buffer it represents)
    'hexylammonium': ('hexylammonium (HA+)', 'CCCCCC[NH3+]', 'hexylammonium acetate (HAA)'),
    'dibutylammonium': ('dibutylammonium (DBA+)', 'CCCC[NH2+]CCCC', 'dibutylammonium acetate (DBAA)'),
    'tributylammonium': ('tributylammonium (TBA+)', 'CCCC[NH+](CCCC)CCCC', 'tributylammonium acetate (TBuAA)'),
    'dipeah': ('N,N-diisopropylethylammonium (DIPEA-H+)', 'CC[NH+](C(C)C)C(C)C', 'DIPEA/HFIP or DIPEA acetate'),
}
PHOSPHATE_O = {'A': (17, 18), 'C': (16, 17), 'G': (13, 14), 'U': (16, 17)}  # as in build_hplc_ternary.py
SILANOL_O = (13, 14)  # c8_silane_proxy, as in the TEAH+ surface seeds


def chain_dihedrals(mol):
    """Heavy-atom dihedrals across C-C bonds whose atoms both carry >= 2 heavy neighbours."""
    heavy = lambda a: [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
    out = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtom(), b.GetEndAtom()
        if (i.GetSymbol(), j.GetSymbol()) != ('C', 'C') or b.IsInRing(): continue
        hi, hj = [n for n in heavy(i) if n != j.GetIdx()], [n for n in heavy(j) if n != i.GetIdx()]
        if hi and hj: out.append((hi[0], i.GetIdx(), j.GetIdx(), hj[0]))
    return out


def build(name, identity, smiles, buffer):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ids, energies = embed(mol, 200, 2026)
    cid = ids[int(np.argmin(energies))]
    conf, dihedrals = mol.GetConformer(cid), chain_dihedrals(mol)
    for d in dihedrals:
        T.SetDihedralDeg(conf, *d, 180.0)
    field = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=cid)
    field.Minimize(maxIts=5000)
    final = [T.GetDihedralDeg(conf, *d) for d in dihedrals]
    if not all(abs(x) >= 150 for x in final): raise RuntimeError(f'{name}: chains did not stay anti: {final}')
    n = next(a for a in mol.GetAtoms() if a.GetSymbol() == 'N')
    donor_h = [a.GetIdx() for a in n.GetNeighbors() if a.GetSymbol() == 'H']
    write(name, mol, cid, {'identity': identity, 'protonation': f'protonated amine, {len(donor_h)} N-H donor(s)',
                           'stereochemistry': 'no stereocenters', 'ion_pair_buffer': buffer,
                           'conformer_selection': {'criterion': 'lowest MMFF94 ETKDG conformer, C-C chain dihedrals set anti and re-minimized',
                                                   'chain_dihedrals_deg': final, 'mmff94_kcal_mol': field.CalcEnergy(), 'embedded': len(ids)},
                           'donor_h': donor_h})
    print(f'{name}: {mol.GetNumAtoms()} atoms, N-H {donor_h}, chain dihedrals {[round(x) for x in final]}')
    return donor_h[0]  # N-H hydrogens of an ammonium are equivalent up to conformation


def place(host, guest, acceptors, donor_h, output, select='clearance'):
    subprocess.run([sys.executable, str(Path(__file__).with_name('build_hplc_complexes.py')), '--host', host, '--guest', guest,
                    '--acceptor', *map(str, acceptors), '--donor-h', str(donor_h), '--select', select, '--output', output], check=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reagents', nargs='+', default=list(REAGENTS), choices=list(REAGENTS))
    p.add_argument('--bases', nargs='+', default=['A'], choices=list(PHOSPHATE_O), help='residues for ion-pair seeds')
    a = p.parse_args()
    for name in a.reagents:
        h = build(name, *REAGENTS[name])
        for base in a.bases:
            place(f'{OUT}/residue_{base}', f'{OUT}/{name}', PHOSPHATE_O[base], h, f'{OUT}/reagents/ionpair/residue_{base}_{name}')
        place(f'{OUT}/c8_silane_proxy', f'{OUT}/{name}', SILANOL_O, h, f'{OUT}/reagents/surface/c8_{name}')
        place(f'{OUT}/c8_silane_proxy', f'{OUT}/{name}', SILANOL_O, h, f'{OUT}/reagents/surface/c8_{name}_contact', 'contact')
        for o in SILANOL_O:  # drop a contact seed that is the same placement as the clearance seed
            same = [json.loads(Path(f'{OUT}/reagents/surface/c8_{name}{s}_{o}_0.json').read_text())['placement']['tilt_deg_azimuth_step_rotation_step']
                    for s in ('', '_contact')]
            if same[0] == same[1]:
                for ext in ('.xyz', '.json'): Path(f'{OUT}/reagents/surface/c8_{name}_contact_{o}_0{ext}').unlink()
                print(f'c8_{name}_contact_{o}_0: identical to the clearance seed, removed')


if __name__ == '__main__':
    main()
