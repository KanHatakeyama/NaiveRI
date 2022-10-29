from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_mol(smiles, atom_index=True):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    if atom_index:
        return mol_with_atom_index(mol)
    else:
        return mol


def mol_with_atom_index(mol):
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',
                                        str(mol.GetAtomWithIdx(idx).GetIdx()))
    return mol


def get_bond_length(conf, i, j):
    return Chem.rdMolTransforms.GetBondLength(conf, i, j)


def get_bond_list(mol):
    bond_list = []
    for bond in mol.GetBonds():
        bond_list.append((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
    return bond_list
