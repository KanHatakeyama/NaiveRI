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
    symbol_id_list = get_symbol_atom_ids(mol)

    bond_list = []
    symbol_bond_list = []
    for bond in mol.GetBonds():
        bid1 = bond.GetBeginAtomIdx()
        bid2 = bond.GetEndAtomIdx()

        # ignore symbol atoms
        symbol_flag = False
        for atom_id in (bid1, bid2):
            if atom_id in symbol_id_list:
                symbol_flag = True

        if symbol_flag:
            symbol_bond_list.append((bid1, bid2))

        bond_list.append((bid1, bid2))

    # add repeating unit bond
    repeat_bond = find_repeat_bond(mol, bond_list)
    bond_list.append(repeat_bond)

    for b in symbol_bond_list:
        bond_list.remove(b)
    return bond_list, repeat_bond


def find_repeat_bond(mol, bond_list):
    symbol_bonds = get_symbol_atom_ids(mol)

    repeat_bond = []
    for bond in bond_list:
        for b in bond:
            if b in symbol_bonds:
                repeat_bond.extend(bond)

    repeat_bond.remove(symbol_bonds[0])
    repeat_bond.remove(symbol_bonds[1])
    return repeat_bond


def get_symbol_atom_ids(mol, unit_symbol="*", inverse_mode=False):
    symbol_list = []
    normal_atom_list = []
    atoms = mol.GetNumAtoms()
    for idx in range(atoms):
        symbol = mol.GetAtomWithIdx(idx).GetSymbol()
        if symbol == unit_symbol:
            symbol_list.append(idx)
        else:
            if inverse_mode:
                normal_atom_list.append(idx)

    if inverse_mode:
        return normal_atom_list

    return symbol_list


def get_atom_ids(mol, unit_symbol="*"):
    return get_symbol_atom_ids(mol, unit_symbol, inverse_mode=True)
