from .atom_dict import atom_van_dict
from .mol_utils import get_bond_length
import numpy as np
from .mol_utils import get_atom_ids, get_bond_list


def calc_mol_volume(mol):
    conf = mol.GetConformers()[0]
    atom_ids = get_atom_ids(mol)

    bond_list, repeat_bond = get_bond_list(mol)

    V = 0

    for target_atom_id in atom_ids:
        dV = calc_dVA(mol, conf, bond_list, target_atom_id, repeat_bond)
        V += dV
        # print("dV",target_atom_id,dV)

    return 0.602*V


def calc_s(mol, conf, bond, repeat_bond=[]):

    atom_id1 = bond[0]
    atom_id2 = bond[1]
    atom1 = mol.GetAtomWithIdx(atom_id1).GetSymbol()
    atom2 = mol.GetAtomWithIdx(atom_id2).GetSymbol()

    Ri = atom_van_dict[atom1]
    Rj = atom_van_dict[atom2]

    if set(repeat_bond) == set((atom_id1, atom_id2)):
        # const length for repeating bond
        di = 1.5
    else:
        di = get_bond_length(conf, atom_id1, atom_id2)

    hi = Ri-(Ri**2+di**2-Rj**2)/(2*di)
    s = 1/3*np.pi*hi**2*(3*Ri-hi)
    return s, Ri


def calc_dVA(mol, conf, bond_list, target_atom_id, repeat_bond):

    sigma = 0
    neighbor_atom_bonds_set = [
        bond for bond in bond_list if target_atom_id in bond]

    # atom order is important!

    neighbor_atom_bonds = []
    for bond in neighbor_atom_bonds_set:
        if target_atom_id != bond[0]:
            bond = bond[1], bond[0]
        neighbor_atom_bonds.append(bond)

    for bond in neighbor_atom_bonds:
        s, Ri = calc_s(mol, conf, bond, repeat_bond)
        #print(s,  bond)
        sigma += s

    dVa = 4/3*np.pi*Ri**3-sigma
    return dVa
