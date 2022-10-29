from .atom_dict import atom_van_dict
from .mol_utils import get_bond_length
import numpy as np


def calc_s(mol, conf, bond, repeat_bond=[]):

    atom_id1 = bond[0]
    atom_id2 = bond[1]
    atom1 = mol.GetAtomWithIdx(atom_id1).GetSymbol()
    atom2 = mol.GetAtomWithIdx(atom_id2).GetSymbol()

    Ri = atom_van_dict[atom1]
    Rj = atom_van_dict[atom2]

    if set(repeat_bond) == set((atom_id1, atom_id2)):
        di = 1.7
    else:
        di = get_bond_length(conf, atom_id1, atom_id2)

    hi = Ri-(Ri**2+di**2-Rj**2)/(2*di)
    s = 1/3*np.pi*hi**2*(3*Ri-hi)
    return s, Ri


def calc_dVA(mol, conf, bond_list, target_atom_id, repeat_bond):

    sigma = 0
    neighbor_atom_bonds = [
        bond for bond in bond_list if target_atom_id in bond]
    for bond in neighbor_atom_bonds:
        s, Ri = calc_s(mol, conf, bond, repeat_bond)
        print(s,  bond)
        sigma += s

    dVa = 4/3*np.pi*Ri**3-sigma
    return dVa
