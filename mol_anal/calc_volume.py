from .atom_dict import atom_van_dict
from .mol_utils import get_bond_length
import numpy as np


def calc_s(mol, conf, bond):

    atom_id1 = bond[0]
    atom_id2 = bond[1]
    atom1 = mol.GetAtomWithIdx(atom_id1).GetSymbol()
    atom2 = mol.GetAtomWithIdx(atom_id2).GetSymbol()

    Ri = atom_van_dict[atom1]
    Rj = atom_van_dict[atom2]
    di = get_bond_length(conf, atom_id1, atom_id2)

    hi = Ri-(Ri**2+di**2-Rj**2)/(2*di)
    s = 1/3*np.pi*hi**2*(3*Ri-hi)
    return s
