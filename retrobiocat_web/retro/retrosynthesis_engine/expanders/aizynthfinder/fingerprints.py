from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np


def get_mol_fingerprint(rd_mol, radius, nbits=None):
    """
    Returns the Morgan fingerprint of the molecule
    """
    if nbits:
        key = (radius, nbits)
    else:
        key = (radius,)

    bitvect = AllChem.GetMorganFingerprintAsBitVect(rd_mol, *key)
    array = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(bitvect, array)
    return array

def get_reaction_fingerprint(reaction, radius, nbits=None):
    """
    Returns a difference fingerprint
    Reaction = ([product1, product2], [substrate1, substrate2]
    """

    products = reaction[0]
    reactants = reaction[1]
    product_fp = sum(get_mol_fingerprint(mol, radius, nbits=nbits) for mol in products)
    reactants_fp = sum(get_mol_fingerprint(mol, radius, nbits=nbits) for mol in reactants)

    return product_fp - reactants_fp  # type: ignore