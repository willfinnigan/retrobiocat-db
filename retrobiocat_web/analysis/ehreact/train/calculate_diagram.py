import numpy as np
from .hasse import extended_hasse
from retrobiocat_web.analysis.ehreact.helpers.rdkit import read_in_smiles_unique, preprocess_seeds



def calculate_diagram_single_reactant(smiles, seed_list, verbose, quiet):
    """
    Computes a Hasse diagram of a list of molecule smiles.

    Parameters
    ----------
    smiles: List[str]
        List of SMILES.
    seed_list: List[str]
        List of SMILES seeds.
    verbose: bool
        Whether to print additional information.
    quiet: bool
        Whether to silence all output.

    Returns
    -------
    d: ehreact.diagram.diagram.Diagram
        The Hasse diagram of the input list of molecules.
    smiles_dict: dict
        A dictionary of the canonicalized input smiles.
    """

    # Create canonical smiles and mol objects:
    smiles, smiles_dict, skipped = read_in_smiles_unique(smiles)

    # Create seeds:
    seeds, rule_dict, num_smiles_seed = preprocess_seeds(seed_list, smiles, smiles_dict)

    if np.sum(num_smiles_seed) < len(smiles):
        raise ValueError("Not all molecules fit a seed. Check given seeds.")
    elif np.sum(num_smiles_seed) > len(smiles):
        raise ValueError("Seeds are not mutually exclusive. Check given seeds.")

    tags_core = {}
    d = extended_hasse(smiles_dict, seeds, rule_dict, tags_core, verbose, quiet)

    return d, smiles_dict

