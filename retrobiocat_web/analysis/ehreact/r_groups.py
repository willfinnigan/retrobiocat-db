from rdkit import Chem
from rdkit.Chem import rdRGroupDecomposition, PandasTools


def r_group_decomposition(list_smis, scaffold_smi):
    """ Do r group decomposition and return resulting dataframe """

    scaffold = Chem.MolFromSmiles(scaffold_smi)

    if scaffold is None:
        print(f'ERROR scaffold {scaffold_smi} is None')
        return None

    # Get mols from product smiles
    mols = []
    for smi in list_smis:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append(mol)

    # Use rdkit RGroupDecompose
    groups, unmatched = rdRGroupDecomposition.RGroupDecompose([scaffold], mols, asSmiles=False, asRows=False)

    # Get a list of the unmatched molecules (they don't match the scaffold)
    unmatched_mols = []
    for i in unmatched:
        unmatched_mols.append(mols[i])

    # Get a list of only the molecules that match the scaffold
    mols_match_scaffold = []
    for i, mol in enumerate(mols):
        if i not in unmatched:
            mols_match_scaffold.append(mol)

    # Use RGroupDecompositionToFrame to get dataframe of R group decomposition
    r_group_df = PandasTools.RGroupDecompositionToFrame(groups, mols_match_scaffold, include_core=True,
                                                        redraw_sidechains=False)

    return r_group_df

def get_core(list_smis, scaffold_smi):
    r_group_df = r_group_decomposition(list_smis, scaffold_smi)
    if len(list_smis) > 0 and r_group_df is not None:
        if 'Core' in r_group_df:
            return r_group_df['Core'][0]

    return Chem.MolFromSmiles(scaffold_smi)

