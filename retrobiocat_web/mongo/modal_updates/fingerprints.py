from typing import List, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.DataStructs import ExplicitBitVect

from retrobiocat_web.mongo.model_queries import specificity_data_query
from retrobiocat_web.mongo.init_db.make_molecule_db import fp_molecules_to_db
from retrobiocat_web.mongo.models.biocatdb_models import Molecule

def get_fingerprints(smis_to_convert: List[str]) -> Tuple[List[ExplicitBitVect], List[str]]:
    mols = []
    converted_smis = []
    for smi in smis_to_convert:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                mols.append(mol)
                converted_smis.append(smi)
        except:
            pass

    fps = [Chem.RDKFingerprint(mol) for mol in mols]
    return fps, converted_smis

def task_update_fingerprints():
    """
    This task should run whenever activity data is added / updated / deleted
    It will update the Molecules document so every smiles in the database has a pre-calculated fingerprint
    """

    # Get a list of smiles in db, and those with fingerprints
    unique_smis = []
    unique_smis += specificity_data_query.get_unique_smiles('substrate_1_smiles')
    unique_smis += specificity_data_query.get_unique_smiles('substrate_2_smiles')
    unique_smis += specificity_data_query.get_unique_smiles('product_1_smiles')
    unique_smis = list(set(unique_smis))

    current_fp_smis = specificity_data_query.get_fp_smiles()

    to_remove = [smi for smi in current_fp_smis if smi not in unique_smis]
    to_add = [smi for smi in unique_smis if smi not in current_fp_smis]

    # Add new fingerprints to database
    fps, smis = get_fingerprints(to_add)
    fp_df = pd.DataFrame({'smiles': smis, 'rdfp_default': fps})
    fp_molecules_to_db(fp_df)

    # Delete old fingerprints
    mols = Molecule.objects(smiles__in=to_remove)
    mols.delete()


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    task_update_fingerprints()
    result = Molecule.objects().as_pymongo()
    print(result)
