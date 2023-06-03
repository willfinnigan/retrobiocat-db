import pandas as pd
from rdkit import DataStructs

from retrobiocat_web.mongo.models.biocatdb_models import Molecule


def load_fp_df_from_mongo():
    mode = 'rdfp_default'

    def convert_to_fp(bitstring):
        return DataStructs.CreateFromBitString(bitstring)

    query_result = Molecule.objects.as_pymongo().only('smiles', mode)
    fp_df = pd.DataFrame(list(query_result))
    fp_df.drop(columns=['_id'], inplace=True, errors='ignore')

    fp_df[mode] = fp_df[mode].apply(convert_to_fp)
    fp_df.rename(columns={mode:'fp'}, inplace=True)

    return fp_df