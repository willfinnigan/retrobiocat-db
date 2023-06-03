from pathlib import Path
from typing import List

import pandas as pd

from retrobiocat_web.mongo.model_queries.fingerprint_query import load_fp_df_from_mongo
from retrobiocat_web.mongo.model_queries.specificity_data_query import query_specificity_data
from retrobiocat_web.similarity.similarity_tools import get_fingerprints_from_fpdf

data_folder = str(Path(__file__).parents[3]) + '/data/retrobiocat'

class MongoSpecificityData():

    def __init__(self):
        self.fp_df = None

    def load_fp_df(self):
        if self.fp_df is None:
            self.fp_df = load_fp_df_from_mongo()
            self.fp_df.set_index('smiles', inplace=True)

    def get_fp_df(self):
        self.load_fp_df()
        return self.fp_df

    def get_single_fp(self, smi: str):
        self.load_fp_df()
        try:
            return self.fp_df.loc[smi]['fp']
        except:
            return None

    def get_fps(self, smis: List[str]):
        self.load_fp_df()
        fps, converted_smis = get_fingerprints_from_fpdf(smis, self.fp_df)
        return fps, converted_smis

    def query_data(self, reaction_name: str = 'All', enzyme_types: List[str] = (), only_reviewed=False) -> pd.DataFrame:
        return query_specificity_data([reaction_name], enzyme_types, only_reviewed=only_reviewed)



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    mongo_data = MongoSpecificityData()
    mongo_data.load_fp_df()
    print(mongo_data.fp_df.head())
    print(mongo_data.get_fps(['CCCO', 'CCC=O']))