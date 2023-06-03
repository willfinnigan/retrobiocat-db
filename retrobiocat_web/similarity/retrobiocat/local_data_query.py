from pathlib import Path
from typing import List

import pandas as pd
from retrobiocat_web.similarity.similarity_tools import make_fp_df, get_fingerprints_from_fpdf

data_folder = str(Path(__file__).parents[2]) + '/data/retrobiocat'

class LocalSpecificityData():
    df = None
    fp_df = None

    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_excel(f"{data_folder}/trial_activity_data.xlsx")

    @classmethod
    def load_fp_df(cls):
        if cls.fp_df is None:
            cls.load_df()
            cls.fp_df = make_fp_df(cls.df, 'product_1_smiles')

    def get_fp_df(self):
        self.load_fp_df()
        return self.fp_df

    def get_single_fp(self, smi: str):
        self.load_df()
        self.load_fp_df()
        try:
            return self.fp_df.loc[smi]['fp']
        except:
            return None

    def get_fps(self, smis: List[str]):
        self.load_df()
        self.load_fp_df()

        fps, converted_smis = get_fingerprints_from_fpdf(smis, self.fp_df)
        return fps, converted_smis

    def query_data(self, reaction_name: str = 'All', enzyme_types: List[str] = ()) -> pd.DataFrame:
        self.load_df()
        df = self.df

        if reaction_name != 'All':
            df = df[df['reaction'] == reaction_name]

        if len(enzyme_types) != 0 and 'All' not in enzyme_types:
            df = df[df['enzyme_type'].isin(enzyme_types)]

        return df




if __name__ == '__main__':
    local_data = LocalSpecificityData()
    #df = local_data.query_data([], [])
    #print(df.head())
    fps, smis = local_data.get_fps(['CCC=O', 'CCCCC=O'])




    # df = local_data.query_data(['Carboxylic acid reduction'], [])
    # print(df.head())
    #
    # df = local_data.query_data(['Reductive amination'], [])
    # print(df.head())
    #
    # df = local_data.query_data(['All'], ['CAR'])
    # print(df.head())
    #
    # from retrobiocat_web.mongo.default_connection import make_default_connection
    # make_default_connection()
    # from retrobiocat_web.mongo.model_queries.specificity_data_query import query_specificity_data
    #
    # result = query_specificity_data(['Carboxylic acid reduction'], [], only_reviewed=True)
    # print(result)