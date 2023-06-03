from pathlib import Path
import pandas as pd
from retrobiocat_web.retro.molecular_similarity_search import make_fingerprints

data_folder = str(Path(__file__).parents[2]) + '/data/retrobiocat'

class LocalSpecificityData():
    df = None
    fp_df = None

    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_excel(f"{data_folder}/trial_activity_data.xlsx")

    @classmethod
    def load_fp_df(cls, fp_settings):
        if cls.fp_df is None:
            cls.load_df()
            cls.fp_df = make_fingerprints.make_fingerprint_df(cls.df, fp_settings[0], fp_settings[1])

    def query_data(self, listReactions, listEnzymes, only_reviewed=False):
        self.load_df()
        df = self.df

        if len(listReactions) != 0 and 'All' not in listReactions:
            df = df[df['reaction'].isin(listReactions)]

        if len(listEnzymes) != 0 and 'All' not in listEnzymes:
            df = df[df['enzyme_type'].isin(listEnzymes)]

        return df




if __name__ == '__main__':
    local_data = LocalSpecificityData()
    #df = local_data.query_data([], [])
    #print(df.head())

    df = local_data.query_data(['Carboxylic acid reduction'], [])
    print(df.head())

    df = local_data.query_data(['Reductive amination'], [])
    print(df.head())

    df = local_data.query_data(['All'], ['CAR'])
    print(df.head())