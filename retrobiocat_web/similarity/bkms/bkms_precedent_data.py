from pathlib import Path
from typing import Sequence

import pandas as pd

data_folder = str(Path(__file__).parents[2]) + '/data/bkms'

class BKMS_Data():
    df = None
    fp_df = None

    @classmethod
    def load_df(cls):
        if cls.df is None:
            cls.df = pd.read_hdf(f"{data_folder}/bkms_metadata.hdf")

    def get_data_by_id(self, query_ids: Sequence[int]=()) -> pd.DataFrame:
        self.load_df()
        p_df = self.df.loc[query_ids]
        return p_df

if __name__ == '__main__':
    bkms_data = BKMS_Data()
    bkms_data.load_df()
    result = bkms_data.get_data_by_id([1,2,3])

