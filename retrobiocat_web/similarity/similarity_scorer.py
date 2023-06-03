from dataclasses import dataclass
from typing import List, Optional, Callable, Tuple
import pandas as pd
from rdkit.DataStructs import ExplicitBitVect

from retrobiocat_web.similarity.similarity_tools import get_fingerprints, get_single_fp, bulk_similarity

Smi = str
DataQueryFunction = Callable[..., pd.DataFrame]
FpQueryFunction = Callable[[List[Smi]], Tuple[List[ExplicitBitVect], List[Smi]]]
RankingFunction = Callable[[pd.DataFrame, int], pd.DataFrame]

@dataclass
class DataColumns():
    product_column: str
    substrate_columns: List[str]
    enzyme_column: str
    id_column: Optional[str]

class SimilarityScorer():

    def __init__(self,
                 columns: DataColumns,
                 data_query_function: callable,
                 fp_query_function: FpQueryFunction = get_fingerprints,
                 ranking_function: Optional[RankingFunction] = None):

        self.columns = columns

        self.data_query_function = data_query_function  # a method for getting dataframe of data
        self.fp_query_function = fp_query_function  # a method for generating fingerprints from a list of smiles
        self.ranking_function = ranking_function  # optionally ranks results

    def score_data(self, target_smi: Smi, topn: int, cutoff: float, **key_atts):
        """
        Finds similar data to target smi and returns it as a dataframe
        :param target_smi: the smiles to compare with
        :param topn: the number of top enzymes to return per similar molecule
        :param cutoff: the similarity cutoff
        """
        data_df = self.data_query_function(**key_atts)  # get the relevant data
        target_fp = get_single_fp(target_smi)  # generate fp for target_smi
        smis = list(data_df[self.columns.product_column].unique())  # the smis to compare target_smi with
        smi_sims = bulk_similarity(target_fp, smis, self.fp_query_function)  # get similarities to smis
        smi_sims = {smi: sim for smi, sim in smi_sims.items() if sim >= cutoff}  # get dict of smi: sim above cutoff

        result_df = data_df[data_df[self.columns.product_column].isin(smi_sims.keys())].copy()  # get df of similiar smis
        result_df['similarity'] = result_df[self.columns.product_column].map(smi_sims)  # add similarity scores
        result_df.sort_values('similarity', ascending=False, inplace=True)  # put most similar first
        if self.ranking_function is not None:
            result_df = self.ranking_function(result_df, topn)  # rank and take only top_n rows per smiles

        return result_df
