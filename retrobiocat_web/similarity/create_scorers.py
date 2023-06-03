from retrobiocat_web.similarity.bkms.bkms_precedent_data import BKMS_Data
from retrobiocat_web.similarity.retrobiocat.local_data_query import LocalSpecificityData
from retrobiocat_web.similarity.retrobiocat.mongo_data_query import MongoSpecificityData
from retrobiocat_web.similarity.retrobiocat.rank_precedents import get_best_enzymes
from retrobiocat_web.similarity.similarity_scorer import DataColumns, SimilarityScorer


def create_retrobiocat_scorer(mode: str) -> SimilarityScorer:
    retrobiocat_data_columns = DataColumns(product_column='product_1_smiles',
                                           substrate_columns=['substrate_1_smiles', 'substrate_2_smiles'],
                                           enzyme_column='enzyme_type',
                                           id_column='_id')

    if mode == 'local':
        local_data = LocalSpecificityData()
        retrobiocat_scorer = SimilarityScorer(columns=retrobiocat_data_columns,
                                              data_query_function=local_data.query_data,
                                              fp_query_function=local_data.get_fps,
                                              ranking_function=get_best_enzymes)
        return retrobiocat_scorer

    if mode == 'mongo':
        mongo_data = MongoSpecificityData()
        retrobiocat_scorer = SimilarityScorer(columns=retrobiocat_data_columns,
                                              data_query_function=mongo_data.query_data,
                                              fp_query_function=mongo_data.get_fps,
                                              ranking_function=get_best_enzymes)
        return retrobiocat_scorer


def create_bkms_scorer(mode: str) -> SimilarityScorer:
    bkms_data_columns = DataColumns(product_column='products',
                                    substrate_columns=[],
                                    enzyme_column='Recommended_Name',
                                    id_column='ID')

    if mode == 'local':
        local_data = BKMS_Data()
        bkms_scorer = SimilarityScorer(columns=bkms_data_columns,
                                       data_query_function=local_data.get_data_by_id)
        return bkms_scorer


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    rbc_scorer = create_retrobiocat_scorer('mongo')
    data = rbc_scorer.score_data('CCCCCCCCC=O', 1, 0.99, reaction_name='All', enzyme_types=['All'])
    #print(data)

    no_product_data = data[data['product_1_smiles'] == ''].iloc[0]

