from retrobiocat_web.similarity.retrobiocat.mongo_data_query import MongoSpecificityData
from retrobiocat_web.similarity.similarity_tools import get_single_fp, get_fingerprints, bulk_similarity
from retrobiocat_web.mongo.default_connection import make_default_connection
make_default_connection()


def test_fps_are_correct():
    smis = ['C=CCC(C(=O)NCc1ccccc1)C(C)O', 'CC(C)C[C@@H](C)N', 'CC(O)CCCCO', 'C=CC(O)CCCCC', 'OC1CCOC1',
            'C[C@H](O)c1cccc(Br)c1', 'N[C@@H](C(=O)N[C@@H]1C(=O)N2C(C(=O)O)=C(Cl)CSC12)c1ccccc1', '',
            'CC1CNC(c2ccccc2)CN1', 'O=C(NCCc1ccccc1)c1ccc2[nH]c3ccccc3c2c1']

    rdkit_fps = [get_single_fp(smi) for smi in smis]
    mongo_data = MongoSpecificityData()
    fps, converted_smis = mongo_data.get_fps(smis)
    assert rdkit_fps[0] == fps[0]

def test_bulk_similarity_empty_string():
    target_smi = 'CCCCCCCCC=O'
    compare_with = ['', '', '']
    target_fp = get_single_fp(target_smi)
    mongo_data = MongoSpecificityData()
    sims = bulk_similarity(target_fp, compare_with, mongo_data.get_fps)
    assert sims == {"": 0}