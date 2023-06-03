import pandas as pd
from retrobiocat_web.mongo.models.retro_tests import TestCascade, TestCascadeResult, TestCascadeRun
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile


def load_tests(excel_path):
    print("Load cascade tests dataframe..")
    if '.xlsx' in excel_path:
        df = pd.read_excel(excel_path)
    else:
        print('Could not load file')
        df = None
    return df

def create_retrobiocat_tests(tests_df, name_prefix='retrobiocat_'):
    for index, row in tests_df.iterrows():
        starting_smiles = [rdkit_smile(smi) for smi in list(row['Substrates'].split(', '))]
        enzymes = list(row['Enzymes'].split(', '))

        new_test = TestCascade(name=f"{name_prefix}{row['Num']}",
                               doi=row['DOI'],
                               target_smiles=rdkit_smile(row['Target Product']),
                               starting_smiles=starting_smiles,
                               enzymes=enzymes,
                               test_type=row['Test type'])
        new_test.save()

    print('retrobiocat cascade tests loaded')

if __name__ == '__main__':
    from pathlib import Path
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    TestCascade.drop_collection()
    TestCascadeResult.drop_collection()
    TestCascadeRun.drop_collection()

    excel_path = str(Path(__file__).parents[3]) + '/scripts/pathway_testing/test_pathways.xlsx'
    df = load_tests(excel_path)
    create_retrobiocat_tests(df)

    print(TestCascade.objects().count())
