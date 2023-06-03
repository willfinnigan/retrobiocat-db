from retrobiocat_web.mongo.models.biocatdb_models import Molecule
from rdkit import Chem
from retrobiocat_web.mongo.model_queries import specificity_data_query


def fp_molecules_to_db(fp_df):
    print("Save fingerprints into molecule objects..")
    for i, row in fp_df.iterrows():
        if len(Molecule.objects(smiles=row['smiles'])) == 0:
            mol = Molecule()
        else:
            mol = Molecule.objects(smiles=row['smiles'])[0]

        mol.smiles = row['smiles']
        mol.mol = Chem.MolFromSmiles(row['smiles']).ToBinary()
        mol['rdfp_default'] = row['rdfp_default'].ToBitString()

        mol.save()
    print("..done")

def get_spec_df_from_mongo():
    print('Load all of specdf from mongo..')
    spec_df = specificity_data_query.query_specificity_data(['All'], ['All'])
    return spec_df




