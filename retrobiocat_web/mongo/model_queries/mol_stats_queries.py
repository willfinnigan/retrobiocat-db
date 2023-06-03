from retrobiocat_web.mongo.models.biocatdb_models import MoleculeStats


def smis_with_descriptors():
    return list(MoleculeStats.objects().distinct('smiles'))


def get_descriptors_for_list_smis(list_smis):
    list_smis = list(set(list_smis))
    if '' in list_smis:
        list_smis.remove('')

    mol_stats_objects = MoleculeStats.objects(smiles__in=list_smis).only('smiles', 'descriptors')

    descriptors_dict = {}
    for mol_stats_obj in mol_stats_objects:
        smi = mol_stats_obj.smiles
        desciptors = mol_stats_obj.descriptors
        descriptors_dict[smi] = desciptors

    return descriptors_dict



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    get_descriptors_for_list_smis(['CCCC=O'])
