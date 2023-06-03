from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.Lipinski import NumHAcceptors, NHOHCount, NOCount, RingCount, NumHDonors

from retrobiocat_web.logging import add_logger
from retrobiocat_web.mongo.model_queries.activity_queries import distinct_product_1_smiles, distinct_substrate_1_smiles, \
    distinct_substrate_2_smiles
from retrobiocat_web.mongo.model_queries.mol_stats_queries import smis_with_descriptors
from retrobiocat_web.mongo.models.biocatdb_models import MoleculeStats


class MolDescriptor_Updater():

    def __init__(self, log_level='WARNING'):
        self.logger = add_logger('MolDescriptor Updater', level=log_level)

    def create_mol_descriptor(self, smi):
        """
        Create a molstats object in the database
        1. get mol (and back to smi to ensure correct format)
        2. optionally get name from pubchem (disable is config for commercial users)
        3. calc descriptors using mordred
        """

        mol = Chem.MolFromSmiles(smi)
        smi = Chem.MolToSmiles(mol)

        descriptors = self._get_descriptors(mol)
        self.logger.debug(f"Descriptors for {smi} = {descriptors}")

        new_mol_stat_obj = MoleculeStats(smiles=smi,
                                         descriptors=descriptors)
        new_mol_stat_obj.save()


    def _get_descriptors(self, mol):
        stats = {'mw': Descriptors.MolWt(mol),
                 'logp': MolLogP(mol),
                 'num_nh_or_oh': NHOHCount(mol),
                 'num_n_or_o': NOCount(mol),
                 'num_h_acceptors': NumHAcceptors(mol),
                 'num_h_donors': NumHDonors(mol),
                 'num_rings': RingCount(mol)
                 }
        return stats



def task_update_molecular_descriptors(log_level='WARNING'):
    # Get a list of smiles in db, and those with descriptors
    all_smis = distinct_product_1_smiles() + distinct_substrate_1_smiles() + distinct_substrate_2_smiles()
    unique_smis = list(set(all_smis))

    descriptor_smis = smis_with_descriptors()

    # Get the smiles which are new
    new_smiles = []
    for smi in unique_smis:
        if smi not in descriptor_smis:
            new_smiles.append(smi)

    # Get smiles where entry no longer needed
    old_smiles = []
    for smi in descriptor_smis:
        if smi not in unique_smis:
            old_smiles.append(smi)

    # Add new smis to database
    updater = MolDescriptor_Updater(log_level=log_level)
    for smi in new_smiles:
        updater.create_mol_descriptor(smi)

    # Delete old entries
    mols = MoleculeStats.objects(smiles__in=old_smiles)
    mols.delete()



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    MoleculeStats.drop_collection()
    task_update_molecular_descriptors(log_level='DEBUG')