from typing import List

from natsort import natsorted
from retrobiocat_web.mongo.models.biocatdb_models import ActivityMol
import mongoengine as db

def get_molecules_in_paper(paper) -> List[ActivityMol]:
    mols = ActivityMol.objects(db.Q(paper=paper)).order_by('name')
    return mols

def get_num_activity_mols_in_paper(paper):
    return ActivityMol.objects(db.Q(paper=paper)).count()

def activity_mol_from_id(mol_id):
    return ActivityMol.objects(id=mol_id).first()

def does_mol_name_already_exist_in_paper(mol_name, paper):
    mols = get_molecules_in_paper(paper)
    for mol in mols:
        if mol.name == mol_name:
            return True
    return False


