import uuid

from rdkit import Chem

from retrobiocat_web.app.retrobiocat.functions import pubchem_funcs
from retrobiocat_web.mongo.model_queries.molecule_queries import get_num_activity_mols_in_paper, \
    does_mol_name_already_exist_in_paper
from retrobiocat_web.mongo.models.biocatdb_models import ActivityMol

def generate_numerical_activity_mol_name(paper):
    num = get_num_activity_mols_in_paper(paper)+1
    new_name = f"new_mol_{num}"
    while does_mol_name_already_exist_in_paper(new_name, paper):
        num += 1
        new_name = f"new_mol_{num}"
    return new_name

def get_chem_name_for_activity_mol(smi):
    compound = pubchem_funcs.get_pubchem_compound_from_smiles(smi)
    if compound != None:
        return compound.iupac_name
    return None

def add_identifier_if_name_already_exists(mol_name, paper):
    if not does_mol_name_already_exist_in_paper(mol_name, paper):
        return mol_name
    num = 1
    new_name = f"{mol_name}_{num}"
    while does_mol_name_already_exist_in_paper(new_name, paper):
        num += 1
        new_name = f"{mol_name}_{num}"
    return new_name


def make_new_activity_molecule(smi, paper, name=None, chem_name=None):
    if name is None:
        name = generate_numerical_activity_mol_name(paper)

    if chem_name is None and smi != "":
        chem_name = get_chem_name_for_activity_mol(smi)

    try:
        mol = Chem.MolFromSmiles(smi)
        smi = Chem.MolToSmiles(mol)
    except Exception as e:
        print('Error parsing smiles with rdkit when generating a new activity mol')
        print(str(e))
        return None

    new_mol = ActivityMol(name=name,
                          chem_name=chem_name,
                          smi=smi,
                          paper=paper)
    new_mol.save()
    return new_mol

def delete_activity_molecule(mol_obj):
    try:
        mol_obj.delete()
        return []
    except Exception as e:
        return [f"error deleting molecule - {e}"]

def update_activity_molecule(mol_obj, mol_name, smi):
    issues = []
    if mol_obj.name != mol_name:
        mol_name = add_identifier_if_name_already_exists(mol_name, mol_obj.paper)
        mol_obj.name = mol_name

    if mol_obj.smi != smi:
        try:
            rd_mol = Chem.MolFromSmiles(smi)
            mol_obj.smi = Chem.MolToSmiles(rd_mol)
        except Exception as e:
            print('Error parsing smiles with rdkit when generating a new activity mol')
            print(str(e))
            return [f'Error parsing smiles with rdkit when generating a new activity mol - {str(e)}']

    mol_obj.chem_name = get_chem_name_for_activity_mol(smi)
    mol_obj.save()
    return issues


def update_activity_molecule_name(mol_obj, mol_name):
    issues = []
    if mol_obj.name != mol_name:
        if does_mol_name_already_exist_in_paper(mol_name, mol_obj.paper):
            issues.append(f'Error - cannot update mol name from {mol_obj.name} to {mol_name} because this name already in paper')
        else:
            mol_obj.name = mol_name
            mol_obj.save()
    return issues
