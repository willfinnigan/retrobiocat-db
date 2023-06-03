from retrobiocat_web.analysis.mol_images import smitosvg_url


def mol_svg_for_mol_table(smi):
    return smitosvg_url(smi, molSize=(100,25))

def get_new_mol_json(mol):
    mol_dict = {'_id': str(mol.id),
                'smi': mol.smi,
                'name': mol.name,
                'chem_name': mol.chem_name,
                'mol': mol_svg_for_mol_table(mol.smi)
                }
    return mol_dict
