from rdkit import Chem
from retrobiocat_web.analysis.drawing import images


def smiles_to_svg(activity_data):
    for i, record in enumerate(activity_data):
        if activity_data[i]["substrate_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_1_smiles"] = url

        if activity_data[i]["substrate_2_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["substrate_2_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["substrate_2_smiles"] = url

        if activity_data[i]["product_1_smiles"] != '':
            mol = Chem.MolFromSmiles(activity_data[i]["product_1_smiles"])
            url = images.moltosvg_url(mol)
            activity_data[i]["product_1_smiles"] = url

    return activity_data