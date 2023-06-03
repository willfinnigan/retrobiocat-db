from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
from retrobiocat_web.app.retrobiocat.functions import get_images
from rdkit import Chem
from rdkit.Chem.BRICS import BRICSDecompose
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile
import re
from retrobiocat_web.retro.network_pathway.network import Network


@bp.route('/_fragment_molecule', methods=['GET', 'POST'])
def fragment_molecule():
    smiles = request.form['smiles']
    smiles = rdkit_smile(smiles)

    if smiles is None or smiles == '':
        result = {'mol_dict': {}}
        return jsonify(result=result)

    mol = Chem.MolFromSmiles(smiles)
    list_smi = list(BRICSDecompose(mol, minFragmentSize=5, keepNonLeafNodes=True))

    list_processed = []
    for smi in list_smi:
        new_smi = re.sub(r"\[(?:[1-9]|[1-9][0-9])\*\]", '*', smi)
        list_processed.append(new_smi)

    mol_dict = {}
    for smi in list_processed:
        mol = Chem.MolFromSmiles(smi)
        img = get_images.moltosvg(mol,molSize=(200,200),kekulize=True)
        mol_dict[smi] = img

    result = {'mol_dict': mol_dict}
    return jsonify(result=result)

@bp.route('/_apply_chemical_steps_molecule', methods=['GET', 'POST'])
def apply_chemical_steps_molecule():
    smiles = request.form['smiles']
    smiles = rdkit_smile(smiles)

    if smiles is None or smiles == '':
        result = {'mol_dict': {}}
        return jsonify(result=result)

    network = Network(target_smiles=smiles)
    network.generate(smiles, 0, calculate_scores=False)
    new_substrate_nodes, new_reaction_nodes = network.add_chemical_step(smiles)

    list_processed = []
    for smi in new_substrate_nodes:
        new_smi = re.sub(r"\[(?:[1-9]|[1-9][0-9])\*\]", '*', smi)
        list_processed.append(new_smi)

    mol_dict = {}
    for smi in list_processed:
        mol = Chem.MolFromSmiles(smi)
        img = get_images.moltosvg(mol,molSize=(200,200),kekulize=True)
        mol_dict[smi] = img

    result = {'mol_dict': mol_dict}
    return jsonify(result=result)