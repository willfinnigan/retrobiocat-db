from rdkit import Chem
from retrobiocat_web.mongo.models.reaction_models import Reaction

def smarts_to_smiles(smarts):
    mol = Chem.MolFromSmarts(smarts)
    return Chem.MolToSmiles(mol)

def extract_seeds(rxn_smarts):
    parts = rxn_smarts.split('>>')
    product = smarts_to_smiles(parts[0])
    substrates = parts[1].split('.')

    substrate_1, substrate_2 = "", ""
    if len(substrates) >= 1:
        substrate_1 = smarts_to_smiles(substrates[0])
    if len(substrates) >= 2:
        substrate_2 = smarts_to_smiles(substrates[1])

    return product, substrate_1, substrate_2


def seeds_from_reaction(rxn_obj: Reaction):
    product_seeds, substrate_1_seeds, substrate_2_seeds = [], [], []

    for smarts in rxn_obj.smarts:
        product, substrate_1, substrate_2 = extract_seeds(smarts)
        product_seeds.append(substrate_1)
        substrate_1_seeds.append(substrate_1)
        substrate_2_seeds.append(substrate_2)

    product_seeds = [s for s in product_seeds if s != ""]
    substrate_1_seeds = [s for s in substrate_1_seeds if s != ""]
    substrate_2_seeds = [s for s in substrate_2_seeds if s != ""]

    return product_seeds, substrate_1_seeds, substrate_2_seeds


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.mongo.model_queries.reaction_queries import reaction_from_name

    name = 'Imine reduction'
    rxn_obj = reaction_from_name(name)

    print(seeds_from_reaction(rxn_obj))

