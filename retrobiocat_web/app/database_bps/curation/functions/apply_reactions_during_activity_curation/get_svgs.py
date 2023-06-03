from rdkit import Chem
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
from retrobiocat_web.retro.visualisation.rdkit_images import moltosvg


def make_svg_dict(products):
    svg_dict = {}
    for list_smis in products:
        for smi in list_smis:
            if smi not in svg_dict:
                mol = Chem.MolFromSmiles(smi)
                svg_dict[smi] = moltosvg(mol)
    return svg_dict

def get_reaction_svg(s1, s2, p):
    if s2 == "":
        smi_rxn = f"{s1}>>{p}"
    else:
        smi_rxn = f"{s1}.{s2}>>{p}"

    svg = smiles_rxn_to_svg(smi_rxn, rxnSize=(500,200))
    return svg

def retro_outcome_svgs(product, outcomes):
    svg_dict = make_svg_dict(outcomes)
    reaction_svgs = []
    for smi_list in outcomes:
        if len(smi_list) == 1:
            reaction_svgs.append(get_reaction_svg(smi_list[0], '', product))
        elif len(smi_list) == 2:
            reaction_svgs.append(get_reaction_svg(smi_list[0], smi_list[1], product))
        else:
            reaction_svgs.append(get_reaction_svg('', '', product))

    return svg_dict, reaction_svgs

def fwd_outcome_svgs(s1, s2, outcomes):
    svg_dict = make_svg_dict(outcomes)
    reaction_svgs = []
    for smi_list in outcomes:
        reaction_svgs.append(get_reaction_svg(s1, s2, smi_list[0]))
    return svg_dict, reaction_svgs

