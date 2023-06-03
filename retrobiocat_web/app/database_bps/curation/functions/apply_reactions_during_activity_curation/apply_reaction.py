from retrobiocat_web.mongo.model_queries.reaction_queries import load_all_reactions_for_single_name
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rule_applicator import RuleApplicator



def retro_rxn(rxn_name, product_smi, combine_enantiomers, use_rdkit):

    rule_applicator = RuleApplicator()
    rule_applicator.config.combine_enantiomers = combine_enantiomers

    rdchiral_rxns, rdkit_rxns, rdchiral_multi, rdkit_multi = load_all_reactions_for_single_name(rxn_name)

    outcomes = []
    if use_rdkit == False:
        outcomes = rule_applicator.apply_rdchiral(product_smi, rdchiral_rxns, multistep_rxns=rdchiral_multi)
    if len(outcomes) == 0:
        outcomes = rule_applicator.apply_rdkit(product_smi, rdkit_rxns, multistep_rxns=rdkit_multi)

    rxn_outcomes = outcomes.get(rxn_name, [])

    no_dupl_outcomes = []
    for prod in rxn_outcomes:
        if prod not in no_dupl_outcomes:
            no_dupl_outcomes.append(prod)
    return no_dupl_outcomes

def reaction_only_takes_one_reactants(rdkit_rxns):
    if len(rdkit_rxns) == 0:
        return False

    first_key = list(rdkit_rxns.keys())[0]
    rxns = rdkit_rxns[first_key]

    if len(rxns) == 0:
        return False

    rxn = rxns[0]

    try:
        num = rxn.GetNumReactantTemplates()
    except:
        num = 0

    if num == 1:
        return True
    return False


def fwd_rxn(substrate_1, substrate_2, rxn_name, combine_enantiomers, use_rdkit):
    rule_applicator = RuleApplicator()
    rule_applicator.config.combine_enantiomers = combine_enantiomers

    rdchiral_rxns, rdkit_rxns, rdchiral_multi, rdkit_multi = load_all_reactions_for_single_name(rxn_name, reverse=True)

    smi = f"{substrate_1}"
    if not reaction_only_takes_one_reactants(rdkit_rxns):
        if substrate_2 != "":
            smi += f".{substrate_2}"

    outcomes = []
    if use_rdkit == False:
        outcomes = rule_applicator.apply_rdchiral(smi, rdchiral_rxns, multistep_rxns=rdchiral_multi)
    if len(outcomes) == 0:
        outcomes = rule_applicator.apply_rdkit(smi, rdkit_rxns, multistep_rxns=rdkit_multi)
    rxn_outcomes = outcomes.get(rxn_name, [])

    no_dupl_outcomes = []
    for prod in rxn_outcomes:
        if prod not in no_dupl_outcomes:
            no_dupl_outcomes.append(prod)
    return no_dupl_outcomes


if __name__ == '__main__':
    from rdkit import Chem
    from rdkit.Chem import rdChemReactions
    rxn = rdChemReactions.ReactionFromSmarts('[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]')

    print(rxn.GetNumReactantTemplates())
