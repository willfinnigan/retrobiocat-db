from retrobiocat_web.mongo.models.reaction_models import Reaction
import mongoengine as db
import time
from rdkit.Chem import rdChemReactions
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rdchiral.initialization import rdchiralReaction

def reaction_names():
    return list(Reaction.objects().distinct('name'))

def reactions_of_type(enzyme_type_abbrevs, names_only=False, include_chemical=True):
    """ Return a list of reaction names which enzyme type can carry out """

    # if All, get all reactions. Otherwise specified enzyme_type abbreviation
    et_q = db.Q(enzyme_types__in=enzyme_type_abbrevs)
    if 'All' in enzyme_type_abbrevs or len(enzyme_type_abbrevs) is None:
        et_q = db.Q()

    ic_q = db.Q()
    if not include_chemical:
        ic_q = db.Q(enzyme_types__nin=['Chemical'])

    # if names_only, just return a list of names.
    if names_only == True:
        reactions_list = Reaction.objects(et_q & ic_q).distinct('name')
        reactions_list.sort()
        return reactions_list

    # otherwise return actual reaction objects
    else:
        reactions = Reaction.objects(et_q & ic_q)
        return reactions




def reaction_from_name(name, get_related=False):
    """Get reaction object from its name """

    if get_related == True:
        q = Reaction.objects(name=name).select_related()
    else:
        q = Reaction.objects(name=name)

    # return None if no match
    if len(q) == 0:
        return None

    return q[0]

def reverse_smarts(smarts):
    m = smarts.find('>>')
    start = smarts[m+2:]
    end = smarts[:m]
    new_smarts = start + '>>' + end
    return new_smarts

def get_reactions(include_experimental=False, include_two_step=False, include_requires_absence_of_water=False, only_reviewed=True):
    if include_experimental:
        q_exp = db.Q()
    else:
        q_exp = db.Q(experimental__ne=True)

    if include_two_step:
        q_two = db.Q()
    else:
        q_two = db.Q(two_step__ne=True)

    if include_requires_absence_of_water:
        q_water = db.Q()
    else:
        q_water = db.Q(requires_absence_of_water__ne=True)

    return Reaction.objects(q_exp & q_two & q_water).select_related()

def load_rxns(query_result, mode='rdchiral', reverse=False):
    rxns = {}
    for rxn in query_result:
        rxns[rxn.name] = []
        for rxn_string in rxn.smarts:
            try:
                if reverse == True:
                    rxn_string = reverse_smarts(rxn_string)
                if mode == 'rdchiral':
                    rxns[rxn.name].append(rdchiralReaction(rxn_string))
                else:
                    rxns[rxn.name].append(rdChemReactions.ReactionFromSmarts(rxn_string))
            except:
                print(f'Error processing rxn for {rxn.name} - mode={mode}, reverse={reverse}')
    return rxns

def get_multi_smarts_from_named_steps(rxn_steps, mode='rdchiral', reverse=False):
    """Rxn steps is a list of lists of reaction objects"""

    processed_steps = []
    processed_steps_strs = []
    for group_steps in rxn_steps:
        group_smas = []
        group_sma_str = []
        for step_rxn in group_steps:
            step = []
            step_str = []
            for rxn_string in step_rxn.smarts:
                if reverse == True:
                    rxn_string = reverse_smarts(rxn_string)
                if mode == 'rdchiral':
                    step.append(rdchiralReaction(rxn_string))
                else:
                    step.append(rdChemReactions.ReactionFromSmarts(rxn_string))
                step_str.append(rxn_string)
            group_smas.append(step)
            group_sma_str.append(step_str)

            if reverse == True:
                group_smas.reverse()
                group_sma_str.reverse()

        processed_steps.append(group_smas)
        processed_steps_strs.append(group_sma_str)

    return processed_steps

def load_multi_step(query_result, mode='rdchiral', reverse=False):
    rxns = {}
    for rxn in query_result:
        if rxn.two_step == True:
            if len(rxn.steps) != 0:
                rxns[rxn.name] = get_multi_smarts_from_named_steps(rxn.steps, mode, reverse=reverse)  # this will become a list of lists of lists: [[[r1a, r1b], [r2a, r2b]], [[r1a, r1b], [r2a, r2b]]]
    return rxns

def load_cofactors(query_result):
    reactionEnzymeCofactorDict = {}

    for rxn in query_result:
        reactionEnzymeCofactorDict[rxn.name] = {}

        for enz in rxn.cofactors:
            cofactor_minus = rxn.cofactors[enz]['cofactors_minus']
            cofactor_plus = rxn.cofactors[enz]['cofactors_plus']
            reactionEnzymeCofactorDict[rxn.name][enz] = [cofactor_plus, cofactor_minus]

    return reactionEnzymeCofactorDict

def load_reactions_and_enzymes(query_result):
    reactions = set()
    enzymes = set()
    reaction_enzyme_map = dict()

    for rxn in query_result:
        reactions.add(rxn.name)
        reaction_enzyme_map[rxn.name] = set()
        for enz in rxn.cofactors:
            enzymes.add(enz)
            reaction_enzyme_map[rxn.name].add(enz)
        reaction_enzyme_map[rxn.name] = list(reaction_enzyme_map[rxn.name])

    return list(reactions), list(enzymes), reaction_enzyme_map

def load_rxn_strings(query_result):
    rxns = {}
    for rxn in query_result:
        rxns[rxn.name] = []
        for rxn_string in rxn.smarts:
            rxns[rxn.name].append(rxn_string)
    return rxns

def load_rules_by_type(query_result):
    rules_by_type = {}
    for rxn in query_result:
        rule_type = rxn.type
        if rule_type not in rules_by_type:
            rules_by_type[rule_type] = []
        rules_by_type[rule_type].append(rxn.name)
    return rules_by_type


def load_all_reactions_for_single_name(reaction_name, reverse=False):
    reaction_query = Reaction.objects(name=reaction_name)
    rdkit_rxns = load_rxns(reaction_query, mode='rdkit', reverse=reverse)
    rdkit_multi = load_multi_step(reaction_query, mode='rdkit', reverse=reverse)
    rdchiral_rxns = load_rxns(reaction_query, reverse=reverse)
    rdchiral_multi = load_multi_step(reaction_query, reverse=reverse)
    return rdchiral_rxns, rdkit_rxns, rdchiral_multi, rdkit_multi


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    t0 = time.time()
    query_result = get_reactions(include_two_step=True)
    rxns = load_rxns(query_result)
    multi_step = load_multi_step(query_result)
    reactionEnzymeCofactorDict = load_cofactors(query_result)
    reactions, enzymes, reaction_enzyme_map = load_reactions_and_enzymes(query_result)
    t1 = time.time()

    print(f"Time to load rules = {round(t1-t0, 3)}")
    #print(reactions)

    multi_step_rxns = load_multi_step(query_result)
    print(multi_step_rxns)