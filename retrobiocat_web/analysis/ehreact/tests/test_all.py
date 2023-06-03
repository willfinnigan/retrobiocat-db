import pytest
from retrobiocat_web.analysis.ehreact.create_hasse_network import hasse_diagrame, hasse_diagram_to_network, \
    get_seed_list, remove_non_matches, ensure_correct_seed_coverage
from retrobiocat_web.analysis.ehreact.get_substrate_summary import summarise
from retrobiocat_web.analysis.data_query.data_query_from_args import data_query_from_args

from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings
from retrobiocat_web.mongo.model_queries.reaction_queries import reactions_of_type

make_default_connection()


args_list = []
for et_abbreviation in all_enzyme_type_strings():
    reactions = reactions_of_type([et_abbreviation], names_only=False, include_chemical=False)
    for reaction in reactions:
        args_list.append((reaction.name, et_abbreviation))
@pytest.mark.parametrize('reaction_name, enzyme_type', args_list)
def test_all_enzymes_and_reactions(reaction_name, enzyme_type):
    dq = data_query_from_args({'reaction': reaction_name,
                               'enzyme_type': enzyme_type})
    smis = dq.unique_smiles()
    seed_list = get_seed_list(dq.reaction_names(), dq.smi_col)
    seed_list = ensure_correct_seed_coverage(seed_list, smis)
    smis, non_matches, ori_smiles = remove_non_matches(seed_list, smis)
    d = hasse_diagrame(smis, seed_list=seed_list)
    graph, start_node = hasse_diagram_to_network(d)
    leafs, not_leaf_children, not_leaf_counts, not_leaf_cores, expandables, parent_core = summarise(start_node,
                                                                                                graph, ori_smiles)
