import pytest
from retrobiocat_web.analysis.ehreact.create_hasse_network import hasse_diagrame, hasse_diagram_to_network, \
    remove_non_matches, ensure_correct_seed_coverage
from retrobiocat_web.analysis.ehreact.get_substrate_summary import summarise, leafs_under_node
from retrobiocat_web.analysis.ehreact.seed_filtering import filter_smiles_with_smarts
from retrobiocat_web.analysis.ehreact.tests import smiles_lists

from retrobiocat_web.mongo.default_connection import make_default_connection

make_default_connection()


@pytest.mark.parametrize('smiles, seeds', [(smiles_lists.car_smiles, smiles_lists.car_seeds),
                                           (smiles_lists.imi_red_smiles, smiles_lists.imi_red_seeds),
                                           (smiles_lists.imi_red_prod_smiles, smiles_lists.imi_red_prod_seeds)
                                           ])
def test_all_smis_captured_in_graph(smiles, seeds):
    smis = filter_smiles_with_smarts(smiles, seeds)
    seeds = ensure_correct_seed_coverage(seeds, smis)
    smis, non_matches, ori_smiles = remove_non_matches(seeds, smis)
    d = hasse_diagrame(smis, seed_list=seeds)
    graph, start_node = hasse_diagram_to_network(d)

    all_leafs = leafs_under_node(start_node, graph)
    #list_ori_smis = [ori_smiles[s] for s in smis]
    print(f"{len(smis)} smiles, {len(all_leafs)} leafs")
    for smi in smis:
        assert smi in graph.nodes

    #assert sorted(all_leafs) == sorted(list_ori_smis)


@pytest.mark.parametrize('smiles, seeds', [(smiles_lists.car_smiles, smiles_lists.car_seeds)])
def test_move_on_one_if_no_leafs_but_has_not_leafs(smiles, seeds):
    smis = filter_smiles_with_smarts(smiles, seeds)
    seeds = ensure_correct_seed_coverage(seeds, smis)
    smis, non_matches, ori_smiles = remove_non_matches(seeds, smis)
    d = hasse_diagrame(smis, seed_list=seeds)
    graph, start_node = hasse_diagram_to_network(d)
    leafs, not_leaf_children_dict, not_leaf_counts, not_leaf_cores, expandables, parent_core = summarise(start_node, graph, ori_smiles)
    assert len(leafs) != 0


@pytest.mark.parametrize('smiles, seeds', [(smiles_lists.imi_red_smiles, smiles_lists.imi_red_seeds),
                                           (smiles_lists.imi_red_prod_smiles, smiles_lists.imi_red_prod_seeds)])
def test_ensure_all_mols_have_a_seed(smiles, seeds):
    smis = filter_smiles_with_smarts(smiles, seeds)
    seeds = ensure_correct_seed_coverage(seeds, smis)
    smis, non_matches, ori_smiles = remove_non_matches(seeds, smis)
    d = hasse_diagrame(smis, seed_list=seeds)
    graph, start_node = hasse_diagram_to_network(d)
    leafs, not_leaf_children_dict, not_leaf_counts, not_leaf_cores, expandables, parent_core = summarise(start_node, graph, ori_smiles)


