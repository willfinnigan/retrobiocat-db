import pytest

from retrobiocat_web.retro.network_pathway import node_analysis
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.mcts.mcts import MCTS
from retrobiocat_web.retro.pathway_search.mcts.tree_node import MCTS_Node

from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine

make_default_connection()


@pytest.mark.parametrize('max_length,num_pathways', [[1,3], [2,6]])
def test_number_of_steps_gives_expected_num_pathways(max_length, num_pathways):
    target_smiles = 'c1ccc([C@@H]2CCCCN2)cc1'
    network = Network(target_smiles=target_smiles)
    mcts = MCTS(network, log_level='DEBUG')
    mcts.config.max_search_time = 60
    mcts.config.maxLength = max_length
    mcts.config.only_solved = False
    mcts.config.biocatalysis_expansion = True
    mcts.config.chemistry_expansion = False
    mcts.config.biosynthesis_expansion_expansion = False
    mcts.run()
    pathways = mcts.get_pathways()

    assert len(pathways) == num_pathways

@pytest.mark.parametrize('max_length,expected_end_node', [[3, "O=CCCCC(=O)c1ccccc1"], [4, "O=C(O)CCCC(=O)c1ccccc1"]])
def test_expected_endnode_is_present(max_length, expected_end_node):
    target_smiles = 'c1ccc([C@@H]2CCCCN2)cc1'
    network = Network(target_smiles=target_smiles)
    mcts = MCTS(network, log_level='WARNING')
    mcts.config.max_search_time = 60
    mcts.config.maxLength = max_length
    mcts.config.only_solved = False
    mcts.config.biocatalysis_expansion = True
    mcts.config.chemistry_expansion = False
    mcts.config.biosynthesis_expansion_expansion = False
    mcts.run()
    pathways = mcts.get_pathways()

    mol_is_present = False
    for pathway in pathways:
        if pathway.end_nodes == [expected_end_node]:
            mol_is_present = True

    assert mol_is_present







