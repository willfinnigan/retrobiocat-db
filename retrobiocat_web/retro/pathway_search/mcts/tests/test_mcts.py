import pytest

import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.retro.network_pathway import node_analysis
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.mcts.mcts import MCTS
from retrobiocat_web.retro.pathway_search.mcts.tree_node import MCTS_Node

from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine

make_default_connection()


class Test_MCTS():

    def test_can_make_root_node(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network)
        assert mcts.root.pathway_nodes == [network.target_smiles]

    def test_can_expand_root_node(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network)
        mcts.root.expand()
        assert len(mcts.root.children) > 0

    def test_selection_of_a_child_executes_option_and_adds_to_graph(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.root.expand()
        child = mcts.root.select(2)
        assert len(network.reaction_nodes) == 1

    def test_can_expand_and_select_multiple_times_manually(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.root.expand()
        child = mcts.root.select(2)

        child.expand()
        new_child = child.select(2)
        assert new_child.depth == 2

    def test_mcts_run_goes_round_multiple_times(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.config.max_search_time = 5
        mcts.run()

    @pytest.mark.parametrize('max_length,search_time', [[1,5], [2,5], [3,10], [4,20]])
    def test_mcts_reaches_but_does_not_go_beyond_max_length(self, max_length, search_time):
        target_smiles = 'c1ccc([C@@H]2CCCCN2)cc1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.config.max_search_time = search_time
        mcts.config.maxLength = max_length
        mcts.run()
        pathways = mcts.get_pathways()

        lengths = []
        for pathway in pathways:
            length = retrobiocat_web.retro.network_pathway.graph_functions.get_pathway_length()
            lengths.append(length)
            if length > mcts.config.maxLength:
                print(pathway.aizynthfinder_tree())

        assert max(lengths) <= mcts.config.maxLength

    def test_aizynthfinder_doesnt_add_duplicate_reactions(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='DEBUG')
        mcts.config.biosynthesis_expansion = False
        mcts.config.biocatalysis_expansion = False
        mcts.config.chemistry_expansion = True
        mcts.config.maxLength = 1
        mcts.config.max_search_time = 10
        mcts.run()

    def test_only_1_pathway_when_max_length_is_zero(self):
        target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.config.biosynthesis_expansion = False
        mcts.config.biocatalysis_expansion = False
        mcts.config.chemistry_expansion = True
        mcts.config.maxLength = 0
        mcts.config.max_search_time = 2
        mcts.run()
        pathways = mcts.get_pathways()

        assert len(pathways) == 1


    def test_retrorules_promiscious_finds_routes(self):
        target_smiles = 'c1ccc([C@@H]2CCCCN2)cc1'
        network = Network(target_smiles=target_smiles)
        mcts = MCTS(network, log_level='WARNING')
        mcts.config.biosynthesis_expansion = True
        mcts.config.biocatalysis_expansion = False
        mcts.config.chemistry_expansion = False
        mcts.config.maxLength = 1
        mcts.config.max_search_time = 20

        mcts.retro_engine.mcts_config.rr_combined_score_threshold = 0
        mcts.retro_engine.mcts_config.rr_score_threshold = 0
        mcts.retro_engine.mcts_config.rr_threshold = 0
        mcts.retro_engine.mcts_config.rr_diameter = 2

        mcts.run()
        pathways = mcts.get_pathways()

        assert len(pathways) > 4







