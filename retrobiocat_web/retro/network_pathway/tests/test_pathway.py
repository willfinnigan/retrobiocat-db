import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway
from retrobiocat_web.mongo.default_connection import make_default_connection
make_default_connection()

from retrobiocat_web.retro.network_pathway.tests.test_pathway_fixtures import test_pathway_1, test_pathway_2



class Test_Pathway():

    def test_pathway_can_be_created(self):
        network = Network(target_smiles='CCC=O')
        pathway = Pathway(['CCC=O'], network)

        assert 'CCC=O' in pathway.sub_graph.nodes

    def test_substrate_in_pathway_has_attribute_for_starting_material(self):
        network = Network(target_smiles='CCC=O')
        pathway = Pathway(['CCC=O'], network)
        assert 'is_starting_material' in pathway.sub_graph.nodes['CCC=O']['attributes']

    def test_end_nodes_are_terminal(self, test_pathway_1):
        pathway = test_pathway_1
        assert pathway.end_nodes == ['O=CCCCC(=O)c1ccccc1']

    def test_can_nake_aizynthfinder_tree(self, test_pathway_1):
        pathway = test_pathway_1
        tree = pathway.aizynthfinder_tree()

        assert tree == {'smiles': 'c1ccc([C@@H]2CCCCN2)cc1', 'depth': 0,
                        'children': [{'smiles': 'c1ccc(C2=NCCCC2)cc1', 'depth': 1,
                                      'children': [{'smiles': 'NCCCCC(=O)c1ccccc1', 'depth': 2,
                                                    'children': [{'smiles': 'O=CCCCC(=O)c1ccccc1', 'depth': 3,
                                                                  'children': []}]}]}]}


    def test_pathway_length_is_correct(self, test_pathway_1):
        assert retrobiocat_web.retro.network_pathway.graph_functions.get_pathway_length() == 3

    def test_long_pathway_length_is_correct(self, test_pathway_2):
        assert retrobiocat_web.retro.network_pathway.graph_functions.get_pathway_length() == 6

    def test_distances_of_end_nodes(self, test_pathway_2):
        distances = test_pathway_2.end_node_distances()
        assert distances == [3, 6, 6]




