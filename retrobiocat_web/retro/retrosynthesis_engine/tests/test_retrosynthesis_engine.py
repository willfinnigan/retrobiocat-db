from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
make_default_connection()

from retrobiocat_web.retro.network_pathway.network import Network

class Test_Retrosynthesis_Engine():

    def test_can_make_a_quick_retrobiocat_network(self):
        network = Network(target_smiles='O=C(O)C(=O)Cc1ccccc1')
        retro = RetrosynthesisEngine(network)
        retro.generate_network(3)
        assert len(network.graph.nodes) > 5

    def test_can_network_and_apply_retrobiocat_single_step(self):
        network = Network(target_smiles='O=C(O)C(=O)Cc1ccccc1')
        retro = RetrosynthesisEngine(network)
        retro.initialise_graph()
        retro.step_retrobiocat(network.target_smiles)
        assert len(network.graph.nodes) > 1

    def test_can_network_and_apply_retrorules_single_step(self):
        network = Network(target_smiles='O=C(O)C(=O)Cc1ccccc1')
        retro = RetrosynthesisEngine(network)
        retro.initialise_graph()
        retro.step_retrorules(network.target_smiles)
        assert len(network.graph.nodes) > 1

    def test_can_network_and_apply_aizynthfinder_single_step(self):
        network = Network(target_smiles='O=C(O)C(=O)Cc1ccccc1')
        retro = RetrosynthesisEngine(network, log_level='DEBUG')
        retro.initialise_graph()
        retro.step_aizynthfinder(network.target_smiles)
        assert len(network.graph.nodes) > 1

    def test_can_network_and_apply_ringbreaker_single_step(self):
        network = Network(target_smiles="c1ccc2c3c([nH]c2c1)CCCCCC3")
        retro = RetrosynthesisEngine(network, log_level='DEBUG')
        retro.initialise_graph()
        retro.step_ringbreaker(network.target_smiles)
        assert len(network.graph.nodes) > 1

    def test_can_network_and_apply_custom_single_step(self):
        network = Network(target_smiles='O=C(O)C(=O)Cc1ccccc1')
        retro = RetrosynthesisEngine(network, log_level='DEBUG')
        retro.initialise_graph()
        product_smis, reaction_names = retro.custom_reaction(network.target_smiles, ['CCC'], 'custom_rxn_1')

        rxn_name = reaction_names[0]

        assert 'CCC' in list(network.graph.successors(rxn_name))
        assert network.target_smiles in list(network.graph.predecessors(rxn_name))

    def test_config_is_same_throughout(self):
        network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
        retro = RetrosynthesisEngine(network, log_level='DEBUG')
        retro.config.source_mol_vendors = ['askcos']
        assert retro.config == retro.custom_expander.mcts_config

    def test_config_changes_carry_through(self):
        network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
        retro = RetrosynthesisEngine(network, log_level='DEBUG')
        retro.config.source_mol_vendors = ['askcos']

        retro.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ["O=C1CCC[C@@H](c2ccccc2)N1"], 'test_rxn')
        assert network.graph.nodes["O=C1CCC[C@@H](c2ccccc2)N1"]['attributes']['is_starting_material'] == 0

    def test_retrorules_settings_change_to_promiscuous_gives_reactions(self):
        network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
        retro = RetrosynthesisEngine(network, log_level='WARNING')
        retro.config.rr_combined_score_threshold = 0
        retro.config.rr_score_threshold = 0
        retro.config.rr_threshold = 0
        retro.config.rr_diameter = 2

        product_smis, reaction_names = retro.step_retrorules(network.target_smiles)
        assert len(reaction_names) > 20


