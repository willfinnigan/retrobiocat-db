import pytest
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.retro.retrosynthesis_engine.graph_control.graph_adder import GraphAdder
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import MoleculeFactory
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_creator import ReactionCreator
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.components.reaction_processor import \
    ReactionProcessor
from retrobiocat_web.retro.retrosynthesis_engine.reaction_processing.reaction_ranker import ReactionRanker

make_default_connection()


@pytest.fixture()
def backwards_reaction():
    network = Network(target_smiles='CCC=O')

    network.graph.add_node('CCC=O', attributes={'node_type': 'substrate'})
    network.graph.add_node('CCCO', attributes={'node_type': 'substrate'})
    network.graph.add_node('reac_1', attributes={'node_type': 'reaction', 'name': 'reac_1'})
    network.graph.add_edge('CCC=O', "reac_1")
    network.graph.add_edge('reac_1', "CCCO")

    output = {'reac_1': [["CCC=O"], ['CCCN']]}
    target = "CCCO"

    return network, output, target

@pytest.fixture()
def network_with_existing_reaction():
    network = Network(target_smiles='CCC=O')

    network.graph.add_node('CCCN', attributes={'node_type': 'substrate'})
    network.graph.add_node('CCC=O', attributes={'node_type': 'substrate'})
    network.graph.add_node('reac_1', attributes={'node_type': 'reaction', 'name': 'reac_1'})
    network.graph.add_edge('CCCN', "reac_1")
    network.graph.add_edge('reac_1', "CCC=O")
    return network


class Test_reaction_processing():

    def test_reactions_are_created(self):
        network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)
        assert len(reactions) == 2

    def test_reactions_are_created_with_metadata(self):
        network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}
        metadata = {'Imine reduction': {'test_metadata': 'hello'}}

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output, metadata=metadata)
        assert reactions[0].metadata['test_metadata'] == 'hello'

    def test_reactions_are_created_with_complexity(self,):
        network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)
        assert reactions[0].complexity_change < 0

    def test_reactions_with_no_products_are_not_created(self,):
        network = Network(target_smiles='CCC=O')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [[]]}
        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)
        assert len(reactions) == 0

    def test_N_can_be_removed(self):
        network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})

        output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}
        reaction_creator = ReactionCreator(network, log_level='DEBUG')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reaction_processor.config.remove_small_mols = True
        reactions = reaction_processor.process(reactions)

        assert len(reactions[0].substrates) == 1

    def test_duplicate_reactions_are_dropped(self, network_with_existing_reaction):
        network = network_with_existing_reaction

        output = {'reac_1': [['CCC=O']]}
        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions('CCCN', output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reactions = reaction_processor.process(reactions)

        assert len(reactions) == 0

    def test_duplicate_reaction_with_different_name_is_not_dropped(self, network_with_existing_reaction):
        network = network_with_existing_reaction

        output = {'a_different_name': [['CCC=O']]}
        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions('CCCN', output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reactions = reaction_processor.process(reactions)

        assert len(reactions) == 1

    def test_duplicate_reaction_with_different_name_is_dropped_if_allowed(self, network_with_existing_reaction):
        network = network_with_existing_reaction

        output = {'a_different_name': [['CCC=O']]}
        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions('CCCN', output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reactions = reaction_processor.process(reactions, duplicates_must_match_name=False)

        assert len(reactions) == 0

    def test_backwards_reactions_are_dropped(self, backwards_reaction):
        network, output, target = backwards_reaction

        reaction_creator = ReactionCreator(network, log_level='DEBUG')
        reactions = reaction_creator.create_new_reactions(target, output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reactions = reaction_processor.process(reactions)

        assert len(reactions) == 1

    def test_allow_backwards_reactions_does_not_drop(self, backwards_reaction):
        network, output, target = backwards_reaction

        reaction_creator = ReactionCreator(network, log_level='DEBUG')
        reactions = reaction_creator.create_new_reactions(target, output)

        reaction_processor = ReactionProcessor(network, log_level='DEBUG')
        reaction_processor.config.allow_backwards = True
        reactions = reaction_processor.process(reactions)

        assert len(reactions) == 2

    def test_reactions_which_are_not_backwards_remain(self, backwards_reaction):
        network, output, target = backwards_reaction

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(target, output)

        reaction_processor = ReactionProcessor(network, log_level='INFO')
        reactions = reaction_processor.process(reactions)

        assert [s.smi for s in reactions[0].substrates] == ['CCCN']

    def test_cyclic_reactions_are_removed(self):
        network = Network(target_smiles='CCC=O')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [['CCC=O']]}

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

        reaction_processor = ReactionProcessor(network, log_level='INFO')
        reactions = reaction_processor.process(reactions)

        assert len(reactions) == 0

    def test_disallowed_products_can_be_removed(self):
        network = Network(target_smiles='[C@H]1(C2=CC=CC=C2)NCCCC1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Imine reduction': [['c1ccc(C2=NCCCC2)cc1', 'N'], ['C1=N[C@H](c2ccccc2)CCC1', 'N']]}
        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

        reaction_processor = ReactionProcessor(network, log_level='INFO')
        reactions = reaction_processor.process(reactions, disallowed_products=['C1=N[C@H](c2ccccc2)CCC1'])

        assert len(reactions) == 1

    def test_can_rank_reactions_by_complexity(self):
        network = Network(target_smiles='O=Cc1ccccc1')
        network.graph.add_node(network.target_smiles, attributes={'node_type': 'substrate'})
        output = {'Bad complexity change': [['[C@H]1(C2=CC=CC=C2)NCCCC1']],
                  'Good complexity change': [['C']],
                  'No complexity change': [['OCc1ccccc1']],
                  }

        reaction_creator = ReactionCreator(network, log_level='INFO')
        reactions = reaction_creator.create_new_reactions(network.target_smiles, output)

        reaction_processor = ReactionProcessor(network, log_level='INFO')
        reactions = reaction_processor.process(reactions)

        reaction_ranker = ReactionRanker(log_level='INFO')
        reactions_ordered = reaction_ranker.rank_by_complexity(reactions)
        assert [r.name for r in reactions_ordered] == ['Good complexity change', 'No complexity change', 'Bad complexity change']





