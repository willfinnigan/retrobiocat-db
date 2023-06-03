import pytest

from retrobiocat_web.mongo.default_connection import make_default_connection
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder_expander import AIZynthfinder_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.custom_expander import CustomRxn_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrobiocat_expander import RetroBioCat_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrorules_expander import RetroRules_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.ringbreaker_expander import RingBreaker_Expander
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
import time

make_default_connection()

class Test_RetroBioCat_Expander():

    def test_complexity_ranking_is_in_correct_order(self):
        network = Network(target_smiles="NCCCCC(=O)c1ccccc1")
        rbc = RetroBioCat_Expander(network)
        reactions = rbc.expand(network.target_smiles)
        assert reactions[0].name == 'Ene reduction'
        assert reactions[1].name == 'Aldehyde amination'


    def test_multistep_rxns_loaded(self):
        network = Network(target_smiles="CCNCc1ccccc1")

        config = RetrosynthesisConfig()
        config.include_two_step = True

        rbc = RetroBioCat_Expander(network, config=config)

        assert rbc.rxn_class.multi_step_rxns != {}

    def test_multistep_reductive_amination_is_applied(self):
        target = "CCNCc1ccccc1"
        network = Network(target_smiles=target)

        config = RetrosynthesisConfig()
        config.include_two_step = True

        retro_engine = RetrosynthesisEngine(network, config=config)

        retro_engine.initialise_graph()
        retro_engine.step_retrobiocat(network.target_smiles)

        assert 'Reductive amination' in [network.graph.nodes[r]['attributes']['name'] for r in network.reaction_nodes]

    def test_can_create_options(self):
        network = Network(target_smiles='CCO')
        expander = RetroBioCat_Expander(network)
        options = expander.get_expansion_options("CCO")
        assert options[0].smi == "CCO"

    def test_options_can_be_executed_to_create_a_reaction(self):
        network = Network(target_smiles='CCO')
        expander = RetroBioCat_Expander(network)
        options = expander.get_expansion_options("CCO")
        reactions = options[0].get_processed_reactions()
        assert reactions[0].reactant.smi == "CCO"

    def test_getting_reactions_from_option_is_very_fast_second_time(self):
        network = Network(target_smiles='CCO')
        expander = RetroBioCat_Expander(network, log_level='WARNING')

        t0 = time.time()
        options = expander.get_expansion_options("CCO")
        reactions = []
        for opt in options:
            reactions.extend(opt.get_processed_reactions())
        t1 = time.time()

        options = expander.get_expansion_options("CCO")
        reactions = []
        for opt in options:
            reactions += opt.get_processed_reactions()
        t2 = time.time()

        time_1 = round(t1 - t0, 4)
        time_2 = round(t2 - t1, 4)

        assert time_2 < 0.01

    def test_that_options_are_in_order_of_complexity_score(self):
        network = Network(target_smiles='NCCCCC(=O)c1ccccc1')
        expander = RetroBioCat_Expander(network)
        options = expander.get_expansion_options("NCCCCC(=O)c1ccccc1")
        assert options[0].score <= options[1].score
        assert options[1].score <= options[2].score
        assert options[2].score <= options[3].score
        assert options[3].score <= options[4].score

class Test_AIZynthfinder_Expander():

    def test_simple_expand(self):
        network = Network(target_smiles="CCO")
        config = RetrosynthesisConfig()
        expander = AIZynthfinder_Expander(network, config=config)
        reactions = expander.expand("CCO")
        assert reactions[0].reactant.smi == 'CCO'

    def test_group_reactions(self):
        network = Network(target_smiles="c1ccc([C@@H]2CCCCN2)cc1")
        config = RetrosynthesisConfig()
        config.aizynth_cutoff_number = 4
        config.aizynth_group_duplicates = True
        expander = AIZynthfinder_Expander(network, config=config)
        reactions = expander.expand("c1ccc([C@@H]2CCCCN2)cc1")
        assert len(reactions) == 2

    def test_can_create_options(self):
        network = Network(target_smiles='CCO')
        expander = AIZynthfinder_Expander(network)
        options = expander.get_expansion_options("CCO")
        assert options[0].smi == "CCO"

    def test_options_can_be_executed_to_create_a_reaction(self):
        network = Network(target_smiles='CCO')
        expander = AIZynthfinder_Expander(network)
        options = expander.get_expansion_options("CCO")
        reactions = options[0].get_processed_reactions()
        assert reactions[0].reactant.smi == "CCO"

    def test_getting_reactions_from_option_is_very_fast_second_time(self):
        network = Network(target_smiles='CCO')
        expander = AIZynthfinder_Expander(network, log_level='WARNING')

        t0 = time.time()
        options = expander.get_expansion_options("CCO")
        reactions = []
        for opt in options:
            reactions.extend(opt.get_processed_reactions())
        t1 = time.time()

        options = expander.get_expansion_options("CCO")
        reactions = []
        for opt in options:
            reactions += opt.get_processed_reactions()
        t2 = time.time()

        time_1 = round(t1 - t0, 4)
        time_2 = round(t2 - t1, 4)

        assert time_2 < 0.01

    def test_that_options_are_in_order_of_policy_score(self):
        network = Network(target_smiles='CCO')
        expander = AIZynthfinder_Expander(network)
        options = expander.get_expansion_options("CCO")
        assert options[0].score >= options[1].score
        assert options[1].score >= options[2].score
        assert options[2].score >= options[3].score
        assert options[3].score >= options[4].score
        assert options[5].score >= options[6].score

class Test_RingBreaker_Expander():

    def test_simple_expand(self):
        target = "c1ccc2c3c([nH]c2c1)CCCCCC3"
        network = Network(target_smiles=target)
        config = RetrosynthesisConfig()
        expander = RingBreaker_Expander(network, config=config)
        reactions = expander.expand(target)
        assert reactions[0].reactant.smi == target

    def test_can_create_options(self):
        target = "c1ccc2c3c([nH]c2c1)CCCCCC3"
        network = Network(target_smiles=target)
        expander = RingBreaker_Expander(network)
        options = expander.get_expansion_options(target)
        assert options[0].smi == target

    def test_options_can_be_executed_to_create_a_reaction(self):
        target = "c1ccc2c3c([nH]c2c1)CCCCCC3"
        network = Network(target_smiles=target)
        expander = RingBreaker_Expander(network)
        options = expander.get_expansion_options(target)
        reactions = options[0].get_processed_reactions()
        assert reactions[0].reactant.smi == target

    def test_getting_reactions_from_option_is_very_fast_second_time(self):
        target = "c1ccc2c3c([nH]c2c1)CCCCCC3"
        network = Network(target_smiles=target )
        expander =  RingBreaker_Expander(network, log_level='WARNING')

        t0 = time.time()
        options = expander.get_expansion_options(target)
        reactions = []
        for opt in options:
            reactions.extend(opt.get_processed_reactions())
        t1 = time.time()

        options = expander.get_expansion_options(target )
        reactions = []
        for opt in options:
            reactions += opt.get_processed_reactions()
        t2 = time.time()

        time_1 = round(t1 - t0, 4)
        time_2 = round(t2 - t1, 4)

        assert time_2 < 0.01

    def test_that_options_are_in_order_of_policy_score(self):
        target = "c1ccc2c3c([nH]c2c1)CCCCCC3"
        network = Network(target_smiles=target)
        expander =  RingBreaker_Expander(network)
        options = expander.get_expansion_options(target)
        assert options[0].score >= options[1].score

class Test_RetroRules_Expander():

    def test_simple_expand(self):
        network = Network(target_smiles="CCO")
        config = RetrosynthesisConfig()
        expander = RetroRules_Expander(network, config=config)
        reactions = expander.expand("CCO")
        assert reactions[0].reactant.smi == 'CCO'

    def test_can_create_options(self):
        network = Network(target_smiles='CCO')
        expander = RetroRules_Expander(network)
        options = expander.get_expansion_options("CCO")
        assert options[0].smi == "CCO"

    def test_options_can_be_executed_to_create_a_reaction(self):
        network = Network(target_smiles='CCCC=O')
        expander = RetroRules_Expander(network)
        options = expander.get_expansion_options("CCCC=O")
        reactions = []
        for opt in options:
            reactions += opt.get_processed_reactions()
        assert reactions[0].reactant.smi == "CCCC=O"

    def test_getting_reactions_from_option_is_very_fast_second_time(self):
        network = Network(target_smiles='CCCC=O')
        expander = RetroRules_Expander(network)

        t0 = time.time()
        options = expander.get_expansion_options("CCCC=O")
        reactions = []
        for opt in options:
            reactions.extend(opt.get_processed_reactions())
        t1 = time.time()

        options = expander.get_expansion_options("CCCC=O")
        reactions = []
        for opt in options:
            reactions += opt.get_processed_reactions()
        t2 = time.time()

        time_1 = round(t1 - t0, 4)
        time_2 = round(t2 - t1, 4)
        assert time_2 < 0.01

    """
    def test_that_options_are_in_order_of_policy_score(self):
        network = Network(target_smiles='CCO')
        expander = RetroRules_Expander(network)
        options = expander.get_expansion_options("CCO")
        assert options[0].score <= options[1].score
        assert options[1].score <= options[2].score
        assert options[2].score <= options[3].score
        assert options[3].score <= options[4].score
        assert options[5].score <= options[6].score
    """
