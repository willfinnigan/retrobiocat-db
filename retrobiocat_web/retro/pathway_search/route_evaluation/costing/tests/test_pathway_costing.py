import pytest

from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway
from retrobiocat_web.retro.pathway_search.route_evaluation.costing.pathway_cost_modal import Pathway_Coster
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine

from retrobiocat_web.mongo.default_connection import make_default_connection
make_default_connection()


@pytest.fixture
def test_pathway_1():
    network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
    retro_engine = RetrosynthesisEngine(network)
    rxns, smis = [], ['c1ccc([C@@H]2CCCCN2)cc1']
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ['c1ccc(C2=NCCCC2)cc1'],
                                                                'Imine reduction')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc(C2=NCCCC2)cc1', ['NCCCCC(=O)c1ccccc1'],
                                                                'Imine formation')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('NCCCCC(=O)c1ccccc1', ['O=CCCCC(=O)c1ccccc1'],
                                                                'Aldehyde amination')
    rxns += reaction_names
    smis += product_smis

    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')

    return pathway

@pytest.fixture
def test_pathway_in_stock():
    network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
    retro_engine = RetrosynthesisEngine(network)

    rxns, smis = [], ['c1ccc([C@@H]2CCCCN2)cc1']
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ['c1ccc(C2=NCCCC2)cc1'],
                                                                'Imine reduction')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc(C2=NCCCC2)cc1', ['NCCCCC(=O)c1ccccc1'],
                                                                'Imine formation')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('NCCCCC(=O)c1ccccc1', ['O=CCCCC(=O)c1ccccc1'],
                                                                'Aldehyde amination')
    rxns += reaction_names
    smis += product_smis
    retro_engine.custom_reaction('O=CCCCC(=O)c1ccccc1', ['O=C(O)CCCC(=O)c1ccccc1'],
                                 'Carboxylic acid reduction')
    rxns += reaction_names
    smis += product_smis

    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')

    return pathway


class Test_Pathway_Cost_Modal():
    """
    Pathway cost model
    - Cost per substrate (either assume, or use askcos)
    - Assume yield per reaction
    - Assume cost per reaction
    """

    def test_that_a_cost_can_be_calculated_for_a_pathway(self):
        network = Network(target_smiles='CCC=O')
        pathway = Pathway(['CCC=O'], network)

        coster = Pathway_Coster()
        cost = coster.cost(pathway)
        print(cost)
        assert cost > 0

    def test_multi_step_pathway_is_costed(self, test_pathway_1):
        pathway = test_pathway_1
        coster = Pathway_Coster()
        cost = coster.cost(pathway)
        print(cost)
        assert cost > 0

    def test_costs_saved_to_pathway(self, test_pathway_1):
        pathway = test_pathway_1
        coster = Pathway_Coster()
        cost = coster.cost(pathway)
        assert pathway.costs[coster.name] == cost

    def test_ordering_of_clusters_by_cost(self, test_pathway_1):
        network = Network(target_smiles='CCC=O')
        test_pathway_2 = Pathway(['CCC=O'], network)
        coster = Pathway_Coster()
        clusters = [[test_pathway_1], [test_pathway_2]]

        ordered_clusters = coster.order_clusters_by_cost(clusters)
        assert ordered_clusters[0][0].costs['default'] < ordered_clusters[1][0].costs['default']