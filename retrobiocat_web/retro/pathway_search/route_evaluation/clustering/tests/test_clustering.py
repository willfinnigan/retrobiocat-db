from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway
from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_clusterer import Route_Clusterer
from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_distancer import Route_Distancer
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
import pytest

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
def test_pathway_2():
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

    product_smis, reaction_names = retro_engine.custom_reaction('O=CCCCC(=O)c1ccccc1', ['O=C(O)CCCC(=O)c1ccccc1'],
                                                                'Carboxylic acid reduction')
    rxns += reaction_names
    smis += product_smis

    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')

    return pathway

@pytest.fixture
def test_pathway_3():
    network = Network(target_smiles='CCCC=O')
    retro_engine = RetrosynthesisEngine(network)
    rxns, smis = [], ['CCCC=O']
    product_smis, reaction_names = retro_engine.custom_reaction('CCCC=O', ['CCCCO'],
                                                                'ALcohol oxidation')
    rxns += reaction_names
    smis += product_smis
    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')

    return pathway

class Test_Get_Distances():

    def test_can_load_model(self):
        distancer = Route_Distancer()
        assert distancer.model is not None

    def test_can_load_trees_with_fingerprints(self, test_pathway_1):
        distancer = Route_Distancer()
        trees = distancer.get_trees([test_pathway_1])
        assert 'fingerprint' in trees[0]['children'][0]

    def test_can_calculate_distances(self, test_pathway_1, test_pathway_2, test_pathway_3):
        distancer = Route_Distancer()
        pathways = [test_pathway_1, test_pathway_2, test_pathway_3]
        distances = distancer.get_distances(pathways)
        assert distances[0][1] < distances[0][2]

    def test_can_cluster(self, test_pathway_1, test_pathway_2, test_pathway_3):
        distancer = Route_Distancer()
        pathways = [test_pathway_1, test_pathway_2, test_pathway_3]
        distances = distancer.get_distances(pathways)
        clusterer = Route_Clusterer()
        labels = list(clusterer.cluster_by_distance(distances, 5))
        assert labels == [0, 0, 1]
