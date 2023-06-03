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
def test_pathway_2():
    network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
    retro_engine = RetrosynthesisEngine(network)
    rxns, smis = [], ['c1ccc([C@@H]2CCCCN2)cc1']
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ['CCOC(=O)N1CCCC[C@H]1c1ccccc1'],
                                                                'Chem_Unassigned*_4')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('CCOC(=O)N1CCCC[C@H]1c1ccccc1', ['CCOC(=O)N1CCC=C[C@H]1c1ccccc1'],
                                                                'Chem_Alkene to alkane*_1')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('CCOC(=O)N1CCC=C[C@H]1c1ccccc1', ['CCOC(=O)Cl', 'C1=C[C@@H](c2ccccc2)NCC1'],
                                                                'Chem_N-acylation to amide*_1')
    rxns += reaction_names
    smis += product_smis

    product_smis, reaction_names = retro_engine.custom_reaction('C1=C[C@@H](c2ccccc2)NCC1',
                                                                ['C=CCCN[C@@H](C=C)c1ccccc1'],
                                                                'Other organometallic C-C bond formation')
    rxns += reaction_names
    smis += product_smis

    product_smis, reaction_names = retro_engine.custom_reaction('C=CCCN[C@@H](C=C)c1ccccc1',
                                                                ['C=CCCN=C(C=C)c1ccccc1'],
                                                                'Imine reduction')
    rxns += reaction_names
    smis += product_smis

    product_smis, reaction_names = retro_engine.custom_reaction('C=CCCN=C(C=C)c1ccccc1',
                                                                ['C=CC(=O)c1ccccc1', 'C=CCCN'],
                                                                'Imine formation')
    rxns += reaction_names
    smis += product_smis

    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')

    return pathway


