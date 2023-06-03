import retrobiocat_web.retro.network_pathway.graph_functions
import retrobiocat_web.retro.network_pathway.rdkit_utils
from retrobiocat_web.logging import add_logger
import networkx as nx

from retrobiocat_web.retro.retrosynthesis_engine.graph_control.graph_adder import GraphAdder
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import MoleculeFactory
from retrobiocat_web.retro.visualisation.visualise import Visualiser


class Network(object):

    def __init__(self, graph=None, target_smiles=None, log_level='WARNING'):
        self.logger = add_logger("Network", level=log_level)

        self.substrate_nodes = []
        self.reaction_nodes = []

        self.target_smiles = target_smiles

        if self.target_smiles is not None:
            self.target_smiles = retrobiocat_web.retro.network_pathway.rdkit_utils.rdkit_smile(target_smiles, warning=True)
            self.logger.debug(f"Target smiles entered as rdkit smiles {self.target_smiles} (original={target_smiles})")

        self.graph = graph
        if self.graph is None:
            self.graph = nx.DiGraph()

            if self.target_smiles is not None:
                self.initialise_graph()

    def initialise_graph(self):
        self.logger.debug(f"Initialising graph")
        target_mol = MoleculeFactory().create_target_mol(self.target_smiles)
        GraphAdder(self).add_starting_mol(target_mol)
        self.logger.debug(f"Graph nodes = {self.graph.nodes}")

    def attributes_dict(self):
        att_dict = {}
        for node in list(self.graph):
            att_dict[node] = self.graph.nodes[node]['attributes']

        return att_dict

    def add_attributes(self, attributes_dict):
        for node in attributes_dict:
            self.graph.nodes[node]['attributes'] = attributes_dict[node]

        self.substrate_nodes = retrobiocat_web.retro.network_pathway.graph_functions.get_substrate_nodes(self.graph)
        self.reaction_nodes = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes(self.graph)

    def get_biocathub_json(self):
        """ Return the network json, but only the information necessary for biocathub """

        # BioCatHub probably want a list of reactions with the relevant info, so thats what we'll give them.
        return_json = []
        for rxn_node in self.reaction_nodes:
            rxn_json = {}
            rxn_json['reaction'] = self.graph.nodes[rxn_node]['attributes']['name']
            selected_enzyme = self.graph.nodes[rxn_node]['attributes']['metadata']['selected_enzyme']
            rxn_json['enzyme'] = selected_enzyme
            rxn_json['cofactors'] = self.graph.nodes[rxn_node]['attributes']['metadata']['enzyme_cofactors'].get(selected_enzyme, [])
            rxn_json['substrates'] = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_substrates(self.graph, rxn_node)
            rxn_json['products'] = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_products(self.graph, rxn_node)
            return_json.append(rxn_json)

        return return_json

if __name__ == '__main__':
    from retrobiocat_web.retro.enzyme_identification import Specificity_Scorer
    from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    target = 'O=C(O)C(=O)Cc1ccccc1'

    network = Network(target_smiles=target)

    retro = RetrosynthesisEngine(network, log_level='INFO')
    retro.generate_network(4, max_nodes=100)

    scorer = Specificity_Scorer(network)
    scorer.score()

    vis = Visualiser()
    nodes, edges = vis.nodes_edges(network.graph)

    biocathub_json = network.get_biocathub_json()


