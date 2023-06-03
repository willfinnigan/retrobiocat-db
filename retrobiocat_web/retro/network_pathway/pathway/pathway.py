import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.network_pathway.pathway import pathway_scoring, pathway_functions



class Pathway(object):

    def __init__(self, list_nodes, network, substrates=None, reactions=None, end_nodes=None,
                 calc_scores=True, attr_dict=None, reaction_enzyme_map=None, log_level='WARNING'):

        self.logger = add_logger("Pathway Obj.", level=log_level)

        # pathway nodes in all forms
        self.list_nodes = list_nodes
        self.network = network
        self.target_smiles = network.target_smiles
        self.sub_graph = pathway_functions.get_pathway_sub_graph(self.network.graph, list_nodes)
        self.substrates = substrates
        self.reactions = reactions
        self.end_nodes = end_nodes
        self._get_nodes_if_none()

        # for loading attributes for pre-made pathway
        self.attributes = {}
        if attr_dict == None:
            self.get_attributes_from_network()
        else:
            self.set_attributes_from_dict(attr_dict)

        # specifies which enzyme for each step
        if reaction_enzyme_map == None:
            self.set_reaction_enzyme_map_from_network()
        else:
            self.reaction_enzyme_map = reaction_enzyme_map

        self.reaction_names = self._get_reaction_names()

        # scores
        self.scores = pathway_scoring.PathwayScores(self, calc_scores=calc_scores)

        self.costs = {}

        # used by grouping for pathway explorer
        self.other_varients = []
        self.other_varients_as_nodes = []

        self.tree = None
        self.pathway_length = -1

    def _get_nodes_if_none(self):
        self.substrates = retrobiocat_web.retro.network_pathway.graph_functions.get_substrate_nodes(self.sub_graph)
        self.reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes(self.sub_graph)
        self.end_nodes = retrobiocat_web.retro.network_pathway.graph_functions.get_nodes_with_no_successors(self.sub_graph)

        self.logger.debug(f"End nodes for pathway are: {self.end_nodes}")

    def __repr__(self):
        """ Returns the pathway nodes as a string - useful for dataframe"""
        return str(self.list_nodes)

    def does_pathway_end_in_spontaneous_chemical_step(self):
        term_rxns = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes_of_list_substrates(self.network.graph, self.end_nodes)
        for rxn in term_rxns:
            if self.network.graph.nodes[rxn]['attributes']['metadata']['selected_enzyme'] == 'Chemical':
                return True
        return False

    def get_attributes_from_network(self):
        """ Gets the attributes of each node from the network graph """
        for node in self.list_nodes:
            self.attributes[node] = self.network.graph.nodes[node]['attributes']

    def set_attributes_from_dict(self, attr_dict):
        """ Loads attributes from the passed dictionary """
        for node in attr_dict:
            self.network.graph.nodes[node]['attributes'] = attr_dict[node]

    def set_reaction_enzyme_map_from_network(self):
        """ Set dict of which enzyme for each step from default options in network """
        self.reaction_enzyme_map = {}
        for node in self.reactions:
            self.reaction_enzyme_map[node] = ''
            if self.network.graph.nodes[node]['attributes']['reaction_type'] == 'retrobiocat':
                self.reaction_enzyme_map[node] = self.network.graph.nodes[node]['attributes']['metadata']['selected_enzyme']

    def _get_reaction_names(self):
        names = []
        for reaction in self.reactions:
            names.append(self.network.graph.nodes[reaction]['attributes']['name'])
        return names

    def get_pathway_length(self):
        #if self.pathway_length == -1:
        #    self.tree = self._get_tree(self.target_smiles)
        #return self.pathway_length
        #self.pathway_length = node_analysis.get_pathway_length(self.target_smiles, self.sub_graph)
        self.pathway_length = retrobiocat_web.retro.network_pathway.graph_functions.get_current_pathway_length(self.sub_graph, self.end_nodes, self.target_smiles)
        return self.pathway_length

    def end_node_distances(self):
        distances = []
        for node in self.end_nodes:
            dist = retrobiocat_web.retro.network_pathway.graph_functions.get_num_reactions_from_target(self.sub_graph, node, self.target_smiles)
            distances.append(dist)
        return distances

    def _calculate_pathway_length(self):
        if self.tree is None:
            self.tree = self._get_tree(self.target_smiles)

    def aizynthfinder_tree(self):
        if self.tree is None:
            self.tree = self._get_tree(self.target_smiles)
        self.logger.debug(f"Tree = {self.tree}")
        return self.tree

    def _get_tree(self, smi, depth=0):
        tree = {'smiles': smi, 'depth': depth, 'children': []}
        child_smis = retrobiocat_web.retro.network_pathway.graph_functions.substrate_precursors(self.sub_graph, smi)
        depth += 1
        for c_smi in child_smis:
            tree['children'].append(self._get_tree(c_smi, depth=depth))

        if self.pathway_length < depth - 1:
            self.pathway_length = depth - 1

        return tree






if __name__ == '__main__':
    from retrobiocat_web.retro.network_pathway.network import Network
    from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()


    network = Network(target_smiles='c1ccc([C@@H]2CCCCN2)cc1')
    retro_engine = RetrosynthesisEngine(network)
    rxns, smis = [], ['c1ccc([C@@H]2CCCCN2)cc1']
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc([C@@H]2CCCCN2)cc1', ['c1ccc(C2=NCCCC2)cc1'], 'Imine reduction')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('c1ccc(C2=NCCCC2)cc1', ['NCCCCC(=O)c1ccccc1'], 'Imine formation')
    rxns += reaction_names
    smis += product_smis
    product_smis, reaction_names = retro_engine.custom_reaction('NCCCCC(=O)c1ccccc1', ['O=CCCCC(=O)c1ccccc1'], 'Aldehyde amination')
    rxns += reaction_names
    smis += product_smis

    pathway = Pathway(rxns + smis,
                      network, log_level='DEBUG')
