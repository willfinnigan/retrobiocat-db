import networkx as nx
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.expanders.aizynthfinder_expander import AIZynthfinder_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.custom_expander import CustomRxn_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrobiocat_expander import RetroBioCat_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.retrorules_expander import RetroRules_Expander
from retrobiocat_web.retro.retrosynthesis_engine.expanders.ringbreaker_expander import RingBreaker_Expander
from retrobiocat_web.retro.retrosynthesis_engine.graph_control.graph_adder import GraphAdder
from retrobiocat_web.retro.retrosynthesis_engine.graph_control.graph_pruner import GraphPruner
from retrobiocat_web.retro.retrosynthesis_engine.mols_and_reactions.molecules_and_reactions import MoleculeFactory


class RetrosynthesisEngine():

    def __init__(self, network, config=None, log_level='WARNING'):
        self.network = network

        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()

        self.logger = add_logger('RetrosynthesisEngine', level=log_level)

        # expanders
        self.retrobiocat_expander = RetroBioCat_Expander(network, config=self.config, log_level=log_level)
        self.retrorules_expander = RetroRules_Expander(network, config=self.config, log_level=log_level)
        self.aizynthfinder_expander = AIZynthfinder_Expander(network, config=self.config, log_level=log_level)
        self.ringbreaker_expander = RingBreaker_Expander(network, config=self.config, log_level=log_level)
        self.custom_expander = CustomRxn_Expander(network, config=self.config, log_level=log_level)

        # graph control
        self.graphAdder = GraphAdder(network, log_level=log_level)
        self.graphPruner = GraphPruner(log_level=log_level)

    def initialise_graph(self):
        self.network.initialise_graph()

    def step_retrobiocat(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return [],[]

        reactions = self.retrobiocat_expander.expand(smi)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def options_retrobiocat(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return []

        return self.retrobiocat_expander.get_expansion_options(smi)

    def step_aizynthfinder(self, smi):

        if self._can_rules_be_applied(smi) == False:
            return [],[]

        reactions = self.aizynthfinder_expander.expand(smi)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def options_aizynthfinder(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return []
        return self.aizynthfinder_expander.get_expansion_options(smi)

    def step_ringbreaker(self, smi):

        if self._can_rules_be_applied(smi) == False:
            return [],[]

        reactions = self.ringbreaker_expander.expand(smi)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def options_ringbreaker(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return []
        return self.ringbreaker_expander.get_expansion_options(smi)

    def step_retrorules(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return [],[]

        reactions = self.retrorules_expander.expand(smi)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def options_retrorules(self, smi):
        if self._can_rules_be_applied(smi) == False:
            return []
        return self.retrorules_expander.get_expansion_options(smi)

    def step_custom_reactions(self, smi, rxns, multi_step_rxns=None):
        if self._can_rules_be_applied(smi) == False:
            return [],[]

        reactions = self.custom_expander.expand_custom_rules(smi, rxns, multi_step_smarts=multi_step_rxns)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def custom_reaction(self, product_smiles, substrate_smiles, reaction_name):
        reactions = self.custom_expander.expand(product_smiles, substrate_smiles, reaction_name)
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def delete_reaction_node(self, reaction_node_to_remove):
        deleted = self.graphPruner.delete_reaction_node(self.network, reaction_node_to_remove)
        return deleted

    def generate_network(self, number_steps, max_nodes=False):
        self.initialise_graph()

        listSmiles = [self.network.target_smiles]
        for i in range(number_steps):
            newListSmiles, newListReactions = [], []
            for smi in listSmiles:
                newSmiles, newReactions = self.step_retrobiocat(smi)
                newListSmiles.extend(newSmiles)
                newListReactions.extend(newReactions)

            self.logger.info('-- Step ' + str(i + 1) + ' --')
            self.logger.info(str(len(newListReactions)) + ' reactions added, ' + str(len(self.network.substrate_nodes)) + ' substrate nodes')

            if max_nodes is not False:
                listSmiles = self.graphPruner.prune(self.network, newListSmiles, max_nodes)
            else:
                listSmiles = newListSmiles

    def add_reactions(self, reactions):
        product_smis, reaction_names = self.graphAdder.add_reactions(reactions)
        return product_smis, reaction_names

    def _can_rules_be_applied(self, smi):
        """ Returns true if the provided smiles is in the graph as a substrate"""
        if smi not in list(self.network.graph.nodes):
            self.logger.error('WARNING SMILE NOT IN GRAPH ' + str(smi))
            return False
        if self.network.graph.nodes[smi]['attributes']['node_type'] != 'substrate':
            self.logger.error('Warning - target SMILES was not a substrate')
            self.network.graph.nodes[smi]['attributes']['node_type'] = 'substrate'
        return True


if __name__ == '__main__':
    pass


# O=C(O)C(=O)Cc1ccccc1
# O=C(O)C(=O)Cc1ccccc1


