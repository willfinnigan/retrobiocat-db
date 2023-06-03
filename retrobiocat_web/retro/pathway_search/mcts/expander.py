import itertools

import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.pathway_search.mcts.config import MCTS_Config


class NetworkExpander():
    """
    A class for the expansion of the retrosynthetic network
    by the application of templates
    """

    def __init__(self, network, retro_engine, config=None, log_level='WARNING'):

        self.config = config
        if self.config is None:
            self.config = MCTS_Config()

        self.network = network
        self.retro_engine = retro_engine
        self.logger = add_logger('MCTS-Network-Expander', level=log_level)

    def options_from_end_nodes(self, endNodes, pathway_nodes):
        """ Take end nodes, and expand in those positions """
        # need more policy options for how expansion should work
        # or are these just the network settings..

        options = []
        subgraph = self.network.graph.subgraph(pathway_nodes)
        for node in endNodes:
            length_from_target = retrobiocat_web.retro.network_pathway.graph_functions.get_num_reactions_from_target(subgraph, node, self.network.target_smiles)
            if length_from_target < self.config.maxLength:  # only expand if there is room to do so
                biocat_options = self._options_for_biocatalysis(node)
                chemistry_options = self._options_for_chemistry(node)
                biosynth_options = self._options_for_biosynthesis(node)
                options = self._interleave([biocat_options, chemistry_options, biosynth_options])

        self.logger.debug(f"Options = {options}")

        return options, [o.score for o in options]  # option scores

    def _interleave(self, list_of_options):
        interleaved = [x for x in itertools.chain(*itertools.zip_longest(*list_of_options)) if x is not None]
        return interleaved

    def _options_for_biocatalysis(self, node):
        if self.config.biocatalysis_expansion:
            return self.retro_engine.options_retrobiocat(node)
        return []

    def _options_for_chemistry(self, node):
        if self.config.chemistry_expansion:
            return self.retro_engine.options_aizynthfinder(node)
        return []

    def _options_for_biosynthesis(self, node):
        if self.config.biosynthesis_expansion:
            return self.retro_engine.options_retrorules(node)
        return []

if __name__ == '__main__':
    list_of_options = [[1,2,3,4,5,6,7,8,9], ['a', 'b', 'c'], ['x', 'y', 'z']]
    interleaved = [x for x in itertools.chain(*itertools.zip_longest(*list_of_options)) if x is not None]

