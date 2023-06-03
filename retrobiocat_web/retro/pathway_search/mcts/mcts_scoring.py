from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.mcts.config import MCTS_Config
from retrobiocat_web.retro.pathway_search.mcts.tree_node import MCTS_Node


class MCTSScorer():
    """ This module will score pathways for use in UCB score"""

    def __init__(self, network: Network, config=None):
        self.config = config
        if self.config is None:
            self.config = MCTS_Config()

        self.network = network

    def score(self, node: MCTS_Node):
        """ A quick score to get things working"""

        # +1 if buyable
        # If not buyable, use complexity to score between -0.1 and -1.
        buyable, not_buyable = node._get_buyable_end_nodes()
        buyable_score = 0
        for node in buyable:
            buyable_score += 1
        for node in not_buyable:
            if self.config.adjust_non_buyable_penalties_with_complexity:
                buyable_score += self._get_non_buyable_penalty_from_complexity(node)
            else:
                buyable_score += -1

        buyable_score = buyable_score / (len(buyable) + len(not_buyable))

        return buyable_score

    def _get_non_buyable_penalty_from_complexity(self, smi):
        min_penalty = self.config.min_penalty_for_non_buyable
        rel_complexity = self.network.graph.nodes[smi]['attributes']['relative_complexity']
        if rel_complexity > self.config.complexity_for_minus_1:
            penalty = -1
        elif rel_complexity < self.config.complexity_for_zero:
            penalty = min_penalty
        else:
            penalty = ((-1 - rel_complexity)/-self.config.complexity_for_zero) + min_penalty

        if penalty < -1:
            penalty = -1
        return penalty







