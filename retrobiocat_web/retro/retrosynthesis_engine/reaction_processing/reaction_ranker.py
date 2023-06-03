from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig


class ReactionRanker(object):

    def __init__(self, config=None, log_level='WARNING'):
        self.config = config
        if self.config is None:
            self.config = RetrosynthesisConfig()
        self.logger = add_logger('ReactionRanker', level=log_level)

    def rank_by_complexity(self, reactions, keep_top=None):
        self.logger.info(f"Sorting {len(reactions)} reactions by complexity")
        complexities = [r.complexity_change for r in reactions]
        reactions_sorted = self._sort_by_score(complexities, reactions, reverse=False)

        if keep_top is not None:
            reactions_sorted = reactions_sorted[:keep_top]
            self.logger.info(f"..kept the top {keep_top} reactions")

        return reactions_sorted

    def rank_by_aizynthfinder_policy_score(self, reactions, keep_top):
        self.logger.info(f"Sorting {len(reactions)} aizynthfinder reactions by their policy scores")
        policy_scores = [r.metadata['policy_probability'] for r in reactions]
        reactions_sorted = self._sort_by_score(policy_scores, reactions, reverse=True)

        if keep_top is not None:
            reactions_sorted = reactions_sorted[:keep_top]
            self.logger.info(f"..kept the top {keep_top} reactions")

        return reactions_sorted

    def _sort_by_score(self, scores, nodes, reverse=False):
        scores_nodes = list(zip(scores, nodes))
        sorted_scores_nodes = sorted(scores_nodes, key=lambda x: x[0], reverse=reverse)
        self.logger.debug(f"{sorted_scores_nodes}")
        nodes_sorted = [node for score, node in sorted_scores_nodes]
        return nodes_sorted
