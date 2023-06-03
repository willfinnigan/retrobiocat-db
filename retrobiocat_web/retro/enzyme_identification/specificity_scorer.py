import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.retro.molecular_similarity_search.molecular_similarity import Similarity_Substrate_Specificity_Searcher
from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config


class Specificity_Scorer():
    similarity_searcher = Similarity_Substrate_Specificity_Searcher()

    def __init__(self, network, config=None):
        self.config = config
        if self.config is None:
            self.config = Specificity_Scorer_Config()

        self.network = network

    def score(self):
        if self.config.run_similarity_search:
            list_reaction_nodes = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes(self.network.graph)
            for node in list_reaction_nodes:
                if self.network.graph.nodes[node]['attributes']['reaction_type'] == 'retrobiocat':
                    if 'specificity_scores' not in self.network.graph.nodes[node]['attributes']['metadata']:
                        self._score_reaction(node)
                        self._select_best_enzyme(node)

    def _score_reaction(self, reaction_node):
        self.network.graph.nodes[reaction_node]['attributes']['metadata']['specificity_scores'] = {}
        self.network.graph.nodes[reaction_node]['attributes']['metadata']['enzyme_info'] = {}

        possible_enzymes = self.network.graph.nodes[reaction_node]['attributes']['metadata']['possible_enzymes']
        reaction_name = self.network.graph.nodes[reaction_node]['attributes']['name']
        product = list(self.network.graph.predecessors(reaction_node))[0]
        subOne, subTwo = None, None
        if not self.config.products_only:
            substrates = list(self.network.graph.successors(reaction_node))
            subOne, subTwo = self._determine_subOne_subTwo(substrates)

        for enz in possible_enzymes:
            score, info, = self.similarity_searcher.scoreReaction(reaction_name, enz, product, subOne, subTwo,
                                                                  sim_cutoff=self.config.similarity_threshold,
                                                                  onlyActive=self.config.only_active,
                                                                  only_reviewed=self.config.only_reviewed)

            self.network.graph.nodes[reaction_node]['attributes']['metadata']['specificity_scores'][enz] = score
            self.network.graph.nodes[reaction_node]['attributes']['metadata']['enzyme_info'][enz] = info


    def _select_best_enzyme(self, reaction_node):
        current_enz = self.network.graph.nodes[reaction_node]['attributes']['metadata']['selected_enzyme']
        current_score = self.network.graph.nodes[reaction_node]['attributes']['metadata']['specificity_scores'][current_enz]
        current_score_neg = True
        possible_enzymes = self.network.graph.nodes[reaction_node]['attributes']['metadata']['possible_enzymes']

        for enz in possible_enzymes:
            score = self.network.graph.nodes[reaction_node]['attributes']['metadata']['specificity_scores'][enz]
            if (score > current_score and score != 0) or (abs(score) > current_score and current_score_neg == True):
                self.network.graph.nodes[reaction_node]['attributes']['metadata']['selected_enzyme'] = enz
                current_score = abs(score)
                if score < 0:
                    current_score_neg = True

    def _determine_subOne_subTwo(self, listSubstrateSmiles):
        subOne = None
        subTwo = None
        for smi in listSubstrateSmiles:
            if self.network.graph.nodes[smi]['attributes']['substrate_num'] == 1:
                subOne = smi
            elif self.network.graph.nodes[smi]['attributes']['substrate_num'] == 2:
                subTwo = smi
        return subOne, subTwo