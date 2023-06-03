from matplotlib import pyplot as plt
from scipy.cluster import hierarchy

import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_clusterer import Route_Clusterer
from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_distancer import Route_Distancer
from retrobiocat_web.retro.pathway_search.route_evaluation.costing.pathway_cost_modal import Pathway_Coster
from retrobiocat_web.retro.pathway_search.mcts.config import MCTS_Config
from retrobiocat_web.retro.pathway_search.mcts.tree_node import MCTS_Node
from retrobiocat_web.retro.pathway_search.mcts.mcts_scoring import MCTSScorer
from retrobiocat_web.retro.pathway_search.mcts.expander import NetworkExpander
import time
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine


class MCTS():
    def __init__(self, network: Network, config=None, retro_config=None, log_level='WARNING'):

        self.stop_expansions = False
        self.network = network

        self.config = config
        if self.config is None:
            self.config = MCTS_Config()

        self.retro_engine = RetrosynthesisEngine(self.network, config=retro_config)

        self.retro_engine.initialise_graph()
        self.target_smiles = self.network.target_smiles

        self.root = MCTS_Node(self, None, [network.target_smiles], [network.target_smiles], [], terminal=False, log_level=log_level)
        self.tree_nodes = [self.root]
        self.terminal_nodes = []
        self.solved_nodes = []
        self.pathways = []
        self.expander = NetworkExpander(self.network, self.retro_engine, config=self.config)
        self.pathway_scorer = MCTSScorer(self.network, config=self.config)

        self.logger = add_logger('MCTS', level=log_level)

    def selection(self, currentTreeNode):
        """ Move through search tree selecting nodes with highest ucb until finding a node with 0 children (ie a leaf)"""

        self.logger.debug("SELECTION")

        if (currentTreeNode is self.root) and self.root.terminal is True:
            # if no expansions possible from the root, return None to end the search
            return None

        while len(currentTreeNode.children) != 0:  # loop through search tree until node has no children
            currentTreeNode = currentTreeNode.select(self.config.exploration)

            if currentTreeNode is None:
                return None

        return currentTreeNode

    def rollout(self, currentTreeNode):
        """ Expand tree until terminal node """

        self.logger.debug("ROLLOUT")

        currentTreeNode.expand()
        while currentTreeNode.terminal is False:
            currentTreeNode = currentTreeNode.select(self.config.exploration)
            if currentTreeNode is None:
                return None

            currentTreeNode.expand()

        return currentTreeNode

    def backpropagate(self, currentTreeNode, score):
        self.logger.debug("BACKPROPAGATE")

        currentTreeNode.value += score
        currentTreeNode.visits += 1
        while currentTreeNode.parent is not None:
            currentTreeNode = currentTreeNode.parent
            currentTreeNode.value += score
            currentTreeNode.visits += 1

        if currentTreeNode.depth != 0:
            print("warning - backpropogation didn't reach root")

    def get_pathways(self):
        """ Get all the pathways in the search tree """
        self.pathways = []

        if self.config.only_solved:
            for node in self.solved_nodes:
                pathway = node.get_pathway()
                self.pathways.append(pathway)
        else:
            for node in self.terminal_nodes:
                pathway = node.get_pathway()
                self.pathways.append(pathway)

        # if no pathways at all, return the root.
        if len(self.pathways) == 0:
            self.pathways = [self.root.get_pathway()]

        return self.pathways

    def cluster_and_cost(self, pathways):

        single_step_pathways = []
        if self.config.dont_cluster_single_step_routes:
            single_step_pathways, pathways = self._get_single_step_pathways(pathways)

        distances = Route_Distancer().get_distances(pathways)
        clusters = Route_Clusterer().group_by_cluster(pathways, distances,
                                                      self.config.distance_threshold)
        clusters += [[p] for p in single_step_pathways]

        pathway_coster = Pathway_Coster(instock_cost=self.config.in_stock_cost,
                                        not_instock_cost=self.config.not_in_stock_cost,
                                        reaction_cost=self.config.reaction_cost,
                                        reaction_yield=self.config.reaction_yield)

        clusters = pathway_coster.order_clusters_by_cost(clusters)
        return clusters

    def _get_single_step_pathways(self, pathways):
        multi_step_pathways, single_step_pathways = [], []
        for pathway in pathways:
            if pathway.pathway_length != 1:
                multi_step_pathways.append(pathway)
            else:
                single_step_pathways.append(pathway)
        return single_step_pathways, multi_step_pathways

    def run(self):
        t0 = time.time()
        elapsed = 0
        iterations = 0
        while (elapsed < self.config.max_search_time) and (len(self.terminal_nodes) < self.config.max_pathways):
            selected_tree_node = self.selection(self.root)  # select a tree node with no children, based on ucb scores
            if selected_tree_node is None:
                self.logger.info("Search exhausted")
                break

            rolled_out_final_node = self.rollout(selected_tree_node)  # rollout is just expansion and selection until terminal.
            if rolled_out_final_node is None:
                self.logger.info("Search exhausted")
                break

            score = self.pathway_scorer.score(rolled_out_final_node)  # score the rolled out pathway
            self.backpropagate(rolled_out_final_node, score)  # backpropagate the score
            elapsed = int(round(time.time() - t0,1))
            iterations += 1

        self.logger.info(f"MCTS COMPLETE - {iterations} iterations in {elapsed} seconds")
        self.logger.info(f"{len(self.terminal_nodes)} terminal nodes, {len(self.solved_nodes)} solved nodes")

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    target_smiles = '[C@H]1(C2=CC=CC=C2)NCCCC1'

    network = Network(target_smiles=target_smiles)
    mcts = MCTS(network, log_level='DEBUG')
    mcts.config.max_search_time = 60
    mcts.config.maxLength = 5
    mcts.run()

    print(f"{len(mcts.tree_nodes)} total nodes")
    print(f"{len(mcts.pathways)} terminal nodes")
    print(f"{len(mcts.solved_nodes)} solved nodes")

    visits = []
    for node in mcts.tree_nodes:
        visits.append(node.visits)
    average_visits = sum(visits)/len(visits)
    most_visits = max(visits)
    print(f"Average visits {average_visits}")
    print(f"Max visits {most_visits}")

    pathways = mcts.get_pathways()

    distances = Route_Distancer().get_distances(pathways)
    clusters = Route_Clusterer().group_by_cluster(pathways, distances, 1)
    clusters = Pathway_Coster().order_clusters_by_cost(clusters)

