import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway
import numpy as np
import copy

class Unevaluated_MCTS_Node():

    def __init__(self, mcts, parent, option, heuristic_choice_score, log_level="WARNING"):
        self.mcts = mcts
        self.option = option
        self.heuristic_choice_score = heuristic_choice_score
        self.parent = parent

        self.log_level = log_level
        self.logger = add_logger(f'Unevaluated MCTS_Node', level=log_level)

    def evaluate(self):
        # execute option creating normal MCTS nodes

        reactions = self.option.get_processed_reactions()
        new_nodes = []
        for reaction in reactions:

            # do not allow choices where any the reaction products already exist in the pathway
            if self._reaction_products_not_already_in_pathway(reaction):

                node = self._create_MCTS_Node_from_reaction(reaction)
                new_nodes.append(node)

        return new_nodes

    def _create_MCTS_Node_from_reaction(self, reaction):
        product_smis, reaction_names = self.mcts.retro_engine.graphAdder.add_reactions([reaction])
        new_end_nodes = self._get_new_end_nodes(reaction.reactant.smi, product_smis)
        new_pathway_nodes = self._get_new_pathway_nodes(product_smis, reaction_names)
        new_reaction_nodes = self._get_new_reaction_nodes(reaction_names)
        node = MCTS_Node(self.mcts, self.parent, new_pathway_nodes, new_end_nodes, new_reaction_nodes,
                         heuristic_choice_score=self.heuristic_choice_score,
                         log_level=self.log_level)
        self.mcts.tree_nodes.append(node)
        return node

    def _reaction_products_not_already_in_pathway(self, reaction):
        if any(elem in self.parent.pathway_nodes for elem in [s.smi for s in reaction.substrates] + [reaction.name]):
            return False
        return True

    def _get_new_end_nodes(self, reactant_smi, new_smis):
        new_end_nodes = copy.copy(self.parent.end_nodes)
        new_end_nodes.remove(reactant_smi)
        new_end_nodes.extend(new_smis)
        return new_end_nodes

    def _get_new_pathway_nodes(self, new_smis, new_reaction_names):
        new_pathway_nodes = copy.copy(self.parent.pathway_nodes)
        new_pathway_nodes.extend(new_reaction_names)
        new_pathway_nodes.extend(new_smis)
        return new_pathway_nodes

    def _get_new_reaction_nodes(self, new_reaction_names):
        new_reaction_nodes = copy.copy(self.parent.reaction_nodes)
        new_reaction_nodes.extend(new_reaction_names)
        return new_reaction_nodes

    def calculate_UCB(self, exploration):
        # as soon as this node is selected it will be evaluated.
        # nodes which have never been selected have a ucb of inf
        return np.inf

class MCTS_Node():

    def __init__(self, mcts, parent, all_pathway_nodes, end_pathway_nodes, reaction_nodes,
                 heuristic_choice_score=0, terminal=False,
                 log_level='WARNING'):

        self.mcts = mcts
        self.parent = parent

        if parent is None:
            self.depth = 0
        else:
            self.depth = self.parent.depth + 1

        self.children = []

        self.pathway_nodes = all_pathway_nodes
        self.end_nodes = end_pathway_nodes
        self.reaction_nodes = reaction_nodes

        # for ucb calculation
        self.visits = 0
        self.value = 0

        self.heuristic_choice_score = heuristic_choice_score

        self.terminal = terminal
        self.solved = False
        self.fully_searched = False

        self.pathway = None

        self.log_level = log_level
        self.logger = add_logger(f'MCTS_Node(d={self.depth})', level=log_level)

    def expand(self):
        """ Expands the search tree by adding child nodes from this one """

        self.logger.debug("----- EXPAND -----")

        # if already expanded, or is terminal, then dont expand
        if len(self.children) > 0 or self.terminal:
            self.logger.debug(f"Not expanding because either Num children>0 ({len(self.children)}), or Terminal={self.terminal}")
            return

        # 1. get buyable and non-buyable end nodes.  only expand non-buyable (will revisit that later)
        buyable, non_buyable = self._get_buyable_end_nodes()
        self.logger.debug(f"{len(self.end_nodes)} end nodes in total, {len(non_buyable)} are non-buyable")

        # 2. if no non-buyable end nodes, then this node is terminal and solved
        if len(non_buyable) == 0:
            self._mark_node_terminal()
            self._mark_node_solved()
            self.logger.debug(f"All end nodes are buyable. Node is terminal and solved")
            return

        # 2. look at which of the non_buyables are at max length
        not_at_max, at_max_length = self._non_buyable_end_nodes_not_at_and_at_max_length(non_buyable)

        # option to prevent expansion if there is a non-buyable node at the max length
        # or do we give up if we reach this scenario
        if len(at_max_length) > 0 and self.mcts.mcts_config.stop_expansion_if_nonbuyable_at_max_length:
            self.logger.debug(f"Non-buyable end node has reach max length.  Stoping expansion of tree nodes here. Node is terminal")
            self._mark_node_terminal()
            return

        # 3. expand the network and get options
        # options should be ordered by some scoring system.  (simplest AiZynth policy score)
        # one option should be to stop. - if config is to stop at buyable, only a stop will be returned
        options, scores = self.mcts.expander.options_from_end_nodes(not_at_max, self.pathway_nodes)
        self.logger.debug(f"Generated {len(options)} options")

        # 4. create child nodes from the options
        for option, score in zip(options, scores):
            child = Unevaluated_MCTS_Node(self.mcts, self, option, score, log_level=self.log_level)
            self.children.append(child)

        if len(self.children) == 0:
            self.logger.debug(f"No children (either because max depth or no appliciable templates.  Node is terminal")
            self._mark_node_terminal()

    def select(self, exploration):

        self.logger.debug("----- SELECT -----")

        child = None
        while child is None:
            self.logger.debug("Child is none - finding most promising")
            self.logger.debug(f"{len([c for c in self.children if type(c) is Unevaluated_MCTS_Node])} unevaluated children out of {len(self.children)}")

            child = self.find_most_promising_child(self.children, exploration)
            if type(child) is Unevaluated_MCTS_Node:
                self.logger.debug("Child is not evaluated - evaluating")
                child = self._get_evaluated_child(child)

            if child is None and self.depth == 0 and self.fully_searched:  # ie is root and has no children left to explore, the search is complete
                self.logger.debug(f"Root is fully searched. Search is complete")
                return None

        self.logger.debug(f"Child pathway = {child.pathway_nodes}")
        self.logger.debug(f"Child end nodes = {child.end_nodes}")

        return child

    def _get_evaluated_child(self, unevaluated_child):
        evaluated_children = unevaluated_child.evaluate()
        self.children.remove(unevaluated_child)
        self.children.extend(evaluated_children)

        if len(evaluated_children) == 0:
            return None
        return evaluated_children[0]

    def calculate_UCB(self, exploration):
        if self.terminal or self.fully_searched:
            return -np.inf

        N = self.parent.visits
        ni = self.visits
        C = exploration
        wi = self.value

        if (ni == 0) or (N == 0):
            return np.inf

        ucb = (wi/ni) + (C*(np.sqrt(np.log(N)/ni)))

        return ucb

    def find_best_heuristic_child(self, children_to_select_from):
        scores = [c.heuristic_choice_score for c in children_to_select_from]
        child_scores = list(zip(children_to_select_from, scores))
        sorted_child_scores = sorted(child_scores, key=lambda x: x[1], reverse=True)
        best_child = sorted_child_scores[0][0]
        return best_child

    def find_most_promising_child(self, children, exploration):
        """Get the child with the highest ucb.  Fallback on heuristic for ties """

        # if no children, this node is terminal, then go back to parent
        if len(children) == 0:
            self.logger.debug('Node has no children left, returning Parent')
            self._mark_node_terminal()
            return self.parent

        nodes_with_max_score = []
        max_score = -np.inf
        children_ucbs = []
        for child in children:
            ucb = child.calculate_UCB(exploration)
            children_ucbs.append(ucb)
            if ucb == max_score:
                nodes_with_max_score.append(child)
            if ucb > max_score:
                nodes_with_max_score = [child]
                max_score = ucb

        self.logger.debug(f"UCB scores for children = {children_ucbs}")

        if max_score == -np.inf:  # if all ucb scores are -np.inf, do not return a child, go back to the parent
            self.fully_searched = True
            self.logger.debug(f"All ucbs are -np.inf, so are all either terminal or have been fully explored...")
            self.logger.debug(f"..this node is therefore also fully explored. Marking as 'fully_searched' and going back to parent.")
            return self.parent

        if len(nodes_with_max_score) == 1:  # if a single best score, use this
            return nodes_with_max_score[0]

        return self.find_best_heuristic_child(nodes_with_max_score)  # if multiple best ucb's, then use a heuristic

    def get_pathway(self):
        """ Get the pathway for a particular tree search node """
        if self.pathway is None:
            self.pathway = Pathway(self.pathway_nodes, self.mcts.network)
        return self.pathway

    def _get_buyable_end_nodes(self):
        buyable = []
        not_buyable = []
        for node in self.end_nodes:
            if self.mcts.network.graph.nodes[node]['attributes']['is_starting_material'] == 1:
                buyable.append(node)
            else:
                not_buyable.append(node)
        return buyable, not_buyable

    def _mark_node_terminal(self):
        self.terminal = True
        self.mcts.terminal_nodes.append(self)

    def _mark_node_solved(self):
        self.solved = True
        self.mcts.solved_nodes.append(self)

    def _non_buyable_end_nodes_not_at_and_at_max_length(self, non_buyables):
        not_at_max_length = []
        at_max_length = []

        for node in non_buyables:
            subgraph = self.mcts.network.graph.subgraph(self.pathway_nodes)
            num = retrobiocat_web.retro.network_pathway.graph_functions.get_num_reactions_from_target(subgraph,
                                                                                                      node,
                                                                                                      self.mcts.network.target_smiles)

            if num >= self.mcts.mcts_config.maxLength:
                at_max_length.append(node)
            else:
                not_at_max_length.append(node)


        return not_at_max_length, at_max_length





