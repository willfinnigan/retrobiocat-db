import retrobiocat_web.retro.network_pathway.graph_functions
import networkx as nx
from retrobiocat_web.logging import add_logger
from retrobiocat_web.retro.network_pathway.node_functions import sort_by_score


class GraphPruner():

    def __init__(self, log_level='WARNING'):
        self.logger = add_logger('GraphPruner', level=log_level)

    def prune(self, network_obj, smiles, max_nodes=False):
        def check_if_nodes_need_to_be_removed(network):
            if (max_nodes != False) and (len(network.substrate_nodes) > max_nodes):
                return True
            return False

        def get_smiles_not_removed(removed_smis, smi):
            """Return list of smiles which are not in list removed"""
            not_removed_smis = []
            for smi in smi:
                if smi not in removed_smis:
                    not_removed_smis.append(smi)
            return not_removed_smis

        if check_if_nodes_need_to_be_removed(network_obj):
            removed_nodes = self._run_prune(network_obj, max_nodes, steps=1)
            not_removed = get_smiles_not_removed(removed_nodes, smiles)
            return not_removed

        else:
            return smiles

    def delete_reaction_node(self, network, node_to_remove):
        if node_to_remove not in network.reaction_nodes:
            self.logger.warning('Can not delete non-reaction node')
            return []

        deleted = [node_to_remove]
        self._remove_node(network, node_to_remove)

        # if any node is now not connected to the target, delete also
        connected_nodes = list(nx.dfs_preorder_nodes(network.graph, source=network.target_smiles))

        while len(connected_nodes) != len(list(network.graph.nodes)):
            for node in list(network.graph.nodes):
                if node not in connected_nodes:
                    deleted.append(node)
                    self._remove_node(network, node)
            connected_nodes = list(nx.dfs_preorder_nodes(network.graph, source=network.target_smiles))

        return deleted

    def delete_terminal_reaction_node(self, network, node_to_remove):
        if node_to_remove not in network.reaction_nodes:
            self.logger.warning('Can not delete non-reaction node')
            return []

        substrate_successors = list(network.graph.successors(node_to_remove))

        for substrate in substrate_successors:
            if len(list(network.graph.successors(substrate))) != 0:
                return []

        self._remove_node(network, node_to_remove)
        deleted = [node_to_remove]

        for node in substrate_successors:
            if len(list(network.graph.predecessors(node))) == 0:
                self._remove_node(network, node)
                deleted.append(node)

        return deleted

    def _remove_node(self, network, node):

        network.graph.remove_node(node)

        if node in network.substrate_nodes:
            network.substrate_nodes.remove(node)
        elif node in network.reaction_nodes:
            network.reaction_nodes.remove(node)

    def _get_nodes_to_prune(self, graph, num_to_remove, on_substrates=False):
        def prune_on_reactions(end_nodes, num):
            end_node_reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes_of_list_substrates(graph, end_nodes)
            end_node_reactions = self._check_other_substrates_of_end_node_reactions(graph, end_node_reactions)
            scores = []
            for reaction in end_node_reactions:
                scores.append(graph.nodes[reaction]['attributes']['change_in_complexity'])

            sorted_reactions = sort_by_score(end_node_reactions, scores)
            return sorted_reactions[0:num]

        def prune_on_substrates(end_nodes, num):
            scores = []
            for substrate in end_nodes:
                scores.append(graph.nodes[substrate]['attributes']['complexity'])
            sorted_substrates = sort_by_score(end_nodes, scores, reverse=True)
            reactions_to_prune = []
            while len(reactions_to_prune) < num:
                reactions = retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes_of_list_substrates(graph, [sorted_substrates.pop(0)])
                reactions_to_prune.extend(reactions)
                reactions_to_prune = self._check_other_substrates_of_end_node_reactions(graph, reactions_to_prune)
            return reactions_to_prune

        all_end_nodes = retrobiocat_web.retro.network_pathway.graph_functions.get_nodes_with_no_successors(graph)
        if on_substrates:
            return prune_on_substrates(all_end_nodes, num_to_remove)
        else:
            return prune_on_reactions(all_end_nodes, num_to_remove)

    @staticmethod
    def _check_other_substrates_of_end_node_reactions(graph, end_node_reactions):
        def are_all_substrates_terminal(reaction_to_check):
            substrate_successors = list(graph.successors(reaction_to_check))
            for substrate in substrate_successors:
                if len(list(graph.successors(substrate))) != 0:
                    return False
            return True

        reactions_to_keep = []
        for reaction in end_node_reactions:
            if are_all_substrates_terminal(reaction):
                reactions_to_keep.append(reaction)

        return reactions_to_keep

    def _run_prune(self, network, max_nodes, steps=5, on_substrates=True):
        self.logger.info('Prune network')
        nodes_removed = []
        while len(network.substrate_nodes) >= max_nodes:
            to_remove = self._get_nodes_to_prune(network.graph, steps, on_substrates=on_substrates)

            for node in to_remove:
                if node in list(network.graph.nodes()):
                    deleted_nodes = self.delete_terminal_reaction_node(network, node)
                    nodes_removed.extend(deleted_nodes)
                    if len(network.substrate_nodes) >= max_nodes:
                        break
                else:
                    self.logger.warning(f'Node not present - could not delete - {node}')

        return nodes_removed
