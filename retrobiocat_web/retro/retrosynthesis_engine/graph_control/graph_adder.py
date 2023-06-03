from retrobiocat_web.logging import add_logger

class GraphAdder():

    def __init__(self, network, log_level='WARNING'):
        self.network = network
        self.logger = add_logger('GraphControl', level=log_level)

    def add_reactions(self, reactions):
        reaction_nodes, substrate_nodes = [], []
        for reaction in reactions:
            new_reac_node, new_substrate_nodes = self._add_single_reaction(reaction)
            reaction_nodes.append(new_reac_node)
            substrate_nodes.extend(new_substrate_nodes)

        return substrate_nodes, reaction_nodes

    def add_starting_mol(self, target_mol):
        self.network.graph.add_node(target_mol.smi, attributes={'name': target_mol.smi,
                                                                'complexity': target_mol.complexity,
                                                                'relative_complexity': 0,
                                                                'target': 1,
                                                                'is_starting_material': target_mol.in_stock,
                                                                'starting_material_info': target_mol.in_stock_info,
                                                                'node_type': 'substrate',
                                                                'node_num': self._get_node_number(target_mol.smi, self.network.graph),
                                                                'substrate_num': 0})

        self._add_to_substrate_nodes(target_mol.smi)

    def _add_single_reaction(self, reaction):
        if reaction.graph_node_name is not None:
            self.logger.debug(f'Reaction {reaction.graph_node_name} already added')
            return reaction.graph_node_name, [s.smi for s in reaction.substrates]

        # otherwise, add reaction to graph
        reaction_node = self._add_reaction_node_to_graph(reaction)
        for i, substrate in enumerate(reaction.substrates):
            self._add_substrate_node_to_graph(substrate, reaction_node, num=i)
        return reaction_node, [s.smi for s in reaction.substrates]

    def _add_reaction_node_to_graph(self, reaction):
        unique_reaction_name = f"{reaction.name} ({reaction.unique_id})"

        self.network.graph.add_node(unique_reaction_name, attributes={'name': reaction.name,
                                                                      'node_type': 'reaction',
                                                                      'reaction_type': reaction.rxn_type,
                                                                      'change_in_complexity': reaction.complexity_change,
                                                                      'metadata': reaction.metadata,
                                                                      'node_num': self._get_node_number(
                                                                          unique_reaction_name, self.network.graph)})

        self.network.graph.add_edge(reaction.reactant.smi, unique_reaction_name)
        self._add_to_reaction_nodes(unique_reaction_name)
        reaction.graph_node_name = unique_reaction_name
        return unique_reaction_name

    def _add_substrate_node_to_graph(self, substrate, source_reaction_name, num=0):

        self.network.graph.add_node(substrate.smi, attributes={'name': substrate.smi,
                                                               'complexity': substrate.complexity,
                                                               'relative_complexity': substrate.relative_complexity,
                                                               'is_starting_material': substrate.in_stock,
                                                               'starting_material_info': substrate.in_stock_info,
                                                               'node_type': 'substrate',
                                                               'node_num': self._get_node_number(substrate, self.network.graph),
                                                               'substrate_num': num+1})

        self.network.graph.add_edge(source_reaction_name, substrate.smi)
        self._add_to_substrate_nodes(substrate.smi)

    def _is_molecule_in_graph(self, smi):
        if smi in self.network.graph.nodes:
            return True
        return False

    def _get_node_number(self, node, graph):
        if node in list(graph.nodes()):
            return graph.nodes[node]['attributes']['node_num']
        else:
            return len(list(graph.nodes())) + 1

    def _add_to_substrate_nodes(self, node):
        if node not in self.network.substrate_nodes:
            self.network.substrate_nodes.append(node)

    def _add_to_reaction_nodes(self, node):
        if node not in self.network.reaction_nodes:
            self.network.reaction_nodes.append(node)


