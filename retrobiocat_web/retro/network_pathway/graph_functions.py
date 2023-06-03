import networkx as nx

def get_nodes_with_no_successors(graph):
    no_succ = []
    for node in list(graph.nodes()):
        if len(list(graph.successors(node))) == 0:
            no_succ.append(node)
    return no_succ


def get_reactions_with_smi_as_reactant(graph, smi):
    return list(graph.successors(smi))


def get_reaction_nodes_of_list_substrates(graph, list_substrates):
    reaction_nodes = []
    for substrate in list_substrates:
        predecessors = list(graph.predecessors(substrate))
        reaction_nodes.extend(predecessors)
    return reaction_nodes


def get_terminal_reaction_nodes(graph):
    """ Return a list of reaction nodes which lead to terminal nodes"""
    no_successors = get_nodes_with_no_successors(graph)
    no_successor_rxns = get_reaction_nodes_of_list_substrates(graph, no_successors)
    return no_successor_rxns


def get_reaction_substrates(graph, reaction):
    """Return a list of the substrates for a given reaction"""
    return list(graph.successors(reaction))


def get_reaction_products(graph, reaction):
    """Return a list of the products for a given reaction"""
    return list(graph.predecessors(reaction))


def get_substrate_nodes(graph):
    substrate_nodes = []
    for node in list(graph):
        if 'node_type' not in graph.nodes[node]['attributes']:
            print('WARNING node_type not available for node: ' + str(node))
        if graph.nodes[node]['attributes']['node_type'] == 'substrate':
            substrate_nodes.append(node)
    return substrate_nodes


def get_reaction_nodes(graph):
    reaction_nodes = []
    for node in list(graph):
        if graph.nodes[node]['attributes']['node_type'] == 'reaction':
            reaction_nodes.append(node)
    return reaction_nodes


def check_substrates_nx_predecessor(graph, smiles_to_check, origin_smile):
    any_are_predecessor = False

    # Get all the substrate nodes which are predecessors of origin_smile
    substrate_predecessors = []
    reaction_predecessors = list(graph.predecessors(origin_smile))
    for reaction_node in reaction_predecessors:
        substrate_predecessors.extend(list(graph.predecessors(reaction_node)))

    # do any of these match the smiles we want to check?
    for smiles in smiles_to_check:
        if smiles in substrate_predecessors:
            any_are_predecessor = True

    return any_are_predecessor


def substrate_precursors(graph, smile):
    substrate_precursors = []
    reaction_precursors = list(graph.successors(smile))

    for reaction in reaction_precursors:
        substrate_precursors.extend(list(graph.successors(reaction)))

    return substrate_precursors


def get_num_reactions_from_target(graph, smile, target):
    """ Returns the number of reaction nodes in shortest path from target to source"""
    num_in_path = nx.shortest_path_length(graph, source=target, target=smile)
    num_reacions_in_path = int(num_in_path/2)
    return num_reacions_in_path


def get_pathway_length(target_smi, subgraph):
    """Returns the length of the pathway as the max distance an end node is from target"""

    def get_depth(smi, depths, d):
        d += 1
        children = substrate_precursors(subgraph, smi)
        depths += [d]
        for c_smi in children:
            new_d, new_depths = get_depth(c_smi, depths, d)
            depths += new_depths
        return d, depths

    d, depths = get_depth(target_smi, [], 0)
    max_depth = max(depths) - 1
    return max_depth


def get_current_pathway_length(graph, endNodes, target):
    currentLength = 0
    for endNode in endNodes:
        length = get_num_reactions_from_target(graph, endNode, target)
        if length > currentLength:
            currentLength = length
    return currentLength


def get_reaction_names(reactions_list, graph):
    """ Returns a list of reaction names (minus uuid) """
    reaction_names = []
    for reaction in reactions_list:
        name = graph.nodes[reaction]['attributes']['name']
        reaction_names.append(name)
    return reaction_names