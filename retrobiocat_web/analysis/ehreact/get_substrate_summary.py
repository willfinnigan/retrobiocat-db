import networkx as nx

from retrobiocat_web.analysis.ehreact.r_groups import get_core


def children_of_node(node, graph):
    return list(graph.successors(node))

def only_leafs(list_of_nodes, graph):
    leafs = []
    for node in list_of_nodes:
        if not node_has_children(node, graph):
            leafs.append(node)
    return leafs

def node_has_children(node, graph):
    if len(list(graph.successors(node))) != 0:
        return True
    return False

def not_leaf_children(list_of_nodes, graph):
    not_leafs = []
    for node in list_of_nodes:
        if node_has_children(node, graph):
            not_leafs.append(node)
    return not_leafs

def nodes_under_node(node, graph):
    successors = nx.nodes(nx.dfs_tree(graph, node))
    return successors

def leafs_under_node(node, graph):
    successors = nodes_under_node(node, graph)  # all the nodes under this one
    leaf_successors = []
    for s_node in successors:
        if not node_has_children(s_node, graph):
            leaf_successors.append(s_node)
    return leaf_successors


def recursively_move_single_child_nodes_along(leafs, not_leafs, graph, max_moves=25):
    single_child_not_leafs = [n for n in not_leafs if len(children_of_node(n, graph)) == 1]
    count = 0
    while len(single_child_not_leafs) != 0:
        if count > max_moves:
            break
        leafs, not_leafs = move_nodes_with_one_child_on_one(leafs, not_leafs, graph)
        single_child_not_leafs = [n for n in not_leafs if len(children_of_node(n, graph)) == 1]
        count += 1
    return leafs, not_leafs

def move_nodes_with_one_child_on_one(leafs, not_leafs, graph):
    new_not_leafs = []
    for node in not_leafs:
        children = children_of_node(node, graph)
        if len(children) == 1:
            only_child = children[0]
            if node_has_children(only_child, graph):
                new_not_leafs.append(only_child)
            else:
                leafs.append(only_child)
        else:
            new_not_leafs.append(node)
    return leafs, new_not_leafs

def remove_not_leafs_with_less_than_min_examples(leafs, not_leafs, graph, min_examples):
    new_not_leafs = []
    for node in not_leafs:
        under_leafs = leafs_under_node(node, graph)
        if len(under_leafs) <= min_examples:
            leafs += under_leafs
        else:
            new_not_leafs.append(node)
    return leafs, new_not_leafs

def get_expandable_not_leafs(graph, not_leaf_children_dict, min_examples):
    expandables = []
    for node, children in not_leaf_children_dict.items():
        leafs, not_leafs = get_leafs_and_not_leafs(node, graph, min_examples)
        if len(not_leafs) != 0:
            expandables.append(node)
    return expandables

def get_leafs_and_not_leafs(start_node, graph, min_examples):

    children = children_of_node(start_node, graph)
    leafs = only_leafs(children, graph)
    not_leafs = not_leaf_children(children, graph)
    leafs, not_leafs = recursively_move_single_child_nodes_along(leafs, not_leafs, graph)
    leafs, not_leafs = remove_not_leafs_with_less_than_min_examples(leafs, not_leafs, graph, min_examples)

    return leafs, not_leafs

def convert_to_ori_smiles(list_smis, ori_smiles):
    return [ori_smiles[s] for s in list_smis]

def summarise(start_node, graph, ori_smiles, min_examples=2):
    if len(graph.nodes) == 0:
        return [], {}, {}, {}, [], ''
    elif len(graph.nodes) == 1:
        node = list(graph.nodes)[0]
        return [node], {}, {}, {}, [], node

    leafs, not_leafs = get_leafs_and_not_leafs(start_node, graph, min_examples)


    count = 0
    while len(leafs) == 0 and len(not_leafs) == 1:
        if count > 5:
            break
        start_node = not_leafs[0]
        leafs, not_leafs = get_leafs_and_not_leafs(start_node, graph, min_examples)
        count += 1

    not_leaf_children_dict = {n: leafs_under_node(n, graph) for n in not_leafs}
    not_leaf_counts = {n: len(v) for n, v in not_leaf_children_dict.items()}
    not_leaf_counts = dict(sorted(not_leaf_counts.items(), key=lambda item: item[1], reverse=True))
    not_leaf_cores = {n: get_core(not_leaf_children_dict[n], n) for n in not_leaf_counts.keys()}
    expandables = get_expandable_not_leafs(graph, not_leaf_children_dict, min_examples)
    parent_core = get_core(leafs, start_node)

    ori_leafs = convert_to_ori_smiles(leafs, ori_smiles)
    ori_not_leaf_children_dict = {n: convert_to_ori_smiles(v, ori_smiles) for n, v in not_leaf_children_dict.items()}

    if '' in ori_leafs:
        ori_leafs.remove('')

    print(ori_leafs)

    return ori_leafs, ori_not_leaf_children_dict, not_leaf_counts, not_leaf_cores, expandables, parent_core


if __name__ == '__main__':
    from retrobiocat_web.analysis.data_query import get_data
    from retrobiocat_web.analysis.ehreact.create_hasse_network import hasse_diagrame, hasse_diagram_to_network, \
        filter_smiles_with_smarts
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    dq = get_data.DataQuery(reaction='Carboxylic acid reduction')
    smis = filter_smiles_with_smarts(dq.unique_smiles(), 'C=O')
    d = hasse_diagrame(smis)
    graph, start_node = hasse_diagram_to_network(d)

    print(len(smis))

    all_leafs, leafs, not_leaf_counts, not_leaf_cores = summarise(start_node, graph)
