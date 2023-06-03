import networkx as nx
from retrobiocat_web.analysis.ehreact.seed_filtering import generate_non_match_seed, \
    are_any_seed_substructures_of_each_other, get_non_matches
from retrobiocat_web.analysis.ehreact.train.calculate_diagram import calculate_diagram_single_reactant
from retrobiocat_web.analysis.ehreact.train.hasse import one_node_hasse
from retrobiocat_web.mongo.model_queries.reaction_queries import reaction_from_name


def get_seed_list(reaction_names, smi_col):
    seed_list = []
    if reaction_names is None:
        return []

    for name in reaction_names:
        rxn_obj = reaction_from_name(name)
        if rxn_obj is not None:
            if smi_col == 'product_1_smiles' and rxn_obj.product_seeds is not None:
                seed_list += rxn_obj.product_seeds
            elif smi_col == 'substrate_1_smiles' and rxn_obj.substrate_1_seeds is not None:
                seed_list += rxn_obj.substrate_1_seeds
            elif smi_col == 'substrate_2_smiles' and rxn_obj.substrate_2_seeds is not None:
                seed_list += rxn_obj.substrate_2_seeds

    return seed_list

def ensure_correct_seed_coverage(seeds, smiles, should_generate_non_match_seed=True):
    if should_generate_non_match_seed == True or len(seeds) == 0:
        non_match_seeds = generate_non_match_seed(smiles, seeds)
        seeds += non_match_seeds
    if len(seeds) > 1:
        seeds = are_any_seed_substructures_of_each_other(seeds)
    return seeds

def remove_non_matches(seeds, smiles):
    non_matches, matches, ori_smiles = get_non_matches(smiles, seeds)
    ori_non_matches = [ori_smiles[s] for s in non_matches]
    return matches, ori_non_matches, ori_smiles

def hasse_diagrame(smiles, seed_list=None):
    if seed_list is None:
        seed_list = []

    if len(smiles) == 0:
        return None

    if len(smiles) == 1:
        return one_node_hasse(smiles[0])

    d, smi_dict = calculate_diagram_single_reactant(smiles, seed_list, False, True)
    return d

def get_children_edges(node):
    """ Recursive function which will return all children and edges from a node """
    children = []
    edges = []
    for edge in node.edges_to_child:
        children.append(str(edge.child_node.key))
        edges.append([str(edge.parent_node.key), str(edge.child_node.key)])

        recursive_children, recursive_edges = get_children_edges(edge.child_node)
        children += recursive_children
        edges += recursive_edges
    return children, edges

def get_real_root(root):
    return root.edges_to_child[0].child_node

def to_networkx(root):
    children, edges = get_children_edges(root)

    g = nx.DiGraph()
    if len(children) == 0:
        g.add_node(root.key, level=0)
        return g

    else:
        for node in children:
            g.add_node(node, level=0)

        for edge in edges:
            g.add_edge(edge[0], edge[1])
        return g

def hasse_diagram_to_network(hasse_diagrame):
    if hasse_diagrame is None:
        graph = nx.DiGraph()
        graph.add_node('')
        return graph, ''

    root = get_real_root(hasse_diagrame.root)
    graph = to_networkx(root)
    start_node = root.key
    return graph, start_node



if __name__ == '__main__':
    from retrobiocat_web.analysis.ehreact.get_substrate_summary import leafs_under_node
    from retrobiocat_web.analysis.data_query.data_query_from_args import data_query_from_args
    from retrobiocat_web.analysis.ehreact.seed_filtering import filter_smiles_with_smarts
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    reaction_name = 'Amination (TPL)'
    enzyme_type = 'TPL'
    dq = data_query_from_args({'reaction': reaction_name,
                               'enzyme_type': enzyme_type,
                               'mol_type': 'substrate_2_smiles'})

    print(dq.unique_smiles())

    seed_list = get_seed_list(dq.reaction_names(), dq.smi_col)
    print(seed_list)
    smis = filter_smiles_with_smarts(dq.unique_smiles(), seed_list)
    seed_list = ensure_correct_seed_coverage(seed_list, smis, should_generate_non_match_seed=True)
    print(seed_list)

    smis, non_matches, ori_smiles = remove_non_matches(seed_list, smis)
    print(f"{len(non_matches)} non matches, {len(smis)} matches")

    d = hasse_diagrame(smis, seed_list=seed_list)
    print(d)
    graph, start_node = hasse_diagram_to_network(d)
    print(f"graph nodes = {graph.nodes}")

    print(f"num graph nodes = {len(graph.nodes)}")
    print(f"num graph edges = {len(graph.edges)}")
    print(f"start node = {start_node}")


    # list_ori_smis = [ori_smiles[s] for s in smis]
    # for smi in list_ori_smis:
    #     if smi not in graph.nodes():
    #         print(f"{smi} not in graph")

    all_leafs = leafs_under_node(start_node, graph)
    print(f"{len(smis)} smiles, {len(all_leafs)} leafs")

    print(all_leafs)
    # all_leafs2 = only_leafs(list_ori_smis, graph)
    # print(f"{len(smis)} smiles, {len(all_leafs2)} leafs")

    #for smi in list_ori_smis:
        #if smi not in all_leafs:
            #print(smi)
            #print(leafs_under_node(smi, graph))

