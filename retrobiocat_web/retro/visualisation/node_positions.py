

def fix_target_position(list_nodes, graph):
    for node_dict, i in enumerate(list_nodes):
        node = node_dict['id']
        if 'target' in graph.nodes[node]['attributes']:
            list_nodes[i]['fixed'] = True
    return list_nodes

