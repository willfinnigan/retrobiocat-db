import json
import os

import networkx as nx
import pandas as pd


def save_to_csv(graph, folder):
    df_graph = nx.to_pandas_edgelist(graph)

    att_dict = {}
    for node in list(graph):
        att_dict[node] = graph.nodes[node]

    df_graph.to_csv(f"{folder}/graph.csv", index=False)

    with open(f'{folder}/attributes.json', 'wb') as outfile:
        outfile.write(json.dumps(att_dict).encode("utf-8"))

def load_graph_from_csv(folder):
    if not os.path.exists(f"{folder}/graph.csv") or not os.path.exists(f"{folder}/attributes.json"):
        print(f"No saved SSN found in {folder}, could not load")
        return False

    df_graph = pd.read_csv(f"{folder}/graph.csv")
    graph = nx.from_pandas_edgelist(df_graph, edge_attr=['weight', 'i'])
    return graph

def load_attrs_from_csv(folder, graph):
    att_dict = json.load(open(f'{folder}/attributes.json'))

    # Nodes with no edges are not in edge list..
    for node in list(att_dict.keys()):
        if node not in list(graph.nodes):
            add_protein_node(graph, node, alignments_made=False)

    nx.set_node_attributes(graph, att_dict)
    return graph

def add_protein_node(graph, node_name, alignments_made=False):
    """ If a protein is not already in the graph, then add it """

    if 'UniRef50' in node_name:
        node_type = 'uniref'
    else:
        node_type = 'biocatdb'

    if node_name not in graph.nodes:
        graph.add_node(node_name, node_type=node_type,
                            alignments_made=alignments_made)
        return 1

    if alignments_made == True:
        graph.nodes[node_name]['alignments_made'] = True

    return 0