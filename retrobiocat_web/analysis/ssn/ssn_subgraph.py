from retrobiocat_web.analysis.ssn.ssn_main import SSN


def keep_top_n_edges(graph, top_n=50):
    for node in graph.nodes:
        node_edges = graph.edges(node, data=True)
        sorted_node_edges = sorted(node_edges, key=lambda x: x[2]['weight'], reverse=True)
        edges_to_keep = sorted_node_edges[0:top_n]
        edges_to_delete = sorted_node_edges[top_n:]


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser

    enzyme_type = 'AAD'
    ssn = SSN(enzyme_type)
    ssn.load()

    keep_top_n_edges(ssn.graph)
