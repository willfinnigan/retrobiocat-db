from pathlib import Path

import networkx as nx
import pandas as pd

from retrobiocat_web.analysis.ssn.ssn_save_load.ssn_save_load_sqlite import save_to_sqlite, load_edgelist_from_sqlite


def test_can_save_a_graph_and_load_the_data():
    filepath = 'test_sqlite.sqlite'
    Path(filepath).unlink(missing_ok=True)

    graph = nx.Graph()
    graph.add_edge('node_1', 'node_2', weight=100, i=0.39)
    graph.add_edge('node_2', 'node_3', weight=110, i=0.5)

    save_to_sqlite(graph, filepath)
    edgelist = load_edgelist_from_sqlite(filepath)

    assert dict(edgelist.iloc[0]) == {'source': 'node_1', 'target': 'node_2', 'weight': 100, 'i': 0.39}
    assert dict(edgelist.iloc[1]) == {'source': 'node_2', 'target': 'node_3', 'weight': 110, 'i': 0.5}



def test_can_save_a_graph_and_load_the_graph_back():
    filepath = 'test_sqlite.sqlite'
    Path(filepath).unlink(missing_ok=True)

    graph = nx.Graph()
    graph.add_edge('node_1', 'node_2', weight=100, i=0.39)
    graph.add_edge('node_2', 'node_3', weight=110, i=0.5)

    save_to_sqlite(graph, filepath)
    edgelist = load_edgelist_from_sqlite(filepath)
    new_graph = nx.from_pandas_edgelist(edgelist, edge_attr=['weight', 'i'])
    assert list(new_graph.edges) == [('node_1', 'node_2'), ('node_2', 'node_3')]


def test_can_alter_sqlite_file():
    filepath = 'test_sqlite.sqlite'
    Path(filepath).unlink(missing_ok=True)

    graph = nx.Graph()
    graph.add_edge('node_1', 'node_2', weight=100, i=0.39)
    graph.add_edge('node_2', 'node_3', weight=110, i=0.5)

    save_to_sqlite(graph, filepath)
    edgelist = load_edgelist_from_sqlite(filepath)
    new_graph = nx.from_pandas_edgelist(edgelist, edge_attr=['weight', 'i'])

    new_graph.add_edge('node_3', 'node_4', weight=120, i=0.6)
    save_to_sqlite(new_graph, filepath)

    edgelist2 = load_edgelist_from_sqlite(filepath)
    new_graph2 = nx.from_pandas_edgelist(edgelist2, edge_attr=['weight', 'i'])

    assert new_graph2.edges == new_graph.edges
