import json
import os
import time
from pathlib import Path

import pandas as pd
import networkx as nx
import sqlite3
import sqlalchemy



def save_to_sqlite(ssn_graph, filepath):
    Path(filepath).unlink(missing_ok=True)
    conn_string = f"sqlite:///{filepath}"
    engine = sqlalchemy.create_engine(conn_string)
    engine.connect()

    df_graph = nx.to_pandas_edgelist(ssn_graph)
    df_graph.to_sql('edgelist', con=engine, if_exists='replace')

    att_dict = {}
    for node in list(ssn_graph):
        att_dict[node] = ssn_graph.nodes[node]

    df_attr = pd.DataFrame.from_dict(att_dict, orient='index')
    df_attr['node'] = df_attr.index
    df_attr.reset_index(drop=True, inplace=True)

    df_attr.to_sql('attrs', con=engine, if_exists='replace')


def load_edgelist_from_sqlite(filepath, alignment_score=None, node_list=None):

    cmd = f"SELECT source, target, weight, i FROM edgelist"

    if alignment_score is not None:
        cmd += f" WHERE weight >= {alignment_score}"

    if node_list is None:
        pass
    else:
        if alignment_score is None:
            cmd += ' WHERE '
        else:
            cmd += ' AND '
        node_list_str = str(node_list).replace('[', '(').replace(']', ')')
        cmd += f"(source IN {node_list_str} OR target in {node_list_str})"

    with sqlite3.connect(filepath) as conn:
        cursor = conn.cursor()
        cursor.execute(cmd)
        result = cursor.fetchall()

    df = pd.DataFrame(result, columns=['source', 'target', 'weight', 'i'])

    return df

def load_attrs_from_sqlite(filepath, node_list=None):
    cmd = f"SELECT node, node_type, alignments_made FROM 'attrs'"
    if node_list is None:
        pass
    else:
        node_list_str = str(node_list).replace('[', '(').replace(']', ')')
        cmd += f" WHERE node IN {node_list_str}"

    with sqlite3.connect(filepath) as conn:
        cursor = conn.cursor()
        cursor.execute(cmd)
        result = cursor.fetchall()

    attr_dict = {}
    for data in result:
        attr_dict[data [0]] = {"node_type": data[1], "alignments_made": bool(data[2])}
    return attr_dict


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    from retrobiocat_web.analysis.ssn.ssn_main import SSN

    enzyme_type = 'TA'
    ssn = SSN(enzyme_type)
    ssn.load()
    #save_to_sqlite(ssn.graph, 'test_ssn.sqlite')
    t0 = time.time()
    filepath = f"{ssn.save_path}/{ssn.enzyme_type}_sqlite_ssn.sqlite"

    #result = load_edgelist_from_sqlite('test_ssn.sqlite', node_list=['VfTA_L56A-I259T'])
    t1 = time.time()
    print(result)
    print(f"Time = {round(t1-t0, 4)} seconds")
