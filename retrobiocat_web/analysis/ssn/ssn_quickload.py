import time
from pathlib import Path
import os
import pandas as pd
from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser
from retrobiocat_web.analysis import analysis_paths

class SSN_quickload(object):

    def __init__(self, enzyme_type, log_level=0):
        self.enzyme_type = enzyme_type
        self.ssn_folder = analysis_paths.ssn_folder
        self.save_path = f'{self.ssn_folder}/{self.enzyme_type}'
        self.ssn = None
        self.log_level = log_level
        self.df = None
        self.vis = SSN_Visualiser(enzyme_type, hidden_edges=False, log_level=log_level)

    def load_df(self):
        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        #self.df = dask.dataframe.read_csv(f"{self.save_path}/graph.csv")
        #self.df = self.df.compute()
        self.df = pd.read_csv(f"{self.save_path}/graph.csv")

        t1 = time.time()

        self.log(f"Loaded df in {round(t1-t0, 1)} seconds")

    def get_edges(self, selected_node, alignment_score):
        t0 = time.time()
        edges = []
        df = self.df[(self.df['source'] == selected_node) | (self.df['target'] == selected_node)]
        df = df[df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            edges.append(self.vis.get_vis_edge(row['target'], row['source'], row['weight']))
        t1 = time.time()

        self.log(f"Retrieved edges for {self.enzyme_type} node at alignment score {alignment_score} in {round(t1-t0, 1)} seconds")

        return edges

    def get_all_edges(self, alignment_score):
        edges = []
        df = self.df[self.df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            edges.append(self.vis.get_vis_edge(row['target'], row['source'], row['weight']))

        return edges


    def get_multiple_edges(self, nodes, alignment_score):
        t0 = time.time()
        edges = []
        df = self.df[(self.df['source'].isin(nodes)) | (self.df['target'].isin(nodes))]
        df = df[df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            edges.append(self.vis.get_vis_edge(row['target'], row['source'], row['weight']))
        t1 = time.time()

        self.log(
            f"Retrieved multiple edges for {len(nodes)} nodes of {self.enzyme_type} node at alignment score {alignment_score} in {round(t1 - t0, 1)} seconds")

        return edges

    def get_connected_nodes(self, list_nodes, alignment_score):
        t0 = time.time()
        connected_nodes = set()

        df = self.df[((self.df['source'].isin(list_nodes)) | (self.df['target'].isin(list_nodes)))]
        df = df[df['weight'] >= alignment_score]

        for index, row in df.iterrows():
            if row['target'] not in list_nodes:
                connected_nodes.add(row['target'])
            if row['source'] not in list_nodes:
                connected_nodes.add(row['source'])

        connected_nodes = list(connected_nodes)
        t1 = time.time()
        self.log(f"Retrieved {len(connected_nodes)} connected nodes for {self.enzyme_type} node at alignment score {alignment_score} in {round(t1 - t0, 1)} seconds")

        return connected_nodes

    def get_all_connected_nodes(self, node, alignment_score, log=True):
        t0 = time.time()

        connected_nodes = [node]
        new_nodes = [node]
        while len(new_nodes) != 0:
            df = self.df[((self.df['source'].isin(new_nodes)) | (self.df['target'].isin(new_nodes)))]
            df = df[df['weight'] >= alignment_score]

            new_nodes = []
            for index, row in df.iterrows():
                if row['target'] not in connected_nodes:
                    connected_nodes.append(row['target'])
                    new_nodes.append(row['target'])
                if row['source'] not in connected_nodes:
                    connected_nodes.append(row['source'])
                    new_nodes.append(row['source'])

        t1 = time.time()
        if log:
            self.log(f"Retrieved all connected nodes to {node} for {self.enzyme_type} node at alignment score {alignment_score} in {round(t1 - t0, 1)} seconds")

        return connected_nodes



    def log(self, msg, level=1):
        if self.log_level >= level:
            print("SSN_quickload: " + msg)


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection

    make_default_connection()

    ql = SSN_quickload('IRED', log_level=1)
    ql.load_df()
    edges = ql.get_multiple_edges(['UniRef50_Q2TW47'], 45)


    """
    print(f"Num edges default = {len(ssn.graph.edges)}")


    """

    # ssn.delete_edges_below_alignment_score(40, min_ident=0.35)
    # print(f"Num edges = {len(ssn.graph.edges)}")

    # ssn.delete_edges_below_alignment_score(45)
    # print(f"Num edges = {len(ssn.graph.edges)}")

    """
    filtered_graph = ssn.get_graph_filtered_edges(40)
    print(f"Num edges alignment 40 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(45)
    print(f"Num edges alignment 45 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(50)
    print(f"Num edges alignment 50 = {len(filtered_graph.edges)}")

    filtered_graph = ssn.get_graph_filtered_edges(55)
    print(f"Num edges alignment 55 = {len(filtered_graph.edges)}")

    """

    # for i in range(20,100,10):
    #   new_graph = ssn.get_graph_filtered_edges(i)
    #  num_edges = len(new_graph.edges)
    # print(f"Score = {i}, num edges = {num_edges}")

    # 1. Set minimum number of nodes for a cluster
    # 2. Move down alignment score.  For each score where there is a different number of scores, visualise this.
