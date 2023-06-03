import time
import json
from retrobiocat_web.analysis.data_query import get_data
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering

from rdkit import DataStructs
import numpy as np

import copy


def find_max_list(list):
    list_len = [len(i) for i in list]
    return max(list_len)

class EnzymeSplitter(object):
    enzyme_col = 'enzyme_name'

    def __init__(self, data_query, max_enzymes=100, log_level=0):

        self.log_level = log_level

        self.max_enzymes = max_enzymes

        self.activity_data = data_query.activity_data
        self.unique_enzymes = get_data.get_unique_enzymes(self.activity_data, self.enzyme_col)

    def find_max_list(self, list):
        list_len = [len(i) for i in list]
        return max(list_len)

    def split(self):

        n_clusters = 2
        t0 = time.time()

        if len(self.unique_enzymes) > self.max_enzymes:
            self.log("Splitting into enzymes with and without sequences")
            embeddings = get_data.get_enzyme_embeddings(self.unique_enzymes)
            seq_enzymes, non_seq_enzymes = self.split_enzymes_into_seq_none_seq(self.unique_enzymes, embeddings)

            while find_max_list(seq_enzymes) > self.max_enzymes:
                self.log(f"Splitting enzymes with sequences into {n_clusters} clusters by agglomerative clustering")
                seq_enzymes = self.split_by_clustering(embeddings, n_clusters)
                n_clusters += 1

            while find_max_list(non_seq_enzymes) > self.max_enzymes:
                self.log(f"Splitting enzymes with-out sequences - now {len(non_seq_enzymes)+1} groups")
                non_seq_enzymes = self.split_nonseq_enzymes(non_seq_enzymes)

            sets_of_unique_enzymes = seq_enzymes + non_seq_enzymes
            lengths_of_lists = [len(enz_list) for enz_list in sets_of_unique_enzymes]
            t1 = time.time()
            self.log(f"Split complete in {round(t1-t0, 2)} seconds - {len(sets_of_unique_enzymes)} groups")
            self.log(f"Lengths are {str(lengths_of_lists)}")

        else:
            self.log(f"No split needed, {len(self.unique_enzymes)} is less than the max {self.max_enzymes} enzymes")
            sets_of_unique_enzymes = [self.unique_enzymes]

        return sets_of_unique_enzymes

    def split_enzymes_into_seq_none_seq(self, enzymes, embeddings_map):
        enz_with_embedding = list(embeddings_map.keys())
        other_enzymes = [enz for enz in enzymes if enz not in enz_with_embedding]
        return [enz_with_embedding], [other_enzymes]

    def split_nonseq_enzymes(self, non_seq_enzymes):
        # find largest list
        large_list_len = 0
        largest_index = None
        for i, list_enz in enumerate(non_seq_enzymes):
            if len(list_enz) > large_list_len:
                large_list_len = len(list_enz)
                largest_index = i

        large_list = non_seq_enzymes.pop(largest_index)
        half_index = int(len(large_list ) / 2)
        non_seq_enzymes.append(large_list[:half_index])
        non_seq_enzymes.append(large_list[half_index:])

        return non_seq_enzymes

    def split_by_clustering(self, embeddings_map, n_clusters):

        enz_with_embeddings = list(embeddings_map.keys())
        embeddings = list(embeddings_map.values())

        ag_clus = AgglomerativeClustering(n_clusters=n_clusters)
        ag_clus.fit(embeddings)
        cluster_labels = ag_clus.labels_

        clustered_enz = {}
        for i, enz in enumerate(enz_with_embeddings):
            if str(cluster_labels[i]) not in clustered_enz:
                clustered_enz[str(cluster_labels[i])] = []
            clustered_enz[str(cluster_labels[i])].append(enz)

        clustered_enz_list = []
        for enz_list in clustered_enz.values():
            clustered_enz_list.append(enz_list)

        return clustered_enz_list

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"EnzymeSplitter ({level}): {msg}")

class MoleculeSplitter(object):

    def __init__(self, data_query, max_molecules=150, smi_col='product_1_smiles', log_level=0):

        self.log_level = log_level

        self.max_molecules = max_molecules
        self.smi_col = smi_col

        self.activity_data = data_query.activity_data
        self.unique_smiles = get_data.get_unique_smiles(self.activity_data, self.smi_col)
        self.fps = get_data.get_substrates_fps(self.unique_smiles)


    def split(self):

        t0 = time.time()

        n_clusters = 2
        molecule_groups = [self.unique_smiles]

        if find_max_list(molecule_groups) > self.max_molecules:
            distance_matrix = self.get_distance_matrix()

            while find_max_list(molecule_groups) > self.max_molecules:
                self.log(f"Splitting molecules into {n_clusters} clusters by agglomerative clustering")
                molecule_groups = self.split_mols_by_clustering(n_clusters, distance_matrix)
                n_clusters += 1

            lengths_of_lists = [len(enz_list) for enz_list in molecule_groups]
            t1 = time.time()
            self.log(f"Split complete in {round(t1-t0, 2)} seconds - {len(molecule_groups)} groups")
            self.log(f"Lengths are {str(lengths_of_lists)}")

        else:
            self.log(f"No split needed, {len(self.unique_smiles)} is less than the max {self.max_molecules} molecules")

        return molecule_groups

    def get_distance_matrix(self):
        size = len(self.unique_smiles)
        hmap = np.empty(shape=(size, size))
        for i, fp_i in enumerate(self.fps):
            similarities = DataStructs.BulkTanimotoSimilarity(fp_i, list(self.fps))
            for j, sim in enumerate(similarities):
                hmap[i, j] = sim

        return hmap

    def split_mols_by_clustering(self, n_clusters, distance_matrix):

        ag_clus = AgglomerativeClustering(n_clusters=n_clusters, linkage='complete', affinity='precomputed')
        ag_clus.fit(distance_matrix)
        cluster_labels = ag_clus.labels_

        clustered_smi = {}
        for i, smi in enumerate(self.unique_smiles):
            if str(cluster_labels[i]) not in clustered_smi:
                clustered_smi[str(cluster_labels[i])] = []
            clustered_smi[str(cluster_labels[i])].append(smi)

        clustered_smi_list = []
        for smi_list in clustered_smi.values():
            clustered_smi_list.append(smi_list)

        return clustered_smi_list

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"MoleculeSplitter ({level}): {msg}")



def get_mol_dqs(dqs, max_molecules, log_level=0):
    new_dqs = []
    for dq in dqs:
        mol_groups = MoleculeSplitter(dq, max_molecules=max_molecules, log_level=log_level).split()

        if len(mol_groups) > 1:
            for list_mols in mol_groups:
                sub_dq = copy.copy(dq)
                sub_dq.smiles = list_mols
                sub_dq.get_data()
                new_dqs.append(sub_dq)

        else:
            new_dqs.append(dq)

    return new_dqs

def get_enz_dqs(dqs, max_enzymes, log_level=0):
    new_dqs = []
    for dq in dqs:
        enz_groups = EnzymeSplitter(dq, max_enzymes=max_enzymes, log_level=log_level).split()

        if len(enz_groups) > 1:
            for list_enz in enz_groups:
                sub_dq = copy.copy(dq)
                sub_dq.enzyme_names = list_enz
                sub_dq.get_data()
                new_dqs.append(sub_dq)

        else:
            new_dqs.append(dq)

    return new_dqs

def split(dq, max_molecules, max_enzymes, log_level=0):
    dqs = [dq]
    dqs = get_mol_dqs(dqs, max_molecules, log_level=log_level)
    dqs = get_enz_dqs(dqs, max_enzymes, log_level=log_level)
    return dqs


if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')
    #make_default_connection(host='192.168.0.88', database='test')

    dq = get_data.DataQuery(enzyme_type='IRED', log_level=1, remove_negative='all')
    dqs = split(dq, 300, 300, log_level=1)

    print(len(dqs))
