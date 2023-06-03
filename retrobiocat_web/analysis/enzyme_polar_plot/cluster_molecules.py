from sklearn import cluster
from scipy.cluster.vq import vq
from rdkit import DataStructs
from sklearn.metrics.pairwise import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import AgglomerativeClustering
import numpy as np
from retrobiocat_web.analysis.data_query import get_data

def organise_groups(unique_smis, fps, labels, kmeans_clusterer, show_log=False):
    label_dict = {}
    for smi, label in zip(unique_smis, labels):
        if label not in label_dict:
            label_dict[label] = []
        label_dict[label].append(smi)

    # Get representative molecules
    centroids = kmeans_clusterer.cluster_centers_
    closest, distances = vq(centroids, fps)

    label_reps = {}
    rep_labels = {}
    sorted_labels = sorted([int(l) for l in label_dict.keys()])
    for l, smi_index in zip(sorted_labels, closest):
        label_reps[str(l)] = unique_smis[smi_index]
        rep_labels[unique_smis[smi_index]] = str(l)

    if show_log == True:
        for label in label_reps:
            print(label)
            print(f"Rep smi: {label_reps[label]}")
            print(f"Group smis: {label_dict[label]}")

    return label_dict, label_reps, rep_labels

def get_clusters(unique_smis, fps, n_clusters, random_seed=0):
    kmeans = cluster.KMeans(n_clusters=n_clusters, random_state=random_seed)
    kmeans.fit(fps)
    labels = [str(l) for l in kmeans.labels_]
    label_dict, label_reps, rep_labels = organise_groups(unique_smis, fps, labels, kmeans)

    return labels, label_dict, label_reps, rep_labels

def check_reps_appear_in_groups(rep_labels, label_dict):
    for smi, label in rep_labels.items():
        group_smis = label_dict[label]
        if smi not in group_smis:
            print(f'Warning - label {label} rep smi {smi} not in group')
            print(group_smis)
            print()



def order_smis(rep_labels):
    rep_smis = list(rep_labels.keys())
    size = len(rep_smis)
    hmap = np.empty(shape=(size, size))
    rep_fps = get_data.get_substrates_fps(rep_smis)
    for i, fp_i in enumerate(rep_fps):
        similarities = DataStructs.BulkTanimotoSimilarity(fp_i, list(rep_fps))
        for j, sim in enumerate(similarities):
            hmap[i, j] = sim

    linked = linkage(hmap, 'single')
    results = dendrogram(linked, no_plot=True, distance_sort='ascending')

    colours = results['color_list']
    sub_numbered = list(map(int, results['ivl']))

    substrates_ordered = []
    for i in sub_numbered:
        substrates_ordered.append(rep_smis[i])

    return substrates_ordered



if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    from retrobiocat_web.analysis.data_query import get_data

    dq = get_data.DataQuery(enzyme_type='CAR', reaction="Carboxylic acid reduction", smi_col='substrate_1_smiles', log_level=1)
    unique_smis = dq.unique_smiles()
    fps = get_data.get_substrates_fps(unique_smis)
    fp_dict = get_data.make_fp_dict(unique_smis, fps)

    labels, label_dict, label_reps, rep_labels = get_clusters(unique_smis, fps, 12, random_seed=0)

    check_reps_appear_in_groups(rep_labels, label_dict)

    substrates_ordered = order_smis(rep_labels)
    print(substrates_ordered)

