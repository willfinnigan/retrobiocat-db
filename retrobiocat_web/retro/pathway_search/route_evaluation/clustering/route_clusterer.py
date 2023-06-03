import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

class Route_Clusterer():

    def cluster_by_distance(self, distances, distance_threshold):
        return self._cluster(distances, None, distance_threshold=distance_threshold)

    def group_by_cluster(self, pathways, distances, distance_threshold):
        if len(pathways) == 0:
            return []
        elif len(pathways) == 1:
            return [[pathways[0]]]
        elif len(pathways) == 2:
            return [[pathways[0]], [pathways[1]]]

        labels = self.cluster_by_distance(distances, distance_threshold)

        groups = {}
        for pathway, label in zip(pathways, labels):
            if label not in groups:
                groups[label] = []
            groups[label].append(pathway)

        list_groups = []
        for label, group in groups.items():
            list_groups.append(group)

        return list_groups

    def linkage_matrix(self, distances, linkage='single', affinity='precomputed', n_clusters=None, distance_threshold=0.1):
        """Produces a linkage matrix for the creation of a dendogram """

        model = AgglomerativeClustering(linkage=linkage, affinity=affinity, n_clusters=n_clusters, distance_threshold=distance_threshold)
        model.fit(distances)
        counts = np.zeros(len(model.distances_))
        matrix = np.column_stack([model.children_, model.distances_, counts])
        return matrix.astype(float)

    def _optimized_cluster(self, distances, max_clusters):
        """
        Optimize the number of cluster based  Silhouette metric.
        Additional arguments to the clustering algorithm
        can be passed in as key-word arguments
        :param max_clusters: the maximum number of clusters to consider
        :return: the cluster index for each observation
        """

        max_score = None
        best_size = None
        optimization_scores = []
        for n_clusters in range(2, min(max_clusters + 1, len(distances))):
            clusters = self._cluster(distances, n_clusters)
            score = silhouette_score(distances, clusters, metric="precomputed")
            optimization_scores.append(score)
            if best_size is None or score > max_score:
                max_score = score
                best_size = n_clusters

        if best_size is None:
            best_size = max_clusters
        return self._cluster(distances, best_size)

    def _cluster(self, distances, n_clusters, linkage='single', affinity="precomputed", distance_threshold=None):
        model = AgglomerativeClustering(linkage=linkage, affinity=affinity, n_clusters=n_clusters, distance_threshold=distance_threshold)
        model.fit(distances)
        return model.labels_



