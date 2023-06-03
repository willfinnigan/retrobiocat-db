import statistics


from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser


class SSN_Cluster_Precalculator(object):

    def __init__(self, ssn, log_level=1):
        self.ssn = ssn

        self.start = self.ssn.min_score
        self.end = 250
        self.step = 5
        self.min_cluster_size = 6
        self.max_cluster_edges = 1500000

        self.log_level = log_level

    def precalulate(self, num=5, current_num_clusters=0):
        precalculated_nodes = {}
        num_at_alignment_score = {}
        identity_at_alignment_score = {}
        count = 0
        for alignment_score in range(self.start, self.end, self.step):
            visualiser = SSN_Visualiser(self.ssn.enzyme_type, log_level=0)

            clusters, graph = visualiser.get_clusters_and_subgraph(self.ssn, alignment_score)
            num_clusters = self._get_num_clusters(clusters)
            max_edges = self._get_max_edges(graph, clusters)
            if max_edges > self.max_cluster_edges:
                self.log(f"Edges ({max_edges}) larger than max_edges ({self.max_cluster_edges}), skipping alignment score {alignment_score}")
                continue

            if (num_clusters > current_num_clusters) and (num_clusters != 0):
                count += 1
                self.log(f"Adding new set of {num_clusters} clusters at alignment score {alignment_score}")

                current_num_clusters = num_clusters

                pos_dict = visualiser.get_cluster_positions(graph, clusters)
                nodes, edges = visualiser.get_nodes_and_edges(graph, pos_dict)
                ident = self._get_ident(graph)

                num_at_alignment_score[str(alignment_score)] = num_clusters
                precalculated_nodes[str(alignment_score)] = nodes
                identity_at_alignment_score[str(alignment_score)] = ident

            if count == num:
                return precalculated_nodes, num_at_alignment_score, identity_at_alignment_score

        return precalculated_nodes, num_at_alignment_score, identity_at_alignment_score

    def _get_ident(self, graph):
        identities = []
        for edge in list(graph.edges(data=True)):
            identities.append(edge[2]['i'])

        if len(identities) >= 2:
            ident = [round(statistics.mean(identities), 2),
                     round(statistics.stdev(identities, statistics.mean(identities)), 2)]
        else:
            ident = [0, 0]

        return ident

    def _get_num_clusters(self, clusters):
        num = 0
        for cluster in clusters:
            if len(cluster) >= self.min_cluster_size:
                num += 1
        return num


    def _get_max_edges(self, graph, clusters) -> int:
        """ Get the maximum number of edges in any one cluster"""
        max_edge = 0
        for cluster in clusters:
            sub_graph = graph.subgraph(cluster)
            edges = sub_graph.number_of_edges()
            if edges > max_edge:
                max_edge = edges

        return max_edge

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"Cluster_Precalc({self.ssn.enzyme_type}): {msg}")


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    from retrobiocat_web.analysis.ssn.ssn_main import SSN

    enzyme_type = 'IRED'
    ssn = SSN(enzyme_type)
    ssn.load_sqlite()

    visualiser = SSN_Visualiser(enzyme_type, log_level=1)
    clusters, graph = visualiser.get_clusters_and_subgraph(ssn, 40)

    total_edges = 0
    total_nodes = 0
    for cluster in clusters:
        sub_graph = graph.subgraph(cluster)
        edges = sub_graph.number_of_edges()
        nodes = len(sub_graph.nodes)
        print(f"Number of edges = {sub_graph.number_of_edges()}")
        print(f"Number of nodes = {len(sub_graph.nodes)}")
        total_edges += edges
        total_nodes += nodes

    print(f"Total nodes = {total_nodes}")
    print(f"Total edges = {total_edges}")
