from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, UniRef50
import networkx as nx
from networkx.drawing.nx_agraph import graphviz_layout
import json
import palettable
from retrobiocat_web.analysis.ssn.ssn_positions import ClusterPositioner
import mongoengine as db
from retrobiocat_web.analysis.ssn import node_density

class SSN_Visualiser(object):

    def __init__(self, enzyme_type, hidden_edges=True, log_level=0):
        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
        self.node_metadata = self._find_uniref_metadata()

        self.edge_colour = {'color': 'lightgrey'}
        self.edge_width = 2
        self.hidden_edges = hidden_edges
        self.border_width = 0.8
        self.border_colour = 'black'
        self.border_width_selected = 5
        self.opacity = 0.8
        self.luminosity = 'bright'
        self.node_colour = f'rgba(5, 5, 168, {self.opacity})'
        self.highlight_background = "rgb(255, 161, 161)"
        self.node_size = 30
        self.node_shape = 'dot'

        self.swiss_prot_colour = 'gold'

        self.log_level = log_level
        self.cluster_positioner = ClusterPositioner()

    def scale_func(self, pos_dict, scale):
        new_pos_dict = {}
        for node, positions in pos_dict.items():
            x = positions[0] * scale
            y = positions[1] * scale
            new_pos_dict[node] = (x, y)
        return new_pos_dict

    def scale_cluster_alt(self, pos_dict, target_median_distance):
        if len(pos_dict) <= 1:
            return pos_dict

        median_distance = node_density.get_median_distance(pos_dict)
        while median_distance > target_median_distance:
            pos_dict = self.scale_func(pos_dict, 0.75)
            median_distance = node_density.get_median_distance(pos_dict)

        while median_distance < target_median_distance:
            pos_dict = self.scale_func(pos_dict, 1.5)
            median_distance = node_density.get_median_distance(pos_dict)

        return pos_dict

    def scale_cluster_alt2(self, pos_dict, target_density, box_distance, graph):
        if len(pos_dict) <= 1:
            return pos_dict

        influencial_node = nx.voterank(graph, number_of_nodes=1)[0]

        density = node_density.calculate_node_density(pos_dict, influencial_node, box_distance=box_distance)
        while density < target_density:
            pos_dict = self.scale_func(pos_dict, 0.75)
            density = node_density.calculate_node_density(pos_dict, influencial_node, box_distance=box_distance)

        while density > target_density:
            pos_dict = self.scale_func(pos_dict, 1.5)
            density = node_density.calculate_node_density(pos_dict, influencial_node, box_distance=box_distance)

        return pos_dict

    def scale_cluster(self, pos_dict, target_density):
        if len(pos_dict) <= 1:
            return pos_dict

        def node_density(pos_dict):
            # get bounding box
            x_values = []
            y_values = []
            for node, positions in pos_dict.items():
                x_values.append(positions[0])
                y_values.append(positions[1])

            box = (max(x_values)-min(x_values), max(y_values)-min(y_values))

            box_size = box[0]*box[1]
            if box_size == 0:
                box_size = 1
            density = len(pos_dict) / box_size
            return density

        density = node_density(pos_dict)
        while density < target_density:
            pos_dict = self.scale_func(pos_dict, 0.75)
            density = node_density(pos_dict)

        while density > target_density:
            pos_dict = self.scale_func(pos_dict, 1.5)
            density = node_density(pos_dict)

        return pos_dict

    def visualise(self, ssn, alignment_score):
        clusters, graph = self.get_clusters_and_subgraph(ssn, alignment_score)
        pos_dict = self.get_cluster_positions(graph, clusters)

        nodes, edges = self.get_nodes_and_edges(graph, pos_dict)

        return nodes, edges

    def get_cluster_positions(self, graph, clusters):

        pos_dict = {}
        for i, cluster in enumerate(clusters):
            self.log(f"Getting layout for cluster {i + 1} of {len(clusters)}")
            sub_graph = graph.subgraph(cluster)
            #scale = 750 + (20 * len(cluster))
            #scale = None

            #if len(cluster) > 200:
            #cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="sfdp")
            #elif len(cluster) > 10:
            #    cluster_positions = nx.nx_pydot.pydot_layout(sub_graph, prog="neato")
            #else:


            if len(cluster) <= 20:
                cluster_positions = nx.spring_layout(sub_graph)
            else:
                cluster_positions = graphviz_layout(sub_graph, prog="sfdp")  #sfdp

            cluster_positions = self.scale_cluster(cluster_positions, 1.0E-04)

            cluster_positions = self.cluster_positioner.position(cluster_positions)
            cluster_positions = self.cluster_positioner.round_positions(cluster_positions, round_to=0)
            pos_dict.update(cluster_positions)

        return pos_dict

    def get_clusters_and_subgraph(self, ssn, alignment_score):
        graph = ssn.get_graph_filtered_edges(alignment_score)
        clusters = list(nx.connected_components(graph))
        clusters.sort(key=len, reverse=True)
        graph = self._add_cluster_node_colours(graph, clusters)

        return clusters, graph

    def _add_cluster_node_colours(self, graph, clusters):
        self.log(f"Colouring clusters.. (opacity={self.opacity})")
        colours = palettable.cartocolors.qualitative.Vivid_10.mpl_colors
        for i, cluster in enumerate(clusters):
            c = colours.pop(0)
            colours.append(c)
            cc = f"rgba({int(c[0] * 255)},{int(c[1] * 255)},{int(c[2] * 255)},{self.opacity})"
            for j, node in enumerate(cluster):
                graph.nodes[node]['colour'] = cc
                graph.nodes[node]['cluster_label'] = i+1

        return graph

    def get_nodes_and_edges(self, graph, pos_dict, full_graph=None):
        nodes = []
        edges = []
        for name in graph.nodes:
            colour = graph.nodes[name].get('colour', None)
            cluster_label = graph.nodes[name].get('cluster_label', None)
            nodes.append(self._get_vis_node(name, pos_dict=pos_dict, colour=colour, cluster_label=cluster_label))

        if full_graph is None:
            for edge in graph.edges:
                weight = graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self.get_vis_edge(edge[0], edge[1], weight))
        else:
            for edge in full_graph.edges:
                weight = full_graph.get_edge_data(edge[0], edge[1], default={'weight': 0})['weight']
                edges.append(self.get_vis_edge(edge[0], edge[1], weight))

        nodes = self.add_swissprot_border(nodes)
        nodes = self.add_rhea_or_pdb(nodes)
        nodes = self.sort_nodes(nodes)

        return nodes, edges

    def _get_vis_node(self, node_name, pos_dict=None, colour=None, cluster_label=None):
        if colour is None:
            colour = self.node_colour
        if cluster_label is None:
            cluster_label = 'ND'

        if 'UniRef50' in node_name:
            node_type = 'uniref'
        else:
            node_type = 'biocatdb'
            colour = 'darkred'

        metadata = self.node_metadata.get(node_name, {})
        protein_name = metadata.get('protein_name', '')
        tax = metadata.get('tax', '')
        if protein_name != '':
            title = f"{protein_name} - {tax}"
        else:
            title = node_name

        node = {'id': node_name,
                'size': self.node_size,
                'borderWidth': self.border_width,
                'borderWidthSelected': self.border_width_selected,
                'color': {'background': colour,
                          'border': self.border_colour,
                          'highlight': {'border': colour,
                                        'background': self.highlight_background}},
                'title': title,
                'cluster_label': cluster_label,
                'label': '',
                'shape': self.node_shape,
                'node_type': node_type,
                'swissprot': 'no',
                'rhea': metadata.get('rhea', []),
                'protein_name': protein_name,
                'pfam_codes': metadata.get('pfam_codes', [])}

        if pos_dict is not None:
            x, y = tuple(pos_dict.get(node_name, (0, 0)))
            node['x'] = x
            node['y'] = y

        return node

    def get_vis_edge(self, edge_one, edge_two, weight):
        # weight = self.graph.get_edge_data(edge_one, edge_two, default={'weight': 0})['weight']
        edge = {'id': f"from {edge_one} to {edge_two}",
                'from': edge_one,
                'to': edge_two,
                'hidden': self.hidden_edges,
                'weight': weight,
                'width': self.edge_width,
                'color': self.edge_colour}
        return edge

    def sort_nodes(self, vis_nodes):
        """ Returns vis_nodes with any nodes marked as node_type='biocatdb' at the front """

        biocatdb_nodes = []
        pdb_nodes = []
        rhea_nodes = []
        swissprot_nodes = []
        other_nodes = []

        for node in vis_nodes:
            if 'biocatdb' in node.get('node_type', ''):
                biocatdb_nodes.append(node)
            else:
                if node['swissprot'] == 'yes':
                    if 'PDB' in node['label']:
                        pdb_nodes.append(node)
                    elif 'Rhea' in node['label']:
                        rhea_nodes.append(node)
                    else:
                        swissprot_nodes.append(node)
                else:
                    other_nodes.append(node)

        return other_nodes + swissprot_nodes + biocatdb_nodes + rhea_nodes + pdb_nodes

    def _find_uniref_metadata(self):
        node_metadata = {}
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).exclude('id', 'enzyme_type', 'sequence',
                                                                             "result_of_blasts_for")

        for seq_obj in unirefs:
            node_metadata[seq_obj.enzyme_name] = json.loads(seq_obj.to_json())

        return node_metadata

    def add_swissprot_border(self, nodes):
        spq = db.Q(sp_annotated=True)
        etq = db.Q(enzyme_type=self.enzyme_type_obj)
        swiss_prot_uniref_ids = list(UniRef50.objects(spq & etq).distinct('enzyme_name'))

        for i, node in enumerate(nodes):
            if node['id'] in swiss_prot_uniref_ids:
                nodes[i]['color']['background'] = self.swiss_prot_colour
                nodes[i]['swissprot'] = 'yes'

        return nodes

    def add_rhea_or_pdb(self, nodes):
        etq = db.Q(enzyme_type=self.enzyme_type_obj)

        pdb_q1 = db.Q(pdbs__ne=None)
        pdb_q2 = db.Q(pdbs__ne=[])
        pdb_ids = list(UniRef50.objects(etq & pdb_q1 & pdb_q2).distinct('enzyme_name'))

        rhea_q1 = db.Q(rhea__ne=None)
        rhea_q2 = db.Q(rhea__ne=[])
        rhea_ids = list(UniRef50.objects(etq & rhea_q1 & rhea_q2).distinct('enzyme_name'))
        for i, node in enumerate(nodes):
            if node['swissprot'] == 'yes':
                if node['id'] in pdb_ids:
                    nodes[i]['label'] += 'PDB'
                if node['id'] in rhea_ids:
                    if nodes[i]['label'] != '':
                        nodes[i]['label'] += ', '
                    nodes[i]['label'] += 'Rhea'
        return nodes

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"SSN_Visualiser: {msg}")

