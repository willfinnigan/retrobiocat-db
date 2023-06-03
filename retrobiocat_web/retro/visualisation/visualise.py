import retrobiocat_web.retro.network_pathway.graph_functions
from retrobiocat_web.retro.visualisation import node_information,\
    edge_information, colour_nodes, colour_edges, node_positions

import jinja2
from IPython.display import HTML
from pathlib import Path

from retrobiocat_web.retro.visualisation import rdkit_images
from retrobiocat_web.retro.visualisation.config import Visualiser_Config


class Visualiser(object):

    def __init__(self, config=None):
        self.config = config
        if self.config is None:
            self.config = Visualiser_Config()

    def html(self, graph, height='750px', width='1000px', options=None):
        nodes, edges = self.nodes_edges(graph)

        if options == None:
            options = {}

        path = str(Path(__file__).parents[0]) + '/visjs_notebook_template.html'
        with open(path) as template_html:
            content = template_html.read()

        template = jinja2.Template(content)
        rendered = template.render(nodes=nodes, edges=edges, height=height, width=width, options=options)
        html_object = HTML(rendered)

        return html_object

    def nodes_edges(self, graph):
        nodes, edges = self._get_nodes_and_edges(graph)
        nodes = self._set_node_positions(nodes, graph)
        nodes = self._add_node_information(nodes, graph)
        nodes = self._add_node_colour(nodes, graph)
        edges = self._add_edge_colour(edges, graph)
        edges = self._add_edge_information(edges, graph)

        return nodes, edges

    def _get_nodes_and_edges(self, graph):
        molSize = self.config.molSize

        nodes = []
        edges = []

        for node in retrobiocat_web.retro.network_pathway.graph_functions.get_substrate_nodes(graph):
            url = rdkit_images.smile_to_svg_url(node, size=molSize)

            nodes.append({"id": node,
                          "shape": "circularImage",
                          "size": 50,
                          "borderWidth": 3,
                          "borderWidthSelected": 6,
                          "color": 'blue',
                          "image": url,
                          "label": "",
                          "title": node,
                          "type": "substrate"})

        for node in retrobiocat_web.retro.network_pathway.graph_functions.get_reaction_nodes(graph):

            nodes.append({"id": node,
                          "size": 10,
                          "borderWidth": 2,
                          "borderWidthSelected": 4,
                          "color": "darkgrey",
                          "label": graph.nodes[node]['attributes']['name'],
                          "title": node,
                          "shape": "dot",
                          "type": "reaction",
                          "reaction_type": graph.nodes[node]['attributes'].get('reaction_type', 'retrobiocat'),
                          'metadata': graph.nodes[node]['attributes'].get('metadata', {})
                          })

        for e in graph.edges:
            edges.append({"id": f"from {str(e[1])} to {str(e[0])}",
                          "from": e[1],
                          "to": e[0],
                          "arrows": "to",
                          "color": "grey"})

        return nodes, edges

    def _add_node_information(self, nodes, graph):
        nodes = node_information.add_substrate_info(nodes, graph)
        nodes = node_information.add_enzyme_info(nodes, graph)
        return nodes

    def _add_edge_information(self, edges, graph):
        if self.config.display_cofactors:
            edges = edge_information.add_cofactors(graph, edges)
        return edges

    def _add_edge_colour(self, edges, graph):
        if self.config.colour_arrows == 'Complexity change':
            edges = colour_edges.colour_edges_by_change_in_complexity(graph, edges)
        return edges

    def _add_node_colour(self, nodes, graph):

        self.nodes = colour_nodes.colour_substrate_nodes_by_building_block(nodes, graph)

        if self.config.colour_reactions == 'Substrate specificity':
            nodes = colour_nodes.colour_reactions_by_activity(nodes, graph, show_negative=self.config.show_negative_enzymes)
        elif self.config.colour_reactions == 'Complexity change':
            nodes = colour_nodes.colour_reactions_by_change_in_complexity(nodes, graph)

        nodes = colour_nodes.colour_target_node(nodes, graph)
        return nodes

    def _set_node_positions(self, nodes, graph):
        if self.config.fix_target == True:
            node_positions.fix_target_position(nodes, graph)
        return nodes


if __name__ == '__main__':
    from retrobiocat_web.retro.network_pathway.network import Network
    network = Network()




