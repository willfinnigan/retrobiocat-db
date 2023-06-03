import networkx as nx
import json
from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request, current_app
from retrobiocat_web.retro.network_pathway.network import Network


@bp.route('/_change_enzyme', methods=['GET', 'POST'])
def change_enzyme():
    selected_node = request.form['selected_node']
    selected_enzyme = request.form['selected_enzyme']

    task_id = request.form['task_id']
    data = json.loads(current_app.redis.get(task_id))
    graph_dict = json.loads(data['graph_dict'])
    attr_dict = json.loads(data['attr_dict'])
    target_smiles = data['target_smiles']
    network_options = json.loads(data['network_options'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.settings = network_options
    network.add_attributes(attr_dict)

    network.calculate_scores()

    network.graph.nodes[selected_node]['attributes']['selected_enzyme'] = selected_enzyme

    data['attr_dict'] = json.dumps(network.attributes_dict())
    current_app.redis.mset({task_id: json.dumps(data)})
    time_to_expire = 15 * 60  # 15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    successors = list(network.graph.successors(selected_node))
    predecessors = list(network.graph.predecessors(selected_node))

    subgraph = network.graph.subgraph([selected_node]+successors+predecessors)
    nodes, edges = network.get_visjs_nodes_and_edges(graph=subgraph)

    result = {'nodes': nodes,
              'edges': edges}

    return jsonify(result=result)