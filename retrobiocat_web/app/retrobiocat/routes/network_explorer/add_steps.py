from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
from rq.job import get_current_job
import networkx as nx
import json
from flask import current_app

from retrobiocat_web.app.retrobiocat.functions.load_save_network import load_network_components_from_redis, \
    update_retrosynthesis_config_for_reaction_settings
from retrobiocat_web.app.retrobiocat.routes.network_explorer.functions import add_new, delete_nodes_and_edges

#ajax call used by network_explorer
@bp.route('/_delete_step', methods=['GET', 'POST'])
def delete_step():
    reaction = request.form['reaction']
    task_id = request.form['task_id']
    redis_data = json.loads(current_app.redis.get(task_id))
    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_data)

    to_delete = retro_engine.delete_reaction_node(reaction)
    nodes = []
    edges = []

    redis_data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    redis_data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(redis_data['nodes'], nodes)
    edges = add_new(redis_data['edges'], edges)
    nodes, edges = delete_nodes_and_edges(to_delete, nodes, edges)
    redis_data['nodes'] = nodes
    redis_data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(redis_data)})
    time_to_expire = 15 * 60  # 15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    result = {'to_delete': to_delete}
    return jsonify(result=result)


#ajax call used by network_explorer
@bp.route('/_step', methods=['GET', 'POST'])
def step():
    clicked_node = request.form['smiles']
    x = request.form['x']
    y = request.form['y']
    task_id = request.form['task_id']
    reaction_type = request.form['reaction_type']

    redis_data = json.loads(current_app.redis.get(task_id))

    # update any settings which have changed
    update_settings = {'max_reactions':  int(request.form['max_reactions']),
                       'aizynth_reaction_filter': request.form['aizynth_reaction_mode'],
                       'rr_diameter': int(request.form['rr_diameter']),
                       'rr_threshold': float(request.form['rr_threshold'])}
    redis_data = update_retrosynthesis_config_for_reaction_settings(update_settings, redis_data)

    # load network data
    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_data)

    new_substrate_nodes, new_reaction_nodes = [], []
    if reaction_type == 'biocatalysis':
        new_substrate_nodes, new_reaction_nodes = retro_engine.step_retrobiocat(clicked_node)
    elif reaction_type == 'biosynthesis':
        new_substrate_nodes, new_reaction_nodes = retro_engine.step_retrorules(clicked_node)
    elif reaction_type == 'chemistry':
        new_substrate_nodes, new_reaction_nodes = retro_engine.step_aizynthfinder(clicked_node)
    elif reaction_type == 'ringbreaker':
        new_substrate_nodes, new_reaction_nodes = retro_engine.step_ringbreaker(clicked_node)

    print(reaction_type)

    scorer.score()

    all_new_nodes = [clicked_node] + new_substrate_nodes + new_reaction_nodes
    subgraph = network.graph.subgraph(all_new_nodes)

    nodes, edges = visualiser.nodes_edges(graph=subgraph)

    for i, node in enumerate(nodes):
        nodes[i].update({'x': x,
                         'y': y})

    result = {'nodes': nodes,
              'edges': edges}

    redis_data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    redis_data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(redis_data['nodes'], nodes)
    edges = add_new(redis_data['edges'], edges)
    nodes, edges = delete_nodes_and_edges([], nodes, edges)
    redis_data['nodes'] = nodes
    redis_data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(redis_data)})
    time_to_expire = 15*60   # 15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    return jsonify(result=result)


#ajax call used by main_site to add custom reaction
@bp.route('/_custom_reaction', methods=['GET', 'POST'])
def custom_reaction():
    product_smiles = str(request.form['product'])
    substrate_smiles = str(request.form['substrate']).split('.')
    reaction_name = str(request.form['name'])
    task_id = request.form['task_id']

    redis_data = json.loads(current_app.redis.get(task_id))

    # load network data
    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_data)

    new_substrate_nodes, new_reaction_nodes = retro_engine.custom_reaction(product_smiles, substrate_smiles, reaction_name)

    all_new_nodes_plus_target = new_substrate_nodes + new_reaction_nodes + [product_smiles]
    subgraph = network.graph.subgraph(all_new_nodes_plus_target)
    nodes, edges = visualiser.nodes_edges(graph=subgraph)

    result = {'nodes': nodes,
              'edges': edges,
              }

    redis_data['graph_dict'] = json.dumps(nx.to_dict_of_lists(network.graph))
    redis_data['attr_dict'] = json.dumps(network.attributes_dict())
    nodes = add_new(redis_data['nodes'], nodes)
    edges = add_new(redis_data['edges'], edges)
    redis_data['nodes'] = nodes
    redis_data['edges'] = edges

    current_app.redis.mset({task_id: json.dumps(redis_data)})
    current_app.redis.expire(task_id, 5*60)

    return jsonify(result=result)