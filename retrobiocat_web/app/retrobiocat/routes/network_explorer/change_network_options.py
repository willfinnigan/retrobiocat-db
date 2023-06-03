from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
import json
from flask import current_app

from retrobiocat_web.app.retrobiocat.functions.load_save_network import update_vis_options, \
    load_network_components_from_redis
from retrobiocat_web.app.retrobiocat.routes.network_explorer.functions import add_new

@bp.route("/_change_network_options", methods=["POST"])
def change_network_options():
    task_id = request.form['task_id']

    redis_data = json.loads(current_app.redis.get(task_id))
    update_dict = {'colour_arrows': request.form['edge_colours']}
    redis_data = update_vis_options(update_dict, redis_data)

    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_data)  # load network data

    nodes, edges = visualiser.nodes_edges(network.graph)

    redis_data['nodes'] = add_new(redis_data['nodes'], nodes)
    redis_data['edges'] = add_new(redis_data['edges'], edges)

    current_app.redis.mset({task_id: json.dumps(redis_data)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    result = {'nodes': nodes,
              'edges': edges,
              }

    return jsonify(result=result)

