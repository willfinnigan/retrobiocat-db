from retrobiocat_web.app.retrobiocat import bp
from flask import render_template
import json
from flask import current_app, redirect
import uuid
from retrobiocat_web.mongo.models.user_saves import Network
from retrobiocat_web.app.app import user_datastore
from flask_security import current_user
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details


# This function is not used and can be deleted
def get_users_saves(all_saves):
    user = user_datastore.get_user(current_user.id)
    user_saves = []
    for save_tuple in all_saves:
        if len(Network.objects(uuid=save_tuple[2])) > 0:
            mongo_network = Network.objects(uuid=save_tuple[2])[0]
            if mongo_network.owner == user:
                user_saves.append(save_tuple)
    return user_saves

def load_from_mongo(id):
    data = False
    if len(Network.objects(uuid=id)) > 0:

        if current_user.is_authenticated:
            user = user_datastore.get_user(current_user.id)
        else:
            user = False

        mongo_network = Network.objects(uuid=id)[0]
        if (mongo_network.public == True) or (mongo_network.owner == user):
            data = mongo_network.data
            if mongo_network.owner != user:
                data['save_id'] = str(uuid.uuid4())
                data['save_links'] = []

    data = patch_old_retrobiocat_attributes_to_new(data)
    return data


def patch_old_retrobiocat_attributes_to_new(mongo_data):
    # create default condig options
    attr_dict = json.loads(mongo_data['attr_dict'])

    if 'retrosynthesis_config' not in mongo_data:
        mongo_data['retrosynthesis_config'] = None
        mongo_data['scorer_config'] = None
        mongo_data['vis_config'] = None

    # copy all attributes to metadata if it doesnt exist
    for node in attr_dict:
        if 'attributes' in attr_dict[node]:
            if 'metadata' not in attr_dict[node]['attributes']:
                attr_dict[node]['attributes']['metadata'] = attr_dict[node]['attributes']

    mongo_data['attr_dict'] = json.dumps(attr_dict)

    return mongo_data

def ensure_nodes_have_location(nodes_json):
    for i, node in enumerate(nodes_json):
        if 'x' not in node:
            nodes_json[i]['x'] = 0
        if 'y' not in node:
            nodes_json[i]['y'] = 0
    return nodes_json

def load_from_redis(id):
    if current_app.redis.exists(id):
        return json.loads(current_app.redis.get(id))

    return None

def save_new_redis(data):
    new_task_id = str(uuid.uuid4())
    current_app.redis.mset({new_task_id: json.dumps(data)})
    current_app.redis.expire(new_task_id, 5*60)
    return new_task_id

@bp.route("/network_explorer/<task_id>/", methods=["GET"])
def network_explorer(task_id):

    # if task in mongo
    save_query = Network.objects(uuid=task_id).first()
    if save_query is not None:
        user = user_datastore.get_user(current_user.id)
        if save_query.public is True or save_query.owner == user:
            data = load_from_mongo(task_id)
            new_task_id = save_new_redis(data)
            return redirect('/network_explorer/' + new_task_id + '/')

    # if result is saved in redis
    result = load_from_redis(task_id)
    if result is not None:
        result = load_from_redis(task_id)
        result['nodes'] = ensure_nodes_have_location(result['nodes'])

        return render_template('network_explorer/network_explorer.html',
                               save_id=result['save_id'],
                               save_name=result['save_name'],
                               save_links=result['save_links'],
                               nodes=result['nodes'],
                               edges=result['edges'],
                               options=result['options'],
                               max_reactions=json.loads(result['retrosynthesis_config'])['max_reactions'],
                               rr_diameter=json.loads(result['retrosynthesis_config'])['rr_diameter'],
                               rr_threshold=json.loads(result['retrosynthesis_config'])['rr_threshold'],
                               retrorules_diameters=[],
                               task_id=task_id)


    # otherwise, job should be in the queue
    queue_name = 'network'
    task = current_app.network_queue.fetch_job(task_id)

    if not task:
        task_id = 'task_not_found'
        task_status = 'task_not_found'
        queue_details = 'Error - task not found'
        task_details = 'Error - task not found'
    else:
        task_id = task.get_id()
        task_status = task.get_status(refresh=True)
        queue_details, task_details = queue_task_details(task_id, queue_name)

    if task_status != 'finished':
        return render_template('queue_loading.html', task_queue=queue_name, task_id=task_id,
                               queue_details=queue_details, task_details=task_details,
                               title='Network explorer', ajax_timer=3000, refresh_timer=30000)
    else:
        return redirect('/network_explorer/' + task_id + '/')



