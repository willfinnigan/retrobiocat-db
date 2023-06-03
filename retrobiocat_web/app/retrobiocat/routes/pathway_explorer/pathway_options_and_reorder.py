from flask import jsonify, request, session
from retrobiocat_web.app.retrobiocat import bp
from rq import get_current_job
import json
from flask import current_app
from rq.job import Job

from retrobiocat_web.app.retrobiocat.functions.load_save_network import load_network_components_from_redis
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.app.retrobiocat.functions import pathway_packagaing, pathway_evaluation
import datetime

from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway


def load_pathways(all_pathways_nodes, all_pathway_scores, network):
    pathways = []
    for i, list_nodes in enumerate(all_pathways_nodes):
        pathway = Pathway(list_nodes, network, calc_scores=False)
        pathway.scores.scores_from_dict(all_pathway_scores[i])
        pathways.append(pathway)
    return pathways

def task_reorder_pathways(weights, pathways_id):
    job = get_current_job()
    set_progress_bar(job, 40, 'started')

    pathway_settings = json.loads(current_app.redis.get(pathways_id + '__pathway_settings'))
    pathway_settings.update({'weight_num_enzymes': weights[0],
                             'weight_complexity': weights[1],
                             'weight_starting': weights[2],
                             'weight_known_enzymes': weights[3],
                             'weight_diversity': weights[4]})
    current_app.redis.mset({f"{pathways_id}__pathway_settings": json.dumps(pathway_settings)})
    current_app.redis.expire(pathways_id, 60 * 60)


    network_data = json.loads(current_app.redis.get(pathways_id + '__network'))

    network, retro_engine, scorer, visualiser = load_network_components_from_redis(network_data)

    all_pathways_nodes, all_scores = json.loads(current_app.redis.get(f"{pathways_id}__all_pathways"))
    pathways = load_pathways(all_pathways_nodes, all_scores, network)

    pathway_evaluator = pathway_evaluation.run_evaluate_pathways(pathways, weights)
    pathway_packagaing.package_evaluated_pathways(pathway_evaluator.df, pathways_id)
    pathway_packagaing.package_visjs_pathways(pathways_id)

    set_progress_bar(job, 90, 'complete')

    return {}


@bp.route("/_reorder_pathways_status/<task_id>", methods=["GET"])
def reorder_pathways_status(task_id):
    task = current_app.pathway_queue.fetch_job(task_id)
    task_id = task.get_id()
    task_status = task.get_status(refresh=True)

    seconds_since_active = 0
    try:
        seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
    except:
        pass
    if seconds_since_active >= 3600:
        seconds_since_active -= 3600

    if seconds_since_active > 180 and task_status != 'finished':
        print('Job no longer active')
        print(task_status)
        task_status = 'failed'
    if seconds_since_active > 600:
        print('Job no longer active')
        print(task_status)
        task_status = 'failed'

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        if task.get_status(refresh=True) == 'finished':

            response_object = {
                "status": "success",
                "data": {
                    "task_id": task_id,
                    "task_status": task_status,
                    "task_progress": progress
                }
            }
        else:
            response_object = {
                "status": "success",
                "data": {
                    "task_id": task_id,
                    "task_status": task_status,
                    "task_progress": progress
                }
            }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

#ajax call used by pathway_explorer
@bp.route('/_reorder_pathways', methods=['GET', 'POST'])
def reorder_pathways():

    weights = [json.loads(request.form['weight_num_enzymes']),
               json.loads(request.form['weight_complexity']),
               json.loads(request.form['weight_starting']),
               json.loads(request.form['weight_known_enzymes']),
               json.loads(request.form['weight_diversity'])]
    task = current_app.pathway_queue.enqueue(task_reorder_pathways, weights, request.form['id'])
    reorder_task_id = task.get_id()

    if 'pathway_reorder_task_id' in session:
        old_task_id = session['pathway_reorder_task_id']
        try:
            old_job = Job.fetch(old_task_id, connection=current_app.redis)
            old_job.delete()
        except:
            pass

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": reorder_task_id,
                "task_status": task.get_status(refresh=True),
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202




