import copy

from flask import request, current_app, render_template, jsonify, redirect, url_for
from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.database_bps.db_analysis.functions.task_queues import get_task_id, get_queue, get_job_registeries
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_summary import load_substrate_summary
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details


def get_hasse_id(args):
    hasse_args = copy.deepcopy(args)
    if 'start_node' in hasse_args:
        hasse_args.pop('start_node')
    job_id = get_task_id(hasse_args)
    hasse_id = 'hasse_diagram_for_' + job_id
    return hasse_id

@bp.route('/substrate_summary', methods=['GET'])
def substrate_summary():

    # get all the arguments, job_ids, and queue information together
    args = request.args.to_dict()
    job_id = get_task_id(args)
    hasse_id = get_hasse_id(args)
    queue_name = 'substrate_summary'
    queue = get_queue(current_app, queue_name)
    finished, started = get_job_registeries(queue)

    # if the job is finished, render the page
    if job_id in list(finished.get_job_ids()):
        task = queue.fetch_job(job_id)
        if task != None:
            return render_template('substrate_summary/substrate_summary.html', result=task.result)

    # if there isn't a job queued or running, then make one
    if job_id not in list(started.get_job_ids()):  # otherwise, if the job is not already queued, add it
        queue.enqueue(load_substrate_summary.task_get_hasse_network_and_summarise, hasse_id, args, job_id=job_id)

    # render the loading page
    queue_details, task_details = queue_task_details(job_id, queue_name)
    return render_template('queue_loading.html', task_queue=queue_name, task_id=job_id,
                           queue_details=queue_details, task_details='',
                           title='Loading substrate summary', ajax_timer=1500, refresh_timer=30000)
