from retrobiocat_web.app.main_site import bp
from flask import render_template, jsonify, current_app, request

FAILED_PROGRESS = {'progressbar_style': 'width: 80%',
                    'progressbar_text': 'failed',
                    'progressbar_classname': 'progress-bar bg-danger'
                    }

QUEUEING_PROGRSS = {'progressbar_style': 'width: 33%',
                    'progressbar_text': 'queuing',
                    'progressbar_classname': 'progress-bar'
                    }


@bp.route("/loading_page/<task_id>/<task_queue>", methods=["GET"])
def loading_page(task_id, task_queue):
    args = request.args.to_dict()
    title = args.get('title', 'Working..')
    task_details = args.get('task_details', '')
    ajax_timer = args.get('ajax_timer', 3000)
    refresh_timer = args.get('refresh_timer', 30000)

    queue_details = "Error - Queue not found"  # default
    queue = current_app.redis_queues_dict.get(task_queue, None)
    if queue:
        queue_details = f"{len(queue.jobs)} jobs in queue"

    return render_template('queue_loading.html',
                           task_id=task_id, task_queue=task_queue,
                           queue_details=queue_details, task_details=task_details,
                           title=title, ajax_timer=ajax_timer, refresh_timer=refresh_timer)


@bp.route("/queue_status/<task_id>/<task_queue>", methods=["GET"])
def queue_status(task_id, task_queue):
    queue = current_app.redis_queues_dict.get(task_queue, None)

    status = 'danger'
    task_status = 'failed'
    progress = FAILED_PROGRESS

    if queue:
        task = queue.fetch_job(task_id)
        if task:
            status = 'success'
            task_id = task.get_id()
            try:
                task_status = task.get_status(refresh=True)
            except:
                task_status = 'unknown'

            seconds_since_active = 0
            if seconds_since_active > 3600:  # hack for clocks going back
                seconds_since_active -= 3600

            if task_status == 'failed':
                progress = FAILED_PROGRESS

            elif seconds_since_active > 180 and task_status != 'finished':
                task_status = 'failed'
                progress = FAILED_PROGRESS
            elif 'progress' in task.meta:
                progress = task.meta['progress']
            else:
                progress = QUEUEING_PROGRSS


    response_object = {
        "status": status,
        "data": {
            "task_id": task_id,
            "task_status": task_status,
            "task_progress": progress}
        }
    return jsonify(response_object), 202
