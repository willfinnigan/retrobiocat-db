@bp.route("/substrate_specificity_form_status/<task_id>", methods=["GET"])
def substrate_specificity_status(task_id):

    print(task_id)
    task = current_app.task_queue.fetch_job(task_id)
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    task_id = task.get_id()
    task_status = task.get_status(refresh=True)
    seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
    if seconds_since_active > 3600:
        seconds_since_active = 0

    if seconds_since_active > 180:
        print('Job no longer active')
        task_status = 'failed'

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task_id,
                "task_status": task_status,
                "task_progress" : progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202
