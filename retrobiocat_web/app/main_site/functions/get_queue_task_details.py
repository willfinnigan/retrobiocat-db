from flask import current_app

def queue_task_details(task_id, queue_name):
    queue_details = "Error - Queue not found"  # default
    task_details = 'Error - Task not found'  # default

    queue = current_app.redis_queues_dict.get(queue_name, None)
    if queue:
        queue_details = f"There are currently {len(queue.jobs)} other jobs in {queue_name} queue"
        task = queue.fetch_job(task_id)
        if task:
            task_details = f'Task: {task_id}'

    return queue_details, task_details

def get_task(queue_name, task_id):
    queue = current_app.redis_queues_dict.get(queue_name, None)
    task = queue.fetch_job(task_id)
    return task

def get_task_status(queue_name, task_id):
    queue = current_app.redis_queues_dict.get(queue_name, None)
    if not queue:
        return 'Error'

    task = queue.fetch_job(task_id)
    if not task:
        return 'Error'

    task_id = task.get_id()
    task_status = task.get_status(refresh=True)
    return task_status
