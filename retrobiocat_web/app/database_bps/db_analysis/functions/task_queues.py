from rq.registry import FinishedJobRegistry, StartedJobRegistry

def get_queue(current_app, queue_name):
    return current_app.redis_queues_dict.get(queue_name, None)

def get_task_id(args):
    job_id = ""
    arg_keys = sorted(list(args.keys()))
    for key in arg_keys:
        value = args[key]
        if job_id != "":
            job_id += '__'
        job_id += f"{key}={value}"
    return job_id

def get_job_registeries(queue):
    finished_reg = FinishedJobRegistry(queue=queue)
    started_reg = StartedJobRegistry(queue=queue)
    return finished_reg, started_reg


def args_to_details(args):
    details = ""

    for key, value in args.items():
        details += f"{key}: {value} \n"

    return details