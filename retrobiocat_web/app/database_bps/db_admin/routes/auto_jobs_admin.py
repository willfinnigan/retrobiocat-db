from retrobiocat_web.app.database_bps.db_admin import bp
from flask import render_template,request, jsonify, current_app
from flask_security import roles_required, current_user


@bp.route('/auto_jobs', methods=['GET'])
@roles_required('admin')
def auto_jobs():
    job_queues = current_app.redis_queues
    queues_dict = {}
    for queue in job_queues:
        queues_dict[queue.name] = len(queue.jobs)

    return render_template('auto_jobs.html', status='', running='', queues=queues_dict, last_ran='')

@bp.route('/_clear_queue', methods=['GET', 'POST'])
@roles_required('admin')
def clear_queue():
    name = request.form['name']
    print(f'Clearing {name}')
    for q in current_app.redis_queues:
        if str(q.name) == str(name):
            q.delete()
            result = {'status': 'success',
                      'msg': f'Queue {name} cleared',
                      'issues': []}
            return jsonify(result=result)

    result = {'status': 'danger',
              'msg': f'Queue {name} not found',
              'issues': []}
    return jsonify(result=result)
