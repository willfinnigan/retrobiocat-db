from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify,  request
from flask import current_app


@bp.route('/_keep_session_open', methods=['GET', 'POST'])
def keep_session_open():
    task_id = request.form['task_id']
    current_app.redis.expire(task_id, 5*60)
    return jsonify({'status': 'success'})