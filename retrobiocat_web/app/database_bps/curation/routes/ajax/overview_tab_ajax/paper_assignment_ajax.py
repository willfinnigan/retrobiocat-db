from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.mongo.modal_updates import paper_CUD, user_CUD
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.models.user_models import user_datastore


@bp.route('/_self_assign', methods=['GET', 'POST'])
@roles_required('contributor')
def self_assign():
    """Assign paper to the user"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])
    user = user_datastore.get_user(current_user.id)

    if check_permission.can_self_assign(user) == False:
        result = {'status': 'danger',
                  'msg': 'Can not assign paper',
                  'issues': ['Users may not have more papers assigned to them than the number they have completed',
                             'Please complete papers already assigned to you before taking on additional papers']}
        return jsonify(result=result)

    if paper.owner == None:
        issues = paper_CUD.new_owner(paper, user)
    else:
        issues = ['Paper already has an owner']

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper is now assigned to you',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not assign paper',
                  'issues': issues}
    return jsonify(result=result)

@bp.route('/_un_self_assign', methods=['GET', 'POST'])
@roles_required('contributor')
def unself_assign():
    """Unassign paper if user is allowed"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])
    user = user_datastore.get_user(current_user.id)

    if check_permission.can_unassign(user, paper) == False:
        result = {'status': 'danger',
                  'msg': 'Can not unassign paper',
                  'issues': ['Papers which have had data added can not be unassigned.',
                             'Please leave a comment in the curation portal, or email admin.']}
        return jsonify(result=result)

    # Update paper assignment
    issues = []
    issues += user_CUD.paper_was_unassigned(user)
    issues += paper_CUD.new_owner(paper, None)

    # response
    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper was unassigned',
                  'issues': []}
    else:
        result = {'status': 'danger',
                 'msg': 'Issue in trying to unassign paper',
                 'issues': issues}
    return jsonify(result=result)