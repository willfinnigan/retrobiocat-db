from flask import flash, request, jsonify
from flask_security import roles_required, current_user

from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.modal_updates import paper_CUD


@bp.route('/_delete_paper', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_paper():
    """ If possible, delete the paper """

    paper = paper_queries.paper_from_id(request.form['paper_id'])

    if not check_permission.check_paper_permission(current_user.id, paper):
        result = {'status': 'danger',
                  'msg': 'You are not the owner of this paper',
                  'issues': ['Assign this paper to yourself in order to delete it']}
        return jsonify(result=result)

    issues = paper_CUD.delete_paper(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper deleted',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Issue deleting paper',
                  'issues': issues}

    return jsonify(result=result)