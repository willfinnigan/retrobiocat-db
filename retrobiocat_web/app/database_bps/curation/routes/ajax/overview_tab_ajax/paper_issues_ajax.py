from distutils.util import strtobool

from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.mongo.modal_updates import paper_CUD
from retrobiocat_web.mongo.model_queries import paper_queries


@bp.route('/_paper_issues', methods=['GET', 'POST'])
@roles_required('contributor')
def paper_issues():
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    issue_status = bool(strtobool(request.form['issues']))

    if not check_permission.check_review_persmission(current_user.id, paper):
        result = {'status': 'danger',
                  'msg': 'No access to add issues to this paper',
                  'issues': [f'No paper permission for {paper.short_citation}']}
        return jsonify(result=result)

    issues = paper_CUD.update_paper_issues_status(paper, issue_status)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Issue status updated',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Error updating issue status',
                  'issues': issues}
    return jsonify(result=result)