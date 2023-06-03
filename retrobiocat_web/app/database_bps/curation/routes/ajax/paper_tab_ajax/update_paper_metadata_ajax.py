import json
from distutils.util import strtobool

from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.mongo.modal_updates import paper_CUD
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.models.user_models import user_datastore


@bp.route('/_save_updated_paper_metadata', methods=['GET', 'POST'])
@roles_required('contributor')
def save_updated_paper_metadata():
    """ Save updates to the metadata for a paper """

    user = user_datastore.get_user(current_user.id)
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    update_dict = json.loads(request.form['update_dict'])

    issues = paper_CUD.update_paper(paper, update_dict, user=user)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper information updated',
                  'issues': issues}
    else:
        result = {'status': 'danger',
                  'msg': 'Problem updating paper data',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_update_paper_importance', methods=['GET', 'POST'])
@roles_required('contributor')
def update_paper_importance():
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    if not check_permission.check_review_persmission(current_user.id, paper):
        result = {'status': 'danger',
                  'msg': 'No access to update importance of this paper',
                  'issues': [f'No paper permission for {paper.short_citation}']}
        return jsonify(result=result)

    importance = bool(strtobool(request.form['importance']))
    issues = paper_CUD.update_paper_importance(paper, importance)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper importance updated',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Error updating paper importance',
                  'issues': issues}
    return jsonify(result=result)