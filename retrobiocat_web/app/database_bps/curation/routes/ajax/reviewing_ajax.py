import distutils
from distutils.util import strtobool

from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.routes.ajax.error_results import no_access_to_edit_paper
from retrobiocat_web.app.database_bps.curation.routes.main_tab_routes.overview_tab import get_and_update_status
from retrobiocat_web.mongo.modal_updates import reviewing_updates
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.models.user_models import user_datastore


@bp.route('/_review_paper_metadata', methods=['GET', 'POST'])
@roles_required('contributor')
def review_paper_metadata():
    """ Mark a paper as reviewed """
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    new_reviewed_status = bool(strtobool(request.form['reviewed']))

    if not check_permission.check_review_persmission(current_user.id, paper):
        return no_access_to_edit_paper()

    issues = reviewing_updates.review_paper_metadata(paper, new_reviewed_status)

    get_and_update_status(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper review status updated',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Error while trying to review paper',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_sequences_ready_for_review', methods=['GET', 'POST'])
@roles_required('contributor')
def sequences_ready_for_review():
    """Update status of sequences ready to review"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])
    user = user_datastore.get_user(current_user.id)

    if request.form['ready'] == 'true':
        ready = True
    else:
        ready = False

    if check_permission.check_seq_curation_permission(user, paper) == False:
        return no_access_to_edit_paper()

    issues = reviewing_updates.sequences_ready_to_review(paper, ready)

    get_and_update_status(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Paper ready to review sequences status updated (ready={ready})',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not update paper ready to review sequences status',
                  'issues': issues}
    return jsonify(result=result)

@bp.route('/_sequences_review', methods=['GET', 'POST'])
@roles_required('contributor')
def sequences_review():
    """Update papers sequences reviewed status"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])
    user = user_datastore.get_user(current_user.id)
    reviewed = bool(distutils.util.strtobool(request.form['reviewed']))

    if reviewed == True:
        if check_permission.check_review_persmission(current_user.id, paper) == False:
            return no_access_to_edit_paper()

    elif reviewed == False:
        if check_permission.has_unreview_access(paper, user) == False:
            result = {'status': 'danger',
                      'msg': 'Can not unreview',
                      'issues': ['Only the owner of the paper can unreview']}
            return jsonify(result=result)

    issues = reviewing_updates.review_paper_sequences(paper, user, reviewed)

    get_and_update_status(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Paper sequence review status updated (reviewed={reviewed})',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not update sequence review status',
                  'issues': issues}
    return jsonify(result=result)


@bp.route('/_activity_ready_for_review', methods=['GET', 'POST'])
@roles_required('contributor')
def activity_ready_for_review():
    """Update status of sequences ready to review"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])

    if check_permission.check_paper_permission(current_user.id, paper) == False:
        return no_access_to_edit_paper()

    ready = bool(distutils.util.strtobool(request.form['ready']))
    issues = reviewing_updates.activity_ready_to_review(paper, ready)

    get_and_update_status(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Paper activity data ready to review status updated (ready={ready})',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not update paper ready to review activity status',
                  'issues': issues}
    return jsonify(result=result)

@bp.route('/_activity_review', methods=['GET', 'POST'])
@roles_required('contributor')
def activity_review():
    """Update papers sequences reviewed status"""

    paper = paper_queries.paper_from_id(request.form['paper_id'])
    user = user_datastore.get_user(current_user.id)
    reviewed = bool(distutils.util.strtobool(request.form['reviewed']))

    if reviewed == True:
        if check_permission.check_review_persmission(current_user.id, paper) == False:
            return no_access_to_edit_paper()
    elif reviewed == False:
        if check_permission.has_unreview_access(paper, user) == False:
            result = {'status': 'danger',
                      'msg': 'Can not unreview',
                      'issues': ['Only the owner of the paper can unreview']}
            return jsonify(result=result)

    issues = reviewing_updates.review_paper_activity(paper, user, reviewed)

    get_and_update_status(paper)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Paper activity review status updated (reviewed={reviewed})',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not update activity review status',
                  'issues': issues}
    return jsonify(result=result)

