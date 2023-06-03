from flask import request, jsonify
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.mongo.model_queries.user_queries import get_user_from_id
from retrobiocat_web.mongo.modal_updates import sequence_CUD, paper_CUD
from retrobiocat_web.mongo.model_queries import paper_queries, sequence_queries, activity_queries

@bp.route('/_admin_set_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_set_owner():
    paper = paper_queries.paper_from_id(request.form['paper_id'], get_related=False)
    new_owner_id = request.form['new_owner_id']
    new_owner = get_user_from_id(new_owner_id)

    issues = paper_CUD.new_owner(paper, new_owner)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Paper owner updated',
                  'issues': []}
    else:
        result = {'status': 'damger',
                  'msg': 'Could not update paper owner',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_admin_activity_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_activity_to_owner():
    paper = paper_queries.paper_from_id(request.form['paper_id'], get_related=False)
    activities = activity_queries.activity_in_paper(paper)

    for activity in activities:
        activity.added_by = paper.owner
        activity.save()

    result = {'status': 'success',
              'msg': 'Activity added by updated',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_admin_unassigned_seqs_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_unassigned_seqs_to_owner():
    paper = paper_queries.paper_from_id(request.form['paper_id'], get_related=False)
    seqs = sequence_queries.seqs_of_paper(paper, unassigned_only=True)

    issues = []
    for seq in seqs:
        issues += sequence_CUD.change_owner(seq, paper.owner)

    result = {'status': 'success',
              'msg': 'Unassigned sequences assigned to paper owner',
              'issues': issues}

    return jsonify(result=result)

@bp.route('/_admin_all_seqs_to_owner', methods=['GET', 'POST'])
@roles_required('admin')
def admin_all_seqs_to_owner():
    paper = paper_queries.paper_from_id(request.form['paper_id'], get_related=False)
    seqs = sequence_queries.seqs_of_paper(paper, unassigned_only=False)

    issues = []
    for seq in seqs:
        issues += sequence_CUD.change_owner(seq, paper.owner)

    result = {'status': 'success',
              'msg': 'All sequences assigned to paper owner',
              'issues': issues}

    return jsonify(result=result)

