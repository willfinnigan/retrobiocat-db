from retrobiocat_web.app.database_bps.curation import bp
from flask import request, jsonify, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.app.app import user_datastore
import json
from distutils.util import strtobool

from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.paper_status import task_update_seq_papers_status
from retrobiocat_web.app.database_bps.curation.functions.seq_obj_to_json import get_seq_table_entry
from retrobiocat_web.mongo.modal_updates import sequence_CUD, paper_CUD
from retrobiocat_web.mongo.model_queries import sequence_queries, paper_queries, activity_queries

""" 
Ajax routes for updating sequence database information 
"""



@bp.route('/_save_edited_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def save_edited_sequence():
    """Update a sequence with new data"""

    seq_obj = sequence_queries.seq_obj_from_name(request.form['original_name'])
    user = user_datastore.get_user(current_user.id)
    issues = []

    if not check_permission.check_seq_permissions(current_user.id, seq_obj):
        issues.append('User does not have access to edit this sequence')

    else:
        # update seq info (doesn't update reviewed or owner)

        # if unreviewed, can do a full update..
        if seq_obj.reviewed == False:
            issues += sequence_CUD.update(seq_obj, request.form, user=user)

            # task to update the status of all papers for this sequence
            job_id = f"update_the_status_of_all_{seq_obj.enzyme_name}_papers"
            current_app.task_queue.enqueue(task_update_seq_papers_status, seq_obj.enzyme_name, job_id=job_id)

        # if reviewed, can only update some fields..
        else:
            issues += sequence_CUD.update_non_essential_only(seq_obj, request.form, user=user)

    if len(issues) != 0:
        status = 'danger'
        msg = 'Issues updating sequence'
    else:
        status = 'success'
        msg = 'Sequence updated'

    seq_table = get_seq_table_entry(seq_obj)
    seq_table['_id'] = str(seq_obj.id)

    result = {'status': status,
              'msg': msg,
              'issues': issues,
              'seq_table_entry': seq_table,
              'seq_id': str(seq_obj.id)}

    return jsonify(result=result)


@bp.route('/_change_sequence_assign', methods=['GET', 'POST'])
@roles_required('contributor')
def change_sequence_assign():
    """Change owner of a sequence"""

    issues = []
    user = user_datastore.get_user(current_user.id)
    seq_obj = sequence_queries.seq_obj_from_name(request.form['enzyme_name'])

    # unassign if user is current owner
    if (seq_obj.owner == user) and (request.form['self_assigned'] == 'false'):
        issues += sequence_CUD.change_owner(seq_obj, None)

    # assign if no owner
    elif (seq_obj.owner == None) and (request.form['self_assigned'] == 'true'):
        print(request.form['self_assigned'])
        print('changing')
        issues += sequence_CUD.change_owner(seq_obj, user)

    # if neither, then user doesn't have access
    else:
        issues += ['No access to change assignment']

    # return issues
    if len(issues) == 0:
        status = 'success'
        msg = 'Sequence assigned'
    else:
        status = 'danger'
        msg = 'Error'

    result = {'status': status,
              'msg': msg,
              'issues': issues}
    return jsonify(result=result)


@bp.route('/_mark_sequence_reviewed', methods=['GET', 'POST'])
@roles_required('contributor')
def mark_sequence_reviewed():
    """Review or unreview a sequence"""

    issues = []
    seq_obj = sequence_queries.seq_obj_from_name(request.form['enzyme_name'])

    if not check_permission.check_seq_review(current_user.id, seq_obj):
        issues.append('User does not have access to review this sequence')

    else:
        # update review status
        if request.form['reviewed'] == 'true':
            issues += sequence_CUD.update_reviewed(seq_obj, True)
        else:
            issues += sequence_CUD.update_reviewed(seq_obj, False)

    # return response
    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Sequence review status updated',
                  'issues': issues}
    else:
        result = {'status': 'danger',
                  'msg': 'Could not change sequence review status',
                  'issues': issues}
    return jsonify(result=result)


@bp.route('/_delete_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_sequence():
    issues = []
    seq_obj = sequence_queries.seq_obj_from_name(request.form['to_delete'])

    if not check_permission.check_seq_permissions(current_user.id, seq_obj):
        issues.append('User does not have access to edit this sequence')
    else:
        issues += sequence_CUD.delete_sequence(seq_obj)

    if len(issues) == 0:
        status = 'success'
        msg = 'Sequence deleted'
    else:
        status = 'danger'
        msg = 'Could not delete sequence'

    result = {'status': status,
              'msg': msg,
              'issues': issues}

    return jsonify(result=result)


@bp.route('/_merge_seq', methods=['GET', 'POST'])
@roles_required('contributor')
def merge_sequences():
    """Merge one sequence entry into another"""

    merge_from = sequence_queries.seq_obj_from_name(request.form['to_merge'])
    merge_to = sequence_queries.seq_obj_from_name(request.form['merge_with'])
    save_other_name = bool(strtobool(str(request.form['save_other_name'])))
    issues = []

    # check for permission
    if not check_permission.check_seq_review(current_user.id, merge_from) or \
            not check_permission.check_seq_review(current_user.id, merge_to):
        issues.append('User does not have access to merge these sequences')
    else:
        # do merge if possible
        issues += sequence_CUD.merge_sequences(merge_to, merge_from, save_other_name)

    # return suitable response
    if len(issues) == 0:
        status = 'success'
        msg = 'Sequences merged'
    else:
        status = 'danger'
        msg = 'Could not merge sequences'

    result = {'status': status,
              'msg': msg,
              'issues': issues}

    return jsonify(result=result)


@bp.route('/_add_new_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def add_new_sequence():
    enzyme_type = request.form['enzyme_type']
    new_name = request.form['new_name']
    user = user_datastore.get_user(current_user.id)
    paper = paper_queries.paper_from_id(request.form['paper_id'])

    issues, seq_obj = sequence_CUD.create_new_sequence({'enzyme_name': new_name, 'enzyme_type': enzyme_type}, user=user,
                                              paper=paper)

    if len(issues) == 0:
        seq_table = get_seq_table_entry(seq_obj)
        seq_table['_id'] = str(seq_obj.id)

        result = {'status': 'success',
                  'msg': f'Sequence {new_name} created and added to paper',
                  'issues': issues,
                  'seq_table_entry': seq_table}

    else:
        result = {'status': 'danger',
                  'msg': 'Problem creating new sequence..',
                  'issues': issues}

    return jsonify(result=result)


@bp.route('/_add_existing_sequence', methods=['GET', 'POST'])
@roles_required('contributor')
def add_existing_sequence():
    existing_name = request.form['existing_name'] # can be the enzyme_name or a other_name
    paper = paper_queries.paper_from_id(request.form['paper_id'])

    seq_obj = sequence_queries.seq_obj_from_name(existing_name, include_other_names=True)

    if seq_obj is None:
        result = {'status': 'danger',
                  'msg': 'Sequence could not be found',
                  'issues': ['No sequence with that name in database']}
        return jsonify(result=result)

    issues = sequence_CUD.add_paper(seq_obj, paper)

    if len(issues) == 0:
        seq_table = get_seq_table_entry(seq_obj)
        seq_table['_id'] = str(seq_obj.id)

        # warn if the match was through an other_name
        if seq_obj.enzyme_name != existing_name:
            result = {'status': 'warning',
                      'msg': f'Sequence {seq_obj.enzyme_name} added to paper',
                      'issues': [f'Sequence has other_name of {existing_name}'],
                      'seq_table_entry': seq_table}
        else:
            result = {'status': 'success',
                      'msg': f'Sequence {existing_name} added to paper',
                      'issues': [],
                      'seq_table_entry': seq_table}
    else:
        result = {'status': 'danger',
                  'msg': 'Problem creating adding sequence to paper..',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_update_other_names', methods=['GET', 'POST'])
@roles_required('contributor')
def update_other_names():
    """Update a single other names entry"""

    enzyme_name = request.form['enzyme_name']
    other_names_dict = json.loads(request.form['other_names_dict'])
    seq_obj = sequence_queries.seq_obj_from_name(enzyme_name, get_related=False)

    if other_names_dict['existing_name'] == '':
        issues = sequence_CUD.create_other_name(seq_obj, other_names_dict)
    else:
        issues = sequence_CUD.update_other_name(seq_obj, other_names_dict)

    if len(issues) != 0:
        result = {'status': 'danger',
                  'msg': f'Issue updating other names for {enzyme_name}',
                  'issues': issues}
    else:
        result = {'status': 'success',
                  'msg': f'Other names for {enzyme_name} updated',
                  'issues': issues,
                  'name': other_names_dict['new_name']}

    return jsonify(result=result)

@bp.route('/_delete_other_names', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_other_names():
    """Delete a single other names entry"""

    enzyme_name = request.form['enzyme_name']
    other_name = request.form['other_name']
    seq_obj = sequence_queries.seq_obj_from_name(enzyme_name, get_related=False)
    issues = sequence_CUD.delete_other_name(seq_obj, other_name)

    if len(issues) != 0:
        result = {'status': 'danger',
                  'msg': f'Issue deleting other name {other_name} for {enzyme_name}',
                  'issues': issues}
    else:
        result = {'status': 'success',
                  'msg': f'Other names for {enzyme_name} deleted',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_save_alt_naming_selection', methods=['GET', 'POST'])
@roles_required('contributor')
def save_alt_naming_selection():
    alt_names_to_set = json.loads(request.form['alt_names_to_set'])
    paper_id = request.form['paper_id']

    issues = sequence_CUD.update_other_names_paper_use(alt_names_to_set, paper_id)

    if len(issues) != 0:
        result = {'status': 'danger',
                  'msg': f'Issue updating paper use of other names',
                  'issues': issues}
    else:
        result = {'status': 'success',
                  'msg': f'Other names for paper have been updated',
                  'issues': issues}

    return jsonify(result=result)

@bp.route('/_remove_seq_from_paper', methods=['GET', 'POST'])
@roles_required('contributor')
def remove_sequence():
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    seq_obj = sequence_queries.seq_obj_from_name(request.form['enzyme_name'])

    issues = sequence_CUD.remove_paper(seq_obj, paper)

    if len(issues) == 0:
        issues += paper_CUD.tag_paper_with_enzyme_types(paper)
        result = {'status': 'success',
                  'msg': 'Sequence removed from paper',
                  'issues': [],
                  'seq_id': str(seq_obj.id)}
    else:
        result = {'status': 'danger',
                  'msg': 'Issue removing sequence from paper',
                  'issues': issues,
                  'seq_id': str(seq_obj.id)}
    return jsonify(result=result)

@bp.route('/_remove_all_seq_from_paper', methods=['GET', 'POST'])
@roles_required('contributor')
def remove_all_seq_from_paper():
    paper = paper_queries.paper_from_id(request.form['paper_id'])

    if activity_queries.activity_in_paper(paper, count_only=True) != 0:
        result = {'status': 'danger',
                  'msg': 'Can not remove sequences - activity data still attached to some sequences',
                  'issues': ['Please remove references to any sequences in the activity section before trying to remove all']}
        return jsonify(result=result)

    else:
        issues = []
        seqs = sequence_queries.seqs_of_paper(paper)
        for seq in seqs:
            issues += sequence_CUD.remove_paper(seq, paper)

    if len(issues) == 0:
        issues += paper_CUD.tag_paper_with_enzyme_types(paper)
        result = {'status': 'success',
                  'msg': 'Sequences removed from paper',
                  'issues': issues}
    else:
        issues += paper_CUD.tag_paper_with_enzyme_types(paper)
        result = {'status': 'success',
                  'msg': 'Issues whilst trying to remove all sequences',
                  'issues': issues}

    return jsonify(result=result)

