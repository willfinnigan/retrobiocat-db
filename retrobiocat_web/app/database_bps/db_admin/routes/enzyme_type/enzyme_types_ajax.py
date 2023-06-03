from flask import request, jsonify
from flask_security import roles_required

from retrobiocat_web.app.database_bps.db_admin import bp
from retrobiocat_web.mongo.modal_updates import enzyme_type_CUD
from retrobiocat_web.mongo.model_queries import enzyme_type_queries

@bp.route('/_load_enzyme_type_data', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def load_enzyme_type_data():
    """Ajax load enzyme type details"""

    enz_type = enzyme_type_queries.enzyme_type_from_name(request.form['enzyme_type'])

    if enz_type.rep_reaction != None:
        rep_reaction = enz_type.rep_reaction.name
    else:
        rep_reaction = ''

    result = {'name': enz_type.enzyme_type,
              'description': enz_type.description,
              'full_name': enz_type.full_name,
              'other_abbreviations': enz_type.other_abbreviations,
              'rep_reaction': rep_reaction}

    return jsonify(result=result)

@bp.route('/_save_enzyme_type_changes', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def save_enzyme_type_changes():
    """Ajax call to update the details of an enzyme type"""

    new_name = request.form['new_name']
    description = request.form['description']
    full_name = request.form['full_name']
    other_abbreviations = request.form['other_abbreviations']
    rep_reaction_name = request.form['rep_reaction']

    enz_type = enzyme_type_queries.enzyme_type_from_name(request.form['original_name'])

    issues = enzyme_type_CUD.update_enzyme_type(enz_type,
                                                new_name, description, full_name, other_abbreviations, rep_reaction_name)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f'Changes to enzyme_type {new_name} saved',
                  'issues': []}
    else:
        result = {'status': 'danger',
                  'msg': f"Issues updating enzyme_type {request.form['original_name']}",
                  'issues': issues}
    return jsonify(result=result)

@bp.route('/_merge_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def merge_enzyme_types():
    """Ajax call to merge two enzyme types"""

    et_to_merge = enzyme_type_queries.enzyme_type_from_name(request.form['to_merge'])
    et_merge_with = enzyme_type_queries.enzyme_type_from_name(request.form['merge_with'])

    issues = enzyme_type_CUD.merge_enzyme_type(et_to_merge, et_merge_with)

    if len(issues) != 0:
        result = {'status': 'danger',
                  'msg': "Issues trying to merge enzyme types",
                  'issues': issues}
    else:
        result = {'status': 'success',
                  'msg': "Enzyme types merged successfully",
                  'issues': issues}

    return jsonify(result=result)


@bp.route('/_delete_enzyme_type', methods=['GET', 'POST'])
@roles_required('enzyme_types_admin')
def delete_enzyme_types():
    """Ajax call to delete an enzyme type"""

    enz_type = enzyme_type_queries.enzyme_type_from_name(request.form['to_delete'])
    issues = enzyme_type_CUD.delete_enzyme_type(enz_type)

    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': f"Enzyme type {request.form['to_delete']} deleted",
                  'issues': issues}
    else:
        result = {'status': 'success',
                  'msg': f"Could not delete Enzyme type {request.form['to_delete']}",
                  'issues': issues}

    return jsonify(result=result)



