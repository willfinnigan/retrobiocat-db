from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax.error_msgs import no_permission, \
    mol_not_found, general_error
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.model_queries.molecule_queries import activity_mol_from_id, get_molecules_in_paper
from retrobiocat_web.mongo.modal_updates import molecules_CUD

@bp.route('/_delete_activity_molecule', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_activity_molecule():
    mol_id = request.form['mol_id']
    paper_id = request.form['paper_id']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    print(mol_id)

    mol_to_delete = activity_mol_from_id(mol_id)
    if mol_to_delete is None:
        return mol_not_found()

    issues = molecules_CUD.delete_activity_molecule(mol_to_delete)
    if len(issues) != 0:
        return general_error(issues)

    result = {'status': 'success',
              'msg': 'Molecule deleted',
              'issues': []}
    return jsonify(result=result)

def should_mol_be_deleted(mode, mol):
    if (mode == 'all') or (mol.smi == "") or (mol.smi == None):
        return True
    return False

@bp.route('/_delete_many_paper_molecules', methods=['GET', 'POST'])
@roles_required('contributor')
def delete_many_paper_molecules():
    mode = request.form['mode']
    paper_id = request.form['paper_id']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    deleted_ids = []
    all_issues = []
    for mol in get_molecules_in_paper(paper):
        mol_id = str(mol.id)
        if should_mol_be_deleted(mode, mol):
            issues = molecules_CUD.delete_activity_molecule(mol)
            if len(issues) == 0:
                deleted_ids.append(mol_id)
            else:
                all_issues += issues

    if len(all_issues) != 0:
        return general_error(all_issues)

    result = {'status': 'success',
              'msg': 'Molecules deleted',
              'issues': [],
              'deleted': deleted_ids}
    return jsonify(result=result)
