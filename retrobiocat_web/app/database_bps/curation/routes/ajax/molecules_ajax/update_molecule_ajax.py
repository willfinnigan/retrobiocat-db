from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax.error_msgs import no_permission, \
    mol_not_found, general_error, molecule_name_already_exists, molecule_name_cannot_be_empty
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax.new_molecule_ajax import generate_update_dict
from retrobiocat_web.mongo.modal_updates import molecules_CUD
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.model_queries.molecule_queries import activity_mol_from_id, \
    does_mol_name_already_exist_in_paper


@bp.route('/_update_activity_molecule', methods=['GET', 'POST'])
@roles_required('contributor')
def update_activity_molecule():
    paper_id = request.form['paper_id']
    smi = request.form['smi']
    mol_name = request.form['mol_name']
    mol_id = request.form['mol_id']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    mol = activity_mol_from_id(mol_id)
    if mol is None:
        return mol_not_found()

    if mol_name == "":
        return molecule_name_cannot_be_empty()

    issues = molecules_CUD.update_activity_molecule(mol, mol_name, smi)
    if len(issues) != 0:
        return general_error(issues)

    result = {'status': 'success',
              'msg': 'Molecule updated',
              'issues': [],
              'update_dict': generate_update_dict(mol)}
    return jsonify(result=result)

@bp.route('/_update_activity_molecule_name_only', methods=['GET', 'POST'])
@roles_required('contributor')
def update_activity_molecule_name_only():
    paper_id = request.form['paper_id']
    mol_name = request.form['mol_name']
    mol_id = request.form['mol_id']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    mol = activity_mol_from_id(mol_id)
    if mol is None:
        return mol_not_found()

    if mol_name != mol.name and does_mol_name_already_exist_in_paper(mol_name, paper):
        return molecule_name_already_exists(mol_name)

    if mol_name == "":
        return molecule_name_cannot_be_empty()

    issues = molecules_CUD.update_activity_molecule_name(mol, mol_name)
    if len(issues) != 0:
        return general_error(issues)

    result = {'status': 'success',
              'msg': 'Molecule updated',
              'issues': [],
              'update_dict': generate_update_dict(mol)}
    return jsonify(result=result)

