from flask import request, jsonify
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.analysis.mol_images import smitosvg_url
from retrobiocat_web.analysis.mol_name_lookup import get_smi_from_name
from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.clean_str_of_spaces import remove_end_spaces_if_str
from retrobiocat_web.app.database_bps.curation.functions.mol_obj_to_json import mol_svg_for_mol_table
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax.error_msgs import no_permission, \
    no_name_found, error_creating_mol
from retrobiocat_web.mongo.modal_updates import molecules_CUD
from retrobiocat_web.mongo.model_queries import paper_queries


def generate_update_dict(new_mol):
    update_dict = {'_id': str(new_mol.id),
                   'smi': new_mol.smi,
                   'mol': mol_svg_for_mol_table(new_mol.smi),
                   'name': new_mol.name,
                   'chem_name': new_mol.chem_name}
    return update_dict

@bp.route('/_activity_mol_from_name', methods=['POST'])
@roles_required('contributor')
def activity_mol_from_name():
    mol_name = request.form['mol_name']
    mol_name = remove_end_spaces_if_str(mol_name)

    paper_id = request.form['paper_id']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    smi = get_smi_from_name(mol_name)
    if smi is None:
        return no_name_found()

    new_mol = molecules_CUD.make_new_activity_molecule(smi, paper, name=None, chem_name=None)
    if new_mol is None:
        return error_creating_mol()

    result = {'status': 'success',
              'msg': 'Molecule added',
              'issues': [],
              'update_dict': generate_update_dict(new_mol)}
    return jsonify(result=result)



@bp.route('/_new_activity_mol', methods=['GET', 'POST'])
@roles_required('contributor')
def new_activity_mol():
    paper_id = request.form['paper_id']
    smi = request.form['smi']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    if not check_permission.check_paper_permission(current_user.id, paper):
        return no_permission()

    new_mol = molecules_CUD.make_new_activity_molecule(smi, paper, name=None, chem_name=None)
    if new_mol is None:
        return error_creating_mol()

    result = {'status': 'success',
              'msg': 'Molecule created',
              'issues': [],
              'update_dict': generate_update_dict(new_mol)}
    return jsonify(result=result)

