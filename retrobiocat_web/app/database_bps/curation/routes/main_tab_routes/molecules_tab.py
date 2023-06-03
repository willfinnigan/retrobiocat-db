from flask import flash, redirect, url_for, render_template
from flask_login import current_user
from flask_security import roles_required
from natsort import natsorted

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.app.database_bps.curation.functions.mol_obj_to_json import get_new_mol_json
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.model_queries.molecule_queries import get_molecules_in_paper
from retrobiocat_web.mongo.models.user_models import user_datastore


def get_molecules_data(paper):
    mol_data = []
    for mol in get_molecules_in_paper(paper):
        new_dict = get_new_mol_json(mol)
        mol_data.append(new_dict)

    ordered_mol_data = natsorted(mol_data, key=lambda x: x['name'])
    return ordered_mol_data

@bp.route('/molecules_tab/<paper_id>', methods=['GET'])
@roles_required('contributor')
def molecules_tab(paper_id):
    paper = paper_queries.paper_from_id(paper_id, get_related=True)
    user = user_datastore.get_user(current_user.id)

    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))

    if not check_permission.check_seq_curation_permission(user, paper):
        flash('No access to edit sequences for this entry', 'fail')
        return redirect(url_for("main_site.home"))

    mol_data = get_molecules_data(paper)
    reviewed_dict = get_reviewed_dict(paper)

    return render_template('curation_molecules/molecules_tab.html',
                           paper_id=paper.id,
                           paper_title=paper.title,
                           reviewed_dict=reviewed_dict,
                           molecules_data=mol_data)

