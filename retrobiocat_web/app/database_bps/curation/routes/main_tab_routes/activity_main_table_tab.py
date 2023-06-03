from flask import render_template, flash, redirect, url_for
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_access_dict import get_access_dict
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.app.database_bps.tables.get_table_data.get_activity_edit_table_data import get_activity_data
from retrobiocat_web.app.database_bps.tables.get_table_data.get_sequence_table_data import enzyme_names_to_paper_names
from retrobiocat_web.mongo.model_queries import paper_queries, sequence_queries
from retrobiocat_web.mongo.model_queries.molecule_queries import get_molecules_in_paper
from retrobiocat_web.mongo.model_queries.reaction_queries import reaction_names, reactions_of_type
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_of_paper, enzyme_type_abbreviations_in_paper
from retrobiocat_web.mongo.models.user_models import user_datastore

def get_paper_molecules(paper):
    mols = get_molecules_in_paper(paper)

    mol_list = []
    for mol in mols:
        mol_list.append((mol.name, mol.smi, mol.svg, str(mol.id)))

    return mol_list

def get_smi_to_name_map(paper):
    smi_to_name_map = {}
    name_to_smi_map = {}
    for mol in get_molecules_in_paper(paper):
        smi_to_name_map[mol.smi] = mol.name
        name_to_smi_map[mol.name] = mol.smi
    return smi_to_name_map, name_to_smi_map


def get_enzyme_type_reactions(paper):
    et_ebbrevs = enzyme_type_abbreviations_in_paper(paper)
    reactions = reactions_of_type(et_ebbrevs, names_only=True, include_chemical=False)
    return reactions

@bp.route('/curation_activity_main/<paper_id>', methods=['GET'])
@roles_required('contributor')
def curation_activity_main(paper_id):
    paper = paper_queries.paper_from_id(paper_id, get_related=True)
    user = user_datastore.get_user(current_user.id)

    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))

    if not check_permission.check_seq_curation_permission(user, paper):
        flash('No access to edit sequences for this entry', 'fail')
        return redirect(url_for("main_site.home"))

    reviewed_dict = get_reviewed_dict(paper)
    access_dict = get_access_dict(paper, user)
    enzyme_alt_names_to_use = sequence_queries.get_alt_naming_for_paper(paper)

    activity_data = get_activity_data(paper, enzyme_alt_names_to_use)
    reactions = get_enzyme_type_reactions(paper)
    if len(reactions) == 0:
        reactions = reaction_names()
    enzyme_names = seqs_of_paper(paper, names_only=True)
    enzyme_names = enzyme_names_to_paper_names(enzyme_names, enzyme_alt_names_to_use)

    paper_molecules = get_paper_molecules(paper)
    smi_to_name_map, name_to_smi_map = get_smi_to_name_map(paper)

    return render_template('curation_activity/curation_activity.html',
                           paper_id=paper.id,
                           paper_title=paper.title,
                           reviewed_dict=reviewed_dict,
                           access_dict=access_dict,
                           activity_data=activity_data,
                           reactions=reactions,
                           enzyme_names=enzyme_names,
                           paper_molecules=paper_molecules,
                           smi_to_name_map=smi_to_name_map,
                           name_to_smi_map=name_to_smi_map)
