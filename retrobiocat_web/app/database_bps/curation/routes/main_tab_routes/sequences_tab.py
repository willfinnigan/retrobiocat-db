from flask import flash, redirect, url_for, render_template
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_access_dict import get_access_dict
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.app.database_bps.tables.get_table_data.get_sequence_table_data import get_sequence_table_data_for_paper
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings
from retrobiocat_web.mongo.models.user_models import user_datastore


@bp.route('/sequence_curation_tab/<paper_id>', methods=['GET'])
@roles_required('contributor')
def sequence_curation_tab(paper_id):
    paper = paper_queries.paper_from_id(paper_id, get_related=True)
    user = user_datastore.get_user(current_user.id)

    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))

    if not check_permission.check_seq_curation_permission(user, paper):
        flash('No access to edit sequences for this entry', 'fail')
        return redirect(url_for("main_site.home"))

    enzyme_data = get_sequence_table_data_for_paper(paper)
    reviewed_dict = get_reviewed_dict(paper)
    access_dict = get_access_dict(paper, user)

    enzyme_types = all_enzyme_type_strings()

    seq_table_options = {'table_height': '60vh',
                         'show_header_filters': True,
                         'lock_enzyme_types': False}

    return render_template('curation_sequences/sequences_tab.html',
                           paper_id=paper.id,
                           paper_title=paper.title,
                           access_dict=access_dict,
                           reviewed_dict=reviewed_dict,
                           seq_data=enzyme_data,
                           enzyme_types=enzyme_types,
                           seq_table_options=seq_table_options)