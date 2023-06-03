from flask import render_template, flash, redirect, url_for, request
from flask_security import roles_required, current_user
from retrobiocat_web.app.app import user_datastore
from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.app.database_bps.tables.get_table_data.get_papers_table_data import process_papers_to_table
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import enzyme_type_from_name
from retrobiocat_web.mongo.model_queries.paper_queries import papers_owned_by_user, papers_tagged_with_type, get_all_papers

@bp.route('/my_papers', methods=['GET', 'POST'])
@roles_required('contributor')
def my_papers():
    user = user_datastore.get_user(current_user.id)
    papers = papers_owned_by_user(user)
    papers_data = process_papers_to_table(papers)

    papers_button_columns = ['edit']
    if user.has_role('paper_adder'):
        papers_button_columns.append('delete')

    paper_table_options = {'table_height': '80vh'}

    return render_template('edit_papers/edit_papers.html',
                           paper_table_options=paper_table_options,
                           papers_data=papers_data,
                           title=f"Papers assigned to {user.first_name} {user.last_name}")

@bp.route('/edit_papers', methods=['GET', 'POST'])
@roles_required('super_contributor')
def edit_papers():
    papers = get_all_papers()
    papers_data = process_papers_to_table(papers)

    paper_table_options = {'table_height': '80vh'}

    return render_template('edit_papers/super_user_edit_papers.html',
                           paper_table_options=paper_table_options,
                           papers_data=papers_data,
                           title="Super contributor access to all papers")

@bp.route('/enz_champ_papers/<enzyme_type>', methods=['GET'])
@roles_required('enzyme_champion')
def enzyme_champion_papers(enzyme_type):
    user = user_datastore.get_user(current_user.id)
    enzyme_type_obj = enzyme_type_from_name(enzyme_type)
    if enzyme_type_obj not in user.enzyme_champion:
        flash('No access', 'danger')
        return redirect(url_for('main_site.home'))

    papers = papers_tagged_with_type(enzyme_type)
    papers_data = process_papers_to_table(papers)

    paper_table_options = {'table_height': '80vh'}

    return render_template('edit_papers/super_user_edit_papers.html',
                           papers_data=papers_data,
                           paper_table_options=paper_table_options,
                           title=f"Enzyme champion for {enzyme_type} papers")

