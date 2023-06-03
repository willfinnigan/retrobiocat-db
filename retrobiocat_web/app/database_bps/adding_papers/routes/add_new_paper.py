from flask import render_template, flash, redirect, url_for, session
from flask_security import roles_required, current_user

from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.list_parse import list_to_string
from retrobiocat_web.app.database_bps.adding_papers import bp
from retrobiocat_web.app.database_bps.adding_papers.forms import PaperInfo
from retrobiocat_web.analysis.uniprot_and_web import crossref_lookup, pubmed

from flask_wtf import FlaskForm
from wtforms import SubmitField, StringField

from retrobiocat_web.mongo.modal_updates import paper_CUD
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.models.user_models import User

class AddPaper(FlaskForm):
    doi = StringField("DOI")
    submit = SubmitField('Go')

@bp.route('/launch_add_paper', methods=['GET', 'POST'])
@roles_required('paper_adder')
def launch_add_paper():
    form = AddPaper()

    if form.validate_on_submit() == True:
        doi = form.data.get('doi', 'failed').replace(' ','')
        list_html_to_remove = ['https://doi.org/', 'http://doi.org/', 'http://dx.doi.org/']
        for to_remove in list_html_to_remove:
            if to_remove in doi:
                doi = doi.replace(to_remove, '')
        session['doi'] = doi
        return redirect(url_for("adding_papers.create_paper"))

    else:
        return render_template('add_paper_workflow/add_paper_start.html', form=form)

@bp.route("/create_paper/", methods=["GET", "POST"])
@roles_required('paper_adder')
def create_paper():
    form = PaperInfo()
    user = User.objects(id=current_user.id)[0]

    doi = session.get('doi', 'failed')
    doi = doi.replace(' ', '').lower()  # clean the doi to a standard format

    if doi == 'failed':
        flash("DOI required", 'error')
        return redirect(url_for("adding_papers.launch_add_paper"))

    # if the paper is already in the database, go to data entry page if allowed
    if paper_queries.check_if_paper_exists(doi) == True:
        paper = paper_queries.paper_from_doi(doi)

        if check_permission.check_paper_permission(user.id, paper):
            flash("Paper already in the database", 'success')
            return redirect(url_for("curation.curation_overview", paper_id=paper.id))

        elif paper.owner == None:
            flash('Paper already in database, self-assign the paper to add data', 'warning')
            return redirect(url_for("adding_papers.launch_add_paper"))

        elif paper.owner != user and not user.has_role('super_contributor'):
            flash("Paper already in the database and is assigned to another user", 'danger')
            return redirect(url_for("adding_papers.launch_add_paper"))
        else:
            flash("error", 'danger')
            return redirect(url_for("adding_papers.launch_add_paper"))

    # if paper not in the database, then lookup fields and create form
    else:
        title, authors_list, journal, date, cite_mini = crossref_lookup.get_metadata_from_crossref(doi)

        if cite_mini == '':
            flash('Failed to get paper from crossref', 'warning')
            title, authors_list, journal, date, cite_mini = pubmed.query_pubmed(doi)

        form.title.data = title
        form.authors.data = list_to_string(authors_list)
        form.journal.data = journal
        form.date.data = date
        form.short_cit.data = cite_mini
        form.doi.data = doi.replace(' ', '').lower()

        can_self_assign = check_permission.can_self_assign(user)
        can_make_high_priority = check_permission.check_has_make_high_priority_access(user)

        if cite_mini == '':
            flash("Paper not found in crossref, pubmed or db, please add manually", 'fail')
        else:
            flash("Paper found, please check information", 'success')

    # unless user clicks submit, validation will be false
    if form.validate_on_submit() == False:
        return render_template('add_paper_workflow/edit_paper_information.html',
                               form=form,
                               can_self_assign=can_self_assign,
                               can_make_high_priority=can_make_high_priority)

    # form is validated, so create new paper
    else:
        # create the new paper
        issues, new_paper = paper_CUD.create_new_paper(form.data, user=user, self_assign=form.self_assign.data)

        # if no issues, paper created successfully
        if len(issues) == 0:
            flash("Paper saved", 'success')
            session.pop('doi')
            if new_paper.owner == user:
                return redirect(url_for("curation.curation_overview", paper_id=new_paper.id))
            else:
                return redirect(url_for("adding_papers.launch_add_paper"))

        else:
            flash(issues, 'danger')
            return render_template('add_paper_workflow/edit_paper_information.html',
                                   form=form,
                                   can_self_assign=can_self_assign,
                                   can_make_high_priority=can_make_high_priority)


