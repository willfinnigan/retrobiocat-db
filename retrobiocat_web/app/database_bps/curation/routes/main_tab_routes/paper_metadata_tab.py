from flask import flash, redirect, url_for, render_template
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_access_dict import get_access_dict
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.app.database_bps.curation.functions.list_parse import list_to_string
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id
from retrobiocat_web.mongo.models.user_models import user_datastore


def get_paper_data(paper, user):
    self_assigned = ''
    other_user = ''

    if check_permission.has_unreview_access(paper, user) == False:
        has_no_unreview_access = 'disabled'
    else:
        has_no_unreview_access = ''

    if paper.owner == user:
        self_assigned = 'checked'
    elif paper.owner is not None and paper.owner != '':
        other_user = 'disabled'

    paper_owner_name = 'None'
    if paper.owner is not None:
        paper_owner_name = f"{paper.owner.first_name} {paper.owner.last_name}, {paper.owner.affiliation}"

    importance = ''
    if paper.high_importance == True:
        importance = 'checked'

    paper_dict = {'short_cit': paper.short_citation,
                  'doi': paper.doi,
                  'date': paper.date,
                  'title': paper.title,
                  'journal': paper.journal,
                  'authors': list_to_string(paper.authors),
                  'tags': list_to_string(paper.tags),
                  'self_assigned': self_assigned,
                  'has_no_unreview_access': has_no_unreview_access,
                  'disable_for_other_user': other_user,
                  'id': paper.id,
                  'importance_checked': importance,
                  'paper_owner_name': paper_owner_name}

    return paper_dict





@bp.route('/edit_paper_metadata/<paper_id>', methods=['GET'])
@roles_required('contributor')
def edit_paper_metadata(paper_id):
    user = user_datastore.get_user(current_user.id)
    paper = paper_from_id(paper_id, get_related=True)

    # check the paper exists and the user has access to view the curation page
    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))
    if not check_permission.check_paper_permission(current_user.id, paper):
        flash('No access to edit this entry', 'fail')
        return redirect(url_for("main_site.home"))

    # dicts with information for rendering the edit paper metadata page
    paper_dict = get_paper_data(paper, user)
    access_dict = get_access_dict(paper, user)
    reviewed_dict = get_reviewed_dict(paper)

    return render_template('curation_paper_metadata/paper_metadata.html',
                           paper_id=paper_id,
                           paper_title=paper.title,
                           paper=paper_dict,
                           access_dict=access_dict,
                           reviewed_dict=reviewed_dict)