"""
Curation starts on the status page, with a tab showing the status, inc. comments and issues checkbox,
and a second tab with links to data curation pages.
The admin tab is also here if the user has this access.
"""
from flask import render_template, flash, redirect, url_for
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import paper_status, check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.app.main_site.functions.comment_functions import get_comments_for_a_paper
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id
from retrobiocat_web.mongo.models.user_models import user_datastore



def get_and_update_status(paper):
    """Get the status of the paper"""
    paper_progress_text, paper_progress = paper_status.paper_metadata_status(paper)
    sequence_progress_text, sequence_progress = paper_status.sequences_status(paper)
    activity_progress_text, activity_progress = paper_status.activity_status(paper)
    status, status_colour = paper_status.get_status(paper_progress, sequence_progress, activity_progress, paper)

    paper.status = status
    paper.save()

    status_dict = {'review_checked': " ",
                   'review_disabled': 'disabled',
                   'reviewed_by': " ",
                   'issues_checked': " ",
                   'issues_disabled': 'disabled',
                   'importance_checked': " ",
                   'importance_disabled': 'disabled',
                   'no_relevant_activity_checked': " ",
                   'no_relevant_activity_disabled': "disabled",
                   'status': status,
                   'status_colour': status_colour,
                   'paper_progress': paper_progress,
                   'paper_progress_text': paper_progress_text,
                   'sequences_progress': sequence_progress,
                   'sequences_progress_text': sequence_progress_text,
                   'activity_progress': activity_progress,
                   'activity_progress_text': activity_progress_text}


    if paper.reviewed == False:
        status_dict['no_relevant_activity_disabled'] = " "

    if check_permission.check_review_persmission(current_user.id, paper):
        status_dict['review_disabled'] = " "
        status_dict['issues_disabled'] = " "
        status_dict['importance_hidden'] = " "
        status_dict['no_relevant_activity_disabled'] = " "

    if paper.reviewed_by != None:
        rb = paper.reviewed_by
        status_dict['reviewed_by'] = f"{rb.first_name} {rb.last_name}, {rb.affiliation}"

    if paper.reviewed == True:
        status_dict['review_checked'] = 'checked'
    if paper.has_issues == True:
        status_dict['issues_checked'] = 'checked'
    if paper.high_importance == True:
        status_dict['importance_checked'] = 'checked'
    if paper.no_relevant_activity == True:
        status_dict['no_relevant_activity_checked'] = 'checked'

    return status_dict

def get_ownership_dict(paper, user):
    """Gets a dict of information about the ownership of the paper for rendering"""
    paper_ownership = {'owner_name': 'None',
                       'self_assign_checked': '',
                       'disable_self_assign': ''}
    if paper.owner is not None:
        paper_ownership['owner_name'] = f"{paper.owner.first_name} {paper.owner.last_name}, {paper.owner.affiliation}"
    if paper.owner == user:
        paper_ownership['self_assign_checked'] = 'checked'
    elif paper.owner is not None and paper.owner != '':
        paper_ownership['disable_self_assign'] = 'disabled'

    return paper_ownership


@bp.route('/curation_overview/<paper_id>', methods=['GET'])
@roles_required('contributor')
def curation_overview(paper_id):
    user = user_datastore.get_user(current_user.id)
    paper = paper_from_id(paper_id, get_related=True)

    # check the paper exists and the user has access to view the curation page
    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))
    if not check_permission.check_paper_permission(current_user.id, paper):
        flash('No access to edit this entry', 'fail')
        return redirect(url_for("main_site.home"))

    # various dicts of information for rendering the curation launch page page
    reviewed_dict = get_reviewed_dict(paper)
    comments = get_comments_for_a_paper(paper, user)
    status_dict = get_and_update_status(paper)
    paper_ownership = get_ownership_dict(paper, user)



    return render_template('curation_overview/curation_overview.html',
                           paper_id=paper.id,
                           paper_title=paper.title,
                           reviewed_dict=reviewed_dict,
                           paper_ownership=paper_ownership,
                           status=status_dict,
                           comments=comments,
                           )