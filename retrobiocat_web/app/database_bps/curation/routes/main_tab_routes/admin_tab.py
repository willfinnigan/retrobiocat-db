from flask import render_template
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id
from retrobiocat_web.mongo.model_queries.user_queries import get_all_contributors
from retrobiocat_web.mongo.models.user_models import user_datastore


def get_admin_dict(paper):
    """Gets a dict of information for admin tab"""
    admin_dict = {}

    all_contributors = get_all_contributors()

    contributor_choices = [(f"{c.first_name} {c.last_name}", str(c.id)) for c in all_contributors]
    admin_dict['contributors'] = [('Unassigned', '')] + contributor_choices

    if paper.owner != None:
        admin_dict['owner'] = str(paper.owner.id)
    else:
        admin_dict['owner'] = ''

    return admin_dict


@bp.route('/curation_admin/<paper_id>', methods=['GET'])
@roles_required('admin')
def curation_admin(paper_id):
    user = user_datastore.get_user(current_user.id)
    paper = paper_from_id(paper_id, get_related=True)
    admin_dict = get_admin_dict(paper)
    reviewed_dict = get_reviewed_dict(paper)

    return render_template('curation_admin/curation_admin.html',
                           reviewed_dict=reviewed_dict,
                           paper_id=paper.id,
                           paper_title=paper.title,
                           admin_dict=admin_dict)
