import json

from flask import flash, redirect, url_for, render_template, session, current_app
from flask_login import current_user
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.curation.functions.get_reviewed_dict import get_reviewed_dict
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.models.user_models import user_datastore

"""
OCR tab
1. Upload images
2. Process
3. New page showing extracted molecules -> buttons to load these.
"""


@bp.route('/molecules_osr_tab/<paper_id>', methods=['GET'])
@roles_required('contributor')
def molecules_osr_tab(paper_id):
    paper = paper_queries.paper_from_id(paper_id, get_related=True)
    user = user_datastore.get_user(current_user.id)

    if paper is None:
        flash('Paper has not been added yet, please add to the database first', 'fail')
        return redirect(url_for("adding_papers.launch_add_paper"))

    if not check_permission.check_seq_curation_permission(user, paper):
        flash('No access to edit sequences for this entry', 'fail')
        return redirect(url_for("main_site.home"))

    reviewed_dict = get_reviewed_dict(paper)

    save_id = f"{paper_id}_osr_smiles"

    save_id = f"{paper_id}_osr_smiles"
    if current_app.redis.exists(save_id):
        redis_content = json.loads(current_app.redis.get(save_id))
    else:
        redis_content = []

    return render_template('curation_molecules/molecules_osr/molecules_osr_tab.html',
                           paper_id=str(paper.id),
                           paper_title=paper.title,
                           redis_content=redis_content,
                           reviewed_dict=reviewed_dict)