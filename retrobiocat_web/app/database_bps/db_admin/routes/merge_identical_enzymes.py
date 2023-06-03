from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.app.database_bps.db_admin import bp
from flask import render_template, redirect, url_for
from retrobiocat_web.analysis.sequence_analysis.find_identical_sequences import find_identical_sequences
from retrobiocat_web.mongo.model_queries import enzyme_type_queries
from flask_security import current_user
from retrobiocat_web.app.app import user_datastore

@bp.route("/identical_enzymes/<enzyme_type>", methods=["GET"])
def identical_enzymes(enzyme_type):

    matches = find_identical_sequences(enzyme_type)

    return render_template('identical_enzymes.html', enzyme_type=enzyme_type, matches=matches)

@bp.route("/merge_identical_enzymes/<enzyme_type>", methods=["GET"])
def merge_identical_enzymes(enzyme_type):

    et_obj = enzyme_type_queries.enzyme_type_from_name(enzyme_type)
    user = user_datastore.get_user(current_user.id)

    if check_permission.check_team_permission(user, et_obj):
        matches = find_identical_sequences(enzyme_type)
        return render_template('merge_identical_enzymes.html', enzyme_type=enzyme_type, matches=matches)

    # redirect if no access
    return redirect(url_for('main_site.home'))