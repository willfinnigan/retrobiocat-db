from flask import render_template, redirect, url_for

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.mongo.models.user_models import User, Role
from retrobiocat_web.mongo.model_queries.activity_queries import num_activity_owned_by_user
from retrobiocat_web.mongo.model_queries.sequence_queries import num_sequences_owned_by_user
from retrobiocat_web.mongo.model_queries.paper_queries import num_papers_owned_by_user

@bp.route("/contributors", methods=["GET"])
def contributors():
    contributor_role = Role.objects(name='contributor')[0]
    contributors = User.objects(roles=contributor_role)

    papers_dict = {}
    sequence_dict = {}
    activity_dict = {}
    for user in contributors:
        username = f"{user.first_name} {user.last_name}, {user.affiliation}"
        papers_dict[username] = num_papers_owned_by_user(user)
        sequence_dict[username] = num_sequences_owned_by_user(user)
        activity_dict[username] = num_activity_owned_by_user(user)

    papers_dict = {k: v for k, v in sorted(papers_dict.items(), key=lambda item: item[1], reverse=True)}
    papers_dict = {k: v for k, v in papers_dict.items() if v != 0}
    sequence_dict = {k: v for k, v in sorted(sequence_dict.items(), key=lambda item: item[1], reverse=True)}
    sequence_dict = {k: v for k, v in sequence_dict.items() if v != 0}
    activity_dict = {k: v for k, v in sorted(activity_dict.items(), key=lambda item: item[1], reverse=True)}
    activity_dict = {k: v for k, v in activity_dict.items() if v != 0}

    return render_template('leaderboard.html',
                           top_papers=papers_dict,
                           top_sequences=sequence_dict,
                           top_activity=activity_dict)


@bp.route("/leaderboard", methods=["GET"])
def leaderboard():
    return redirect(url_for('curation.contributors'))
