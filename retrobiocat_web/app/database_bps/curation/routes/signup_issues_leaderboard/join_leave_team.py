from flask_login import current_user
from flask_security import roles_required
from flask import redirect, url_for, flash

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import enzyme_type_from_name
from retrobiocat_web.mongo.models.user_models import user_datastore


@bp.route('/join_team/<enzyme_type>', methods=['GET', 'POST'])
@roles_required('contributor')
def join_team(enzyme_type):
    user = user_datastore.get_user(current_user.id)
    enzyme_team_role = user_datastore.find_role('enzyme_teams')
    if not user.has_role(enzyme_team_role):
        user_datastore.add_role_to_user(user, enzyme_team_role)

    enzyme_type_object = enzyme_type_from_name(enzyme_type)
    if enzyme_type_object not in user.enzyme_teams:
        user.enzyme_teams.append(enzyme_type_object)
        user.save()

    flash(f'You joined the {enzyme_type} team', 'success')
    return redirect(url_for('db_analysis.enzyme_homepage', enzyme_type=enzyme_type))

@bp.route('/leave_team/<enzyme_type>', methods=['GET', 'POST'])
@roles_required('contributor')
def leave_team(enzyme_type):
    user = user_datastore.get_user(current_user.id)
    enzyme_team_role = user_datastore.find_role('enzyme_teams')
    enzyme_champ_role = user_datastore.find_role('enzyme_champion')
    enzyme_type_object = enzyme_type_from_name(enzyme_type)

    if enzyme_type_object in user.enzyme_teams:
        user.enzyme_teams.remove(enzyme_type_object)
        user.save()

    if enzyme_type_object in user.enzyme_champion:
        user.enzyme_enzyme_champion.remove(enzyme_type_object)
        user.save()

    if len(user.enzyme_teams) == 0 and user.has_role(enzyme_team_role):
        user_datastore.remove_role_from_user(user, enzyme_team_role)

    if len(user.enzyme_champion) == 0 and user.has_role(enzyme_champ_role):
        user_datastore.remove_role_from_user(user, enzyme_champ_role)

    flash(f'You left the {enzyme_type} team', 'success')
    return redirect(url_for('db_analysis.enzyme_homepage', enzyme_type=enzyme_type))