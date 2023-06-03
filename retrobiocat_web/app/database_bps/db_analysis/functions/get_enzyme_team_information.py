from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Paper
from retrobiocat_web.mongo.models.user_models import User, Role
import mongoengine as db

def get_team(enzyme_type_obj):

    def get_info(user):
        user_string = f"{user.first_name} {user.last_name}"
        affiliation_string = f"{user.affiliation}"
        num_papers = Paper.objects(db.Q(tags=enzyme_type_obj.enzyme_type) & db.Q(owner=user)).count()
        return [user_string, affiliation_string, num_papers]

    team_info = {'enzyme_champions': [],
                 'team_members': []}

    champions = User.objects(enzyme_champion=enzyme_type_obj)
    team_members = User.objects(db.Q(enzyme_teams=enzyme_type_obj) & db.Q(enzyme_champion__ne=enzyme_type_obj))

    for user in team_members:
        team_info['team_members'].append(get_info(user))
    for user in champions:
        team_info['enzyme_champions'].append(get_info(user))

    return team_info


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    et = "CAR"
    et_obj = EnzymeType.objects(enzyme_type=et)[0]
    team_info = get_team(et_obj)

    print(team_info)

