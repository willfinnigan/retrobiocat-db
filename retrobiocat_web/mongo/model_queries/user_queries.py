from retrobiocat_web.mongo.models.user_models import User


def get_all_contributors():
    return User.contributors()


def get_user_from_id(user_id):
    if user_id == '':
        return None
    return User.objects(id=user_id).first()