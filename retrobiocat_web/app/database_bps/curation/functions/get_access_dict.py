from retrobiocat_web.app.database_bps.curation.functions import check_permission


def get_access_dict(paper, user):
    access_dict = {'review_disabled': review_disabled(paper, user),
                   'unreview_disabled': unreview_disabled(paper, user),
                   'importance_disabled': importance_checkbox_disabled(paper, user)}
    return access_dict

def unreview_disabled(paper, user):
    if check_permission.has_unreview_access(paper, user) == False:
        return 'disabled'
    return ''

def review_disabled(paper, user):
    if check_permission.check_review_persmission(user.id, paper):
        return ''
    return 'disabled'

def importance_checkbox_disabled(paper, user):
    if check_permission.check_review_persmission(user.id, paper):
        return ''
    return 'disabled'