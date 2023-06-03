from retrobiocat_web.mongo.model_queries.activity_queries import activity_in_paper
from retrobiocat_web.mongo.model_queries.paper_queries import num_papers_owned_by_user_which_need_data, \
    num_papers_owned_by_user_which_are_complete
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_of_paper
from retrobiocat_web.mongo.model_queries.user_queries import get_user_from_id


def check_paper_permission(user_id, paper):
    user = get_user_from_id(user_id)

    if paper is None:
        return False
    if user.has_role('super_contributor'):
        return True
    if paper.owner == user:
        return True
    if user.has_role('enzyme_champion'):
        champ_types = [e.enzyme_type for e in user.enzyme_champion]
        if any(i in champ_types for i in paper.tags):
            return True
    return False

def check_review_persmission(user_id, paper):
    user = get_user_from_id(user_id)
    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_champion'):
        champ_types = [e.enzyme_type for e in user.enzyme_champion]
        if any(i in champ_types for i in paper.tags):
            return True

    return False

def check_has_make_high_priority_access(user):
    if user.has_role('super_contributor'):
        return True

    if user.has_role('paper_added'):
        return True

    return False


def check_seq_permissions(user_id, seq):
    user = get_user_from_id(user_id)

    # no editing reviewed sequences
    if seq.reviewed == True:
        return False

    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_champion'):
        team_types = [e.enzyme_type for e in user.enzyme_champion]
        if seq.enzyme_type in team_types:
            return True

    if seq.owner == user:
        return True

    return False

def check_seq_partial_permissions(user_id, seq_obj):
    """Checks if user has permission for partial edits to a sequence"""

    user = get_user_from_id(user_id)

    if seq_obj.reviewed == True:
        return False

    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_teams'):
        team_types = [e.enzyme_type for e in user.enzyme_teams]
        if seq_obj.enzyme_type in team_types:
            return True

    if seq_obj.owner == user:
        return True

    return False



def check_seq_review(user_id, seq_obj):
    """If user has permission to review this sequence, return True"""

    user = get_user_from_id(user_id)
    if user.has_role('super_contributor'):
        return True
    if user.has_role('enzyme_champion'):
        champ_types = [e.enzyme_type for e in user.enzyme_champion]
        if seq_obj.enzyme_type in champ_types:
            return True

    return False

def check_team_permission(user, enzyme_type_obj):
    try:
        if user.has_role('super_contributor'):
            return True

        if user.has_role('enzyme_champion'):
            if enzyme_type_obj in user.enzyme_champion:
                return True
    except:
        return False

    return False

def can_self_assign(user):
    num_papers_need_data = num_papers_owned_by_user_which_need_data(user)
    num_papers_with_data = num_papers_owned_by_user_which_are_complete(user)

    if num_papers_need_data > num_papers_with_data:
        return False
    else:
        return True

def can_unassign(user, paper):
    """
    Can the user unassign this paper?
    Only if it contains no data or they are an enzyme champion
    """

    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_champion'):
        team_types = [e.enzyme_type for e in user.enzyme_champion]
        if any(item in paper.tags for item in team_types):
            return True

    if seqs_of_paper(paper, count_only=True) == 0 and activity_in_paper(paper, count_only=True) == 0:
        return False
    else:
        return True

def has_unreview_access(paper, user):
    """Determine if the user is allowed to unreview"""

    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_champion'):
        team_types = [e.enzyme_type for e in user.enzyme_champion]
        if any(item in paper.tags for item in team_types):
            return True

    if paper.owner == user:
        return True

    return False

def check_seq_curation_permission(user, paper):
    """Determines if user has access to curate sequences for this paper"""
    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_teams'):
        team_types = [e.enzyme_type for e in user.enzyme_teams]
        if any(item in paper.tags for item in team_types):
            return True

    if paper.owner == user:
        return True

    return False

def check_type_seq_curation_permission(user, enzyme_type):
    """Determines if user has access to curate sequences for this paper"""
    if user.has_role('super_contributor'):
        return True

    if user.has_role('enzyme_teams'):
        team_types = [e.enzyme_type for e in user.enzyme_teams]
        if enzyme_type in team_types:
            return True

    return False
