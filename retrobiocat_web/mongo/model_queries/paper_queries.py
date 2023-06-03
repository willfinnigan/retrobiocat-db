from retrobiocat_web.mongo.models.biocatdb_models import Paper
import mongoengine as db

def get_all_papers():
    return Paper.objects()

def get_papers_that_need_data(enzyme_type=None):
    e_type_query = db.Q()
    if enzyme_type is not None:
        e_type_query = db.Q(tags=enzyme_type)

    q_no_user = db.Q(owner=None)
    q_no_data = db.Q(status__nin=['Awaiting review', 'Complete'])

    return Paper.objects(q_no_user & q_no_data & e_type_query).order_by('-status').select_related()

def get_high_importance_papers_with_no_owner(tag=None):
    hi_q = db.Q(high_importance=True)
    assigned_q = db.Q(owner=None)

    tag_q = db.Q()
    if tag is not None:
        tag_q = db.Q(tags=tag)

    return Paper.objects(hi_q & assigned_q & tag_q).order_by('-status').select_related()


def paper_from_id(paper_id, get_related=False):
    """ Get paper from its id """
    if get_related == False:
        paper = Paper.objects(id=paper_id).first()
    else:
        paper = Paper.objects(id=paper_id).first().select_related()
    return paper

def get_multiple_papers_from_ids(list_paper_ids, get_related=False):
    """ Get list of papers from a list of ids"""
    if get_related == False:
        papers = Paper.objects(id__in=list_paper_ids)
    else:
        papers = Paper.objects(id__in=list_paper_ids).select_related()
    return papers


def check_if_paper_exists(doi):
    """Checks if there is a paper with this doi already in the database"""

    num_papers = Paper.objects(doi__iexact=doi).count()
    if num_papers == 0:
        return False

    return True

def paper_from_doi(doi):
    """Return paper for a given doi"""

    return Paper.objects(doi__iexact=doi).first()

def get_papers_where_seq_is_also_reviewed(seq_with_related):
    """Return papers where sequence has been reviewed"""
    papers = []
    for paper in seq_with_related.papers:
        if paper.seq_reviewed == True:
            papers.append(paper)
    return papers

def papers_owned_by_user(user):
    return Paper.objects(owner=user).order_by('-status')

def papers_tagged_with_type(enzyme_type_abbrev):
    return Paper.objects(tags=enzyme_type_abbrev).order_by('-status')

def num_papers_owned_by_user(user):
    return Paper.objects(owner=user).count()

def num_papers_owned_by_user_which_need_data(user):
    q_user = db.Q(owner=user)
    q_no_data = db.Q(status__nin=['Complete - Awaiting review', 'Complete'])
    return Paper.objects(q_user & q_no_data).count()

def num_papers_owned_by_user_which_are_complete(user):
    q_user = db.Q(owner=user)
    q_has_data = db.Q(status__in=['Complete - Awaiting review', 'Complete'])
    return Paper.objects(q_user & q_has_data).count()



