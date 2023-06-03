import time
from typing import List, Optional

from retrobiocat_web.mongo.models.retro3_models import Project, SearchResult

def get_project_for_user(user) -> List[Project]:
    return Project.objects(owner=user)

def get_project(project_id) -> Project:
    return Project.objects(id=project_id).first().select_related()

def get_search(search_id) -> SearchResult:
    return SearchResult.objects(id=search_id).first()

def get_pathway_from_search_obj(search_id, pathway_set, pathway_index):

    if pathway_set == 'solved':
        try:
            return SearchResult.objects(id=search_id).only('solved_pathways').fields(slice__solved_pathways=[pathway_index, 1]).first().solved_pathways[0]['reaction_ids']
        except Exception as e:
            print(str(e))
            return None
    else:
        try:
            return SearchResult.objects(id=search_id).only('all_pathways').fields(slice__all_pathways=[pathway_index, 1]).first().all_pathways[0]['reaction_ids']
        except Exception as e:
            print(str(e))
            return None

def get_recent_search_started(project_obj: Project) -> Optional[SearchResult]:
    searches = project_obj.search_results
    if len(searches) == 0:
        return None
    order_by_data = sorted(searches, key=lambda x: x.datetime, reverse=True)
    return order_by_data[0]

def get_recent_search_completed(project_obj: Project) -> Optional[SearchResult]:
    searches = project_obj.search_results
    searches = [s for s in searches if hasattr(s, 'datetime_completed')]
    complete_searches = [s for s in searches if s.datetime_completed is not None]
    if len(complete_searches) == 0:
        return None
    order_by_data = sorted(complete_searches, key=lambda x: x.datetime_completed, reverse=True)
    return order_by_data[0]

def get_pathway_hash(reaction_ids):
    return str(hash(str(sorted(reaction_ids))))

def is_pathway_saved_in_project(project, reaction_ids):
    reaction_ids = sorted(reaction_ids)
    if reaction_ids in project.saved_pathways:
        print(f"Pathway is saved")
        return True
    return False


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    t0 = time.time()
    pathway = get_pathway('634d2d962c29029676750259', 'all', 200)
    print(time.time()-t0)
    print(pathway)

    t0 = time.time()
    search = get_search('634d2d962c29029676750259')
    pathway = search.all_pathways[200]['reaction_ids']
    print(time.time() - t0)
    print(pathway)