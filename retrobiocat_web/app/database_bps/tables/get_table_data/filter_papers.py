from retrobiocat_web.mongo.model_queries.activity_queries import activities_of_reaction
from retrobiocat_web.mongo.model_queries.sequence_queries import seq_obj_from_name


def filter_papers_by_enzyme_name(papers, enzyme_name):
    seq = seq_obj_from_name(enzyme_name, get_related=True)
    seq_papers = [paper for paper in seq.papers]
    filtered_papers = [paper for paper in papers if paper in seq_papers]
    return filtered_papers

def filter_papers_by_reaction(papers, reaction):
    activities = activities_of_reaction(reaction, select_related=True)
    act_papers = [act.paper for act in activities]
    act_papers = list(set(act_papers))
    filtered_papers = [paper for paper in papers if paper in act_papers]
    return filtered_papers