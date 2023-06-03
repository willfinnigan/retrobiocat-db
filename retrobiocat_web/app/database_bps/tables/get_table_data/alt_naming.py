from retrobiocat_web.mongo.model_queries import sequence_queries
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id


def get_alt_naming(args):
    enzyme_alt_names_to_use = {}
    if 'paper_id' in args:
        paper = paper_from_id(args['paper_id'])
        enzyme_alt_names_to_use = sequence_queries.get_alt_naming_for_paper(paper)
    return enzyme_alt_names_to_use