from retrobiocat_web.mongo.model_queries.activity_queries import activity_in_paper
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_of_paper


def create_short_paper_dict_for_table(paper):
    """ Returns a dictionary in the format the paper tabulator table """

    paper_dict = {'_id': str(paper.id),
                  'short_citation': paper.short_citation,
                  'title': paper.title,
                  'doi': paper.doi}

    return paper_dict

def create_paper_dict_for_table(paper, include_counts=False):
    """ Returns a dictionary in the format the paper tabulator table """

    paper_dict = {'_id': str(paper.id),
                  'short_citation': paper.short_citation,
                  'title': paper.title,
                  'doi': paper.doi,
                  'status': paper.status,
                  'tags': paper.tags,
                  'owner': paper.owner_text()}

    if include_counts:
        paper_dict.update({'num_activity': activity_in_paper(paper, count_only=True),
                           'num_seqs': seqs_of_paper(paper, count_only=True)})

    return paper_dict

def process_papers_to_table(list_paper_objs, include_counts=False, short=False):
    """ Create data for tabulator by creating a list of dicts """

    table_data = []

    for paper in list_paper_objs:
        if short:
            paper_dict = create_short_paper_dict_for_table(paper)
        else:
            paper_dict = create_paper_dict_for_table(paper, include_counts=include_counts)
        table_data.append(paper_dict)

    return table_data


