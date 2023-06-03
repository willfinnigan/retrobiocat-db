from retrobiocat_web.mongo.models.biocatdb_models import Paper
from retrobiocat_web.mongo.model_queries import sequence_queries, enzyme_type_queries, activity_queries, paper_queries


def create_new_paper(data_dict, user=None, self_assign=False):

    # clean the doi to a standard format
    doi = data_dict['doi'].replace(' ', '').lower()

    # check if paper is already in the database
    if paper_queries.check_if_paper_exists(doi) == True:
        return ['Paper with this doi already exists'], None

    try:
        # if self assigning, make the current user the owner
        if self_assign != False and user is not None:
            owner = user
        else:
            owner = None

        # create the new paper
        new_paper = Paper(doi=doi,
                          short_citation=data_dict['short_cit'],
                          title=data_dict['title'],
                          html='https://doi.org/' + data_dict['doi'],
                          journal=data_dict['journal'],
                          date=data_dict['date'],
                          tags=data_dict['tags'].split(', '),
                          authors=data_dict['authors'].split(', '),
                          high_importance=data_dict['high_importance'],
                          owner=owner,
                          added_by=user,
                          status='Data required')
        new_paper.save()

        print(f'Added new paper entry with doi {doi}')

        return [], new_paper

    except Exception as e:
        return ['Error while creating new paper', str(e)], None

def new_owner(paper, user):
    """ Update the owner of a sequence"""

    try:
        paper.owner = user
        paper.save()
        return []
    except Exception as e:
        return [f"Couldn't update user for sequence: {paper.short_citation}", str(e)]

def tag_paper_with_enzyme_types(paper):
    """Update the tags for a paper to include enzyme types for attached sequences"""

    try:
        # this will remove all enzyme_type tags to then add back only those that are present
        all_types = enzyme_type_queries.all_enzyme_type_strings()
        for enz_type in all_types:
            if enz_type in paper.tags:
                paper.tags.remove(enz_type)

        # adds enzyme_type tags for enzymes that are in the paper
        seqs = sequence_queries.seqs_of_paper(paper)
        for seq in seqs:
            if (seq.enzyme_type not in paper.tags) and (seq.enzyme_type.lower() != 'chemical'):
                paper.tags.append(seq.enzyme_type)

        paper.save()
        return []
    except Exception as e:
        return [f'Error updating enzyme type tags for {paper.short_citation}', str(e)]

def update_paper(paper, data_dict, user=None):
    if paper.metadata_reviewed == True:
        return ["Cannot update paper while the metadata is marked as 'Reviewed'"]

    issues = []
    try:
        paper.short_citation = data_dict['short_cit']
        paper.doi = data_dict['doi'].replace(' ', '')
        paper.html = 'https://doi.org/' + data_dict['doi']
        if data_dict['date'] != "":
            paper.date = data_dict['date']
        paper.title = data_dict['title']
        paper.journal = data_dict['journal']
        paper.authors = data_dict['authors'].split(', ')
        paper.tags = data_dict['tags'].split(', ')

        if user not in paper.edits_by and user is not None:
            paper.edits_by.append(user)

        paper.save()

        if paper.seq_reviewed:
            issues += tag_paper_with_enzyme_types(paper)
        return issues
    except Exception as e:
        return [f'Error while updating paper {paper.short_citation}', str(e)]

def delete_paper(paper):

    issues = []

    # Do not delete papers with sequences attached
    if sequence_queries.seqs_of_paper(paper, count_only=True) != 0:
        issues.append('Please remove any sequences from paper before deleting')

    # Do not delete papers with activity data attached
    if activity_queries.activity_in_paper(paper, count_only=True) != 0:
        issues.append('Please remove any activity data from paper before deleting')

    # If no issues, delete the paper
    if len(issues) == 0:
        try:
            paper.delete()
        except Exception as e:
            issues += [f'Error while trying to delete paper {paper.short_citation}', str(e)]

    return issues

def update_paper_issues_status(paper, issue_status):
    try:
        paper.has_issues = issue_status
        paper.save()
        return []
    except Exception as e:
        return [f"Error updating issue status for {paper.short_citation}", str(e)]

def update_paper_importance(paper, important):
    try:
        paper.high_importance = important
        paper.save()
        return []
    except Exception as e:
        return [f"Error updating high importance flag for {paper.short_citation}", str(e)]

