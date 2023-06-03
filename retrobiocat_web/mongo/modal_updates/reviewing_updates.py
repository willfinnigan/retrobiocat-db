
from retrobiocat_web.mongo.model_queries import sequence_queries, enzyme_type_queries, activity_queries, paper_queries
from retrobiocat_web.mongo.modal_updates import sequence_CUD, bioinformatics_CUD, activity_CUD, paper_CUD


def review_paper_metadata(paper, reviewed):
    try:
        paper.metadata_reviewed = reviewed
        paper.save()
        return []
    except Exception as e:
        return [f"Error updating paper metadata reviewed status for {paper.short_citation}", str(e)]


def sequences_ready_to_review(paper, ready):
    if ready == False and paper.seq_reviewed == True:
        return ['Sequences have already been reviewed']

    try:
        paper.seq_review_ready = ready
        paper.save()
        return []
    except Exception as e:
        return [f"Error updating sequences ready to review status for {paper.short_citation}", str(e)]


def review_paper_sequences(paper, user, reviewed):
    """Update a papers sequences reviewed fields, along with attached sequences as appropriate"""

    sequences = sequence_queries.seqs_of_paper(paper, get_related=True)
    issues = []

    if reviewed == True:
        print('reviewing..')

        # first check all sequences are ok to review
        names = sequence_queries.seq_names_in_paper_not_ok_to_review(paper)
        for name in names:
            issues.append(f'{name} needs a protein sequence (or be marked as unavailable)')

        # if no issues, them update the sequences (check off because bulk checked above)
        if len(issues) == 0:
            print('updating seqs..')
            for seq in sequences:
                issues += sequence_CUD.update_reviewed(seq, True, check=False, update_bio=False)

        # mark bioinformatics for updating for enzyme types in paper
        if len(issues) == 0:
            enzyme_type_names = list(enzyme_type_queries.enzyme_type_names_in_paper(paper))
            for enzyme_type in enzyme_type_names:
                bioinformatics_CUD.enzyme_type_needs_blasting(enzyme_type)
                bioinformatics_CUD.SSN_needs_updating(enzyme_type)

    # if reviewed is False, un-review any seqs with no other reviewed papers
    elif reviewed == False:
        print('unreviewing..')
        for seq in sequences:
            papers = paper_queries.get_papers_where_seq_is_also_reviewed(seq)
            if len(papers) == 1:
                issues += sequence_CUD.update_reviewed(seq, False)

    else:
        return ['Error - reviewed is not a bool']

    if len(issues) != 0:
        return issues

    print('updating paper..')
    paper.seq_reviewed = reviewed
    paper.seq_review_ready = reviewed
    if reviewed == True:
        paper.seq_reviewed_by = user
    else:
        paper.seq_reviewed_by = None
    paper.save()

    if len(issues) == 0 and reviewed is True:
        paper_CUD.tag_paper_with_enzyme_types(paper)

    if paper.seq_reviewed != reviewed:
        return ['Error - could not update sequence review status']

    return []


def activity_ready_to_review(paper, ready):
    if ready == False and paper.activity_reviewed == True:
        return ['Activity has already been reviewed']

    try:
        paper.activity_review_ready = ready
        paper.save()
        return []
    except Exception as e:
        return [f"Error updating activity ready to review status for {paper.short_citation}", str(e)]


def review_paper_activity(paper, user, reviewed):
    """Update a papers sequences reviewed fields, along with attached sequences as appropriate"""

    activities = activity_queries.activity_in_paper(paper)

    issues = []

    # update all the activity rows to reviewed
    for act_obj in activities:
        issues += activity_CUD.update_reviewed(act_obj, reviewed)
        if len(issues) != 0:
            break

    # update the paper activity_reviewed status
    if len(issues) == 0:
        print('updating paper..')
        try:
            paper.activity_reviewed = reviewed
            paper.activity_review_ready = reviewed
            if reviewed == False:
                paper.activity_reviewed_by = None
            else:
                paper.activity_reviewed_by = user
            paper.save()
        except Exception as e:
            return [f"Error updating paper sequence reviewed status for {paper.short_citation}", str(e)]

    return issues
