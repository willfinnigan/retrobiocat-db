
def get_reviewed_dict(paper):
    reviewed_dict = {'paper_metadata_reviewed': paper_metadata_reviewed(paper),
                     'seq_review_ready': seqs_review_ready(paper),
                     'seq_reviewed': seqs_reviewed(paper),
                     'seq_reviewed_by': seqs_reviewed_by(paper),
                     'activity_review_ready': activity_review_ready(paper),
                     'activity_reviewed': activity_reviewed(paper),
                     'activity_reviewed_by': activity_reviewed_by(paper)}
    return reviewed_dict



def paper_metadata_reviewed(paper):
    if paper.metadata_reviewed == True:
        return 'checked'
    return ''


def seqs_review_ready(paper):
    if paper.seq_review_ready == True:
        return 'checked'
    return ''

def seqs_reviewed(paper):
    if paper.seq_reviewed == True:
        return 'checked'
    return ''

def seqs_reviewed_by(paper):
    if paper.seq_reviewed_by != None:
        return f"Reviewed by {paper.seq_reviewed_by.first_name} {paper.seq_reviewed_by.last_name}, {paper.seq_reviewed_by.affiliation}"
    return ''

def activity_review_ready(paper):
    if paper.activity_review_ready == True:
        return 'checked'
    return ''

def activity_reviewed(paper):
    if paper.activity_reviewed == True:
        return 'checked'
    return ''

def activity_reviewed_by(paper):
    if paper.activity_reviewed_by != None:
        return f"Reviewed by {paper.activity_reviewed_by.first_name} {paper.activity_reviewed_by.last_name}, {paper.activity_reviewed_by.affiliation}"
    return ''



