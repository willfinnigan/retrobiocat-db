from retrobiocat_web.mongo.model_queries import sequence_queries, activity_queries
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_of_paper


def task_update_seq_papers_status(enzyme_name):
    seq = sequence_queries.seq_obj_from_name(enzyme_name, get_related=True)
    for paper in seq.papers:
        tag_paper_with_enzyme_types(paper)
        paper_progress_text, paper_progress = paper_metadata_status(paper)
        sequence_progress_text, sequence_progress = sequences_status(paper)
        activity_progress_text, activity_progress = activity_status(paper)
        status, status_colour = get_status(paper_progress, sequence_progress, activity_progress, paper)

        paper.status = status
        paper.save()


def tag_paper_with_enzyme_types(paper):
    all_types = all_enzyme_type_strings()
    for enz_type in all_types:
        if enz_type in paper.tags:
            paper.tags.remove(enz_type)

    seqs = seqs_of_paper(paper)
    for seq in seqs:
        if (seq.enzyme_type not in paper.tags) and (seq.enzyme_type.lower() != 'chemical'):
            paper.tags.append(seq.enzyme_type)

    paper.save()


def paper_metadata_status(paper):
    if paper.metadata_reviewed == True:
        progress = 100

    else:
        progress = 90
        required = []
        if paper.short_citation is None or paper.short_citation == '':
            progress -= 10
            required.append('Short citation')
        if paper.title is None or paper.title == '':
            progress -= 10
            required.append('Title')
        if paper.authors is None or paper.authors == [] or paper.authors == ['']:
            progress -= 10
            required.append('Authors')
        if paper.journal is None or paper.journal == '':
            progress -= 10
            required.append('Journal')
        if paper.date is None:
            progress -= 10
            required.append('Date')
        if paper.tags is None or paper.tags == [] or paper.tags == ['']:
            progress -= 10
            required.append('Tags')

    if progress == 100:
        progress_text = 'Complete'
    elif progress == 90:
        progress_text = 'Ready for review'
    else:
        progress_text = f'Required: {required}'

    return progress_text, str(progress)+"%"

def sequences_status(paper):
    if paper.seq_reviewed == True:
        progress = 100
        progress_text = 'Complete'
        return progress_text, str(progress) + '%'

    elif paper.seq_review_ready == True:
        progress = 90
        progress_text = 'Ready for review'
        return progress_text, str(progress) + '%'

    else:
        seqs = sequence_queries.seqs_of_paper(paper)

        if len(seqs) != 0:
            progress = 75
            progress_text = 'Sequences added, data entry ongoing..'

            for seq in seqs:
                if (seq.sequence is None or seq.sequence == '') and (seq.sequence_unavailable is False):
                    progress = 50
                    progress_text = 'Sequences which have been added require a protein sequence, or have protein sequence marked as unavailable'
        else:
            progress_text = 'Sequences required'
            progress = 5

        return progress_text, str(progress) + '%'

def activity_status(paper):
    if paper.activity_reviewed == True:
        progress = 100
        progress_text = 'Complete'
        return progress_text, str(progress) + '%'

    elif paper.activity_review_ready == True:
        progress = 90
        progress_text = 'Ready for review'
        return progress_text, str(progress) + '%'

    else:
        num_act = activity_queries.activity_in_paper(paper, count_only=True)
        if num_act != 0:
            progress = 50
            progress_text = "Activity data curation started, in progress"
        else:
            progress = 5
            progress_text = "Activity data required"

        return progress_text, str(progress) + '%'

def get_status(paper_progress, sequence_progress, activity_progress, paper):
    paper_progress = int(paper_progress[:-1])
    sequence_progress = int(sequence_progress[:-1])
    activity_progress = int(activity_progress[:-1])

    if paper_progress == 100 and sequence_progress == 100 and activity_progress == 100:
        status, colour = 'Complete', 'green'
    elif paper_progress <= 90 and sequence_progress == 100 and activity_progress == 100:
        status, colour = 'Complete but requires paper metadata review', 'green'
    elif sequence_progress == 90 and activity_progress == 90:
        status, colour = 'Awaiting review', 'green'

    else:
        if activity_progress == 90:
            status, colour = 'Activity data ready for review', 'orange'
        elif sequence_progress == 90:
            status, colour = 'Enzymes ready for review', 'orange'
        elif sequence_progress == 100 and activity_progress >= 5:
            status, colour = 'Sequences finished, activity data curation required', 'orange'
        elif activity_progress >= 10 or sequence_progress >= 10:
            status, colour = 'Data curation started', 'red'
        else:
            status, colour = 'Data required', 'red'

    if paper.has_issues == True:
        status, colour = 'Issues need to be resolved', 'red'

    return status, colour

def update_status(paper):
    paper_progress_text, paper_progress = paper_metadata_status(paper)
    sequence_progress_text, sequence_progress = sequences_status(paper)
    activity_progress_text, activity_progress = activity_status(paper)
    status, status_colour = get_status(paper_progress, sequence_progress, activity_progress, paper)

    paper.status = status
    paper.save()
