import numpy as np
from rq import get_current_job

from retrobiocat_web.analysis.drawing.activity_data_smiles_to_svg import smiles_to_svg
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.mongo.functions import process_get_activity_data
from retrobiocat_web.similarity.create_scorers import create_retrobiocat_scorer


def task_get_spec_data(enzyme_type, reaction, product, similarity_cutoff, num_choices):
    job = get_current_job()
    set_progress_bar(job, 40, 'started')

    scorer = create_retrobiocat_scorer('mongo')
    activity_df = scorer.score_data(product, num_choices, similarity_cutoff, reaction_name=reaction, enzyme_types=[enzyme_type])

    set_progress_bar(job, 80, 'data retrieved')

    if activity_df is None:
        print('Activity df is none')
        return []

    if len(activity_df.index) == 0:
        print('Len activity df index is 0')
        return []

    processed_activity_df = activity_df[process_get_activity_data.COLUMNS]
    processed_activity_df['similarity'] = activity_df['similarity']
    processed_activity_df = processed_activity_df.round(2)
    processed_activity_df.replace(np.nan, '', inplace=True)
    processed_activity_df.replace(True, 'True', inplace=True)
    processed_activity_df.replace(False, 'False', inplace=True)

    activity_data = processed_activity_df.to_dict(orient='records')
    activity_data = process_get_activity_data.process_activity_data(activity_data)
    activity_data = smiles_to_svg(activity_data)
    set_progress_bar(job, 90, 'data processed')

    return activity_data
