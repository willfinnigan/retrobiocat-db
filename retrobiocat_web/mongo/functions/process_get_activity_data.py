import numpy as np
import copy

COLUMNS = ['reaction',
           'enzyme_type',
           'enzyme_name',
           'short_citation',
           'html_doi',
           'cascade_num',
           'substrate_1_smiles',
           'substrate_2_smiles',
           'product_1_smiles',
           'temperature',
           'ph',
           'solvent',
           'other_conditions',
           'notes',
           'reaction_vol',
           'formulation',
           'biocat_conc',
           'kcat',
           'km',
           'mw',
           'substrate_1_conc',
           'substrate_2_conc',
           'specific_activity',
           'conversion',
           'conversion_time',
           'categorical',
           'binary',
           'selectivity',
           'paper',
           '_id']

MONGO_COLS = copy.copy(COLUMNS)
MONGO_COLS.remove('_id')
MONGO_COLS.append('id')


def process_activity_data(activity_data):
    for i, record in enumerate(activity_data):
        activity_data[i]['paper'] = str(activity_data[i]['paper'])
        if 'id' in activity_data[i]:
            activity_data[i]['_id'] = str(activity_data[i]['id'])
        else:
            activity_data[i]['_id'] = str(activity_data[i]['_id'])

        if 'added_by' in activity_data[i]:
            activity_data[i]['added_by'] = str(activity_data[i]['added_by'])

        # if 'conversion' in activity_data[i]:
        #     activity_data[i]['conversion'] = float(activity_data[i]['conversion'])
        # if 'conversion_time' in activity_data[i]:
        #     activity_data[i]['conversion_time'] = float(activity_data[i]['conversion_time'])
        # if 'specific_activity' in activity_data[i]:
        #     activity_data[i]['specific_activity'] = float(activity_data[i]['specific_activity'])

        for key in activity_data[i]:
            if activity_data[i][key] is True:
                activity_data[i][key] = "True"
            if activity_data[i][key] is False:
                activity_data[i][key] = "False"
            if activity_data[i][key] == np.nan:
                activity_data[i][key] = ""
            if type(activity_data[i][key]) == float:
                activity_data[i][key] = round(activity_data[i][key], 2)

    return activity_data
