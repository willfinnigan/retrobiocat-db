from retrobiocat_web.app.database_bps.curation.functions.free_text_standardisation import standardise_forumation
from retrobiocat_web.mongo.functions.save_activity_data import activity_functions
from retrobiocat_web.mongo.models.biocatdb_models import Activity
from retrobiocat_web.mongo.model_queries import sequence_queries
import mongoengine as db


def update_or_create_activity(data_dict, paper, user=None):

    if paper.activity_reviewed == True:
        return ['Can not update any rows for the activity data of a paper if its been marked as activity reviewed']

    issues = []

    # if _id in database, update that entry
    activity = None
    if data_dict.get('_id', '') != '':
        query = Activity.objects(id=data_dict['_id'])
        if len(query) != 0:
            activity = query[0]

    # Otherwise make a new entry
    if activity is None:
        activity = Activity()
        activity.paper = paper  # link paper to activity
        activity.added_by = user  # new data is added by the current user
        activity.save()

    # sequence must be in database
    seq_obj = None
    if data_dict.get('enzyme_name', None) == None:
        issues.append('enzyme_name not present in data_dict')
    else:
        seq_obj = sequence_queries.seq_obj_from_name(data_dict.get('enzyme_name'))
        if seq_obj is None:
            issues.append(f"No sequence in database with name {data_dict.get('enzyme_name')}")

    # data_dict must have a reaction
    if 'reaction' not in data_dict:
        issues.append('No reaction is specified')


    # add or update activity data
    if len(issues) == 0:
        try:
            activity.enzyme_type = seq_obj.enzyme_type
            activity.enzyme_name = seq_obj.enzyme_name
            activity.reaction = data_dict['reaction']
            activity.short_citation = paper.short_citation
            activity.html_doi = paper.html
            activity.paper = paper
            activity.substrate_1_smiles = str(data_dict.get('substrate_1_smiles', ''))
            activity.substrate_2_smiles = str(data_dict.get('substrate_2_smiles', ''))
            activity.product_1_smiles = str(data_dict.get('product_1_smiles', ''))
            activity.temperature = activity_functions.check_is_float(data_dict.get('temperature', None))
            activity.ph = activity_functions.check_is_float(data_dict.get('ph', None))

            activity.solvent = str(data_dict.get('solvent', ''))
            activity.other_conditions = str(data_dict.get('other_conditions', ''))
            activity.notes = str(data_dict.get('notes', ''))
            activity.reaction_vol = str(data_dict.get('reaction_vol', ''))
            activity.formulation = standardise_forumation(str(data_dict.get('formulation', None)))
            activity.biocat_conc = str(data_dict.get('biocat_conc', ''))
            activity.kcat = activity_functions.check_is_float(data_dict.get('kcat', None))
            activity.km = activity_functions.check_is_float(data_dict.get('km', None))
            activity.mw = activity_functions.check_is_float(data_dict.get('mw', None))

            activity.substrate_1_conc = str(data_dict.get('substrate_1_conc', ''))
            activity.substrate_2_conc = str(data_dict.get('substrate_2_conc', ''))

            activity.specific_activity = activity_functions.check_is_float(data_dict.get('specific_activity', None))
            activity.conversion = activity_functions.check_is_float(data_dict.get('conversion', None))
            activity.conversion_time = activity_functions.check_is_float(data_dict.get('conversion_time', None))
            activity.categorical = str(data_dict.get('categorical', None))
            if data_dict.get('binary', None) == 1:
                binary = True
                if activity.categorical == 'None':  # fix issue where None is being converted to "None"
                    activity.categorical = None
            elif data_dict.get('binary', None) == 0:
                binary = False
                activity.categorical = "None"
            else:
                binary = None
            activity.binary = binary
            activity.selectivity = str(data_dict.get('selectivity', ''))
            activity.auto_generated = False

            # add current user to edits_by
            if user not in paper.edits_by and user is not None:
                paper.edits_by.append(user)
                paper.save()
                activity.edits_by.append(user)

            activity.save()
            return []
        except Exception as e:
            issues += ['Error while trying to add or update activity data', str(e)]
    return issues

def delete_other_data(ids_to_keep, paper):
    """ Delete any activity from a paper which isn't in the list of id's to keep"""

    if paper.activity_reviewed == True:
        return ['Can not update any rows for the activity data of a paper if its been marked as activity reviewed']

    try:
        p_Q = db.Q(paper=paper)
        data_id_Q = db.Q(id__nin=ids_to_keep)
        Activity.objects(p_Q & data_id_Q).delete()
        return []
    except Exception as e:
        return [f'Error deleting activity data in {paper.short_citation}', str(e)]

def update_reviewed(act_obj, reviewed):
    try:
        act_obj.reviewed = reviewed
        act_obj.save()
        return []
    except Exception as e:
        return ['Error trying to update activity data review status', str(e)]

def delete_activity(act_obj):
    act_obj.delete()