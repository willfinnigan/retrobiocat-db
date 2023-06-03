import uuid

from retrobiocat_web.analysis.reaction_summary.ph_histogram import pH_Processor
from retrobiocat_web.analysis.reaction_summary.temperature_histogram import Temperature_Processor
from retrobiocat_web.analysis.uniprot_and_web.crossref_lookup import get_metadata_from_crossref
from retrobiocat_web.analysis.uniprot_and_web.pubmed import query_pubmed
from retrobiocat_web.app.database_bps.curation.functions import paper_status
from retrobiocat_web.app.database_bps.curation.functions.free_text_standardisation import standardise_forumation
from retrobiocat_web.app.database_bps.curation.functions.paper_status import tag_paper_with_enzyme_types
from retrobiocat_web.app.database_bps.db_admin import bp
from flask import render_template, jsonify, current_app
from flask_security import roles_required

from retrobiocat_web.mongo.modal_updates.molecular_descriptors import task_update_molecular_descriptors
from retrobiocat_web.mongo.model_queries.sequence_queries import num_seqs_for_ssn_for_enzyme_type
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence, Tag, EnzymeType, SSN_record
import time
from retrobiocat_web.mongo.functions.mongo_dump import execute_mongo_dump
from pathlib import Path
from rq.registry import ScheduledJobRegistry

from retrobiocat_web.analysis.ssn.ssn_main import SSN
from rdkit import Chem
from retrobiocat_web.mongo.modal_updates.fingerprints import task_update_fingerprints
from retrobiocat_web.analysis.uniprot_and_web.update_uniref_details import task_get_uniref_info
from retrobiocat_web.mongo.modal_updates import paper_CUD, sequence_CUD, activity_CUD


@bp.route('/other_admin_functions', methods=['GET', 'POST'])
@roles_required('admin')
def other_admin_functions():
    return render_template('init_db/other_admin_functions.html')


@bp.route('/_delete_sequences_no_paper', methods=['GET', 'POST'])
@roles_required('admin')
def delete_sequences_no_paper():
    current_app.db_queue.enqueue(task_delete_sequences_no_paper)
    return jsonify(result={'status': 'success',
                           'msg': f'Job added to queue to delete sequences',
                           'issues': []})


def task_delete_sequences_no_paper():
    print('Deleting sequences with no papers')
    count = 0
    seqs = Sequence.objects()
    for seq in seqs:
        if len(seq.papers) == 0:
            seq.delete()
            count += 1

    print(f"Deleted {count} sequences")


@bp.route('/_assign_papers', methods=['GET', 'POST'])
@roles_required('admin')
def secret_assign_papers():
    current_app.db_queue.enqueue(task_assign_papers)
    current_app.db_queue.enqueue(task_get_paper_metadata)
    current_app.db_queue.enqueue(biocatdb_init_complete)

    result = {'status': 'success',
              'msg': 'assigning papers',
              'issues': []}

    return jsonify(result=result)


def biocatdb_init_complete():
    print("SUBSTRATE SPECIFICITY DATABASE INITIALISATION COMLPETE")


def task_assign_papers():
    users = User.objects()
    papers = Paper.objects()

    for paper in papers:
        paper_status.update_status(paper)

    for user in users:
        usernames = get_usernames(user)
        for paper in papers:
            if paper.added_by is None or paper.added_by == '':
                activities = Activity.objects(paper=paper)
                for activity in activities:
                    if does_username_match(usernames, activity.added_by_string):
                        activity.added_by = user
                        activity.save()

                        if paper.added_by is None or paper.added_by == '':
                            paper.added_by = user
                            paper.owner = user
                            paper.save()


def get_usernames(user):
    usernames = []
    if len(User.objects(last_name=user.last_name)) == 1 and len(user.last_name) > 2:
        usernames.append(str(user.last_name).lower())
    if len(User.objects(first_name=user.first_name)) == 1 and len(user.first_name) > 2:
        usernames.append(str(user.first_name).lower())

    usernames.append(f"{user.first_name} {user.last_name}".lower())
    usernames.append(f"{user.first_name} {user.last_name}, {user.affiliation}".lower())

    for name in usernames:
        if len(name) > 3:
            usernames.remove(name)

    return usernames


def does_username_match(usernames, added_by_str):
    added_by_str = str(added_by_str).lower()
    if added_by_str != 'nan':
        for name in usernames:
            if name != 'tom':
                if name in added_by_str or name == added_by_str:
                    return True
    return False


def task_add_sequence_data(df):
    users = User.objects()

    for i, row in df.iterrows():
        seq_query = Sequence.objects(enzyme_name=row['enzyme_name'])
        if len(seq_query) != 0:
            seq = seq_query[0]
            if row['sequence'] != '' and row['sequence'] is not None:
                seq.sequence = str(row['sequence'])

            for user in users:
                usernames = get_usernames(user)
                if does_username_match(usernames, row['added_by']):
                    seq.added_by = user
                    seq.owner = user

            seq.save()


def task_get_paper_metadata():
    papers = Paper.objects()
    for paper in papers:
        tag_paper_with_enzyme_types(paper)

        if paper.authors is None or paper.authors == [''] or paper.authors == []:
            title, authors_list, journal, date, cite_mini = get_metadata_from_crossref(paper.doi)
            if cite_mini == '':
                title, authors_list, journal, date, cite_mini = query_pubmed(paper.doi)
                time.sleep(3)

            if cite_mini != '':
                paper.title = title
                paper.authors = authors_list
                paper.journal = journal
                paper.date = date
                paper.short_citation = cite_mini
                paper.save()

                paper_status.update_status(paper)


@bp.route('/_orphan_enzymes', methods=['GET', 'POST'])
@roles_required('admin')
def orphan_enzymes():
    current_app.db_queue.enqueue(task_search_for_orphan_enzymes)

    result = {'status': 'success',
              'msg': 'search for orphan enzymes',
              'issues': []}

    return jsonify(result=result)


def task_search_for_orphan_enzymes():
    activity_enzyme_names = list(set(Activity.objects().distinct('enzyme_name')))
    for name in activity_enzyme_names:
        if len(Sequence.objects(enzyme_name=name)) == 0:
            enzyme_type = Activity.objects(enzyme_name=name)[0].enzyme_type
            new_seq = Sequence(enzyme_name=name,
                               enzyme_type=enzyme_type)
            new_seq.save()
            print(f"found orphan enzyme, added sequence entry for {name} - {enzyme_type}")


@bp.route('/_find_tags', methods=['GET', 'POST'])
@roles_required('admin')
def find_tags():
    seqs = Sequence.objects()
    n_tags = Tag.objects(n_term=True).distinct('seq')
    n_tags = sorted(n_tags, key=len, reverse=True)
    c_tags = Tag.objects(c_term=True).distinct('seq')
    c_tags = sorted(c_tags, key=len, reverse=True)

    print(n_tags)

    for seq in seqs:
        for n_tag in n_tags:
            if n_tag == seq.sequence[0:len(n_tag)]:
                seq.n_tag = n_tag
                seq.sequence = seq.sequence[len(n_tag):]
                if seq.sequence[0] != 'M':
                    seq.sequence = 'M' + seq.sequence
                print(f"Found N term: {n_tag}")
                print(f"Removed from seq: {seq.sequence}")

        for c_tag in c_tags:
            if c_tag == seq.sequence[-len(c_tag):]:
                seq.c_tag = c_tag
                seq.sequence = seq.sequence[:-len(c_tag):]
                print(f"Found C term: {c_tag}")
                print(f"Removed from seq: {seq.sequence}")

        seq.save()

    result = {'status': 'success',
              'msg': 'Searching for tags',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_convert_to_pdb_schema', methods=['GET', 'POST'])
@roles_required('admin')
def convert_to_pdb_schema():
    seqs = Sequence.objects()
    for seq in seqs:
        if seq.structure == True:
            seq.pdb = seq.accession
            seq.accession = ''
            seq.structure = None
            seq.save()

    result = {'status': 'success',
              'msg': 'Coverting to pdb schema',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_clear_all_redis_jobs', methods=['GET', 'POST'])
@roles_required('admin')
def clear_all_redis_jobs():

    registry = ScheduledJobRegistry(queue=current_app.auto_jobs)
    for job_id in registry.get_job_ids():
        registry.remove(job_id, delete_job=True)

    for queue in current_app.redis_queues:
        queue.delete()

    result = {'status': 'success',
              'msg': 'Cleared all redis jobs',
              'issues': []}

    return jsonify(result=result)


def task_ensure_correct_sequence_naming():
    seqs = Sequence.objects()
    for seq in seqs:
        seq.update_name(seq.enzyme_name)

@bp.route('/_ensure_correct_sequence_naming', methods=['GET', 'POST'])
@roles_required('admin')
def ensure_correct_sequence_naming():
    current_app.db_queue.enqueue(task_ensure_correct_sequence_naming)

    result = {'status': 'success',
              'msg': f'Removing invalid chars sequences',
              'issues': []}

    return jsonify(result=result)



@bp.route('/_update_reviewed_status', methods=['GET', 'POST'])
@roles_required('admin')
def update_reviewed_status():
    current_app.db_queue.enqueue(task_update_review_status)

    result = {'status': 'success',
              'msg': f'Updating all review status based on original paper.reviewed',
              'issues': []}

    return jsonify(result=result)

def task_update_review_status():
    papers = Paper.objects(reviewed=True)

    for paper in papers:
        user = paper.reviewed_by
        paper_CUD.review_paper_metadata(paper, True)
        paper_CUD.review_paper_sequences(paper, user, True)
        paper_CUD.review_paper_activity(paper, user, True)
        paper_status.update_status(paper)







@bp.route('/_set_all_seqs_to_reblast', methods=['GET', 'POST'])
@roles_required('admin')
def set_all_seqs_to_reblast():
    seqs = Sequence.objects()
    for seq in seqs:
        seq.blast = None
        seq.save()

    for enz_type_obj in EnzymeType.objects():
        enz_type_obj.bioinformatics_status = 'Queued for update'
        enz_type_obj.save()

    result = {'status': 'success',
              'msg': f'Bioinformatics status reset',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_clear_empty_ssns', methods=['GET', 'POST'])
@roles_required('admin')
def clear_empty_ssns():
    ssn_records = SSN_record.objects().select_related()

    for ssn_r in ssn_records:
        enzyme_type = ssn_r.enzyme_type
        if num_seqs_for_ssn_for_enzyme_type(enzyme_type.enzyme_type) == 0:
            ssn_r.delete()

    result = {'status': 'success',
              'msg': f'Empty SSNs removed',
              'issues': []}
    return jsonify(result=result)


@bp.route('/_clear_autoprocessed_activity_data', methods=['GET', 'POST'])
@roles_required('admin')
def clear_autoprocessed_activity_data():
    activity_data = Activity.objects(auto_generated=True)

    for act in activity_data:
        act.delete()

    result = {'status': 'success',
              'msg': f'Autoprocessed data deleted',
              'issues': []}

    return jsonify(result=result)




@bp.route('/_clear_ssn_position_info', methods=['GET', 'POST'])
@roles_required('admin')
def clear_ssn_position_info():
    ssn_records = SSN_record.objects().select_related()

    for ssn_r in ssn_records:
        enzyme_type = ssn_r.enzyme_type.enzyme_type
        ssn = SSN(enzyme_type)
        ssn.db_object = ssn_r
        ssn.clear_position_information()
        ssn_r.status = 'To update'
        ssn_r.precalc_status = 'To update'
        ssn_r.save()

    result = {'status': 'success',
              'msg': f'Cleared SSN position info',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_ensure_all_rdkit', methods=['GET', 'POST'])
@roles_required('admin')
def ensure_all_rdkit():
    current_app.db_queue.enqueue(task_ensure_all_rdkit_smiles)

    result = {'status': 'success',
              'msg': f'Converting all smiles to rdkit smiles',
              'issues': ['This will take some time..']}

    return jsonify(result=result)


def task_ensure_all_rdkit_smiles():
    activity = Activity.objects()

    for act in activity:
        need_to_save = False
        try:
            if act.substrate_1_smiles != '' or act.substrate_1_smiles is not None:
                mol = Chem.MolFromSmiles(act.substrate_1_smiles)
                smi = Chem.MolToSmiles(mol)
                if smi != act.substrate_1_smiles:
                    act.substrate_1_smiles = smi
                    need_to_save = True
        except:
            print(f"Error converting substrate_1_smiles - {act.substrate_1_smiles}")

        try:
            if act.substrate_2_smiles != '' or act.substrate_2_smiles is not None:
                mol = Chem.MolFromSmiles(act.substrate_2_smiles)
                smi = Chem.MolToSmiles(mol)
                if smi != act.substrate_2_smiles:
                    act.substrate_2_smiles = smi
                    need_to_save = True
        except:
            print(f"Error converting substrate_2_smiles - {act.substrate_2_smiles}")

        try:
            if act.product_1_smiles != '' or act.product_1_smiles is not None:
                mol = Chem.MolFromSmiles(act.product_1_smiles)
                smi = Chem.MolToSmiles(mol)
                if smi != act.product_1_smiles:
                    act.product_1_smiles = smi
                    need_to_save = True
        except:
            print(f"Error converting product_1_smiles - {act.product_1_smiles}")

        if need_to_save == True:
            act.save()

    current_app.task_queue.enqueue(task_update_fingerprints, job_id=f"update_molecule_fingerprints_{str(uuid.uuid4())}")
    current_app.task_queue.enqueue(task_update_molecular_descriptors, job_id=f"update_molecule_descriptors_{str(uuid.uuid4())}")

def task_ensure_positive_categorical_data_has_positive_binary():
    activity = Activity.objects()
    for act in activity:
        if act.categorical == "High" or act.categorical == "Medium" or act.categorical == "Low":
            if act.binary != True:
                act.binary = True
                act.save()

@bp.route('/_ensure_all_cat_binary_match', methods=['GET', 'POST'])
@roles_required('admin')
def ensure_all_cat_binary_match():
    current_app.db_queue.enqueue(task_ensure_positive_categorical_data_has_positive_binary)

    result = {'status': 'success',
              'msg': f'Ensuring all binary data match categorical',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_force_all_uniref_update', methods=['GET', 'POST'])
@roles_required('admin')
def force_all_uniref_update():
    enzyme_types = EnzymeType.objects()

    for et in enzyme_types:
        job_id = f"force_uniref_update_get_uniref_info_for_{et.enzyme_type}_enzymes"
        current_app.db_queue.enqueue(task_get_uniref_info, et.enzyme_type, True, job_id=job_id)

    result = {'status': 'success',
              'msg': f'Updating all UniRef info',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_force_no_duplicate_other_names', methods=['GET', 'POST'])
@roles_required('admin')
def force_no_duplicate_other_names():
    """Ensure other_names field in sequences has no duplicates"""

    seqs = Sequence.objects()

    for seq in seqs:
        for name in seq.other_names:
            if Sequence.objects(enzyme_name=name).count() + Sequence.objects(other_names=name).count() > 1:
                seq.other_names.remove(name)
                seq.save()

    result = {'status': 'success',
              'msg': f'Updated other names field',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_clean_seq_list', methods=['GET', 'POST'])
@roles_required('admin')
def clean_seq_lists():
    """Remove random blank entries in list fields for sequences"""
    def remove_empty_space(seq_list):
        if ' ' in seq_list:
            seq_list.remove(' ')
        if '' in seq_list:
            seq_list.remove('')
        return seq_list

    def remove_duplicates(seq_list):
        seq_list = list(set(seq_list))
        return seq_list

    seqs = Sequence.objects()

    for seq_obj in seqs:
        seq_obj.other_names = remove_empty_space(seq_obj.other_names)
        seq_obj.other_names = remove_duplicates(seq_obj.other_names)
        seq_obj.other_identifiers = remove_empty_space(seq_obj.other_identifiers)
        seq_obj.other_identifiers = remove_duplicates(seq_obj.other_identifiers)
        seq_obj.save()

    result = {'status': 'success',
              'msg': f'Updated sequence list fields',
              'issues': []}

    return jsonify(result=result)


@bp.route('/_clear_all_embeddings', methods=['GET', 'POST'])
@roles_required('admin')
def clear_all_embeddings():
    """Removes all embeddings from sequences"""

    seqs = Sequence.objects()

    for seq in seqs:
        seq.unirep = None
        seq.seqvec = None
        seq.save()

    result = {'status': 'success',
              'msg': f'Deleted all embeddings',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_convert_other_names_to_data_dicts', methods=['GET', 'POST'])
@roles_required('admin')
def convert_other_names_to_data_dicts():
    seqs = Sequence.objects()

    for seq in seqs:
        for name in seq.other_names:
            data_dict = {'existing_name': name}
            sequence_CUD.update_other_names(seq, data_dict)

    result = {'status': 'success',
              'msg': f'Converted other names from lists to embedded docs',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_mongo_dump', methods=['GET', 'POST'])
@roles_required('admin')
def mongo_dump():
    execute_mongo_dump()

    result = {'status': 'success',
              'msg': f'Mongo dump command initiated',
              'issues': []}

    return jsonify(result=result)

from retrobiocat_web.mongo.models.user_models import user_datastore

@bp.route('/_delete_all_users', methods=['GET', 'POST'])
@roles_required('admin')
def delete_all_users():
    """ Deletes all accounts except the admin account"""

    admin_email = current_app.config['ADMIN_EMAIL']
    users = User.objects(email__ne=admin_email)  # users which dont have the admin email address

    for user in users:
        user.delete()

    result = {'status': 'success',
              'msg': f'All user accounts apart from default admin account are not deleted',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_delete_seq_act_with_no_paper', methods=['GET', 'POST'])
@roles_required('admin')
def delete_seq_act_with_no_paper():
    """ Deletes any sequence or activity data which has no paper"""

    # activity with no papers
    acts = Activity.objects(paper=None)
    for act_obj in acts:
        activity_CUD.delete_activity(act_obj)

    # sequences with no papers
    seqs = Sequence.objects(papers__size=0)
    for seq_obj in seqs:
        sequence_CUD.delete_sequence(seq_obj)

    result = {'status': 'success',
              'msg': f'Deleted {len(acts)} activities, and {len(seqs)} sequences',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_force_ph_and_temp_into_float', methods=['GET', 'POST'])
@roles_required('admin')
def force_ph_and_temp_into_float():
    acts = Activity.objects()
    print(len(acts))
    temp_processor = Temperature_Processor()
    ph_processor = pH_Processor()
    for activity in acts:
        if activity.temperature is not None:
            try:
                float_temp = temp_processor.parse_single_temp(activity.temperature)
                if isinstance(float_temp, float):
                    activity.temperature = float_temp + 0.1
                    activity.temperature -= 0.1
                else:
                    activity.temperature = None

            except Exception as e:
                activity.temperature = None

        if activity.ph is not None:
            try:
                float_ph = ph_processor.parse_single_ph(activity.ph)
                if isinstance(float_ph, float):
                    activity.ph = float_ph + 0.1
                    activity.ph -= 0.1
                else:
                    activity.ph = None

            except Exception as e:
                print(e)
                activity.ph = None

        activity.save()


    print('Done')

    result = {'status': 'success',
              'msg': f'Updated ph and temperature fields to be floats',
              'issues': []}

    return jsonify(result=result)

@bp.route('/_force_formulations_to_standard_format', methods=['GET', 'POST'])
@roles_required('admin')
def force_formulations_to_standard_format():
    acts = Activity.objects()
    for act in acts:
        if act.formulation is not None:
            act.formulation = standardise_forumation(act.formulation)
            act.save()

    print('Done')
    result = {'status': 'success',
              'msg': f'Updated formulation fields to standard format',
              'issues': []}

    return jsonify(result=result)






if __name__ == '__main__':
    data_folder = str(Path(__file__).parents[4]) + '/retro/data/buyability'
    print(data_folder)