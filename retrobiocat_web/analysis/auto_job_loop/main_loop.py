import time

from retrobiocat_web.analysis.auto_job_loop.bicatdb_sequence_jobs import task_check_and_generate_biocatdb_embeddings, \
    task_add_uniref_codes
from retrobiocat_web.analysis.auto_job_loop.blast_sequences import complete_blast_for_enzyme_type_if_required
from retrobiocat_web.analysis.auto_job_loop.check_uniref import do_unirefs_need_updating_using_random_checks, full_uniref_check, ensure_all_unirefs_have_blast_source
from retrobiocat_web.analysis.auto_job_loop.update_ssn import update_ssn_if_required

from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.analysis.enzyme_type_info_scores import enzyme_scores


def enzyme_type_loop(enzyme_type):

    # 1. do unirefs need updating?
    ensure_all_unirefs_have_blast_source(enzyme_type)
    if do_unirefs_need_updating_using_random_checks(enzyme_type, num_to_check=25):
        full_uniref_check(enzyme_type)

    # 2. do sequences need blasting?
    complete_blast_for_enzyme_type_if_required(enzyme_type)

    # 3. does ssn need updating?
    update_ssn_if_required(enzyme_type)


def main_loop(queue, min_time=60*60):
    t0 = time.time()

    for enzyme_type in EnzymeType.objects(enzyme_type='ADH'):
        print(f"Running bioinformatics jobs for {enzyme_type.enzyme_type}")
        enzyme_type_loop(enzyme_type)

    queue.enqueue(enzyme_scores.update_scores)
    queue.enqueue(task_check_and_generate_biocatdb_embeddings)
    queue.enqueue(task_add_uniref_codes)

    time_running = time.time() - t0
    if time_running < min_time:
        print(f"Time running {time_running} is less than min time {min_time}, sleeping for remainder")
        time.sleep(min_time - time_running)






