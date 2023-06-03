from rq.registry import StartedJobRegistry
from flask import current_app
import time
from retrobiocat_web.analysis.uniprot_and_web import embl_restfull

def get_jobs_in_queue():
    queued_blast_jobs = len(current_app.blast_queue.jobs)
    started_blast_jobs = len(StartedJobRegistry('blast', connection=current_app.redis).get_job_ids())

    queued_process_blast_jobs = len(current_app.process_blasts_queue.jobs)
    started_process_blast_jobs = len(StartedJobRegistry('process_blasts', connection=current_app.redis).get_job_ids())

    queued_alignment_jobs = len(current_app.alignment_queue.jobs)
    started_alignment_jobs = len(StartedJobRegistry('alignment', connection=current_app.redis).get_job_ids())

    num_jobs = queued_blast_jobs + started_blast_jobs \
           + queued_process_blast_jobs + started_process_blast_jobs \
           + queued_alignment_jobs + started_alignment_jobs

    return num_jobs

def do_blasts_need_doing_for_enzyme_type(enzyme_type):
    # check all the seqs of enzyme type and update blast status
    embl_restfull.update_blast_status(enzyme_type.enzyme_type)

    if enzyme_type.bioinformatics_status != 'Complete':
        return True
    else:
        return False


def complete_blast_for_enzyme_type_if_required(enzyme_type):
    # don't blast if not allowed
    if current_app.config['ALLOW_UNIPROT'] == False:
        return

    print(f"Running blasts for {enzyme_type}")
    if do_blasts_need_doing_for_enzyme_type(enzyme_type) == True:
        embl_restfull.set_blast_jobs(enzyme_type.enzyme_type)

