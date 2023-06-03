from retrobiocat_web.analysis.auto_job_loop.blast_sequences import get_jobs_in_queue
from retrobiocat_web.analysis.uniprot.uniref_data_extraction import check_id_match
from retrobiocat_web.analysis.uniprot.xml_retrieval import UniRef_XML_Retriever
from retrobiocat_web.mongo.model_queries.sequence_queries import num_seqs_for_ssn_for_enzyme_type
from retrobiocat_web.mongo.models.biocatdb_models import SSN_record, UniRef50, Sequence, EnzymeType
from retrobiocat_web.mongo.models.user_models import User

from flask import current_app
import mongoengine as db
import random
import time
from datetime import datetime
from dateutil.relativedelta import relativedelta

from retrobiocat_web.analysis.ssn import ssn_tasks
from retrobiocat_web.mongo.models.other_models import AutoJobStatus
from retrobiocat_web.analysis.sequence_analysis import get_unirep
from retrobiocat_web.analysis.enzyme_type_info_scores import enzyme_scores
from retrobiocat_web.analysis.uniprot_and_web.update_uniref_details import task_get_uniref_info
from retrobiocat_web.analysis.uniprot_and_web.get_uniref_codes import blast_for_accession
from retrobiocat_web.analysis.uniprot_and_web import embl_restfull

def is_ssn_creation_disabled():
    if current_app.config['ALLOW_SSN_CREATION'] == False:  # if SSN creation is disabled
        ssn_records = SSN_record.objects().select_related()
        for ssn_r in ssn_records:
            if ssn_r.status != 'SSN creation is disabled':
                ssn_r.status = 'SSN creation is disabled'
                ssn_r.save()
        return True

    else:
        return False

def can_ssn_be_updated(ssn_record):
    if ssn_record.status != 'Complete' and ssn_record.enzyme_type.bioinformatics_status == 'Complete':
        return True
    else:
        return False

def create_update_ssn_job(ssn_record):
    enzyme_type = ssn_record.enzyme_type.enzyme_type
    job_name = f"{enzyme_type}_expand_ssn"
    current_app.alignment_queue.enqueue(ssn_tasks.task_expand_ssn, enzyme_type, job_id=job_name)
    print(f'Queued SSN job for {enzyme_type}')

def mark_ssn_for_update(enzyme_type_obj):
    ssn_query = SSN_record.objects(enzyme_type=enzyme_type_obj)
    if len(ssn_query) != 0:
        ssn_record = SSN_record.objects(enzyme_type=enzyme_type_obj)[0]
        ssn_record.status = 'Queued for update'
        ssn_record.precalc_status = 'Queued for update'
        ssn_record.save()


def check_should_make_a_new_ssn_record(enzyme_type):
    if enzyme_type.bioinformatics_status != 'Complete':
        return False

    if enzyme_type in SSN_record.objects().distinct('enzyme_type'):
        print(f"Error - SSN for {enzyme_type.enzyme_type} already exists")
        return False

    num_unirefs = UniRef50.objects(enzyme_type=enzyme_type).count()
    num_biocatdb_seqs = num_seqs_for_ssn_for_enzyme_type(enzyme_type.enzyme_type)

    if num_unirefs != 0 and num_biocatdb_seqs != 0:
        print(f"No SSN for {enzyme_type.enzyme_type}, but blasts are complete and sequences present..  creating SSN record.")
        return True

    return False


def get_or_create_ssn_record(enzyme_type):
    ssn_record = SSN_record.objects(enzyme_type=enzyme_type).first()
    if ssn_record is not None:
        return ssn_record

    if check_should_make_a_new_ssn_record(enzyme_type):
        return SSN_record(enzyme_type=enzyme_type)
    else:
        return None

def update_ssn_if_required(enzyme_type):
    if is_ssn_creation_disabled():
        return

    ssn_record = get_or_create_ssn_record(enzyme_type)
    if ssn_record is None:
        return

    if can_ssn_be_updated(ssn_record):
        create_update_ssn_job(ssn_record)  # then create/expand the ssn

