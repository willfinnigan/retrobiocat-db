from flask import current_app
import mongoengine as db

from retrobiocat_web.analysis.sequence_analysis import get_unirep
from retrobiocat_web.analysis.uniprot_and_web.get_uniref_codes import blast_for_accession
from retrobiocat_web.mongo.models.biocatdb_models import Sequence


def task_add_uniref_codes():
    """Adds any uniref codes to biocatdb sequences which have no identifier"""

    # don't check uniref if not allowed
    if current_app.config['ALLOW_UNIPROT'] == False:
        return

    print('add unirep codes job')
    if len(current_app.blast_queue.jobs) + len(current_app.process_blasts_queue.jobs) + len(current_app.alignment_queue.jobs) == 0:
        seq_q1 = db.Q(sequence__ne=None)
        seq_q2 = db.Q(sequence__ne='')
        acc_q1 = db.Q(accession=None)
        acc_q2 = db.Q(accession='')

        seqs_to_blast = Sequence.objects(seq_q1 & seq_q2 & (acc_q1 | acc_q2))
        for seq in seqs_to_blast:
            job_id = f"blast_for_accession_for_{seq.enzyme_name}"
            current_app.blast_queue.enqueue(blast_for_accession, seq, job_id=job_id)


def task_check_and_generate_biocatdb_embeddings():
    """Generates embeddings for biocat sequences"""
    print('generate embeddings job')

    num_embeddings_jobs = len(current_app.embeddings_queue.jobs)
    batch_size = 10
    if num_embeddings_jobs == 0:
        seq_names = Sequence.objects(db.Q(sequence__ne="") &
                                     db.Q(sequence__ne=None) &
                                     db.Q(sequence_unavailable__ne=True) &
                                     db.Q(unirep=None)).distinct('enzyme_name')

        chunked_seqs = [seq_names[i:i + batch_size] for i in range(0, len(seq_names), batch_size)]
        for seqs in chunked_seqs:
            current_app.embeddings_queue.enqueue(get_unirep.generate_and_save_embedding, seqs)

    else:
        print(f"Number of embeddings jobs = {num_embeddings_jobs}")