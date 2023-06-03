from retrobiocat_web.mongo.models.biocatdb_models import Sequence
import mongoengine as db
from retrobiocat_web.analysis.uniprot_and_web.embl_restfull import BlastRunner, BlastParser


def run_uniref100_blast(sequence):
    br = BlastRunner()
    br.db = 'uniref100'
    br.alignments = 50

    blast_record = br.run(sequence)
    return blast_record

def get_match_from_blast_record(blast_record, seq_obj):
    query_length = len(seq_obj.sequence)
    match = None
    for alignment in blast_record.alignments:
        identifier = alignment.hit_id.replace('UR100:UniRef100_', '')
        if len(alignment.hsps) == 1:
            hsp = alignment.hsps[0]
            coverage = hsp.align_length / query_length
            identity = hsp.identities / hsp.align_length
            if coverage == 1 and identity == 1:
                match = identifier
                break
    return match

def get_sequence(enzyme_name):
    name_q = db.Q(enzyme_name=enzyme_name)
    seq_q1 = db.Q(sequence__ne=None)
    seq_q2 = db.Q(sequence__ne='')
    acc_q1 = db.Q(accession=None)
    acc_q2 = db.Q(accession='')

    seq_q = Sequence.objects(name_q & seq_q1 & seq_q2 & (acc_q1 | acc_q2))
    if len(seq_q) == 0:
        return None

    else:
        return seq_q[0]

def blast_for_accession(seq_obj):

    if seq_obj.sequence == None or seq_obj.sequence == "":
        return False

    if seq_obj == None:
        return False

    blast_record = run_uniref100_blast(seq_obj.sequence)
    match = get_match_from_blast_record(blast_record, seq_obj)

    if match != None:
        seq_obj.accession = match
        seq_obj.save()
        return True

    else:
        seq_obj.accession = '-'
        seq_obj.save()
        return False


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    seq = Sequence.objects(enzyme_type='AAD')[1]

    match = blast_for_accession(seq)

    print(match)


