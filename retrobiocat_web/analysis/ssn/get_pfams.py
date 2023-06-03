from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50, SSN_record
from retrobiocat_web.analysis.uniprot_and_web import rhea_api
import mongoengine as db
import time

def get_pfams(enzyme_type):
    et_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    unirefs = UniRef50.objects(enzyme_type=et_obj).distinct('pfams')
    pfams = {}
    for ur_pfam in unirefs:
        pfams.update(ur_pfam)
    return pfams

def count_pfams(enzyme_type, pfams):
    et_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    q_et = db.Q(enzyme_type=et_obj)
    pfam_counts = {}
    for pfam in pfams:
        q_pfam = db.Q(pfam_codes=pfam)
        pfam_counts[pfam] = UniRef50.objects(q_et & q_pfam).count()
    return pfam_counts

def order_pfams_by_count(pfams, pfam_counts):
    pfam_counts = {k: v for k, v in sorted(pfam_counts.items(), key=lambda item: item[1], reverse=True)}
    pfams = {k: pfams[k] for k in pfam_counts.keys()}
    return pfams, pfam_counts



if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    t0 = time.time()
    pfams = get_pfams('CAR')
    t1 = time.time()
    print(pfams)
    print(f"{round(t1-t0, 2)} seconds")

    t2 = time.time()
    counts = count_pfams('IRED', pfams)
    t3 = time.time()
    print(counts)
    print(f"{round(t3 - t2, 2)} seconds")
