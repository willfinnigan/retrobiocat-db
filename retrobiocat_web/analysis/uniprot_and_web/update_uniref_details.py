from retrobiocat_web.analysis.uniprot.uniprot_data_extraction import get_pfams, get_pdbs, get_rhea, is_swissprot
from retrobiocat_web.analysis.uniprot.uniref_data_extraction import get_uniref_members, get_cluster_name, \
    strip_uniref_name, get_uniprot_id
from retrobiocat_web.analysis.uniprot.xml_retrieval import UniProt_XML_Retriever, UniRef_XML_Retriever
from retrobiocat_web.mongo.models.biocatdb_models import UniRef50, EnzymeType
import mongoengine as db

def update_uniref_obj(uniref):
    uniref_xml_dict = UniRef_XML_Retriever().get_xml(uniref.enzyme_name)

    cluster_name = get_cluster_name(uniref_xml_dict)
    rep_name = strip_uniref_name(cluster_name)

    uniref.cluster_name = cluster_name
    uniref.rep_name = rep_name

    uni90, uni100, uniprot, uniprot_kb = get_uniref_members(uniref_xml_dict)
    uni100_codes = []
    for u in uni100:
        code = u.replace('UniRef100_', '')
        uni100_codes.append(code)

    uniref.uni100 = uni100_codes
    uniref.num_uni90 = len(uni90)
    uniref.num_uni100 = len(uni100)
    uniref.num_uniprot = len(list(uniprot.keys()))

    uniprot_id = get_uniprot_id(uniref.cluster_name)
    if uniprot_id[0:2] == 'UP':
        uniprot_id = ""

    if uniprot_id != "":
        uniprot_xml_dict = UniProt_XML_Retriever().get_xml(uniprot_id)
        pfams = get_pfams(uniprot_xml_dict)
        pdbs = get_pdbs(uniprot_xml_dict)
        rhea = get_rhea(uniprot_xml_dict)
        sp = is_swissprot(uniprot_xml_dict)
    else:
        pfams = {}
        pdbs = []
        rhea = []
        sp = False

    uniref.pfams = pfams
    uniref.pdbs = pdbs
    uniref.rhea = rhea
    uniref.pfam_codes = list(pfams.keys())
    uniref.sp_annotated = sp
    uniref.save()

def task_get_uniref_info(enzyme_type, update_all):
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    etq = db.Q(enzyme_type=enzyme_type_obj)

    if update_all == True:
        rq = db.Q()
    else:
        rq = db.Q(rep_name=None)

    unirefs = UniRef50.objects(etq & rq)
    for uniref in unirefs:
        update_uniref_obj(uniref)

