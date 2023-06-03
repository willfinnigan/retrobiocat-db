from retrobiocat_web.analysis.auto_job_loop.update_ssn import mark_ssn_for_update
from retrobiocat_web.analysis.uniprot.xml_retrieval import UniRef_XML_Retriever
from retrobiocat_web.mongo.models.biocatdb_models import UniRef50
from flask import current_app
import random
import time
import mongoengine as db


def does_xml_exist(uniref):
    name = uniref.enzyme_name
    uniref_xml_dict = UniRef_XML_Retriever().get_xml(name)
    if uniref_xml_dict is not False:
        return True
    return False

def do_unirefs_need_updating_using_random_checks(enzyme_type_obj, num_to_check=30):
    # don't check uniref if not allowed
    if current_app.config['ALLOW_UNIPROT'] == False:
        return False

    unirefs = UniRef50.objects(enzyme_type=enzyme_type_obj)

    if len(unirefs) == 0:
        return False

    for i in range(num_to_check):
        rand_uniref = random.choice(unirefs)
        if does_xml_exist(rand_uniref):
            time.sleep(0.5)
        else:
            time.sleep(0.5)
            print(f'Identified mismatches with online uniref entries for {enzyme_type_obj.enzyme_type}')
            return True

    print(f'Uniref checks ok')
    return False

def full_uniref_check(enzyme_type_obj):
    # don't check uniref if not allowed
    if current_app.config['ALLOW_UNIPROT'] == False:
        return

    unirefs = UniRef50.objects(enzyme_type=enzyme_type_obj).select_related()
    if len(unirefs) != 0:
        for ur in unirefs:
            if does_xml_exist(ur):
                time.sleep(0.5)
            else:
                print(f"{ur.enzyme_name} doesnt have an xml, deleting..")
                for seq in ur.result_of_blasts_for:
                    try:
                        seq.blast = None
                        seq.save()
                    except:
                        print(f'Could not update blast status for source of {ur.enzyme_name}')
                ur.delete()

    mark_ssn_for_update(enzyme_type_obj)

    enzyme_type_obj.bioinformatics_status = 'Queued for update'
    enzyme_type_obj.save()

    print(f"Removed UniRef50 entries which no longer exist for {enzyme_type_obj.enzyme_type}")

def ensure_all_unirefs_have_blast_source(enzyme_type):
    print('Checking for unirefs with no blast source..')
    uniref_query = UniRef50.objects(db.Q(result_of_blasts_for__size=0) & db.Q(enzyme_type=enzyme_type))
    for uniref in uniref_query:
        print(f"Deleting {uniref.enzyme_name}")
        uniref.delete()

