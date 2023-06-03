from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, SSN_record
from retrobiocat_web.analysis.ssn import ssn_main

def enzyme_type_needs_blasting(enzyme_type):
    """Mark that the sequences for enzyme_type need blasting"""

    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    enz_type_obj.bioinformatics_status = 'Queued for update'
    enz_type_obj.save()

def SSN_needs_updating(enzyme_type):
    """SSN for enzyme_type needs updating"""

    enz_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    ssn_q = SSN_record.objects(enzyme_type=enz_type_obj)
    if len(ssn_q) != 0:
        ssn_record = ssn_q[0]
        ssn_record.status = 'Queued for update'
        ssn_record.save()

def rename_ssn(old_enzyme_type, new_enzyme_type):
    """Renames the directory for an SSN"""

    ssn = ssn_main.SSN(old_enzyme_type)
    ssn.rename_dir(new_enzyme_type)  #  maybe this code should be here rather than a ssn method?


