from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Sequence
from retrobiocat_web.mongo.model_queries import reaction_queries

def enzyme_types_for_reaction(reaction_name):
    reaction = reaction_queries.reaction_from_name(reaction_name, get_related=False)
    enzyme_types = reaction.enzyme_types
    return enzyme_types

def get_enzyme_type_full_name(enzyme_type):
    """Return the full name of the given enzyme_type abbreviation"""
    return EnzymeType.objects(enzyme_type=enzyme_type)[0].full_name

def get_enzyme_type_full_name_dict():
    et_objects = EnzymeType.objects()
    et_dict = {}
    for et in et_objects:
        et_dict[et.enzyme_type] = et.full_name
    return et_dict

def get_enzyme_types_with_full_names(list_of_enzyme_types):
    full_name_dict = {}
    enzyme_type_objects = EnzymeType.objects(enzyme_type__in=list_of_enzyme_types)
    for enz_type_obj in enzyme_type_objects:
        full_name_dict[enz_type_obj.enzyme_type] = enz_type_obj.full_name
    return full_name_dict

def all_enzyme_type_strings():
    """Return a list of all the enzyme type abbreviations"""
    return EnzymeType.objects().distinct('enzyme_type')

def enzyme_type_from_name(enzyme_type, get_related=False):
    """Get enzyme type object from the enzyme type abbreviation"""

    if get_related == True:
        return EnzymeType.objects(enzyme_type=enzyme_type).select_related().first()
    else:
        return EnzymeType.objects(enzyme_type=enzyme_type).first()

def enzyme_type_names_in_paper(paper):
    return Sequence.objects(papers=paper).distinct('enzyme_type')

def enzyme_types_dict(names_only=True):
    if names_only==True:
        et_list = EnzymeType.objects().only('enzyme_type', 'full_name').as_pymongo()
        return [{x: et[x] for x in ['enzyme_type', 'full_name']} for et in et_list]

    return EnzymeType.objects().as_pymongo()

def enzyme_types_ordered_by_database_score():
    return EnzymeType.objects().order_by('-database_score').select_related()