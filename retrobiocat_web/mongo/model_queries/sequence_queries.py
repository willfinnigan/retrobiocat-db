from retrobiocat_web.mongo.models.biocatdb_models import Activity, Sequence
import json
import mongoengine as db

def seq_obj_from_name(enzyme_name, get_related=False, include_other_names=True) -> Sequence:
    """ Get sequence object from enzyme_name (or optionally other_names) """

    seq = Sequence.objects(enzyme_name=enzyme_name).first()

    # if include_other_names, then search for a match here also
    if seq is None and include_other_names == True:
        seq = Sequence.objects(other_names_data__name=enzyme_name).first()

    if get_related == True and seq is not None:
        seq.select_related()

    return seq

def seqs_of_type(enzyme_type):
    """ Return a list of sequence names for enzyme type"""
    if enzyme_type == 'All':
        query = db.Q()
    else:
        query = db.Q(enzyme_type=enzyme_type)

    sequences = Sequence.objects(query).distinct('enzyme_name')
    sequences.sort()

    return sequences

def seqs_of_type_and_or_names(enzyme_type, enzyme_names, as_pymongo=False):
    if enzyme_type is None:
        enzyme_type_query = db.Q()
    else:
        enzyme_type_query = db.Q(enzyme_type=enzyme_type)

    if enzyme_names is None:
        enzyme_names_query = db.Q()
    else:
        enzyme_names_query = db.Q(enzyme_name__in=enzyme_names)

    if as_pymongo:
        return Sequence.objects(enzyme_type_query & enzyme_names_query).as_pymongo()
    return Sequence.objects(enzyme_type_query & enzyme_names_query)

def seqs_from_list_of_names(enzyme_names):
    seqs = Sequence.objects(db.Q(enzyme_name__in=enzyme_names))
    return seqs

def seqs_of_type_with_other_names(enzyme_type):
    """ Return a list of sequence names for enzyme type, but also includes all the other_names too"""

    if enzyme_type == 'All':
        enzyme_names = Sequence.objects().distinct('enzyme_name')
        other_names = Sequence.objects().distinct('other_names')
        sequences = enzyme_names + other_names
        sequences.sort()
    else:
        enzyme_names = Sequence.objects(enzyme_type=enzyme_type).distinct('enzyme_name')
        other_names = Sequence.objects(enzyme_type=enzyme_type).distinct('other_names')
        sequences = enzyme_names + other_names
        sequences.sort()

    return sequences

def does_seq_exist(enzyme_name):
    if Sequence.objects(enzyme_name=enzyme_name).count() == 0:
        return False
    return True

def seqs_of_paper(paper, unassigned_only=False, count_only=False, get_related=False, names_only=False):

    q_paper = db.Q(papers=paper)

    if unassigned_only == True:
        q_owner = db.Q(owner=None)
    else:
        q_owner  = db.Q()

    if count_only == True:
        return Sequence.objects(q_paper & q_owner).count()

    if names_only == True:
        return list(Sequence.objects(q_paper & q_owner).distinct('enzyme_name'))

    if get_related == True:
        return Sequence.objects(q_paper & q_owner).select_related()
    else:
        return Sequence.objects(q_paper & q_owner)


def enzyme_type_abbreviations_in_paper(paper):
    ets = Sequence.objects(papers=paper).distinct('enzyme_type')
    return ets

def seqs_from_paper(paper, select_related=False):
    if select_related:
        return Sequence.objects(papers=paper).order_by('enzyme_type').select_related()
    else:
        return Sequence.objects(papers=paper).order_by('enzyme_type')

def seq_objs_of_enzyme_type(et_abbrev, select_related=False):

    if et_abbrev == 'All':
        query = db.Q()
    else:
        query = db.Q(enzyme_type=et_abbrev)

    if select_related:
        return Sequence.objects(query).order_by('enzyme_type').select_related()
    else:
        return Sequence.objects(query).order_by('enzyme_type')


def seq_names_in_paper_not_ok_to_review(paper):
    names = Sequence.objects(db.Q(papers=paper) &
                           (db.Q(sequence='') | db.Q(sequence=None)) &
                           db.Q(sequence_unavailable__ne=True)).distinct('enzyme_name')
    return names

def unique_products(seq):
    smis = Activity.objects(enzyme_name=seq.enzyme_name).distinct('product_1_smiles')
    return smis

def get_embeddings(list_seq_objs):
    """Return a list of embeddings from a list of seq objects"""

    list_embeddings = []
    for seq in list_seq_objs:
        if seq.unirep is not None:
            embedding = json.loads(seq.unirep)
            list_embeddings.append(embedding)

    return list_embeddings

def seq_objs_from_seq_names(list_seq_names):
    """Return a list of sequence objects from a list of sequence names"""

    return Sequence.objects(enzyme_name__in=list_seq_names)

def seq_obj_other_names(seq_obj):
    """ Return a list of other names for a seq object """
    other_names = []
    for data_dict in seq_obj.other_names_data:
        other_names.append(data_dict['name'])
    return other_names

def seq_obj_other_names_data(seq_obj):
    """ Return a dict of other names data for a seq object """
    other_names_data = []
    for data_dict in seq_obj.other_names_data:
       data = {'name': data_dict['name'],
               'n_tag': data_dict['n_tag'],
               'c_tag': data_dict['c_tag'],
               'notes': data_dict['notes']}

       data = {k:v for k,v in data.items() if v is not None}
       other_names_data.append(data)
    return other_names_data

def num_sequences_owned_by_user(user):
    return Sequence.objects(owner=user).count()

def get_alt_naming_for_paper(paper):
    seqs = seqs_of_paper(paper)
    alt_names = {}
    for seq in seqs:
        for other_name_data in seq.other_names_data:
            if paper in other_name_data.papers:
                alt_names[seq.enzyme_name] = other_name_data.name
    return alt_names

def num_seqs_for_ssn_for_enzyme_type(enzyme_type_str):
    q_et = db.Q(enzyme_type=enzyme_type_str)
    q_rev = db.Q(reviewed=True)
    q_bio_ignore = db.Q(bioinformatics_ignore__ne=True)
    q_seq = db.Q(sequence__ne=None)
    q_seq2 = db.Q(sequence__ne="")
    q_seq_av = db.Q(sequence_unavailable=False)

    return Sequence.objects(q_et & q_rev & q_bio_ignore & q_seq & q_seq2 & q_seq_av).count()

def all_seqs(select_related):
    if select_related:
        return Sequence.objects().order_by('enzyme_type').select_related()
    else:
        return Sequence.objects().order_by('enzyme_type')

def sequences_with_protein_sequence():
    seqs = Sequence.objects(db.Q(sequence__ne="") &
                     db.Q(sequence__ne=None) &
                     db.Q(sequence_unavailable__ne=True))
    return seqs

if __name__ == '__main__':
    import mongoengine as db
    from retrobiocat_web.mongo.default_connection import make_default_connection
    from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType
    import yaml

    make_default_connection()

    num = num_seqs_for_ssn_for_enzyme_type('ADH(FAD)')
    print(num)


