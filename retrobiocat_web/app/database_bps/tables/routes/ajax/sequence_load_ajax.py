from flask import request, jsonify
from flask_security import current_user
from retrobiocat_web.app.app import user_datastore

from retrobiocat_web.analysis.uniprot_and_web import lookup_accession
from retrobiocat_web.analysis.sequence_analysis import find_identical_sequences
from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.mongo.model_queries import reaction_queries, sequence_queries, enzyme_type_queries, paper_queries

""" 
Ajax routes for loading sequence information 
"""

@bp.route('/_sequences_of_same_type', methods=['GET', 'POST'])
def get_sequences_of_same_type():
    seq_obj = sequence_queries.seq_obj_from_name(request.form['enzyme_name'])
    sequences = sequence_queries.seqs_of_type(seq_obj.enzyme_type)
    seq_dict = {}
    for seq in sequences:
        seq_dict[seq] = seq

    result = {'sequences': seq_dict}
    return jsonify(result=result)

@bp.route('/_sequences_of_type', methods=['GET', 'POST'])
def get_sequences_of_type():
    sequences = sequence_queries.seqs_of_type(request.form['enzyme_type'])
    seq_dict = {}
    for seq in sequences:
        seq_dict[seq] = seq
    result = {'sequences': seq_dict}
    return jsonify(result=result)

@bp.route('/_sequences_of_type_with_other_names', methods=['GET', 'POST'])
def get_sequences_of_type_with_other_names():
    sequences = sequence_queries.seqs_of_type_with_other_names(request.form['enzyme_type'])
    seq_dict = {}
    for seq in sequences:
        seq_dict[seq] = seq
    result = {'sequences': seq_dict}
    return jsonify(result=result)


@bp.route('/_reactions_of_type', methods=['GET', 'POST'])
def get_reactions_of_type():
    enzyme_type = request.form['enzyme_type']
    reaction_list = reaction_queries.reactions_of_type([enzyme_type], names_only=False)

    # convert to dict for easier creation of a select
    reaction_dict = {}
    for name in reaction_list:
        reaction_dict[name] = name

    result = {'reactions': reaction_dict}

    return jsonify(result=result)

@bp.route('/_load_sequence_data', methods=['GET', 'POST'])
def load_sequence_data():
    """Returns a dict with all the sequence data for a given sequence name"""

    enzyme_name = request.form['enzyme_name']

    if enzyme_name == '':
        return jsonify(result={})

    seq_obj = sequence_queries.seq_obj_from_name(enzyme_name, get_related=True)

    sequences = sequence_queries.seqs_of_type(seq_obj.enzyme_type)
    seq_dict = {}
    for seq in sequences:
        seq_dict[seq] = seq

    can_edit = 'none'
    self_assigned = False
    other_user = False
    can_review = False
    if current_user.is_authenticated:
        user = user_datastore.get_user(current_user.id)
        if check_permission.check_seq_permissions(current_user.id, seq_obj):
            can_edit = 'full'
        elif check_permission.check_seq_partial_permissions(current_user.id, seq_obj):
            can_edit = 'partial'

        if seq_obj.owner == user:
            self_assigned = True
        else:
            if seq_obj.owner != '' and seq_obj.owner is not None:
                other_user = True

        if check_permission.check_seq_review(current_user.id, seq_obj):
            can_review = True

    if seq_obj.owner is None:
        owner = ''
    else:
        owner = f"{seq_obj.owner.first_name} {seq_obj.owner.last_name}, {seq_obj.owner.affiliation}"

    other_names = ''
    list_other_names = [data['name'] for data in seq_obj.other_names_data]
    for i, name in enumerate(list_other_names):
        other_names += name
        if (len(list_other_names) > 1) and (i < len(list_other_names)-1):
            other_names += ', '

    other_names_data = []
    for data in seq_obj.other_names_data:
        data_dict = {'name': data.name,
                     'n_tag': data.n_tag,
                     'c_tag': data.c_tag,
                     'notes': data.notes}
        data_dict = {k: v for k, v in data_dict.items() if v != '' and v is not None}

        #if len(data_dict.keys()) > 1:  # this code to only send alt names with some data
        other_names_data.append(data_dict)

    other_identifiers = ''
    for i, ident in enumerate(seq_obj.other_identifiers):
        other_identifiers += ident
        if (len(seq_obj.other_identifiers) > 1) and (i < len(seq_obj.other_identifiers)-1):
            other_identifiers += ', '

    enzyme_type_full = enzyme_type_queries.get_enzyme_type_full_name(seq_obj.enzyme_type)

    if seq_obj.n_tag is None:
        seq_obj.n_tag = ''
    if seq_obj.c_tag is None:
        seq_obj.c_tag = ''
    if seq_obj.pdb is None:
        seq_obj.pdb = ''

    result = {'enzyme_type': seq_obj.enzyme_type,
              'enzyme_name': seq_obj.enzyme_name,
              'sequence': seq_obj.sequence,
              'sequence_unavailable': seq_obj.sequence_unavailable,
              'n_tag': seq_obj.n_tag,
              'c_tag': seq_obj.c_tag,
              'accession': seq_obj.accession,
              'other_identifiers': other_identifiers,
              'pdb': seq_obj.pdb,
              'mutant_of': seq_obj.mutant_of,
              'sequences': seq_dict,
              'notes': seq_obj.notes,
              'bioinformatics_ignore': seq_obj.bioinformatics_ignore,
              'can_edit': can_edit,
              'can_review': can_review,
              'self_assigned': self_assigned,
              'owner_is_another_user': other_user,
              'other_names': other_names,
              'other_names_list': list_other_names,
              'other_names_data': other_names_data,
              'owner': owner,
              'enzyme_type_full': enzyme_type_full,
              'seq_reviewed': seq_obj.reviewed}

    return jsonify(result=result)

@bp.route('/_load_sequence_papers', methods=['GET', 'POST'])
def load_sequence_papers():
    enzyme_name = request.form['enzyme_name']
    seq = sequence_queries.seq_obj_from_name(enzyme_name, get_related=True)

    papers_list = []
    for paper in seq.papers:
        paper_dict = {}
        paper_dict['_id'] = str(paper.id)
        paper_dict['short_citation'] = paper.short_citation
        paper_dict['doi'] = paper.doi
        paper_dict['title'] = paper.title

        papers_list.append(paper_dict)

    result = {'papers': papers_list}

    return jsonify(result=result)

@bp.route('/_get_sequence_from_uniprot',methods=['GET', 'POST'])
def get_sequence_from_uniprot():
    accession = request.form['accession']

    seq, loaded_from = lookup_accession.get_sequence(accession)

    if seq != "":
        result = {'status': 'success',
                  'msg': f'Sequence loaded from {loaded_from}',
                  'issues': [],
                  'seq': seq}

    else:
        result = {'status': 'danger',
                  'msg': 'Sequence not found',
                  'issues': [],
                  'seq': seq}

    return jsonify(result=result)

@bp.route('/_does_sequence_match',methods=['POST'])
def check_sequence_match():
    sequence_to_check = request.form['sequence']
    enzyme_type = request.form['enzyme_type']
    match = find_identical_sequences.does_sequence_match_check(enzyme_type, sequence_to_check)
    return jsonify(result={'match': match})

@bp.route('/_load_sequence_other_names_data',methods=['POST'])
def load_sequence_other_names_data():
    enzyme_name = request.form['enzyme_name']
    existing_name = request.form['existing_name']
    seq = sequence_queries.seq_obj_from_name(enzyme_name, get_related=False)
    other_names_data = sequence_queries.seq_obj_other_names_data(seq)

    # get only specified name (or blank dict if not present already)
    single_other_name = {'name': '', 'n_tag': '', 'c_tag': '', 'notes': ''}
    for data in other_names_data:
        if data.get('name', None) == existing_name:
            single_other_name.update(data)

    return jsonify(result={'other_names_data': single_other_name})


@bp.route('/_get_possible_alternative_naming_for_paper', methods=['GET', 'POST'])
def get_possible_alternative_naming_for_paper():
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    sequences = sequence_queries.seqs_of_paper(paper, get_related=True)

    alt_names = {}  # 'name': [[alt_names], selected]
    for seq in sequences:
        all_names = [seq.enzyme_name]
        selected = None
        for other_name_data in seq.other_names_data:
            all_names.append(other_name_data.name)
            if paper in other_name_data.papers:
                selected = other_name_data.name
        if selected is None:
            selected = seq.enzyme_name

        alt_names[seq.enzyme_name] = [all_names, selected]

    result = {'alt_names': alt_names}
    return jsonify(result=result)

