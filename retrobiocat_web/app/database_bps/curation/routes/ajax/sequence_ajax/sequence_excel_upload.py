from retrobiocat_web.app.database_bps.curation import bp
from flask import request, jsonify
from flask_security import roles_required, current_user

from retrobiocat_web.app.database_bps.curation.functions.clean_str_of_spaces import remove_end_spaces_if_str
from retrobiocat_web.app.database_bps.curation.functions.seq_obj_to_json import get_seq_table_entry
from retrobiocat_web.app.app import user_datastore
from werkzeug.utils import secure_filename
import os
import pandas as pd
import numpy as np

from retrobiocat_web.app.database_bps.curation.routes.ajax.error_results import not_a_post_request_error, \
    not_an_excel_file
from retrobiocat_web.mongo.functions.save_sequence_data.save_sequence_functions import sequence_check

from retrobiocat_web.analysis.uniprot_and_web import lookup_accession
from retrobiocat_web.mongo.modal_updates.sequence_CUD import update_sequence

from retrobiocat_web.mongo.model_queries import sequence_queries
from retrobiocat_web.mongo.modal_updates import sequence_CUD
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id

INVALID_NAME_CHARS = [".", "(", ")", "'", "/"]

@bp.route('/_upload_sequence_excel', methods=['GET', 'POST'])
@roles_required('contributor')
def upload_sequence_excel_existing():
    """Route for uploading excel file"""

    if request.method != 'POST':
        return not_a_post_request_error()

    excel_file = request.files['file_seq']
    paper = paper_from_id(request.form['paper_id_field'])
    existing_or_new = request.form['new_or_existing']
    filename = secure_filename(excel_file.filename)

    if filename[-5:] != '.xlsx':
        return not_an_excel_file()

    excel_file.save(filename)
    df = pd.read_excel(filename)

    # limit number of rows to 400 if user is not super contributor
    if not current_user.has_role('super_contributor') and len(df.index) > 400:
        result = {'status': 'danger',
                  'msg': 'Can not load more than 400 rows',
                  'issues': ['Please email your excel to an admin for addition'],
                  'new_seq_msgs': [],
                  'new_seqs_data': []}
        return jsonify(result=result)

    data_list = process_uploaded_excel(df)
    os.remove(filename)

    if existing_or_new == 'new':
        save_issues, seqs_added, new_seq_msgs = excel_create_new_seqs(data_list, paper)
    elif existing_or_new == 'existing':
        save_issues, seqs_added, new_seq_msgs = excel_add_existing_sequences(data_list, paper)
    else:
        result = {'status': 'danger',
                  'msg': 'Error processing file',
                  'issues': ["Mode must be either 'new' or 'existing'"]}
        return jsonify(result=result)

    new_seqs_data = {}

    for seq in seqs_added:
        new_seqs_data[str(seq.id)] = get_seq_table_entry(seq)

    if len(save_issues) == 0:
        result = {'status': 'success',
                  'msg': 'Sequences saved and added to paper',
                  'issues': [],
                  'new_seq_msgs': new_seq_msgs,
                  'new_seqs_data': new_seqs_data}
    else:
        result = {'status': 'warning',
                  'msg': 'Sequences saved with some issues:',
                  'issues': save_issues,
                  'new_seq_msgs': new_seq_msgs,
                  'new_seqs_data': new_seqs_data}

    return jsonify(result=result)



def process_uploaded_excel(df):
    """Process the uploaded excel into a data list for adding sequences"""

    cols = ['enzyme_type', 'enzyme_name', 'other_names', 'sequence',
            'sequence_unavailable', 'accession', 'structure',
            'mutant_of', 'notes']

    cols_to_keep = [c for c in cols if c in list(df.columns)]
    df = df[cols_to_keep]
    df.replace(np.nan, '', inplace=True)

    data_list = df.to_dict(orient='records')

    for i, data in enumerate(data_list):
        for key, value in data_list[i].items():
            data_list[i][key] = remove_end_spaces_if_str(value)  # removes any trailing end spaces from any of the cells

        # process to correct format
        if 'sequence_unavailable' in data_list[i]:
            if data_list[i]['sequence_unavailable'] == '':
                data_list[i]['sequence_unavailable'] = 'False'
        if 'structure' in data_list[i]:
            if data_list[i]['structure'] == '':
                data_list[i]['structure'] = 'False'
        if 'sequence' in data_list[i]:
            data_list[i]['sequence'] = data_list[i]['sequence'].replace('\n', '')
            data_list[i]['sequence'] = data_list[i]['sequence'].replace(' ', '')

    return data_list

def attempt_accession_lookup_if_needed(seq_obj):
    if seq_obj is None:
        return None
    if seq_obj.sequence != '' and seq_obj.sequence is not None:
        return None
    if seq_obj.accession == '' or seq_obj.accession is None:
        return None

    seq_str, loaded_from = lookup_accession.get_sequence(seq_obj.accession)

    if seq_str != '':
        update_sequence(seq_obj, seq_str)

    return None

def excel_create_new_seqs(data_list, paper):
    """
    Take data list after uploading excel and try to add new sequences
    Will only create new sequences.  Any existing names will be ignored.
    If no sequence but accension is provided, try to lookup sequence
    """

    user = user_datastore.get_user(current_user.id)
    issues = []
    new_seq_msgs = []
    enzyme_types = all_enzyme_type_strings()
    seqs_added = []

    for seq_dict in data_list:
        if seq_dict.get('enzyme_name', '') == '':
            issues.append(f"Sequence must have a name")
        elif any(l in seq_dict.get('enzyme_name', '') for l in INVALID_NAME_CHARS):
            issues.append(f"Sequence cannot contain the following characters: {INVALID_NAME_CHARS}")
        else:

            # if name exists already, cant create..
            if sequence_queries.seq_obj_from_name(seq_dict['enzyme_name'], include_other_names=True) is not None:
                issues.append(f"A sequence already exists with the name {seq_dict['enzyme_name']}")

            # otherwise create new sequence
            else:

                # sequence must have an enzyme type
                if seq_dict.get('enzyme_type', '') not in enzyme_types:
                    print(f"Enzyme type {seq_dict.get('enzyme_type', '')} does not exist")
                    issues.append(f"Enzyme type {seq_dict.get('enzyme_type', '')} does not exist")

                # any protein sequence must use correct characters
                elif sequence_check(seq_dict.get('sequence', '')) == False:
                    print(f"Amino acid sequence for {seq_dict['enzyme_name']} uses incorrect amino acid characters")
                    issues.append(
                        f"Amino acid sequence for {seq_dict['enzyme_name']} uses incorrect amino acid characters")
                else:
                    print('Creating new sequence..')

                    seq_creation_issues, seq = sequence_CUD.create_new_sequence(seq_dict, user=user, paper=paper)
                    issues += seq_creation_issues

                    attempt_accession_lookup_if_needed(seq)

                    if seq is not None:
                        seqs_added.append(seq)
                        new_seq_msgs.append(f"{seq.enzyme_name} - {seq.enzyme_type} created and added to paper")

    return issues, seqs_added, new_seq_msgs

def clear_empty_values(seq_dict):
    keys = list(seq_dict.keys())
    for key in keys:
        if seq_dict[key] == '':
            seq_dict.pop(key)
    return seq_dict

def excel_add_existing_sequences(data_list, paper):
    """
    Add enzymes which match from an excel upload
    """

    user = user_datastore.get_user(current_user.id)
    issues = []
    new_seq_msgs = []
    seqs_added = []

    for seq_dict in data_list:
        already_added = False
        already_owned = False

        # if theres no sequence but there is an accession number, try to look up the sequence
        if seq_dict.get('accession', '') != '' and seq_dict.get('sequence', '') == '':
            accession = seq_dict.get('accession', '')
            seq, loaded_from = lookup_accession.get_sequence(accession)
            if seq != "":
                seq_dict['sequence'] = seq

        # check enzyme name is ok
        if seq_dict.get('enzyme_name','') == '':
            issues.append(f"Sequence must have a name")
        elif any(l in seq_dict.get('enzyme_name', '') for l in INVALID_NAME_CHARS):
            issues.append(f"Sequence cannot contain the following characters: {INVALID_NAME_CHARS}")
        else:
            seq = sequence_queries.seq_obj_from_name(seq_dict['enzyme_name'], include_other_names=True)

            # if no sequence found, can't add
            if seq is None:
                issues.append(f"No sequence with name {seq_dict['enzyme_name']} found")
            else:

                # add seq to paper if not already there
                if paper not in seq.papers:
                    seq.papers.append(paper)
                else:
                    already_added = True

                # does the user own this sequence? if so can update..
                if seq.owner == user or seq.owner is None:
                    if seq.reviewed == False:
                        seq_dict = clear_empty_values(seq_dict)
                        sequence_CUD.update(seq, seq_dict, user)
                else:
                    already_owned = True

                issue = f'{seq.enzyme_name} found'
                if already_owned == True:
                    issue += f' - is owned by another user so no sequence data was updated'
                elif seq.reviewed == False:
                    issue += f' - not reviewed, so <b> sequence data was updated</b>'
                elif seq.reviewed == True:
                    issue += f' - is marked reviewed, so cannot update any fields</b>'

                if already_added == True:
                    issue += f" - sequence was already previously added to this paper."
                else:
                    issue += f" - sequence added to paper"

                issues.append(issue)

                seq.save()
                seqs_added.append(seq)

    return issues, seqs_added, new_seq_msgs
