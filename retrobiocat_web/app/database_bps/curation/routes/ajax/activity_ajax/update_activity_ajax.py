import uuid

from flask import request, jsonify, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.app.app import user_datastore
import json

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import check_permission
from retrobiocat_web.mongo.functions.save_activity_data.process_activity_data import ActivityDataProcessor
from retrobiocat_web.mongo.modal_updates.molecular_descriptors import task_update_molecular_descriptors
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.modal_updates import activity_CUD
from retrobiocat_web.mongo.modal_updates import fingerprints
from retrobiocat_web.mongo.model_queries.sequence_queries import seqs_from_list_of_names
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def filter_for_edited_and_new_data(table_data, edited_rows):
    edited_and_new = []
    for data_dict in table_data:
        if '_id' not in data_dict:
            edited_and_new.append(data_dict)
        elif data_dict['_id'] == '' or data_dict['_id'] in edited_rows:
            edited_and_new.append(data_dict)
    return edited_and_new


def get_sequences(data):
   enzymes = []
   for row in data:
       enzyme_name = row.get('enzyme_name', None)
       if enzyme_name is not None and enzyme_name not in enzymes:
           enzymes.append(enzyme_name)

   seqs = seqs_from_list_of_names(enzymes)
   return seqs


def is_mw_within_acceptable_range(mw):
    if mw < 4000:
        return False
    elif mw > 500000:
        return False
    else:
        return True

def get_mw_dict(data):
    mw_dict = {}
    seqs = get_sequences(data)
    for seq in seqs:
        if seq.sequence is not None and seq.sequence != "":
            try:
                prot_analysis = ProteinAnalysis(seq.sequence)
                mw = round(prot_analysis.molecular_weight(), 2)
                if is_mw_within_acceptable_range(mw):
                    mw_dict[seq.enzyme_name] = mw
            except Exception as e:
                print(f"Error trying to calculate MW for protein {seq.enzyme_name}")
                print(str(e))
    return mw_dict

def add_mws(data):
    mw_dict = get_mw_dict(data)
    print(mw_dict)
    for i, row in enumerate(data):
        if 'km' in row and 'kcat' in row:
            enzyme_name = row.get('enzyme_name', None)
            mw = mw_dict.get(enzyme_name, '')
            data[i]['mw'] = mw
    return data

@bp.route('/_save_activity_data', methods=['GET', 'POST'])
@roles_required('contributor')
def save_activity_data():
    """Update activity data from tabulator"""

    user = user_datastore.get_user(current_user.id)
    paper = paper_queries.paper_from_id(request.form['paper_id'])

    table_data = json.loads(request.form['data'])
    edited_rows = json.loads(request.form['edited_rows'])  # a list of none empty ids which have been edited

    # check the user is actually allowed to edit this data
    if not check_permission.check_paper_permission(current_user.id, paper):
        result = {'status': 'danger',
                  'msg': 'No access to modify this paper',
                  'issues': ['Paper not assigned to user, and not a super_contributor']}
        return jsonify(result=result)

    # first delete ids no long in tabulator table
    ids = [d['_id'] for d in table_data if '_id' in d]
    ids = [i for i in ids if i != '']  # remove blanks
    issues = activity_CUD.delete_other_data(ids, paper)  # deletes any activity for paper not given in list of ids

    # if no issues with deleting, then process new data
    processed_data = []
    if len(issues) == 0:
        edited_and_new = filter_for_edited_and_new_data(table_data, edited_rows)  # only update edited or new rows
        edited_and_new = add_mws(edited_and_new)
        processed_data, issues = ActivityDataProcessor().process_data(edited_and_new, paper)

    # if processing ok, then processed to save the new data
    if len(issues) == 0:
        for data_dict in processed_data:
            update_create_issues = activity_CUD.update_or_create_activity(data_dict, paper, user=user)
            issues += update_create_issues

    # if new data was saved ok, then queue fingerprints update and mol descriptors
    if len(issues) == 0:
        current_app.task_queue.enqueue(fingerprints.task_update_fingerprints, job_id=f"update_molecule_fingerprints_{str(uuid.uuid4())}")
        current_app.task_queue.enqueue(task_update_molecular_descriptors, job_id=f"update_molecule_descriptors_{str(uuid.uuid4())}")

    # return response
    if len(issues) == 0:
        result = {'status': 'success',
                  'msg': 'Activity data updated',
                  'issues': [],
                  'data': processed_data}

    else:
        result = {'status': 'danger',
                  'msg': 'Issues encountered while trying to save activity data..',
                  'issues': issues}
    return jsonify(result=result)



if __name__ == '__main__':
    prot_analysis = ProteinAnalysis('MSKHIGIFGLGAMGTALAAKYLEHGYKTSVWNRTTAKAIPLVEQGAKLASTISEGVNANDLIIICLLNNQVVEDALRDALQTLPSKTIVNLTNGTPNQARKLADFVTSHGARYIHGGIMAVPTMIGSPHAVLLYSGESLELFQSIESHLSLLGMSKYLGTDAGSASLHDLALLSGMYGLFSGFLHAVALIKSGQDTSTTATGLLPLLTPWLSAMTGYLSSIAKQIDDGDYATQGSNLGMQLAGVENIIRAGEEQRVSSQMILPIKALIEQAVGEGHGGEDLSALIEYFKVGKNVD')
    mw = prot_analysis.molecular_weight()
    print(mw)






