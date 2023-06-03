from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import request, jsonify
from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.mongo.models.user_models import User
from retrobiocat_web.app.app import csrf
from retrobiocat_web.mongo.model_queries import enzyme_type_queries, sequence_queries
from retrobiocat_web.retro.network_pathway.network import Network


def check_for_access(user_email, api_key):
    if user_email == None:
        return False, 'No user'
    user = User.objects(email=user_email).first()
    if user == None:
        return False, 'No user'
    if not user.has_role('full_data_access'):
        return False, 'No data access role'
    if user.api_key == None:
        return False, 'No api key'
    if user.api_key != api_key:
        return False, 'Api key incorrect'

    if user.has_role('full_data_access') and (user.api_key == api_key):
        return True, 'OK'

    return False, 'Other issue'

@bp.route("/api/enzyme_types", methods=["POST"])
@csrf.exempt
def enzyme_types():
    """Returns a list of dicts for enzyme_types """

    data = request.json
    user_email = data.get('user_email', None)
    api_key = data.get('api_key', None)
    names_only = data.get('names_only', True)

    access, msg = check_for_access(user_email, api_key)
    if access == False:
        return jsonify(failed=msg)

    et_dict = enzyme_type_queries.enzyme_types_dict(names_only=names_only)

    return jsonify(et_dict)

@bp.route("/api/data_api", methods=["POST"])
@csrf.exempt
def data_api():
    """Returns activity data for a given query"""
    data = request.json
    user_email = data.get('user_email', None)
    api_key = data.get('api_key', None)

    access, msg = check_for_access(user_email, api_key)
    if access == False:
        return jsonify(failed=msg)

    if 'query_json' not in data:
        return jsonify(failed='No query supplied')

    query_dict = data['query_json']
    dq = get_data.DataQuery(enzyme_type=query_dict.get('enzyme_type', None),
                            reaction=query_dict.get('reaction', None),
                            paper_id=query_dict.get('paper_id', None),
                            only_reviewed=query_dict.get('only_reviewed', True),
                            smi_col=query_dict.get('smi_col', 'product_1_smiles'),
                            smiles=query_dict.get('smiles', None),
                            remove_negative=query_dict.get('remove_negative', None))
    return jsonify(dq.activity_data)


@bp.route("/api/enzyme_sequences_api", methods=["POST"])
@csrf.exempt
def enzyme_sequences_api():
    """Returns activity data for a given query"""
    data = request.json
    user_email = data.get('user_email', None)
    api_key = data.get('api_key', None)

    access, msg = check_for_access(user_email, api_key)
    if access == False:
        return jsonify(failed=msg)

    if 'query_json' not in data:
        return jsonify(failed='No query supplied')

    query_dict = data['query_json']
    enzyme_type = query_dict.get('enzyme_type', None)
    enzyme_names = query_dict.get('enzyme_names', None)

    seq_data = sequence_queries.seqs_of_type_and_or_names(enzyme_type, enzyme_names, as_pymongo=True)
    return jsonify(seq_data)

@bp.route("/api/rule_api", methods=["POST"])
@csrf.exempt
def rule_api():
    """Returns result"""
    data = request.json
    user_email = data.get('user_email', None)
    api_key = data.get('api_key', None)

    access, msg = check_for_access(user_email, api_key)
    if access == False:
        return jsonify(failed=msg)

    if 'query_json' not in data:
        return jsonify(failed='No query reaction supplied')

    query_dict = data['query_json']

    reaction = query_dict.get('reaction', None)
    product = query_dict.get('product', '')
    options = query_dict.get('options', {})

    network = Network(target_smiles=product, number_steps=1, include_two_step=True)
    if options is not None:
        network.update_settings(options)

    # if reaction is specified, then limit reactions to only this one
    if reaction is not None:
        network.rxns = {reaction: network.rxns.get(reaction, [])}
        network.multi_step_rxns = {reaction: network.multi_step_rxns.get(reaction, [])}

    network.generate(product, 1)

    return_json = network.get_biocathub_json()
    return jsonify(return_json)



if __name__ == "__main__":
    import requests

    url = "http://0.0.0.0:5000/api/data_api"
    query_json = {'enzyme_type': 'CAR'}
    data = {
        "user_email": 'admin1@email.com',
        "api_key": "wifheoifh23jb34pgbp3i4ijbg3434urh34fyb32iuewiubp4g4g3bbovcnwo43fu4t43ibg43ijb34gibj",
        "query_json": query_json,
        'names_only': True
    }
    response = requests.post(url, json=data)
    #print(response.text)

    url = "http://0.0.0.0:5000/api/enzyme_types"
    response = requests.post(url, json=data)
    #print(response.text)

    url = "http://0.0.0.0:5000/api/rule_api"
    data.update({'query_json': {'product': 'CCCC=O'}})
    response = requests.post(url, json=data)
    print(response.text)





