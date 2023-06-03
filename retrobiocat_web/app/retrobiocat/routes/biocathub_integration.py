from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request, current_app
import networkx as nx
import requests
import json

from retrobiocat_web.app.retrobiocat.functions.load_save_network import load_network_components_from_redis
from retrobiocat_web.retro.network_pathway.network import Network

@bp.route('/_get_biocathub_network_json', methods=['GET', 'POST'])
def get_biocathub_network_json():
    """ Ajax request to load json for biocathub """

    if current_app.config['ALLOW_BIOCATHUB_INTEGRATION'] == False:
        result = {'status': 'danger',
                  'msg': 'Biocathub integration is blocked',
                  'issues': [],
                  'biocathub_json': {}}

        return jsonify(result=result)

    network_id = request.form['network_id']

    # load the network
    redis_data = json.loads(current_app.redis.get(network_id))
    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_data)

    # get biocathub_json
    biocathub_json = network.get_biocathub_json()

    # return the result
    result = {'status': 'success',
              'msg': 'Loaded network json ok',
              'issues': [],
              'biocathub_json': biocathub_json}

    return jsonify(result=result)


@bp.route('/_send_data_to_retrobiohub', methods=['GET', 'POST'])
def send_data_to_retrobiohub_endpoint():
    biocathub_json = json.loads(request.form['biocathub_json'])
    url = 'https://retrobiohub.org'
    print('sending json now..')
    print(biocathub_json)

    status_code, id = retrobiohub_endpoint(url, biocathub_json)
    if status_code == 200:
        print('success')
        status = 'success'
        msg = f"Success - data sent - click link to launch biocathub"
        redirect_url = url + f'/rbh/?name={id}'
    else:
        status = 'danger'
        msg = f'Error - status code {status_code}'
        redirect_url = ""

    return jsonify(result={'status': status,
                            'status_code': status_code,
                            'msg': msg,
                            'issues': [],
                            'redirect_url': redirect_url})


def retrobiohub_endpoint(url, data):
    """
    Sends data to the /retrobiohub post endpoint.
    This should save the data in minimongo and return a json with an id,
    which should be used in the subsequent redirect
    """

    endpoint = url + '/retrobiohub/'
    response = requests.post(endpoint, json=data)
    status_code = response.status_code
    if status_code != 200:
        return status_code, ""
    else:
        response_json = json.loads(json.loads(response.text))['id']
        return status_code, response_json


if __name__ == '__main__':
    #test_url = 'http://127.0.0.1:5000'
    test_url = 'https://retrobiohub.org'

    test_data = [{"substrates": ["CC(N)c1ccccc1"],
                  "products": ["CC(=O)c1ccccc1"],
                  "cofactors": [["Amine"], ["Carbonyl"]],
                  "enzyme": "TA",
                  "reaction": "Secondary amine deamination",
                  }]



    # top one currently causes error
    test_data2 = [{'cofactors': [['NADP'], ['NADPH']], 'enzyme': 'IRED', 'products': ['c1ccc([C@@H]2CCCCN2)cc1'], 'reaction': 'Imine reduction', 'substrates': ['C1=N[C@H](c2ccccc2)CCC1']},
                  {'cofactors': [['NADP'], ['NADPH']], 'enzyme': 'IRED', 'products': ['c1ccc([C@@H]2CCCCN2)cc1'], 'reaction': 'Imine reduction', 'substrates': ['c1ccc(C2=NCCCC2)cc1']}]

    status_code, id = retrobiohub_endpoint(test_url, test_data2)
    print(status_code)
    print(id)

    id_url = test_url + '/rbh/?name=test_id'
    response = requests.get(id_url)
    print(response.text)


