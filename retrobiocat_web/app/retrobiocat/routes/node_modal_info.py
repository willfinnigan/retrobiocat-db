from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg, get_rxn_smiles
from retrobiocat_web.app.retrobiocat import bp
from flask import jsonify, request
import json
from flask import current_app
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.analysis.drawing import images
from retrobiocat_web.app.retrobiocat.functions import pubchem_funcs


def format_enzyme_info(info_dict):
    cols_to_ignore = ['smiles_reaction', 'paper_id', 'activity_id']

    if info_dict == False:
        return ''

    formatted = ''

    for key, value in info_dict.items():
        if key not in cols_to_ignore:

            if key == 'DOI':
                formatted += f"{key}: <a href='{value}' target='_blank'>{value}</a> <br>"
            else:
                formatted += (str(key) + ': ' + str(value) + '<br>')

        # formatted += '<hr>'

    return formatted

@bp.route('/_get_top_biocatdb_hits', methods=['GET', 'POST'])
def get_top_biocatdb_hits():
    reaction_node = request.form['reaction_node']
    network_id = request.form['network_id']
    substrates = json.loads(request.form['parents'])
    products = json.loads(request.form['children'])
    label = request.form['label']
    enzyme = request.form['enzyme']

    reaction_smiles = f"{substrates[0]}"
    if len(substrates) > 1:
        reaction_smiles += f".{substrates[1]}"
    reaction_smiles += f">>{products[0]}"
    query_reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=(600, 100))

    data = json.loads(current_app.redis.get(network_id))
    attr_dict = json.loads(data['attr_dict'])
    print(attr_dict)

    if enzyme == 'selected_enzyme':
        enzyme = attr_dict[reaction_node]['metadata']['selected_enzyme']
    node_info = attr_dict[reaction_node]['metadata']['enzyme_info'][enzyme]

    if node_info is False:
        print('No similar reactions found')
        result = {'node_info': '',
                  'product_keys': [],
                  'query_reaction_svg': query_reaction_svg,
                  'reaction_name': label,
                  'enzyme_name': enzyme}
        return jsonify(result=result)

    else:

        for product_key in node_info:
            node_info[product_key]['formatted_info'] = format_enzyme_info(node_info[product_key])
        

        for product_key in node_info:
            try:
                node_info[product_key]['reaction_svg'] = smiles_rxn_to_svg(node_info[product_key]['smiles_reaction'], rxnSize=(400,75))
            except Exception as e:
                print(str(e))
                node_info[product_key]['reaction_svg'] = ''


        result = {'node_info': node_info,
                  'product_keys': sorted(list(node_info.keys()), reverse=True),
                  'query_reaction_svg': query_reaction_svg,
                  'reaction_name': label,
                  'enzyme_name': enzyme}

        return jsonify(result=result)

@bp.route('/_get_reaction_svg', methods=['GET', 'POST'])
def get_reaction_svg():
    substrates = json.loads(request.form['parents'])
    products = json.loads(request.form['children'])
    label = request.form['label']

    reaction_smiles = get_rxn_smiles(substrates, products)
    query_reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=(600, 100))

    result = {'query_reaction_svg': query_reaction_svg,
              'reaction_name': label}
    return jsonify(result=result)


@bp.route('/_get_possible_enzymes', methods=['GET', 'POST'])
def get_possible_enzymes():
    reaction_node = request.form['reaction_node']
    network_id = request.form['network_id']

    data = json.loads(current_app.redis.get(network_id))
    attr_dict = json.loads(data['attr_dict'])


    if  attr_dict[reaction_node]['reaction_type'] != 'retrobiocat':
        result = {'possible_enzymes': [],
                  'selected_enzyme': ''}
        return jsonify(result=result)

    enzyme = attr_dict[reaction_node]['metadata']['selected_enzyme']
    possible_enzymes = attr_dict[reaction_node]['metadata']['possible_enzymes']

    if '' in possible_enzymes:
        possible_enzymes.remove('')

    choices = []
    for enz in possible_enzymes:
        enz_full = EnzymeType.objects(enzyme_type=enz)[0].full_name
        choices.append((f"{enz}", f"{enz} - {enz_full}"))

    result = {'possible_enzymes': choices,
              'selected_enzyme': enzyme}

    return jsonify(result=result)

@bp.route('/_get_buyable_info', methods=['GET', 'POST'])
def get_buyable_info():
    substrate = request.form['substrate']
    network_id = request.form['network_id']
    data = json.loads(current_app.redis.get(network_id))
    attr_dict = json.loads(data['attr_dict'])
    available = attr_dict[substrate]['is_starting_material']
    buyable_info = attr_dict[substrate]['starting_material_info']

    return jsonify(result={'buyable_info': buyable_info,
                           'available': available})

@bp.route('/_get_substrate_img', methods=['GET', 'POST'])
def get_substrate_img():
    substrate = request.form['substrate']
    img = images.smitosvg_url(substrate)
    result = {'substrate_image': img}

    return jsonify(result=result)

@bp.route('/_get_pubchem_cid', methods=['GET', 'POST'])
def get_pubchem_cid():
    substrate = request.form['substrate']

    compound = None
    if current_app.config['ALLOW_PUBCHEM_LOOKUP'] == True:
        compound = pubchem_funcs.get_pubchem_compound_from_smiles(substrate)
    if compound != None:
        name = compound.iupac_name
        cid = compound.cid
    else:
        name = 'Not found'
        cid = ''

    result = {'pubchem_cid': cid,
              'substrate_name': name}

    return jsonify(result=result)


