import json

from flask import request, jsonify
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.apply_reaction import \
    retro_rxn, fwd_rxn
from retrobiocat_web.app.database_bps.curation.functions.reaction_name_autodetection import ReactionName_AutoDetector
from retrobiocat_web.mongo.model_queries.reaction_queries import reactions_of_type
from retrobiocat_web.mongo.model_queries.sequence_queries import seq_obj_from_name


def get_reaction_name_from_product(product_smi, possible_reactions):
    product_smi = product_smi.replace('@', '')
    for reaction in possible_reactions:
        outcomes = retro_rxn(reaction.name, product_smi, True, True)
        if len(outcomes) != 0:
            return reaction.name

    return ''

def get_reaction_name_from_substrates(substrate_1, substrate_2, possible_reactions):
    substrate_1 = substrate_1.replace('@', '')
    substrate_2 = substrate_2.replace('@', '')

    for reaction in possible_reactions:
        outcomes = fwd_rxn(substrate_1, substrate_2, reaction.name, True, True)
        if len(outcomes) != 0:
            return reaction.name

    return ''

def get_possible_reactions_from_enzyme_name(enzyme_name):
    seq_obj = seq_obj_from_name(enzyme_name, include_other_names=True)
    reactions = reactions_of_type([seq_obj.enzyme_type])
    return reactions

def get_reaction_name_from_smiles(row_data):
    if row_data.get('reaction', '') != '':
        return ''

    enzyme_name = row_data.get('enzyme_name', '')
    possible_reactions = get_possible_reactions_from_enzyme_name(enzyme_name)

    reaction_name = ''
    if row_data.get('product_1_smiles', '') != '':
        product = row_data['product_1_smiles']
        reaction_name = get_reaction_name_from_product(product, possible_reactions)
    elif row_data.get('substrate_1_smiles', '') != '':
        substrate_1 = row_data.get('substrate_1_smiles', '')
        substrate_2 = row_data.get('substrate_2_smiles', '')
        reaction_name = get_reaction_name_from_substrates(substrate_1, substrate_2, possible_reactions)

    return reaction_name

def combine_rows_by_smis_and_enzyme_types(rows):
    rows_by_product = {}
    rows_by_substrates = {}
    for i, row_data in enumerate(rows):
        enzyme_type = row_data.get('enzyme_type', None)
        product = row_data.get('product_1_smiles', None)
        substrate_1 = row_data.get('substrate_1_smiles', None)
        substrate_2 = row_data.get('substrate_2_smiles', None)

        if product is not None:
            product_enzyme_type = (product, enzyme_type)
            if product_enzyme_type not in rows_by_product:
                rows_by_product[product_enzyme_type] = []
            rows_by_product[product_enzyme_type].append(row_data)

        elif substrate_1 is not None:
            substrates_enzyme_type = (substrate_1, substrate_2, enzyme_type)
            if substrates_enzyme_type not in rows_by_substrates:
                rows_by_substrates[substrates_enzyme_type] = []
            rows_by_substrates[substrates_enzyme_type].append(row_data)

    return rows_by_product, rows_by_substrates

def run_reaction_name_autodetection(rows):
    rows_by_product, rows_by_substrates = combine_rows_by_smis(rows)

    for p_et_tuple, rows in rows_by_product.items():
        product = p_et_tuple[0]
        enzyme_type = p_et_tuple[1]
        reaction_name = get_reaction_name_from_product(product, possible_reactions)


def add_enzyme_types(rows):
    enzyme_names = []
    for i, row_data in enumerate(rows):
        name = row_data.get('enzyme_name', None)
        if name is not None:
            enzyme_names.append(name)

    enzyme_types = {}
    for name in enzyme_names:
        seq_obj = seq_obj_from_name(name, include_other_names=True)
        if seq_obj is not None:
            enzyme_types[name] = seq_obj.enzyme_type

    for i, row_data in enumerate(rows):
        name = row_data.get('enzyme_name', None)
        if name in enzyme_types:
            rows[i]['enzyme_type'] = enzyme_types[name]
    return rows


def get_just_index_and_name(rows):
    new_rows = []
    for i, row_data in enumerate(rows):
        n = row_data.get('n', None)
        reaction = row_data.get('reaction', None)
        if n is None or reaction is None:
            pass
        else:
            new_rows.append({"n": n, 'reaction': reaction})
    return new_rows

@bp.route('/_autodetect_reaction_names', methods=['GET', 'POST'])
@roles_required('contributor')
def autodetect_reaction_names():
    rows = json.loads(request.form['rows'])

    new_rows = ReactionName_AutoDetector().run(rows)
    reaction_rows = get_just_index_and_name(new_rows)

    result = {'status': 'success',
              'msg': 'Auto-detect applied',
              'issues': [],
              'new_rows': reaction_rows}
    return jsonify(result=result)

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    row_data = {'enzyme_name': 'mpCAR', 'substrate_1_smiles': 'CCC(=O)O'}
    name = get_reaction_name_from_smiles(row_data)
    print(name)


