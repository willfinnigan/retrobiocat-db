import json

from flask import request, jsonify
from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.apply_reaction import \
    fwd_rxn, retro_rxn
from flask_security import roles_required

from retrobiocat_web.app.database_bps.curation.functions.table_reaction_rules_apply import ActivityTable_Rxn_Apply


def products_but_no_substrates(product, substrate_1):
    if (product != '') and (substrate_1 == ''):
        return True
    return False

def substrates_but_no_product(product, substrate_1):
    if product == '' and substrate_1 != '':
        return True
    return False

def retro_result_to_table_item(outcomes, item):
    if len(outcomes) == 1:
        if len(outcomes[0]) == 1:
            item['substrate_1_smiles'] = outcomes[0][0]
        elif len(outcomes[0]) == 2:
            item['substrate_1_smiles'] = outcomes[0][0]
            item['substrate_2_smiles'] = outcomes[0][1]
    return item

def fwd_result_to_table_item(outcomes, item):
    if len(outcomes) == 1:
        if len(outcomes[0]) == 1:
            item['product_1_smiles'] = outcomes[0][0]
    return item

@bp.route('/_apply_reaction_rules_table', methods=['GET', 'POST'])
@roles_required('contributor')
def apply_reaction_rules_table():
    rows = json.loads(request.form['rows'])

    new_rows = ActivityTable_Rxn_Apply().run(rows)

    result = {'status': 'success',
              'msg': 'Reaction rules applied - please check the resulting substrate and product smiles',
              'new_rows': new_rows,
              'issues':
                  ['Reaction rules applied in the forward direction are currently not tested, please check these.']}
    return jsonify(result=result)