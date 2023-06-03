
from flask import request, jsonify

from retrobiocat_web.app.database_bps.curation import bp

from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.apply_reaction import \
    fwd_rxn, retro_rxn
from retrobiocat_web.app.database_bps.curation.functions.apply_reactions_during_activity_curation.get_svgs import \
    retro_outcome_svgs, fwd_outcome_svgs
from retrobiocat_web.mongo.models.reaction_models import Reaction
from flask_security import roles_required
from distutils.util import strtobool


def reaction_not_found_error(reaction_name):
    result = {'status': 'danger',
              'msg': 'No reaction on RetroBioCat with this name',
              'issues': [f"{reaction_name} was not found in RetroBioCat"]}
    return jsonify(result=result)

def no_products_error():
    result = {'status': 'danger',
              'msg': 'No products',
              'issues': []}
    return jsonify(result=result)

def reaction_applied_result(outcomes, reaction_svgs, svg_dict):
    result = {'status': 'success',
              'msg': '',
              'products': outcomes,
              'reaction_svgs': reaction_svgs,
              'svg_dict': svg_dict,
              'issues': []}
    return jsonify(result=result)

@bp.route('/_apply_reaction_fwd', methods=['GET', 'POST'])
@roles_required('contributor')
def apply_reaction_fwd():
    # get reaction
    rxn_name = request.form['reaction']
    if Reaction.objects(name=rxn_name).first() is None:
        reaction_not_found_error(rxn_name)

    # get substrates
    substrate_1 = request.form['s1']
    substrate_2 = request.form['s2']
    if bool(strtobool(request.form['no_stereo'])):
        substrate_1 = substrate_1.replace('@', '')
        substrate_2 = substrate_2.replace('@', '')

    # apply reaction
    combine_enantiomers = bool(strtobool(request.form['combine_enantiomers']))
    use_rdkit = bool(strtobool(request.form['use_rdkit']))
    outcomes = fwd_rxn(substrate_1, substrate_2, rxn_name, combine_enantiomers, use_rdkit)
    if len(outcomes) == 0:
        return no_products_error()

    # return result
    svg_dict, reaction_svgs = fwd_outcome_svgs(substrate_1, substrate_2, outcomes)
    return reaction_applied_result(outcomes, reaction_svgs, svg_dict)

@bp.route('/_apply_reaction_rev', methods=['GET', 'POST'])
@roles_required('contributor')
def apply_reaction_rev():

    # get reaction
    rxn_name = request.form['reaction']
    if Reaction.objects(name=rxn_name).first() is None:
        reaction_not_found_error(rxn_name)

    # get product
    product = request.form['p']
    if bool(strtobool(request.form['no_stereo'])):
        product = product.replace('@', '')

    # apply reaction
    combine_enantiomers = bool(strtobool(request.form['combine_enantiomers']))
    use_rdkit = bool(strtobool(request.form['use_rdkit']))
    outcomes = retro_rxn(rxn_name, product, combine_enantiomers, use_rdkit)
    if len(outcomes) == 0:
        return no_products_error()

    # return result
    svg_dict, reaction_svgs = retro_outcome_svgs(product, outcomes)
    return reaction_applied_result(outcomes, reaction_svgs, svg_dict)


