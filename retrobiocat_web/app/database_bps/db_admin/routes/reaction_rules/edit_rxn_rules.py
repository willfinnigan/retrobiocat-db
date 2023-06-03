from retrobiocat_web.app.database_bps.db_admin import bp
from flask import render_template, jsonify, request
import yaml

from retrobiocat_web.app.database_bps.db_admin.functions.reaction_rule_functions import yaml_conversion
from retrobiocat_web.app.database_bps.db_admin.functions.reaction_rule_functions.reaction_tests import ReactionTester
from retrobiocat_web.app.retrobiocat.functions.get_images import rxntosvg
from flask_security import roles_required
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType
from retrobiocat_web.retro.network_pathway.network import Network

from retrobiocat_web.mongo.models.biocatdb_models import Activity
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile

from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
from retrobiocat_web.retro.retrosynthesis_engine.rule_application.rdchiral.initialization import rdchiralReaction
from retrobiocat_web.retro.visualisation.visualise import Visualiser


def get_reactions():
    reactions = list(Reaction.objects().distinct('name'))
    reactions.sort()
    return reactions

def get_enzymes():
    enzymes = list(EnzymeType.objects().distinct('enzyme_type'))
    enzymes.sort()
    return enzymes

def get_rxn_smarts_yaml(rxn_yaml_dict):
    list_smarts = rxn_yaml_dict['smarts']
    return yaml.dump(list_smarts)

def get_seed_yaml(rxn_yaml_dict, seed_name):
    seeds = rxn_yaml_dict[seed_name]
    yaml_str = yaml.dump(seeds)
    if yaml_str == '[]\n':
        yaml_str = ""
    return yaml_str

def get_rxn_steps_yaml(rxn_yaml_dict):
    list_steps_names = rxn_yaml_dict['steps']
    return yaml.dump(list_steps_names)

def get_enzyme_cofactor_yaml(rxn_yaml_dict):
    enzymes = rxn_yaml_dict['enzymes']
    return yaml.dump(enzymes)

def get_reaction_yaml_dict(reaction_name):
    reaction = Reaction.objects(name=reaction_name).select_related()[0]
    yaml_dict = yaml_conversion.reaction_to_yaml_dict({}, reaction)
    return yaml_dict

def get_positive_tests(rxn_yaml_dict):
    pos = rxn_yaml_dict['positive_tests']
    return yaml.dump(pos)

def get_negative_tests(rxn_yaml_dict):
    neg = rxn_yaml_dict['negative_tests']
    return yaml.dump(neg)

def get_steps_objs(groups_list_rxn_names):
    groups_objs = []
    for list_rxn_names in groups_list_rxn_names:
        rxn_objs = []
        for name in list_rxn_names:
            rxn_q = Reaction.objects(name=name)
            if len(rxn_q) == 1:
                rxn_objs.append(rxn_q[0])
        groups_objs.append(rxn_objs)
    return groups_objs

@bp.route('/rule_editor')
@roles_required('rxn_rules_admin')
def rule_editor():
    args = request.args.to_dict()

    go_to_reaction = args.get('reaction', '')

    reactions = get_reactions()
    enzymes = get_enzymes()
    return render_template('rxn_rule_editor/new_rxn_rules_page.html', reactions=reactions, enzymes=enzymes,
                           go_to_reaction=go_to_reaction)


@bp.route("/_load_rule", methods=["POST"])
@roles_required('rxn_rules_admin')
def load_rule():
    reaction_name = request.form['selected_rule']

    if reaction_name == 'Empty template':
        rxn_name = ""
        rxn_smarts = ""
        rxn_enzyme_cofactor = ""
        example_rxn_smiles = ""
        positive_tests = ""
        negative_tests = ""
        reaction_type = ""
        experimental = True
        two_step = False
        rxn_steps = ""
        requires_absence_of_water = False
        product_seeds = ""
        substrate_1_seeds = ""
        substrate_2_seeds = ""
    else:
        yaml_dict = get_reaction_yaml_dict(reaction_name)

        rxn_name = reaction_name
        rxn_smarts = get_rxn_smarts_yaml(yaml_dict[rxn_name])
        rxn_steps = get_rxn_steps_yaml(yaml_dict[rxn_name])
        rxn_enzyme_cofactor = get_enzyme_cofactor_yaml(yaml_dict[rxn_name])
        positive_tests = get_positive_tests(yaml_dict[rxn_name])
        negative_tests = get_negative_tests(yaml_dict[rxn_name])
        reaction_type = yaml_dict[rxn_name]['type']
        rxn = Reaction.objects(name=reaction_name)[0]
        experimental = rxn.experimental
        example_rxn_smiles = rxn.example_rxn_string
        two_step = rxn.two_step
        requires_absence_of_water = rxn.requires_absence_of_water
        product_seeds = get_seed_yaml(yaml_dict[rxn_name], 'product_seeds')
        substrate_1_seeds = get_seed_yaml(yaml_dict[rxn_name], 'substrate_1_seeds')
        substrate_2_seeds = get_seed_yaml(yaml_dict[rxn_name], 'substrate_2_seeds')

    result = {"rxn_name": rxn_name,
              "rxn_smarts": rxn_smarts,
              "rxn_steps": rxn_steps,
              "example_rxn_smiles": example_rxn_smiles,
              "rxn_enzyme_cofactor": rxn_enzyme_cofactor,
              "positive_tests": positive_tests,
              "negative_tests": negative_tests,
              "reaction_type": reaction_type,
              "experimental": experimental,
              "two_step": two_step,
              "requires_absence_of_water": requires_absence_of_water,
              "product_seeds": product_seeds,
              "substrate_1_seeds": substrate_1_seeds,
              'substrate_2_seeds': substrate_2_seeds}

    print(result)

    return jsonify(result=result)

@bp.route("/_visualise_smarts", methods=["POST"])
def visualise_smarts():
    smarts_yaml = request.form['smarts_yaml']
    try:
        smarts_list = yaml.safe_load(smarts_yaml)
    except:
        result = {'status': 'fail',
                  'msg': 'Could not load yaml'}
        return jsonify(result=result)

    try:
        rxn_list = []
        for sma in smarts_list:
            rxn_list.append(rdchiralReaction(sma))
        list_imgs = rxntosvg(rxn_list, rxnSize=(450, 150))
    except:
        result = {'status': 'fail',
                  'msg': 'Could not load rxn imgs'}
        return jsonify(result=result)

    result = {'status': 'success',
              'msg': '',
              'list_imgs': list_imgs}

    return jsonify(result=result)


@bp.route("/_save_reaction", methods=["POST"])
@roles_required('rxn_rules_admin')
def save_reaction():
    rxn_selection = request.form['rxn_selection']
    rxn_name = request.form['rxn_name']
    smarts_yaml = request.form['smarts_yaml']
    steps_yaml = request.form['steps_yaml']
    cofactors = request.form['cofactors']
    example_rxn_smiles = request.form['example_rxn_smiles']
    positive_tests = request.form['positive_tests']
    negative_tests = request.form['negative_tests']
    rxn_type = request.form['rxn_type']
    experimental = request.form['experimental']
    two_step = request.form['two_step']
    requires_absence_of_water = request.form['requires_absence_of_water']
    product_seeds = request.form['product_seeds']
    substrate_1_seeds = request.form['substrate_1_seeds']
    substrate_2_seeds = request.form['substrate_2_seeds']

    if experimental == 'false':
        experimental = False
    elif experimental == 'true':
        experimental = True

    if two_step == 'false':
        two_step = False
    elif two_step == 'true':
        two_step = True

    if requires_absence_of_water == 'false':
        requires_absence_of_water = False
    elif requires_absence_of_water == 'true':
        requires_absence_of_water = True

    if rxn_selection == 'Empty template':
        if len(Reaction.objects(name=rxn_name)) == 0 and rxn_name != '':
            reaction = Reaction(name=rxn_name)
        else:
            result = {'status': 'failed',
                      'msg': "Could not create new reaction - name already exists",
                      'issues' : []}
            return jsonify(result=result)
    else:
        reaction = Reaction.objects(name=rxn_selection)[0]

    if len(Reaction.objects(name=rxn_name)) == 0:
        reaction.name = rxn_name
    reaction.smarts = yaml.load(smarts_yaml, Loader=yaml.FullLoader)
    reaction.enzyme_types = list(yaml.load(cofactors, Loader=yaml.FullLoader).keys())
    reaction.cofactors = yaml.load(cofactors, Loader=yaml.FullLoader)
    reaction.positive_tests = yaml.load(positive_tests, Loader=yaml.FullLoader)
    reaction.negative_tests = yaml.load(negative_tests, Loader=yaml.FullLoader)
    reaction.type = rxn_type
    reaction.experimental = experimental
    reaction.two_step = two_step
    reaction.requires_absence_of_water = requires_absence_of_water
    reaction.example_rxn_string = example_rxn_smiles
    reaction.product_seeds = yaml.load(product_seeds, Loader=yaml.FullLoader)
    reaction.substrate_1_seeds = yaml.load(substrate_1_seeds, Loader=yaml.FullLoader)
    reaction.substrate_2_seeds = yaml.load(substrate_2_seeds, Loader=yaml.FullLoader)
    print(reaction.product_seeds)
    print(product_seeds)

    step_groups = yaml.load(steps_yaml, Loader=yaml.FullLoader)
    reaction.steps = get_steps_objs(step_groups)

    reaction.save()

    refresh = 'False'
    if rxn_selection != reaction.name:
        refresh = 'True'
        activities = Activity.objects(reaction=rxn_selection)
        for act in activities:
            act.reaction = reaction.name
            act.save()

    print("reaction saved")
    result = {'status': 'info',
              'msg': "Reaction saved",
              'issues': [],
              'refresh': refresh}
    return jsonify(result=result)

@bp.route("/_test_reaction", methods=["POST"])
@roles_required('rxn_rules_admin')
def test_reaction():
    rxn_selection = request.form['rxn_selection']
    rxn_name = request.form['rxn_name']
    smarts_yaml = request.form['smarts_yaml']
    steps_yaml = request.form['steps_yaml']
    example_rxn_smiles = request.form['example_rxn_smiles']
    cofactors = request.form['cofactors']
    positive_tests = request.form['positive_tests']
    negative_tests = request.form['negative_tests']
    rxn_type = request.form['rxn_type']
    experimental = request.form['experimental']
    two_step = request.form['experimental']
    requires_absence_of_water = request.form['requires_absence_of_water']
    product_seeds = request.form['product_seeds']
    substrate_1_seeds = request.form['substrate_1_seeds']
    substrate_2_seeds = request.form['substrate_2_seeds']

    tester = ReactionTester()
    tester.run(rxn_selection, rxn_name, smarts_yaml, cofactors,
               positive_tests, negative_tests, rxn_type, experimental, two_step,
               requires_absence_of_water, example_rxn_smiles, steps_yaml,
               product_seeds, substrate_1_seeds, substrate_2_seeds)

    result = {'status': tester.state,
              'msg': tester.get_msg(),
              'issues': tester.issues}
    return jsonify(result=result)

@bp.route("/_test_product_against_rules", methods=["POST"])
def test_product_against_rules():
    target_smiles = request.form['target_smiles']
    smarts = request.form['smarts']
    steps = request.form['steps']
    mode = request.form['mode']

    target_smiles = rdkit_smile(target_smiles)

    try:
        smarts_list = yaml.load(smarts, Loader=yaml.FullLoader)
        steps_list = yaml.load(steps, Loader=yaml.FullLoader)
    except:
        print('test_product_against_rules failed creating smarts list')
        return jsonify(result={'status': 'fail'})

    # get multi_step
    try:
        multi_step_smarts = []
        for group in steps_list:
            group_smas = []
            for step in group:
                step_smas = []
                reaction = Reaction.objects(name=step)[0]

                for sma in reaction.smarts:
                    step_smas.append(sma)

                group_smas.append(step_smas)
            multi_step_smarts.append(group_smas)

    except:
        print('test_product_against_rules failed creating multi_step smarts list')
        return jsonify(result={'status': 'fail'})

    try:
        network = Network(target_smiles=target_smiles)
        retro = RetrosynthesisEngine(network)
        retro.initialise_graph()
        retro.config.rbc_rule_mode = mode

        if request.form['combine_enantiomers'] == 'true':
            retro.config.combine_enantiomers = True
        else:
            retro.config.combine_enantiomers = False

        retro.step_custom_reactions(target_smiles, smarts_list, multi_step_smarts)

        visualiser = Visualiser()
        nodes, edges = visualiser.nodes_edges(network.graph)

        result = {'status': 'success',
                  'nodes': nodes,
                  'edges': edges}
    except Exception as e:
        print(e)
        print('test_product_against_rules failed creating network')
        result = {'status': 'fail'}

    return jsonify(result=result)


