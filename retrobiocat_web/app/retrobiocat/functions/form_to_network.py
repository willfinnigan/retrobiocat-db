from flask import current_app

from retrobiocat_web.retro.network_pathway.network import Network
import json
import networkx as nx

def initialise_network_with_form_data(form_data):
    network = Network(target_smiles=form_data['target_smiles'])
    update_network_options_from_form_data(form_data, None)

    return network

def update_retrosynthesis_config(form_data, config):
    # rule application / processing
    config.retrosynthesis_config.combine_enantiomers = bool(form_data['combine_enantiomers'])
    config.retrosynthesis_config.remove_small_mols = bool(form_data['remove_small'])

    # retrobiocat
    config.retrosynthesis_config.include_experimental = bool(form_data['include_experimental'])
    config.retrosynthesis_config.include_two_step = bool(form_data['include_two_step'])
    config.retrosynthesis_config.include_requires_absence_of_water = bool(form_data['include_requires_absence_of_water'])


def update_network_options_from_form_data(form_data, network):
    update_settings = ['similarity_score_threshold',
                       'colour_reactions',
                       "source_mols_vendors",
                       'source_mols_db_mode',
                       "colour_arrows",
                       "show_negative_enzymes",
                       "only_postitive_enzyme_data",
                       'only_reviewed_activity_data',
                       'specificity_score_substrates']

    if 'max_nodes' in form_data:
        form_data['max_nodes'] = int(form_data['max_nodes'])

    if 'only_reviewed_activity_data' in form_data:
        form_data['only_reviewed_activity_data'] = bool(form_data['only_reviewed_activity_data'])

    update_dict = {}
    for setting_name in update_settings:
        if setting_name in list(form_data.keys()):
            update_dict[setting_name] = form_data[setting_name]
    network.settings.update(update_dict)

def save_network(network, job_id):
    network_data = {'graph_dict': json.dumps(nx.to_dict_of_lists(network.graph)),
                    'target_smiles': str(network.target_smiles),
                    'network_options': json.dumps(network.settings),
                    'attr_dict': json.dumps(network.attributes_dict())}
    current_app.redis.mset({f"{job_id}__network": json.dumps(network_data)})
    current_app.redis.expire(f"{job_id}__network", 60 * 60)