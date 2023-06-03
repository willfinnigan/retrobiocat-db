import json

import networkx as nx
from flask import current_app

from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config
from retrobiocat_web.retro.enzyme_identification.specificity_scorer import Specificity_Scorer
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
from retrobiocat_web.retro.visualisation.config import Visualiser_Config
from retrobiocat_web.retro.visualisation.visualise import Visualiser


def load_network_components_from_redis(redis_data):
    graph_dict = json.loads(redis_data['graph_dict'])
    attr_dict = json.loads(redis_data['attr_dict'])
    target_smiles = redis_data['target_smiles']
    retrosynthesis_config_dict = json.loads(redis_data['retrosynthesis_config'])
    scorer_config_dict = json.loads(redis_data['scorer_config'])
    vis_config_dict = json.loads(redis_data['vis_config'])

    graph = nx.from_dict_of_lists(graph_dict, create_using=nx.DiGraph)
    network = Network(graph=graph, target_smiles=target_smiles)
    network.add_attributes(attr_dict)

    retro_config = RetrosynthesisConfig().update_from_dict(retrosynthesis_config_dict)
    retro_engine = RetrosynthesisEngine(network, config=retro_config)

    specificity_config = Specificity_Scorer_Config().update_from_dict(scorer_config_dict)
    scorer = Specificity_Scorer(network, config=specificity_config)
    scorer.score()

    vis_config = Visualiser_Config().update_from_dict(vis_config_dict)
    visualiser = Visualiser(config=vis_config)

    return network, retro_engine, scorer, visualiser

def update_retrosynthesis_config_for_reaction_settings(update_dict, redis_data):
    retrosynthesis_config_dict = json.loads(redis_data['retrosynthesis_config'])
    settings_updated = False
    for key, value in update_dict.items():
        if retrosynthesis_config_dict.get(key, '') != value:
            retrosynthesis_config_dict[key] = value
            settings_updated = True

    if settings_updated:
        redis_data['retrosynthesis_config'] = json.dumps(retrosynthesis_config_dict)

    return redis_data

def update_vis_options(update_dict, redis_data):
    vis_config_dict = json.loads(redis_data['vis_config'])
    settings_updated = False
    for key, value in update_dict.items():
        if vis_config_dict.get(key, '') != value:
            vis_config_dict[key] = value
            settings_updated = True

    if settings_updated:
        redis_data['vis_config'] = json.dumps(vis_config_dict)

    return redis_data

def save_network(job_id, network, retro_config, scorer_config, vis_config):
    network_data = {'graph_dict': json.dumps(nx.to_dict_of_lists(network.graph)),
                    'target_smiles': str(network.target_smiles),
                    'retrosynthesis_config': json.dumps(retro_config.to_dict()),
                    'scorer_config': json.dumps(scorer_config.to_dict()),
                    'vis_config': json.dumps(vis_config.to_dict()),
                    'attr_dict': json.dumps(network.attributes_dict())}
    current_app.redis.mset({f"{job_id}__network": json.dumps(network_data)})
    current_app.redis.expire(f"{job_id}__network", 60 * 60)








