import json
from flask import current_app

from retrobiocat_web.app.retrobiocat.functions.load_save_network import load_network_components_from_redis
from retrobiocat_web.retro.network_pathway.pathway.pathway import Pathway


def package_evaluated_pathways(pathway_evaluator_df, task_id):
    evaluated_pathways = []
    for index, row in pathway_evaluator_df.iterrows():
        pathway_data = [row['Pathway'].list_nodes]
        for pathway in row['Pathway'].other_varients:
            pathway_data.append(pathway.list_nodes)
        evaluated_pathways.append(pathway_data)

    # a list of lists of pathways
    current_app.redis.mset({f"{task_id}__evaluated_pathways": json.dumps(evaluated_pathways)})
    current_app.redis.expire(f"{task_id}__evaluated_pathways", 60 * 60)


def package_clustered_pathways(clustered_pathways, task_id):
    to_save = []
    for cluster in clustered_pathways:
        save_group = [p.list_nodes for p in cluster]
        to_save.append(save_group)

    current_app.redis.mset({f"{task_id}__evaluated_pathways": json.dumps(to_save)})
    current_app.redis.expire(f"{task_id}__evaluated_pathways", 60 * 60)

def package_all_pathways(task_id, pathways):
    all_pathways = []
    all_scores = []
    for pathway in pathways:
        all_pathways.append(pathway.list_nodes)
        all_scores.append(pathway.scores.scores_dict())

    current_app.redis.mset({f"{task_id}__all_pathways": json.dumps((all_pathways, all_scores))})
    current_app.redis.expire(f"{task_id}__all_pathways", 60 * 60)

def package_visjs_pathways(task_id, max_vis=100):

    redis_network_data = json.loads(current_app.redis.get(task_id + '__network'))
    network, retro_engine, scorer, visualiser = load_network_components_from_redis(redis_network_data)
    evaluated_pathways = json.loads(current_app.redis.get(f"{task_id}__evaluated_pathways"))   # a list of lists of pathways

    for i, pathway_varients in enumerate(evaluated_pathways):
        if i > max_vis:
            break
        pathway_vis_js_data = []
        max_var = len(pathway_varients)
        for nodes in pathway_varients:
            pathway = Pathway(nodes, network, calc_scores=False)
            nodes, edges = visualiser.nodes_edges(pathway.sub_graph)
            pathway_vis_js_data.append((nodes, edges, max_var))
        current_app.redis.mset({f"{task_id}__{i+1}": json.dumps(pathway_vis_js_data)})
        current_app.redis.expire(f"{task_id}__{i+1}", 60 * 60)

def get_visjs_pathway(task_id, pathway_id, varient):
    network_redis_data = json.loads(current_app.redis.get(task_id + '__network'))

    network, retro_engine, scorer, visualiser = load_network_components_from_redis(network_redis_data)

    pathway_redis_data = json.loads(current_app.redis.get(task_id + f'__{pathway_id}'))
    pathway_nodes = pathway_redis_data[varient-1]

    pathway = Pathway(pathway_nodes, network, calc_scores=False)

    nodes, edges = visualiser.nodes_edges(pathway.sub_graph)
    max_var = len(pathway_redis_data)

    return nodes, edges, max_var