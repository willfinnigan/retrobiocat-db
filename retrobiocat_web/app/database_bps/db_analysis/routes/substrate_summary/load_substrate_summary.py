import json

import networkx as nx
from flask import current_app
from rq import get_current_job
from retrobiocat_web.analysis.ehreact.create_hasse_network import hasse_diagrame, hasse_diagram_to_network, \
    get_seed_list, remove_non_matches, ensure_correct_seed_coverage
from retrobiocat_web.analysis.ehreact.get_substrate_summary import summarise, leafs_under_node
from retrobiocat_web.analysis.ehreact.seed_filtering import filter_smiles_with_smarts
from retrobiocat_web.analysis.data_query.data_query_from_args import data_query_from_args
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar




def task_get_hasse_network_and_summarise(hasse_id, args):
    job = get_current_job()
    set_progress_bar(job, 40, 'started')
    if not current_app.redis.exists(hasse_id):
        dq = data_query_from_args(args)

        seed_list = get_seed_list(dq.reaction_names(), dq.smi_col)
        smis = filter_smiles_with_smarts(dq.unique_smiles(), seed_list)
        seed_list = ensure_correct_seed_coverage(seed_list, smis)
        set_progress_bar(job, 60, 'unique smiles retrieved')

        smis, non_matches, ori_smiles = remove_non_matches(seed_list, smis)

        d = hasse_diagrame(smis, seed_list=seed_list)
        set_progress_bar(job, 70, 'hasse diagram made')

        graph, start_node = hasse_diagram_to_network(d)
        set_progress_bar(job, 75, 'hasse network made')

        dict_of_lists = nx.to_dict_of_lists(graph)

        hasse_data = json.dumps(dict_of_lists)
        data = json.dumps({'hasse_data': hasse_data,
                           'start_node': start_node,
                           'non_matches': json.dumps(non_matches),
                           'ori_smiles': json.dumps(ori_smiles)})
        current_app.redis.mset({hasse_id: data})
        current_app.redis.expire(hasse_id, 5 * 60)  # 5 minutes
    else:
        graph, start_node, non_matches, ori_smiles = load_hasse_network_from_redis(hasse_id)

    set_progress_bar(job, 80, 'hasse diagram ready')
    if 'start_node' in args:
        start_node = args.pop('start_node')
        non_matches = []

    set_progress_bar(job, 85, 'summarising data..')
    leafs, not_leaf_children, not_leaf_counts, not_leaf_cores, expandables, parent_core = summarise(start_node, graph, ori_smiles)
    set_progress_bar(job, 90, 'summary ready')

    if '' in leafs:
        leafs.remove('')

    return {'leafs': leafs,
            'non_matches': non_matches,
            'expandable_not_leafs': expandables,
            'user_start_node': start_node,
            'parent_core': parent_core,
            'parent_smi': start_node,
            'not_leaf_counts': not_leaf_counts,
            'not_leaf_cores': not_leaf_cores,
            'not_leaf_children': not_leaf_children,
            'title': 'Substrate summary',
            'args': args}


def load_hasse_network_from_redis(hasse_id):
    data = json.loads(current_app.redis.get(hasse_id))
    start_node = data['start_node']
    non_matches = json.loads(data['non_matches'])
    ori_smiles = json.loads(data['ori_smiles'])
    graph = nx.from_dict_of_lists(json.loads(data['hasse_data']), create_using=nx.DiGraph)
    current_app.redis.expire(hasse_id, 60 * 60)
    return graph, start_node, non_matches, ori_smiles


def get_examples(graph, start_node):
    leafs = leafs_under_node(start_node, graph)
    return leafs