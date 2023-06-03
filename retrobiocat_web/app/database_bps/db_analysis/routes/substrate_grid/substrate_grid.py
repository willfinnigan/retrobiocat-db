import copy
import json

import numpy as np
from flask import current_app, redirect, url_for, request, render_template
from rdkit import Chem, DataStructs
from rq.registry import FinishedJobRegistry, StartedJobRegistry
from scipy.cluster.hierarchy import linkage, dendrogram

from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.analysis.data_query.data_query_from_args import data_query_from_args
from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.database_bps.db_analysis.functions.task_queues import get_queue, get_task_id
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.mongo.models.biocatdb_models import Paper

QUENE_NAME = 'db'

def order_by_agglomerative_clustering(smis):
    fps = get_data.get_substrates_fps(smis)

    if len(smis) <= 1:
        return smis

    size = len(smis)
    hmap = np.empty(shape=(size, size))
    for i, fp_i in enumerate(fps):
        similarities = DataStructs.BulkTanimotoSimilarity(fp_i, list(fps))
        for j, sim in enumerate(similarities):
            hmap[i, j] = sim

    linked = linkage(hmap, 'single')
    results = dendrogram(linked, no_plot=True, distance_sort='ascending')
    sub_numbered = list(map(int, results['ivl']))
    substrates_ordered = [smis[i] for i in sub_numbered]

    return substrates_ordered

def task_make_substrate_grid(args):
    dq = data_query_from_args(args)
    mol_dict = dq.get_mol_stats()
    ordered_smis = order_by_agglomerative_clustering(list(mol_dict.keys()))
    return ordered_smis, mol_dict

def get_title(args):
    title = f"Molecules for "

    reaction = args.get('reaction', None)
    enzyme_type = args.get('enzyme_type', None)
    paper_id = args.get('paper_id', None)
    enzyme_names = args.get('enzyme_names', None)
    mol_type = args.get('mol_type', None)

    if reaction != None:
        title += f"{reaction} reactions "
    if enzyme_type != None:
        title += f"{enzyme_type} enzymes "
    if paper_id != None:
        paper_q = Paper.objects(id=paper_id)
        if len(paper_q) != 0:
            paper = paper_q[0]
            title += f"in {paper.short_citation} "
    if enzyme_names != None:
        title += f" selected enzymes"

    if mol_type == 'product_1_smiles':
        title += " - showing reaction products"
    elif mol_type == 'substrate_1_smiles':
        title += " - showing substrate 1"
    elif mol_type == 'substrate_2_smiles':
        title += " - showing substrate 2"

    return title

@bp.route("/substrate_grid/", methods=["GET"])
def substrate_grid():
    args = request.args.to_dict()
    mol_type = args.get('mol_type', None)
    if mol_type is None or mol_type == '':
        mol_type = args.get('smi_col', 'product_1_smiles')
    args['mol_type'] = mol_type

    job_id = get_task_id(args)
    queue = get_queue(current_app, QUENE_NAME)
    registry = FinishedJobRegistry(queue=queue)
    started_reg = StartedJobRegistry(queue=queue)
    if job_id in list(registry.get_job_ids()):
        task = queue.fetch_job(job_id)
        if task != None:
            ordered_smis, mols = task.result
            svgs = {s: Chem.MolFromSmiles(s) for s in ordered_smis}
            title = get_title(args)
            return render_template('substrate_grid/substrate_grid.html', title=title,
                                   mols=mols, svgs=svgs, smis=ordered_smis, mol_type=mol_type)

    elif job_id not in list(started_reg.get_job_ids()):
        queue.enqueue(task_make_substrate_grid, args, job_id=job_id)

    queue_details, task_details = queue_task_details(job_id, QUENE_NAME)
    return render_template('queue_loading.html', task_queue=QUENE_NAME, task_id=job_id,
                           queue_details=queue_details, task_details='',
                           title='Loading substrates', ajax_timer=1500, refresh_timer=30000)


