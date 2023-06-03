from flask import render_template,  request, redirect, url_for
from flask import current_app
from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.database_bps.db_analysis.functions.task_queues import get_queue, get_job_registeries
from retrobiocat_web.app.database_bps.db_analysis.routes.substrate_similarity_search.similarity_search_task import task_get_spec_data
from retrobiocat_web.app.database_bps.tables.get_table_data.alt_naming import get_alt_naming
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.retro.network_pathway.rdkit_utils import rdkit_smile
from rq.registry import FinishedJobRegistry, StartedJobRegistry
from distutils.util import strtobool

def get_task_id(enzyme_type, reaction, product, similarity_cutoff, num_choices):
    job_id = f"enzyme_types={enzyme_type}__reaction={reaction}__product={product}__" \
             f"similarity_cuttoff={similarity_cutoff}__num_choices={num_choices}"

    return job_id

def title_from_args(enzyme_type, reaction, product, similarity_cuttoff):
    title = "Substrate specificity query"
    if enzyme_type != 'All' and enzyme_type != '':
        title += f' {enzyme_type} enzymes'
    if reaction != 'All' and reaction != '':
        title += f' {reaction} reaction'
    if product != '':
        title += f' similarity search for {product} with a cutoff of {similarity_cuttoff}'
    return title

@bp.route("/substrate_specificity/", methods=["GET"])
def substrate_specificity():

    args = request.args.to_dict()

    enzyme_type = args.get('enzyme_type', '')
    reaction = args.get('reaction', '')
    target_smiles = args.get('target_smiles', '')

    if target_smiles != '':
        target_smiles = rdkit_smile(target_smiles)

    similarity_cuttoff = float(args.get('similarity', 0.5))
    num_choices = int(args.get('num_choices', 10))

    job_id = get_task_id(enzyme_type, reaction, target_smiles, similarity_cuttoff, num_choices)

    queue_name = 'tasks'
    queue = get_queue(current_app, queue_name)
    finished, started = get_job_registeries(queue)

    # if the job is finished, render the page
    if job_id in list(finished.get_job_ids()):
        task = queue.fetch_job(job_id)
        if task != None:
            activity_data = task.result
            title = title_from_args(enzyme_type, reaction, target_smiles, similarity_cuttoff)
            enzyme_alt_names_to_use = get_alt_naming(args)

            return render_template('substrate_similarity_search/table_result_specificity.html',
                                    substrate_specificity_data=activity_data, title=title,
                                    alt_names=enzyme_alt_names_to_use)

    # if there isn't a job queued or running, then make one
    if job_id not in list(started.get_job_ids()):  # otherwise, if the job is not already queued, add it
        queue.enqueue(task_get_spec_data, enzyme_type, reaction, target_smiles, similarity_cuttoff, num_choices, job_id=job_id)

    # render the loading page
    queue_details, task_details = queue_task_details(job_id, queue_name)
    return render_template('queue_loading.html', task_queue=queue_name, task_id=job_id,
                           queue_details=queue_details, task_details='',
                           title="Substrate specificity query", ajax_timer=2000, refresh_timer=30000)




