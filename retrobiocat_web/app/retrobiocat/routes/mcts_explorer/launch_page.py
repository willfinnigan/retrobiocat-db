from flask import render_template, redirect, url_for
from retrobiocat_web.app.retrobiocat import bp
from flask import current_app

from retrobiocat_web.app.retrobiocat.functions.load_save_network import save_network
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.app.retrobiocat.new_forms.enzyme_specificity_form import Specificity_Scorer_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.mcts_form import MCTSExploreForm
from retrobiocat_web.app.retrobiocat.new_forms.retrosynthesis_form import Retrosynthesis_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.vis_form import Visualiser_Config_Form
from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config
from retrobiocat_web.retro.enzyme_identification.specificity_scorer import Specificity_Scorer
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.mcts.config import MCTS_Config
from retrobiocat_web.retro.pathway_search.mcts.mcts import MCTS
import json
from rq import get_current_job
from retrobiocat_web.app.retrobiocat.functions import pathway_packagaing
from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.visualisation.config import Visualiser_Config


@bp.route('/mcts_explorer_form', methods=['GET', 'POST'])
def mcts_explorer_form():
    form = MCTSExploreForm()
    retro_config_form = Retrosynthesis_Config_Form()
    spec_config_form = Specificity_Scorer_Config_Form()
    vis_config_form = Visualiser_Config_Form()

    if form.validate_on_submit() and retro_config_form.validate_on_submit() \
            and spec_config_form.validate_on_submit() and vis_config_form.validate_on_submit():

        if current_app.config['PRODUCTION'] == False:
            log_level = 'DEBUG'
        else:
            log_level = 'WARNING'


        task = current_app.pathway_queue.enqueue(initial_search_job,
                                                 form.form_to_config_dict(),
                                                 retro_config_form.form_to_config_dict(),
                                                 spec_config_form.form_to_config_dict(),
                                                 vis_config_form.form_to_config_dict(),
                                                 log_level
                                                 )
        task_id = task.get_id()
        return redirect(url_for('retrobiocat.mcts_explorer', task_id=task_id))

    return render_template('mcts_explorer/mcts_explorer_form.html',
                           form=form,
                           retro_config_form=retro_config_form,
                           spec_config_form=spec_config_form,
                           vis_config_form=vis_config_form)

def initial_search_job(form_data, retro_config_attrs, spec_config_attrs, vis_config_attrs, log_level):
    job = get_current_job()

    network = Network(target_smiles=form_data['target_smiles'])

    retro_config = RetrosynthesisConfig().update_from_dict(retro_config_attrs)
    mcts_config = MCTS_Config().update_from_dict(form_data)
    vis_config = Visualiser_Config().update_from_dict(vis_config_attrs)
    scorer_config = Specificity_Scorer_Config().update_from_dict(spec_config_attrs)

    mcts = MCTS(network, config=mcts_config, retro_config=retro_config, log_level=log_level)

    set_progress_bar(job, 50, f'running ({mcts.config.max_search_time}s)')
    mcts.run()
    set_progress_bar(job, 75, f'MCTS complete, adding enzyme information')

    scorer = Specificity_Scorer(mcts.network, config=scorer_config)
    scorer.score()

    set_progress_bar(job, 75, f'MCTS and scoring complete, generating pathways')

    pathways = mcts.get_pathways()
    print(f"NUM PATHWAYS = {len(pathways)}")

    clusters = mcts.cluster_and_cost(pathways)

    save_network(job.id, mcts.network, retro_config, scorer_config, vis_config)
    pathway_packagaing.package_all_pathways(job.id, pathways)
    pathway_packagaing.package_clustered_pathways(clusters, job.id)
    pathway_packagaing.package_visjs_pathways(job.id)

    set_progress_bar(job, 90, f'Pathways scored')

    pathway_settings = {'weight_num_enzymes': 1,
                        'weight_complexity': 1,
                        'weight_starting': 1,
                        'weight_known_enzymes': 1,
                        'weight_diversity': 1,
                        'options': {},
                        'num_pathways': len(clusters)}


    current_app.redis.mset({f"{job.id}__pathway_settings": json.dumps(pathway_settings)})
    current_app.redis.expire(job.id, 60 * 60)


@bp.route('/mcts_explorer/<task_id>/', methods=['GET'])
def mcts_explorer(task_id):
    queue_name = 'pathway'

    # if task exists but is not finished, go to loading screen
    task = current_app.pathway_queue.fetch_job(task_id)

    if not task:
        task_id = 'task_not_found'
        task_status = 'task_not_found'
        queue_details = 'Error - task not found'
        task_details = 'Error - task not found'
    else:
        task_id = task.get_id()
        task_status = task.get_status(refresh=True)
        queue_details, task_details = queue_task_details(task_id, queue_name)

    if task_status != 'finished':
        return render_template('queue_loading.html', task_queue=queue_name, task_id=task_id,
                               queue_details=queue_details, task_details=task_details,
                               title='Loading substrate summary', ajax_timer=3000, refresh_timer=30000)

    # nodes, edges, max_varient = get_visjs_pathway(task_id, 1, 1)
    pathway_settings = json.loads(current_app.redis.get(task_id + '__pathway_settings'))
    num_pathways = pathway_settings['num_pathways']
    pathway_data = json.loads(current_app.redis.get(f"{task_id}__1"))
    nodes, edges, max_varient = pathway_data[0]

    return render_template('pathway_explorer/pathway_explorer.html',
                           nodes=nodes,
                           edges=edges,
                           max_varient=max_varient,
                           options=pathway_settings['options'],
                           weight_complexity=pathway_settings['weight_complexity'],
                           weight_num_enzymes=pathway_settings['weight_num_enzymes'],
                           weight_starting=pathway_settings['weight_starting'],
                           weight_known_enzymes=pathway_settings['weight_known_enzymes'],
                           weight_diversity=pathway_settings['weight_diversity'],
                           num_pathways=num_pathways,
                           task_id=task_id)

