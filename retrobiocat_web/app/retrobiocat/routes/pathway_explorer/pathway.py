from flask import render_template, jsonify, request, redirect, url_for
from retrobiocat_web.app.retrobiocat import bp
from rq import get_current_job
import json
from flask import current_app

from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.app.retrobiocat.functions import pathway_packagaing, pathway_evaluation
from retrobiocat_web.app.retrobiocat.functions.load_save_network import save_network
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.app.retrobiocat.new_forms.enzyme_specificity_form import Specificity_Scorer_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.pathway_explorer_form import PathwayExploreForm
from retrobiocat_web.app.retrobiocat.new_forms.retrosynthesis_form import Retrosynthesis_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.vis_form import Visualiser_Config_Form
from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config
from retrobiocat_web.retro.enzyme_identification.specificity_scorer import Specificity_Scorer
from retrobiocat_web.retro.network_pathway.network import Network
from retrobiocat_web.retro.pathway_search.pathway_generation.best_first_search import BFS
from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
from retrobiocat_web.retro.visualisation.config import Visualiser_Config
from retrobiocat_web.retro.visualisation.visualise import Visualiser


@bp.route('/pathway_explorer_form', methods=['GET', 'POST'])
def pathway_explorer_form():
    form = PathwayExploreForm()
    retro_config_form = Retrosynthesis_Config_Form()
    spec_config_form = Specificity_Scorer_Config_Form()
    vis_config_form = Visualiser_Config_Form()

    if form.validate_on_submit() and retro_config_form.validate_on_submit() \
            and spec_config_form.validate_on_submit() and vis_config_form.validate_on_submit():

        task = current_app.pathway_queue.enqueue(task_get_pathways,
                                                 form.data,
                                                 retro_config_form.form_to_config_dict(),
                                                 spec_config_form.form_to_config_dict(),
                                                 vis_config_form.form_to_config_dict())
        task_id = task.get_id()
        return redirect(url_for('retrobiocat.pathway_explorer', task_id=task_id))

    return render_template('pathway_explorer_form/pathway_explorer_form.html',
                           form=form,
                           retro_config_form=retro_config_form,
                           spec_config_form=spec_config_form,
                           vis_config_form=vis_config_form)

@bp.route("/pathway_explorer/<task_id>/", methods=["GET"])
def pathway_explorer(task_id):
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
                               title='Pathway explorer', ajax_timer=3000, refresh_timer=30000)

    #nodes, edges, max_varient = get_visjs_pathway(task_id, 1, 1)
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

#ajax call used by pathway_explorer
@bp.route('/_next_pathway', methods=['GET', 'POST'])
def next_pathway():
    pathway_num = int(request.form['pathway_num'])
    varient_num = int(request.form['varient_num'])
    task_id = request.form['task_id']

    pathway_data = json.loads(current_app.redis.get(f"{task_id}__{pathway_num}"))
    nodes, edges, max_varient = pathway_data[varient_num-1]
    current_app.redis.expire(f"{task_id}__network", 60 * 60)

    result = {"nodes": nodes,
              "edges": edges,
              "max_varient": max_varient}

    jsonified_result = jsonify(result=result) #json.dumps({'result': result})#

    return jsonified_result

def task_get_pathways(form_data, retro_config_attrs, spec_config_attrs, vis_config_attrs):
    job = get_current_job()
    set_progress_bar(job, 40, 'started')

    network = Network(target_smiles=form_data['target_smiles'])

    retro_config = RetrosynthesisConfig().update_from_dict(retro_config_attrs)
    retro_engine = RetrosynthesisEngine(network, config=retro_config)

    retro_engine.generate_network(form_data['number_steps'], max_nodes=400)
    set_progress_bar(job, 60, 'network generated')

    scorer_config = Specificity_Scorer_Config().update_from_dict(spec_config_attrs)
    scorer = Specificity_Scorer(network, config=scorer_config)
    scorer.score()
    set_progress_bar(job, 70, 'network scored')

    vis_config = Visualiser_Config().update_from_dict(vis_config_attrs)
    vis = Visualiser(config=vis_config)

    save_network(job.id, network, retro_config, scorer_config, vis_config)

    bfs = BFS(network=network,
              max_pathways=form_data['max_pathways'],
              max_pathway_length=form_data['number_steps'],
              min_weight=float(form_data['min_weight']),
              print_log=not current_app.config['PRODUCTION'])
    bfs.run()
    pathways = bfs.get_pathways()

    set_progress_bar(job, 80, 'pathways generated')

    pathway_packagaing.package_all_pathways(job.id, pathways)

    pathway_evaluator = pathway_evaluation.run_evaluate_pathways(pathways, [form_data['weight_num_enzymes'],
                                                     form_data['weight_complexity'],
                                                     form_data['weight_starting'],
                                                     form_data['weight_known_enzymes'],
                                                     form_data['weight_diversity']])

    pathway_packagaing.package_evaluated_pathways(pathway_evaluator.df, job.id)
    pathway_packagaing.package_visjs_pathways(job.id)

    set_progress_bar(job, 90, 'pathways scored')

    options = {}

    pathway_settings = {'weight_num_enzymes': form_data['weight_num_enzymes'],
                        'weight_complexity': form_data['weight_complexity'],
                        'weight_starting': form_data['weight_starting'],
                        'weight_known_enzymes': form_data['weight_known_enzymes'],
                        'weight_diversity': form_data['weight_diversity'],
                        'options': options,
                        'num_pathways': len(pathway_evaluator.df.index)}
    current_app.redis.mset({f"{job.id}__pathway_settings": json.dumps(pathway_settings)})
    current_app.redis.expire(job.id, 60 * 60)



