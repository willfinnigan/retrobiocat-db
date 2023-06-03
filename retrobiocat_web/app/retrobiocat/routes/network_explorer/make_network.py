from retrobiocat_web.app.retrobiocat import bp
from flask import render_template, redirect, url_for
from rq import get_current_job
import networkx as nx
import json
from flask import current_app
import uuid

from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.app.retrobiocat.new_forms.enzyme_specificity_form import Specificity_Scorer_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.network_explorer_form import NetworkExploreForm
from retrobiocat_web.app.retrobiocat.new_forms.retrosynthesis_form import Retrosynthesis_Config_Form
from retrobiocat_web.app.retrobiocat.new_forms.vis_form import Visualiser_Config_Form
from retrobiocat_web.retro.enzyme_identification.config import Specificity_Scorer_Config
from retrobiocat_web.retro.enzyme_identification.specificity_scorer import Specificity_Scorer
from retrobiocat_web.retro.network_pathway.network import Network
import datetime

from retrobiocat_web.retro.retrosynthesis_engine.config import RetrosynthesisConfig
from retrobiocat_web.retro.retrosynthesis_engine.retrosynthesis_engine import RetrosynthesisEngine
from retrobiocat_web.retro.visualisation.config import Visualiser_Config
from retrobiocat_web.retro.visualisation.visualise import Visualiser


@bp.route('/network_explorer_form', methods=['GET', 'POST'])
def network_explorer_form():
    form = NetworkExploreForm()
    retro_config_form = Retrosynthesis_Config_Form()
    spec_config_form = Specificity_Scorer_Config_Form()
    vis_config_form = Visualiser_Config_Form()
    task_id = ''


    if form.validate_on_submit() \
            and retro_config_form.validate_on_submit() \
            and spec_config_form.validate_on_submit() \
            and vis_config_form.validate_on_submit():


        # if no steps, make the network immedietly..
        if form.data['number_steps'] == 0:
            task_id = str(uuid.uuid4())
            make_immediate_network(task_id,
                                   form.data,
                                   retro_config_form.form_to_config_dict(),
                                   spec_config_form.form_to_config_dict(),
                                   vis_config_form.form_to_config_dict())

        # if steps required, create a job
        else:
            task = current_app.network_queue.enqueue(task_make_network,
                                                     form.data,
                                                     retro_config_form.form_to_config_dict(),
                                                     spec_config_form.form_to_config_dict(),
                                                     vis_config_form.form_to_config_dict())
            task_id = task.get_id()

        return redirect(url_for('retrobiocat.network_explorer', task_id=task_id))

    return render_template('network_explorer_form/network_explorer_form.html',
                           form=form,
                           retro_config_form=retro_config_form,
                           spec_config_form=spec_config_form,
                           vis_config_form=vis_config_form,
                           task_id=task_id)


def make_immediate_network(task_id, form_data, retro_config_attrs, spec_config_attrs, vis_config_attrs):
    retro_config = RetrosynthesisConfig().update_from_dict(retro_config_attrs)
    specificity_config = Specificity_Scorer_Config().update_from_dict(spec_config_attrs)

    network = Network(target_smiles=form_data['target_smiles'])

    vis_config = Visualiser_Config().update_from_dict(vis_config_attrs)
    vis = Visualiser(config=vis_config)
    nodes, edges = vis.nodes_edges(network.graph)

    default_network_name = 'Network for ' + str(network.target_smiles)

    result = {'save_id': str(uuid.uuid4()),
              'save_links': [],
              'save_name': default_network_name,
              'nodes': nodes,
              'edges': edges,
              'options': json.dumps({}),
              'graph_dict': json.dumps(nx.to_dict_of_lists(network.graph)),
              'target_smiles': str(network.target_smiles),
              'retrosynthesis_config': json.dumps(retro_config.to_dict()),
              'scorer_config': json.dumps(specificity_config.to_dict()),
              'vis_config': json.dumps(vis_config.to_dict()),
              'attr_dict': json.dumps(network.attributes_dict())}

    current_app.redis.mset({task_id: json.dumps(result)})
    time_to_expire = 15 * 60  # 15 mins * 60 seconds
    current_app.redis.expire(task_id, time_to_expire)

    return result

def task_make_network(form_data, retro_config_attrs, spec_config_attrs, vis_config_attrs):
    job = get_current_job()

    set_progress_bar(job, 40, 'started')

    network = Network(target_smiles=form_data['target_smiles'])

    retro_config = RetrosynthesisConfig().update_from_dict(retro_config_attrs)
    retro_engine = RetrosynthesisEngine(network, config=retro_config)
    retro_engine.generate_network(form_data['number_steps'], max_nodes=int(form_data['max_initial_nodes']))

    set_progress_bar(job, 60, 'network_generated')

    specificity_config = Specificity_Scorer_Config().update_from_dict(spec_config_attrs)
    scorer = Specificity_Scorer(network, config=specificity_config)
    scorer.score()

    set_progress_bar(job, 80, 'scores calculated')

    vis_config = Visualiser_Config().update_from_dict(vis_config_attrs)
    vis = Visualiser(config=vis_config)
    nodes, edges = vis.nodes_edges(network.graph)

    default_network_name = 'Network for ' + str(network.target_smiles)

    result = {'save_id':str(uuid.uuid4()),
              'save_links' : [],
              'save_name' : default_network_name,
              'nodes':nodes,
              'edges':edges,
              'options':json.dumps({}),
              'graph_dict':json.dumps(nx.to_dict_of_lists(network.graph)),
              'target_smiles':str(network.target_smiles),
              'retrosynthesis_config': json.dumps(retro_config.to_dict()),
              'scorer_config': json.dumps(specificity_config.to_dict()),
              'vis_config': json.dumps(vis_config.to_dict()),
              'attr_dict': json.dumps(network.attributes_dict())}

    current_app.redis.mset({job.id: json.dumps(result)})
    time_to_expire = 15*60   #15 mins * 60 seconds
    current_app.redis.expire(job.id, time_to_expire)

    return result


def is_job_still_active(task):
    try:
        seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()

        # hack to deal with time savings
        if seconds_since_active > 3600:
            seconds_since_active -= 3600

    except:
        print("Error determining seconds since job was active")
        return True

    task_status = task.get_status(refresh=True)
    if seconds_since_active > 300 and task_status != 'finished':
        print('Job no longer active')
        return False
    if seconds_since_active > 900:
        return False

