from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg, get_rxn_smiles
from flask_security import roles_required, current_user
import datetime
from retrobiocat_web.mongo.models.biocatdb_models import Activity
import mongoengine as db
from rq.job import Job
from rq import get_current_job
from rq.registry import FinishedJobRegistry
from retrobiocat_web.app.database_bps.db_analysis.forms import ScopeForm
from retrobiocat_web.mongo.models.biocatdb_models import Paper
import json
from retrobiocat_web.analysis.data_query import get_data, split_data
from retrobiocat_web.analysis.data_vis.sequence_chemical_space_viewer import SpaceViewer
from bokeh.embed import components
from distutils.util import strtobool

@bp.route('/scope_form', methods=['GET', 'POST'])
def scope_form():
    form = ScopeForm()
    form.set_choices()

    if form.validate_on_submit() == True:
        reaction = form.data['reaction']
        if reaction == '-':
            reaction = None

        enzyme_type = form.data['enzyme_type']
        if enzyme_type == '-':
            enzyme_type = None

        only_reviewed = form.data['only_reviewed']
        mol_type = form.data['molecule']
        remove_negative = form.data['remove_negative']

        return redirect(url_for('db_analysis.scope', reaction=reaction, enzyme_type=enzyme_type,
                                only_reviewed=only_reviewed, mol_type=mol_type, remove_negative=remove_negative))

    return render_template('scope_viewer/scope_viewer_form.html', form=form)

@bp.route("/scope_status/<task_id>", methods=["GET"])
def scope_status(task_id):
    task = current_app.scope_queue.fetch_job(task_id)
    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    task_id = task.get_id()
    task_status = task.get_status(refresh=True)
    seconds_since_active = 0
    try:
        seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
    except:
        pass
    # if seconds_since_active > 180 and task_status != 'finished':
    # print('Job no longer active')
    # print(task_status)
    # task_status = 'failed'
    # if seconds_since_active > 600:
    # print('Job no longer active')
    # print(task_status)
    # task_status = 'failed'

    progress = 'queuing'
    if 'progress' in task.meta:
        progress = task.meta['progress']

    if task:
        response_object = {
            "status": "success",
            "data": {
                "task_id": task_id,
                "task_status": task_status,
                "task_progress": progress
            },
        }
    else:
        response_object = {"status": "error"}

    return jsonify(response_object), 202

def task_make_scope_viewer(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type, remove_negative):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    title = f"Chemical and Sequence space for "
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

    dq = get_data.DataQuery(enzyme_type=enzyme_type,
                            reaction=reaction,
                            only_reviewed=only_reviewed,
                            enzyme_names=enzyme_names,
                            paper_id=paper_id,
                            smi_col=mol_type,
                            remove_negative=remove_negative)

    sv = SpaceViewer(dq)

    layout = sv.create_figure()

    job = get_current_job()
    job.meta['progress'] = 'scope viewer ready'
    job.save_meta()

    script, div = components(layout)

    all_mols = dq.unique_smiles()
    all_seqs = dq.unique_enzymes()

    return {'script': script,
            'all_mols': all_mols,
            'all_seqs': all_seqs,
            'div': div,
            'title': title}

def get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type, remove_negative):
    job_id = f"reaction={reaction}__enzyme_type={enzyme_type}__paperid={paper_id}__" \
             f"enzyme_names={enzyme_names}__moltype={mol_type}__onlyreviewed={only_reviewed}__remove_negative={remove_negative}"
    return job_id

# http://0.0.0.0:5000/heatmap/reaction=Aldehyde%20oxidation__enzyme_type=ADH__paperid=None__enzyme_names=None__moltype=substrate_1_smiles__onlyreviewed=True/
# "localhost:5000/heatmap?enzyme_type=CAR&paper_id=24934239"
@bp.route("/scope/", methods=["GET"])
def scope():
    # options are:  Reaction, Enzyme Type, Paper.  Also specify mol type, but default to product.
    # Must have at least one
    # Also params for splitting the data up

    def check_args(args):
        required_args = ['enzyme_type', 'reaction', 'paper_id', 'enzyme_names']
        allowed_args = required_args + ['only_reviewed', 'mol_type', 'remove_negative']
        has_required = False
        for arg in args:
            if str(arg) not in allowed_args:
                print(f'Error arg - {arg} -  not allowed')
                return False
            elif str(arg) in required_args:
                has_required = True

        if has_required == False:
            print("error need a filter")
            return False

        return True

    args = request.args.to_dict()

    if check_args(args) == False:
        pass # return error msg

    reaction = args.get('reaction', None)
    enzyme_type = args.get('enzyme_type', None)
    enzyme_names = args.get('enzyme_names', None)
    paper_id = args.get('paper_id', None)
    only_reviewed = bool(strtobool(args.get('only_reviewed', 'True')))
    mol_type = args.get('mol_type', 'product_1_smiles')
    remove_negative = args.get('remove_negative', 'False')

    job_id = get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type, remove_negative)

    registry = FinishedJobRegistry(queue=current_app.scope_queue)
    if job_id in list(registry.get_job_ids()):
        task = current_app.scope_queue.fetch_job(job_id)
        if task != None:
            script = task.result['script']
            div = task.result['div']
            title = task.result['title']
            all_mols = task.result['all_mols']
            all_seqs = task.result['all_seqs']
            return render_template('scope_viewer/scope.html', script=script, div=div, title=title,
                                   all_mols=all_mols, all_seqs=all_seqs, args=args,
                                   only_reviewed=only_reviewed, remove_negative=remove_negative,
                                   mol_type=mol_type, paper_id=paper_id, reaction=reaction, enzyme_type=enzyme_type)

    current_app.scope_queue.enqueue(task_make_scope_viewer,
                                      reaction, enzyme_type,
                                      enzyme_names, paper_id,
                                      only_reviewed, mol_type,
                                      remove_negative,
                                      job_id=job_id)

    return render_template('scope_viewer/scope_loading.html', task_id=job_id)
