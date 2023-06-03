from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, request, jsonify, session, current_app, redirect, url_for
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, UniRef50, SSN_record, Sequence
from retrobiocat_web.mongo.models.reaction_models import Reaction
from rq import get_current_job
from retrobiocat_web.analysis.ssn.ssn_main import SSN
from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser
from retrobiocat_web.analysis.ssn.get_pfams import get_pfams, count_pfams, order_pfams_by_count
from retrobiocat_web.app.database_bps.db_analysis.forms import SSN_Form
from retrobiocat_web.analysis.uniprot_and_web import rhea_api
import json
from distutils.util import strtobool
import datetime
from retrobiocat_web.analysis.ssn.ssn_tasks import load_precalculated_vis, load_all_edges
from distutils.util import strtobool
from rq.registry import FinishedJobRegistry, StartedJobRegistry


def get_task_id(enzyme_type, alignment_score, only_biocatdb):
    job_id = f"enzyme_type={enzyme_type}__alignment_score={alignment_score}__only_biocatdb={only_biocatdb}"
    return job_id

# "localhost:5000/ssn?enzyme_type=CAR&alignment_score=40"
@bp.route('/ssn/', methods=['GET'])
def ssn():
    def check_args(args):
        required_args = ['enzyme_type', 'alignment_score']
        allowed_args = required_args + ['only_biocatdb', 'hide_mutants']
        has_required = False
        for arg in args:
            if str(arg) not in allowed_args:
                print(f'Error arg - {arg} -  not allowed')
                return False
            elif str(arg) in required_args:
                has_required = True

        if has_required == False:
            print("missing an argument")
            return False

        return True

    args = request.args.to_dict()

    if check_args(args) == False:
        pass # return error msg

    enzyme_type = args.get('enzyme_type', None)
    alignment_score = int(args.get('alignment_score', 40))
    only_biocatdb = bool(strtobool(args.get('only_biocatdb', 'False')))
    hide_mutants = bool(strtobool(args.get('hide_mutants', 'False')))
    job_id = get_task_id(enzyme_type, alignment_score, only_biocatdb)

    queue = current_app.network_queue
    registry = FinishedJobRegistry(queue=queue)
    started_reg = StartedJobRegistry(queue=queue)
    if job_id in list(registry.get_job_ids()):
        task = queue.fetch_job(job_id)
        if task != None:
            result = task.result
            node_one = result['nodes'][0]
            start_pos = {'x': node_one['x'], 'y': node_one['y']}
            reactions = [c for c in
                         ['-'] + (list(Reaction.objects(enzyme_types=result['enzyme_type']).distinct('name')))]
            return render_template('ssn/ssn.html',
                                   nodes=result['nodes'],
                                   edges=result['edges'],
                                   alignment_score=result['alignment_score'],
                                   pfams=result['pfams'],
                                   pfam_counts=result['pfam_counts'],
                                   rhea=result['rhea'],
                                   rhea_imgs=result['rhea_imgs'],
                                   start_pos=start_pos,
                                   reactions=reactions,
                                   enzyme_type=result['enzyme_type'],
                                   hide_mutants=result['hide_mutants'],
                                   only_biocatdb=result['only_biocatdb'])
        else:
            print("Error ssn task = None")

    elif job_id not in list(started_reg.get_job_ids()):
        queue.enqueue(task_get_ssn, enzyme_type, alignment_score,
                                          hide_mutants, only_biocatdb, job_id=job_id)

    return render_template('ssn/ssn_loading.html', task_id=job_id)


def task_get_ssn(enzyme_type, score, hide_mutants, only_biocatdb):
    job = get_current_job()
    job.meta['progress'] = 'started'
    job.save_meta()

    ssn = SSN(enzyme_type)

    if only_biocatdb is False and str(score) in ssn.db_object.num_at_alignment_score:
        job.meta['progress'] = 'loading pre-positioned ssn'
        job.save_meta()

        nodes = load_precalculated_vis(enzyme_type, str(score))
        edges = [] # sending edges to browser can be slow if network is big
        #edges = load_all_edges(enzyme_type, score)
        #nodes = ssn.db_object.precalculated_vis[str(score)]
        # Need to filter mutants here
    else:
        job.meta['progress'] = 'positioning ssn (slow for large networks)'
        job.save_meta()

        ssn.load(include_mutants=not hide_mutants, only_biocatdb=only_biocatdb)
        vis = SSN_Visualiser(enzyme_type, log_level=1)
        nodes, edges = vis.visualise(ssn, score)
        edges = [] # sending edges to browser can be slow if network is big

    job.meta['progress'] = 'getting pfam and rhea'
    job.save_meta()

    pfams = get_pfams(enzyme_type)
    pfam_counts = count_pfams(enzyme_type, pfams)
    pfams, pfam_counts = order_pfams_by_count(pfams, pfam_counts)

    rheas = rhea_api.get_rheas_for_enzyme_type(enzyme_type)
    rhea_equations, rhea_imgs = rhea_api.get_rhea_dicts(rheas, molSize=(125,90))

    job.meta['progress'] = 'ssn ready'
    job.save_meta()

    result = {'nodes': nodes,
              'edges': edges,
              'pfams': pfams,
              'pfam_counts': pfam_counts,
              'rhea': rhea_equations,
              'rhea_imgs': rhea_imgs,
              'alignment_score': score,
              'enzyme_type': enzyme_type,
              'hide_mutants': hide_mutants,
              'only_biocatdb': only_biocatdb}

    return result

@bp.route("/ssn_status/<task_id>", methods=["GET"])
def ssn_status(task_id):
    task = current_app.network_queue.fetch_job(task_id)

    task_id = task.get_id()
    task_status = task.get_status(refresh=True)
    seconds_since_active = 0
    try:
        print(datetime.datetime.now())
        print(task.last_heartbeat)
        seconds_since_active = (datetime.datetime.now() - task.last_heartbeat).total_seconds()
        print(seconds_since_active)
    except:
        pass
    #if seconds_since_active > 180 and task_status != 'finished':
        #print('Job no longer active')
        #print(task_status)
        #task_status = 'failed'
    #if seconds_since_active > 600:
        #print('Job no longer active')
        #print(task_status)
        #task_status = 'failed'

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

@bp.route('/ssn_form', methods=['GET', 'POST'])
def ssn_form():
    form = SSN_Form()
    form.set_choices()

    task_id = ''

    if 'ssn_task_id' in session:
        old_task_id = session['ssn_task_id']
    else:
        old_task_id = None

    if form.validate_on_submit() == True:
        enzyme_type = form.data['enzyme_type']
        alignment_score = form.data['alignment_score']
        hide_mutants = form.data['hide_mutants']
        only_biocatdb = form.data['only_biocatdb']

        return redirect(url_for('.ssn', enzyme_type=enzyme_type, alignment_score=alignment_score,
                                only_biocatdb=only_biocatdb, hide_mutants=hide_mutants))

    return render_template('ssn/ssn_form.html', form=form, task_id=task_id)

@bp.route("/_ssn_object", methods=["POST"])
def ssn_object():
    enzyme_type = request.form['enzyme_type']
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]
    ssn_obj_q = SSN_record.objects(enzyme_type=enzyme_type_obj)
    #print(f"{len(list(ssn_obj_q))} ssn objects for {enzyme_type}")
    ssn_obj = ssn_obj_q[0]

    num_biocatdb = Sequence.objects(enzyme_type=enzyme_type).count()
    num_uniref = UniRef50.objects(enzyme_type=enzyme_type_obj).count()

    precalc_choices = {}
    for score in ssn_obj.num_at_alignment_score:
        clusters = ssn_obj.num_at_alignment_score[score]
        idt = ssn_obj.identity_at_alignment_score[score]

        choice_text = f"{score}, {clusters} clusters, avg identity {idt[0]} ± {idt[1]}"
        precalc_choices[score] = choice_text

    result = {'status': ssn_obj.status,
              'num_biocatdb': num_biocatdb,
              'num_uniref': num_uniref,
              'precalculated': precalc_choices}
    return jsonify(result=result)



