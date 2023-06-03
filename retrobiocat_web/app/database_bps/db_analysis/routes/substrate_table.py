from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, flash, redirect, url_for, request, jsonify, session, current_app
from flask_security import roles_required, current_user
from distutils.util import strtobool
from rq.registry import FinishedJobRegistry, StartedJobRegistry
from retrobiocat_web.analysis.substrate_table import make_substrate_table
from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.mongo.models.biocatdb_models import Paper

def get_queue(current_app):
    return current_app.heatmap_queue

def get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type):

    job_id = f"reaction={reaction}__enzyme_type={enzyme_type}__paperid={paper_id}__" \
             f"enzyme_names={enzyme_names}__moltype={mol_type}__onlyreviewed={only_reviewed}"
    return job_id

@bp.route('/substrate_table', methods=['GET'])
def substrate_table():

    def check_args(args):
        required_args = ['enzyme_type']
        allowed_args = required_args + ['only_reviewed', 'reaction', 'paper_id', 'enzyme_names', 'mol_type']
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
    mol_type = args.get('mol_type', 'substrate_1_smiles')

    job_id = get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type)

    queue = get_queue(current_app)

    registry = FinishedJobRegistry(queue=queue)
    started_reg = StartedJobRegistry(queue=queue)

    if job_id in list(registry.get_job_ids()):
        task = queue.fetch_job(job_id)
        if task != None:
            svg_cats = task.result['svg_cats']
            smi_stats = task.result['smi_stats']
            title = task.result['title']
            return render_template('substrate_table/substrate_table.html',
                                   svg_cats=svg_cats, smi_stats=smi_stats,
                                   title=title, mol_type=mol_type, enzyme_type=enzyme_type)

    elif job_id not in list(started_reg.get_job_ids()):
        queue.enqueue(task_make_substrate_table,
                      reaction, enzyme_type,
                      enzyme_names, paper_id,
                      only_reviewed, mol_type,
                      job_id=job_id)

    return render_template('substrate_table/substrate_table_loading.html', task_id=job_id)

@bp.route("/substrate_table_status/<task_id>", methods=["GET"])
def substrate_table_status(task_id):
    queue = get_queue(current_app)

    task = queue.fetch_job(task_id)

    task_id = task.get_id()
    task_status = task.get_status(refresh=True)

    print(task_id)
    print(task_status)

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

def task_make_substrate_table(reaction, enzyme_type,
                              enzyme_names, paper_id,
                              only_reviewed, mol_type):

    dq = get_data.DataQuery(reaction=reaction, enzyme_type=enzyme_type,
                            enzyme_names=enzyme_names, paper_id=paper_id,
                            only_reviewed=only_reviewed, smi_col=mol_type,
                            remove_negative='All')

    svg_cats, smi_stats = make_substrate_table.get_substrate_table(dq.activity_data, smi_col=mol_type)

    title = f""
    if mol_type == 'product_1_smiles':
        title += "Products for "
    elif mol_type == 'substrate_1_smiles':
        title += 'Substrates for '
    elif mol_type == 'substrate_2_smiles':
        title += 'Second substrates for '

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

    return {'svg_cats': svg_cats,
            'smi_stats': smi_stats,
            'title': title}


