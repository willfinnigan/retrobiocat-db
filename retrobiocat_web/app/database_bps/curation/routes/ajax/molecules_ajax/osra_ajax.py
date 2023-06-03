import json

from flask import request, current_app, jsonify, session, flash
from flask_security import roles_required
from rq import get_current_job
from werkzeug.utils import secure_filename

from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions import structure_recognition
from retrobiocat_web.app.database_bps.curation.functions.mol_obj_to_json import get_new_mol_json
from retrobiocat_web.app.database_bps.curation.routes.ajax.molecules_ajax.error_msgs import not_post_request, \
    error_processing_mol_images
from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.mongo.modal_updates import molecules_CUD
from retrobiocat_web.mongo.model_queries import paper_queries
from retrobiocat_web.mongo.model_queries.molecule_queries import get_num_activity_mols_in_paper, get_molecules_in_paper
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id

def get_saved_files(items):
    issues, saved_files = [], []
    for name, file in items:
        filename = secure_filename(file.filename)
        saved = structure_recognition.save_uploaded_image(file, filename)
        if saved == True:
            saved_files.append(filename)
        else:
            issues.append(f'error saving {filename} - ensure correct filetype')
    return issues, saved_files

@bp.route('/_upload_molecule_images',methods=['GET', 'POST'])
@roles_required('contributor')
def upload_molecule_images():
    if request.method != 'POST':
        return not_post_request()

    issues, saved_files = get_saved_files(request.files.items())
    paper_id = request.form['mol_paper_id_field']
    paper = paper_queries.paper_from_id(paper_id, get_related=True)

    task = current_app.osra_queue.enqueue(task_process_images, saved_files, paper_id)
    task_id = task.get_id()

    if len(issues) != 0:
        return error_processing_mol_images(issues)

    result = {'status': 'success',
              'msg': f'Uploaded files {len(saved_files)} - now processing',
              'issues': [],
              'task_id': task_id}
    return jsonify(result=result)

def task_process_images(list_filenames, paper_id):
    osra_url = current_app.config['OSRA_API_HOST']
    paper = paper_from_id(paper_id)

    paper_mols = get_molecules_in_paper(paper)
    paper_smis = [mol.smi for mol in paper_mols]

    job = get_current_job()

    set_progress_bar(job, 40, 'started')

    list_filenames.reverse()
    progress_increment = int(50 / len(list_filenames))
    identified_smiles = []
    for i, filename in enumerate(list_filenames):
        new_smis = structure_recognition.process_image(filename, osra_url, log=True)
        #new_smis = ['CCCC=O', 'CCCCN']  # used for testing without actually doing osr
        identified_smiles += new_smis
        set_progress_bar(job, 40 + (progress_increment * i + 1), f'processed {i + 1} of {len(list_filenames)}')

    already_added, not_added = [], []
    for smi in identified_smiles:
        if smi in paper_smis:
            already_added.append(smi)
        else:
            not_added.append(smi)

    # save to redis
    save_id = f"{paper_id}_osr_smiles"
    if current_app.redis.exists(save_id):
        redis_content = json.loads(current_app.redis.get(save_id))
    else:
        redis_content = []

    redis_content += not_added

    current_app.redis.mset({save_id: json.dumps(redis_content)})

    return already_added, not_added


def remove_smi_from_redis(paper_id, smi):
    save_id = f"{paper_id}_osr_smiles"
    if current_app.redis.exists(save_id):
        redis_content = json.loads(current_app.redis.get(save_id))
        if smi in redis_content:
            redis_content.remove(smi)
            current_app.redis.mset({save_id: json.dumps(redis_content)})


@bp.route('/_remove_osr_molecule',methods=['GET', 'POST'])
@roles_required('contributor')
def remove_osr_molecule():
    paper_id = request.form['paper_id']
    smi = request.form['smi']
    remove_smi_from_redis(paper_id, smi)

    result = {'status': 'success',
              'msg': 'Molecule deleted',
              'issues': []}
    return jsonify(result=result)

@bp.route('/_add_smi_to_paper',methods=['GET', 'POST'])
@roles_required('contributor')
def add_smi_to_paper():
    paper_id = request.form['paper_id']
    smi = request.form['smi']
    name = request.form['name']
    paper = paper_from_id(paper_id)

    new_mol = molecules_CUD.make_new_activity_molecule(smi, paper, name=name, chem_name=None)
    if new_mol is not None:
        remove_smi_from_redis(paper_id, smi)
        result = {'status': 'success',
                  'msg': f'Molecule {name} - {smi} saved to paper',
                  'issues': []}
        return jsonify(result=result)
    else:
        result = {'status': 'success',
                  'msg': f'Problem creating molecule {name} - {smi}',
                  'issues': []}
        return jsonify(result=result)