from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template, redirect, url_for, request, jsonify, current_app

from retrobiocat_web.app.main_site.functions.get_queue_task_details import queue_task_details
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg, get_rxn_smiles

from retrobiocat_web.app.main_site.functions.progress_bar import set_progress_bar
from retrobiocat_web.mongo.models.biocatdb_models import Activity
from rq import get_current_job
from rq.registry import FinishedJobRegistry, StartedJobRegistry
from retrobiocat_web.app.database_bps.db_analysis.forms import HeatmapForm
from retrobiocat_web.mongo.models.biocatdb_models import Paper
import json
from distutils.util import strtobool

from retrobiocat_web.analysis.data_query import get_data, split_data

from retrobiocat_web.analysis.data_vis.heatmap import HeatMapper
from bokeh.embed import components

@bp.route('/heatmap_form', methods=['GET', 'POST'])
def heatmap_form():
    form = HeatmapForm()
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
        max_molecules = form.data['max_molecules']
        max_enzymes = form.data['max_enzymes']
        remove_negative = form.data['remove_negative']

        return redirect(url_for('.heatmap', reaction=reaction, enzyme_type=enzyme_type,
                                only_reviewed=only_reviewed, mol_type=mol_type,
                                max_molecules=max_molecules, max_enzymes=max_enzymes,
                                remove_negative=remove_negative))

    return render_template('heatmap/heatmap_form.html', form=form)


def task_make_heatmap(reaction, enzyme_type, enzyme_names,
                      paper_id, only_reviewed, mol_type,
                      max_molecules, max_enzymes, remove_negative, group, mols,
                      mol_list_id):

    job = get_current_job()
    set_progress_bar(job, 40, 'started')

    title = f"Activity data for "
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

    if group != -10:
        title += f" - group {group}"

    if remove_negative == 'all':
        title += " (No negative data)"

    if mol_list_id is not None:
        mols = json.loads(current_app.redis.get(mol_list_id))

    initial_dq = get_data.DataQuery(enzyme_type=enzyme_type,
                            reaction=reaction,
                            only_reviewed=only_reviewed,
                            enzyme_names=enzyme_names,
                            paper_id=paper_id,
                            smi_col=mol_type,
                            remove_negative=remove_negative,
                            smiles=mols)

    set_progress_bar(job, 60, 'data retrieved')

    dqs = split_data.split(initial_dq, max_molecules, max_enzymes, log_level=0)

    if len(dqs) != 1:
        if group >= 1 and group <= len(dqs):
            dq = dqs[group-1]
            hm = HeatMapper(dq, log_level=0, smi_col=mol_type)
            heat_map = hm.create_heatmap()
            script, div = components(heat_map)
            return {'script': script,
                    'div': div,
                    'title': title}

        else:
            return {'groups': len(dqs),
                    'possible_reactions': initial_dq.reactions(),
                    'title': title}

    else:
        set_progress_bar(job, 80, 'heatmap ready')
        hm = HeatMapper(initial_dq, log_level=0, smi_col=mol_type)
        heat_map = hm.create_heatmap()
        script, div = components(heat_map)
        return {'script': script,
                'div': div,
                'title': title}

def get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type, max_molecules, max_enzymes, remove_negative, group, mols, mol_list_id):
    if type(mols) == str:
        mols = mols.replace('/', 'fs')
        mols = mols.replace('#', 'hash')

    job_id = f"reaction={reaction}__enzyme_type={enzyme_type}__paperid={paper_id}__" \
             f"enzyme_names={enzyme_names}__moltype={mol_type}__onlyreviewed={only_reviewed}" \
             f"max_molecules={max_molecules}__max_enzymes={max_enzymes}__remove_negative={remove_negative}__group={group}__mols={mols}__mol_list_id={mol_list_id}"
    return job_id

# http://0.0.0.0:5000/heatmap/reaction=Aldehyde%20oxidation__enzyme_type=ADH__paperid=None__enzyme_names=None__moltype=substrate_1_smiles__onlyreviewed=True/
# "localhost:5000/heatmap?enzyme_type=CAR&paper_id=24934239"
@bp.route("/heatmap/", methods=["GET"])
def heatmap():
    # options are:  Reaction, Enzyme Type, Paper.  Also specify mol type, but default to product.
    # Must have at least one
    # Also params for splitting the data up

    def check_args(args):
        required_args = ['enzyme_type', 'reaction', 'paper_id', 'enzyme_names']
        allowed_args = required_args + ['only_reviewed', 'mol_type', 'max_molecules', 'max_enzymes', 'group', 'remove_negative', 'mols']
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
    max_molecules = int(args.get('max_molecules', '1000'))
    max_enzymes = int(args.get('max_enzymes', '200'))
    group = int(args.get('group', '-10'))
    remove_negative = args.get('remove_negative', 'False')
    mols = args.get('mols', None)
    mol_list_id = args.get('mol_list_id', None)

    job_id = get_task_id(reaction, enzyme_type, enzyme_names, paper_id, only_reviewed, mol_type, max_molecules, max_enzymes, remove_negative, group, mols, mol_list_id)

    registry = FinishedJobRegistry(queue=current_app.heatmap_queue)
    started_reg = StartedJobRegistry(queue=current_app.heatmap_queue)
    if job_id in list(registry.get_job_ids()):
        task = current_app.heatmap_queue.fetch_job(job_id)
        if task != None:
            if 'groups' in task.result:
                buttons = [i + 1 for i in range(task.result['groups'])]
                possible_reactions = task.result['possible_reactions']

                return render_template('heatmap/heatmap_select.html', title=task.result['title'],
                                       reaction=reaction, enzyme_type=enzyme_type, enzyme_names=enzyme_names,
                                       paper_id=paper_id, only_reviewed=only_reviewed, mol_type=mol_type,
                                       max_molecules=max_molecules, max_enzymes=max_enzymes, remove_negative=remove_negative,
                                       buttons=buttons, possible_reactions=possible_reactions, mols=mols, mol_list_id=mol_list_id)
            else:
                script = task.result['script']
                div = task.result['div']
                title = task.result['title']
                return render_template('heatmap/heatmap.html', script=script, div=div,
                                       title=title)

    elif job_id not in list(started_reg.get_job_ids()):
        current_app.heatmap_queue.enqueue(task_make_heatmap,
                                          reaction, enzyme_type,
                                          enzyme_names, paper_id,
                                          only_reviewed, mol_type,
                                          max_molecules, max_enzymes,
                                          remove_negative, group, mols,
                                          mol_list_id,
                                          job_id=job_id)

    queue_details, task_details = queue_task_details(job_id, 'heatmap')
    return render_template('queue_loading.html', task_queue='heatmap', task_id=job_id,
                           queue_details=queue_details, task_details='',
                           title='Loading substrate summary', ajax_timer=1500, refresh_timer=30000)


def get_activity_info(act):
    cat_colours = {'High': 'green',
                   'Medium': 'goldenrod',
                   'Low': 'orange',
                   'None': 'red'}

    activity_info = ""

    activity_info += f"<p> {act.reaction} - "
    if act.binary == True:
        activity_info += """<span style="color: green">Active</span>"""
    else:
        activity_info += """<span style="color: red">Inactive</span>"""

    colour = 'darkgrey'
    if act.categorical != None:
        colour = cat_colours.get(act.categorical, 'darkgrey')
        activity_info += f"""<span style="color: {colour}">, {act.categorical}</span>"""

    if act.conversion != None:
        activity_info += f"""<span style="color: {colour}">, {act.conversion} % conversion over {act.conversion_time} hours</span>"""

    if act.specific_activity != None:
        activity_info += f"""<span style="color: {colour}">, {round(act.specific_activity,2)} umols/min/mg</span>"""

    activity_info += "</p>"
    activity_info += f"<p>{act.short_citation}</p>"
    if act.added_by != None:
        activity_info += f"<p>Curated by {act.added_by.first_name} {act.added_by.last_name} <small>({act.added_by.affiliation})</small></p>"

    return activity_info

@bp.route("/_heatmap_modal_info", methods=["POST"])
def heatmap_modal_info():
    act_ids = request.form.getlist('act_ids[]') #json.loads(request.form['act_ids'])

    activity_data = Activity.objects(id__in=act_ids).select_related()

    reactions = []
    for act in activity_data:
        reaction_smiles = get_rxn_smiles([act.substrate_1_smiles, act.substrate_2_smiles],
                                         [act.product_1_smiles])
        reaction_svg = smiles_rxn_to_svg(reaction_smiles, rxnSize=(500, 100))

        return_dict = {}
        return_dict['reaction_svg'] = reaction_svg
        return_dict['reaction_info'] = get_activity_info(act)
        return_dict['enzyme_name'] = act.enzyme_name
        return_dict['activity_id'] = str(act.id)
        return_dict['paper_id'] = str(act.paper.id)
        reactions.append(return_dict)

    enzyme_name = activity_data[0].enzyme_name
    enzyme_type = activity_data[0].enzyme_type

    return jsonify(result={'reactions': reactions,
                           'enzyme_name': enzyme_name,
                           'enzyme_type': enzyme_type})


if __name__ == "__main__":
    args = {}
    print(int(args.get('group', '-10')))
