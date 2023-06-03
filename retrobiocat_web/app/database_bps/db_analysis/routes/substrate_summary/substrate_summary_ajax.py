import json
import uuid

from flask import request, current_app, render_template, jsonify
from rdkit import Chem

from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.analysis.data_query.data_query_from_args import data_query_from_args


# def active_inactive_counts(smi, mol_type):
#     acts = activity_for_mol(smi, mol_type)
#
#     num_active, num_inactive = 0, 0
#     active_enzymes, active_papers = 0, 0
#     inactive_enzymes, inactive_papers = 0, 0
#
#     actives, inactives = [], []
#     for act in acts:
#         if act.binary == True:
#             actives.append(act)
#         else:
#             inactives.append(act)
#
#     active_papers, active_enzymes = [], []
#     for act in actives:
#         if act.paper not in active_papers:
#             active_papers.append(act)
#         if act.enzyme_name not in active_enzymes:
#             active_enzymes.append(act)
#
#     inactive_papers, inactive_enzymes = [], []
#     for act in inactives:
#         if act.paper not in inactive_papers:
#             inactive_papers.append(act)
#         if act.enzyme_name not in inactive_enzymes:
#             inactive_enzymes.append(act)


@bp.route('/_load_substrate_summary_examples', methods=['POST'])
def load_substrate_summary_examples():
    leafs = json.loads(request.form['leafs'])
    args = json.loads(request.form['args'])
    args['smiles'] = leafs

    dq = data_query_from_args(args)
    mol_stats = dq.get_mol_stats()

    examples = [{'smi': leaf, 'mol': Chem.MolFromSmiles(leaf)} for leaf in leafs]
    html = render_template('substrate_summary/modals/examples_jinja_template.html',
        examples=examples, mol_type=args['mol_type'], mol_stats=mol_stats)
    return jsonify({'html': html})


@bp.route('/_send_heatmap_mol_list', methods=['POST'])
def send_heatmap_mol_list():
    """ Gets a list of smis and saves these to redis for loading into a heatmap """
    mol_list_id = str(uuid.uuid4())
    current_app.redis.mset({mol_list_id: request.form['leafs']})
    current_app.redis.expire(mol_list_id, 120)
    return jsonify({'mol_list_uuid': mol_list_id})

