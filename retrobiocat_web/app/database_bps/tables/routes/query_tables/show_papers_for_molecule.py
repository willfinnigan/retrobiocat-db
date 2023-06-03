from flask import request, render_template, flash, redirect, url_for

from retrobiocat_web.app.database_bps.tables import bp
import mongoengine as db

from retrobiocat_web.app.database_bps.tables.get_table_data.get_papers_table_data import process_papers_to_table
from retrobiocat_web.mongo.models.biocatdb_models import Activity



def safe_smi(smi):
    smi = smi.replace('[hash]', '#')
    smi = smi.replace('fs', '/')
    return smi

# localhost:5000/papers_for_molecule?product_1_smiles='CCCC=O'
@bp.route("/papers_for_molecule", methods=["GET"])
def papers_for_molecule():
    args = request.args.to_dict()
    title = "Papers"

    if 'product_1_smiles' not in args and 'substrate_1_smiles' not in args and 'substrate_2_smiles' not in args:
        flash('No molecules specified')
        return redirect(url_for("main_site.home"))

    revQ = db.Q()
    if 'reviewed' in args:
        revQ = db.Q(reviewed=True)

    prod_Q = db.Q()
    if 'product_1_smiles' in args:
        prod_Q = db.Q(product_1_smiles=safe_smi(args['product_1_smiles']))
        title += f" containing product_1_smiles {args['product_1_smiles']}"

    sub1_Q = db.Q()
    if 'substrate_1_smiles' in args:
        sub1_Q = db.Q(substrate_1_smiles=safe_smi(args['substrate_1_smiles']))
        title += f" containing substrate_1_smiles {args['substrate_1_smiles']}"

    sub2_Q = db.Q()
    if 'substrate_2_smiles' in args:
        sub2_Q = db.Q(substrate_2_smiles=safe_smi(args['substrate_2_smiles']))
        title += f" containing substrate_1_smiles {args['substrate_2_smiles']}"

    et_q = db.Q()
    if 'enzyme_type' in args:
        if args['enzyme_type'] != 'All':
            et_q = db.Q(enzyme_type=args['enzyme_type'])
            title += f" for enzyme_type {args['enzyme_type']}"

    r_q = db.Q()
    if 'reaction' in args:
        if args['reaction'] != 'All':
            r_q = db.Q(reaction=args['reaction'])
            title += f" for reaction {args['reaction']}"

    acts = Activity.objects(revQ & prod_Q & sub1_Q & sub2_Q & et_q & r_q).only('paper').select_related()
    papers = []
    for act in acts:
        if act.paper not in papers:
            papers.append(act.paper)
    table_data = process_papers_to_table(papers)

    paper_table_options = {'table_height': '80vh'}

    return render_template('show_papers/show_papers.html',
                           papers_data=table_data,
                           paper_table_options=paper_table_options,
                           title=title)

