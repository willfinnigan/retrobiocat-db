import json

from flask import request, render_template

from retrobiocat_web.analysis.drawing.activity_data_smiles_to_svg import smiles_to_svg
from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.app.database_bps.tables.get_table_data.alt_naming import get_alt_naming
from retrobiocat_web.app.database_bps.tables.get_table_data.get_show_activity_table_data import activity_data_from_args
from retrobiocat_web.mongo.model_queries import sequence_queries
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id


def title_from_args(args):
    title = 'Activity data'

    if 'enzyme_type' in args:
        title += f" - for {args['enzyme_type']} enzymes"

    if 'enzyme_name' in args:
        title += f" - for {args['enzyme_name']}"
    elif 'enzyme_names' in args:
        title += f" - for {len(args.get('enzyme_names', []))} selected enzymes"

    if 'reaction' in args:
        title += f" - for {args['reaction']} reactions"

    if 'paper_id' in args:
        paper = paper_from_id(args['paper_id'])
        title += f" - in {paper.short_citation}"

    if 'formulation' in args:
        title += f" - with biocatalyst forumation {args['formulation']}"

    if 'solvent' in args:
        title += f" - with solvent {args['solvent']}"

    if 'temp_left' in args:
        title += f"- above {args['temp_left']} deg C. "
    if 'temp_right' in args:
        title += f"- below {args['temp_right']} deg C."

    if 'ph_left' in args:
        title += f"- above pH {args['ph_left']}"
    if 'ph_right' in args:
        title += f"- below pH {args['ph_right']}"

    if args.get('include_negative', True) == False:
        title += " - excluding negative data"

    if args.get('include_auto_generated', False) == True:
        title += " - including automatically generated entries"

    if 'smiles' in args:
        #smiles = json.loads(args.get('smiles', []))
        smiles = args.get('smiles', [])
        smiles = [smiles]
        print(smiles)

        if 'smi_col' in args:
            smi_col = args.get('smi_col', 'product_1_smiles')
        else:
            smi_col = args.get('mol_type', 'product_1_smiles')

        if len(smiles) == 1:
            title += f" - for {smi_col} {smiles[0]}"
        else:
            title += f" - for {len(smiles)} selected {smi_col} smiles"

    if args.get('reviewed', True) == False:
        title += " (including not reviewed)"

    return title



@bp.route("/show_activity/", methods=["GET"])
def show_activity():
    args = request.args.to_dict()  # eg reviewed, enzyme_type, enzyme_name, reaction
    activity_data = activity_data_from_args(args)
    activity_data = smiles_to_svg(activity_data)
    enzyme_alt_names_to_use = get_alt_naming(args)
    title = title_from_args(args)

    return render_template('show_activity/show_activity.html', substrate_specificity_data=activity_data, title=title, alt_names=enzyme_alt_names_to_use)