from flask import render_template, jsonify, session, request, redirect, url_for
import mongoengine as db
from flask_login import current_user

from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.app.database_bps.tables.get_table_data.alt_naming import get_alt_naming
from retrobiocat_web.app.database_bps.tables.get_table_data.get_sequence_table_data import \
    parse_sequence_data_to_tabulator
from retrobiocat_web.mongo.model_queries.paper_queries import paper_from_id
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Paper
from retrobiocat_web.mongo.models.biocatdb_models import Sequence


def seqs_from_args(args):
    revQ = db.Q()
    if 'reviewed' in args:
        revQ = db.Q(reviewed=True)

    enzyme_type_query = db.Q()
    if 'enzyme_type' in args:
        enzyme_type_query = db.Q(enzyme_type=args['enzyme_type'])

    paper_query = db.Q()
    if 'paper_id' in args:
        paper_query = db.Q(papers=args['paper_id'])

    return Sequence.objects(enzyme_type_query & paper_query & revQ)

def title_from_args(args):
    title = "Enzyme sequences"
    if args.get('reviewed', True) == False:
        title += " (including not reviewed)"
    if 'enzyme_type' in args:
        title += f" for {args['enzyme_type']} enzymes"
    if 'paper_id' in args:
        paper = paper_from_id(args['paper_id'])
        title += f" in {paper.short_citation}"

    return title


@bp.route("/sequences", methods=["GET"])
def show_sequences():

    args = request.args.to_dict()
    title = title_from_args(args)
    sequences = seqs_from_args(args)
    enzyme_alt_names_to_use = get_alt_naming(args)

    trim_sequences = True
    can_copy = False
    if current_user.has_role('admin'):
        trim_sequences = False
        can_copy = True

    table_data = parse_sequence_data_to_tabulator(sequences, add_activity_counts=True, trim_sequences=trim_sequences)

    enzyme_types = sorted(list(EnzymeType.objects().distinct("enzyme_type")))

    seq_table_options = {'table_height': '80vh',
                         'show_header_filters': True,
                         'lock_enzyme_types': True}

    return render_template('show_sequences/show_sequences.html',
                           seq_data=table_data,
                           seq_table_options=seq_table_options,
                           enzyme_types=enzyme_types,
                           alt_names=enzyme_alt_names_to_use,
                           can_copy=can_copy,
                           title=title)


