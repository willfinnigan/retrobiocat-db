from flask import request, render_template, flash, redirect, url_for

from retrobiocat_web.app.database_bps.tables import bp
import mongoengine as db

from retrobiocat_web.app.database_bps.tables.get_table_data.filter_papers import filter_papers_by_enzyme_name, \
    filter_papers_by_reaction
from retrobiocat_web.app.database_bps.tables.get_table_data.get_papers_table_data import process_papers_to_table
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity


def papers_from_args(args):
    revQ = db.Q()
    if 'not_reviewed' in args:
        revQ = db.Q(reviewed=False)

    enzyme_type_query = db.Q()
    if 'enzyme_type' in args:
        enzyme_type_query = db.Q(tags=args['enzyme_type'])

    papers = Paper.objects(enzyme_type_query & revQ)

    if 'enzyme_name' in args:
        papers = filter_papers_by_enzyme_name(papers, args['enzyme_name'])

    if 'reaction' in args:
        papers = filter_papers_by_reaction(papers, args['reaction'])

    return papers

def title_from_args(args):
    title = "Papers"
    if 'enzyme_type' in args:
        title += f" - featuring {args['enzyme_type']} enzymes"
    if 'enzyme_name' in args:
        title += f" - featuring {args['enzyme_name']}"
    if 'reaction' in args:
        title += f" - featuring a {args['reaction']} reaction"
    return title

@bp.route("/papers", methods=["GET"])
def show_papers():

    args = request.args.to_dict()  # reviewed, enzyme_type, enzyme_name, reaction
    papers = papers_from_args(args)
    table_data = process_papers_to_table(papers, include_counts=True)
    title = title_from_args(args)

    paper_table_options = {'table_height': '80vh'}

    return render_template('show_papers/show_papers.html',
                           papers_data=table_data,
                           paper_table_options=paper_table_options,
                           title=title)


