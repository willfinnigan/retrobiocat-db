import uuid

from flask import request, render_template
from flask_security import roles_required

from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.app.database_bps.tables.get_table_data.get_papers_table_data import process_papers_to_table
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import all_enzyme_type_strings, \
    get_enzyme_types_with_full_names
from retrobiocat_web.mongo.model_queries.paper_queries import get_papers_that_need_data, \
    get_high_importance_papers_with_no_owner


@bp.route('/papers_need_data', methods=['GET', 'POST'])
@roles_required('contributor')
def papers_that_need_data():
    args = request.args.to_dict()
    enzyme_type = args.get('enzyme_type', None)

    title = "Papers that require curating"
    if 'enzyme_type' in args:
        title += f" - {args['enzyme_type']}"

    papers = get_papers_that_need_data(enzyme_type=enzyme_type)
    papers_data = process_papers_to_table(papers)

    paper_table_options = {'table_height': '80vh'}

    return render_template('edit_papers/self_assign_papers.html',
                           papers_data=papers_data,
                           paper_table_options=paper_table_options,
                           title=title)


def get_tags_from_list_papers(papers):
    tags = []
    for paper in papers:
        for tag in paper.tags:
            if tag not in tags:
                tags.append(str(tag))
    tags = sorted(tags)
    return tags

def get_table_data_per_tag(tags):
    data_by_tag = {}
    for tag in tags:
        hi_papers_of_tag = get_high_importance_papers_with_no_owner(tag=tag)
        papers_data = process_papers_to_table(hi_papers_of_tag, short=True)
        data_by_tag[tag] = papers_data
    return data_by_tag

def variable_name_safe_tag(tag):
    safe_tag = tag.replace('-', '')
    return safe_tag

@bp.route('/priority_papers', methods=['GET'])
def high_importance_papers():
    hi_papers = get_high_importance_papers_with_no_owner()
    tags = get_tags_from_list_papers(hi_papers)
    tag_ids = {tag: variable_name_safe_tag(tag) for tag in tags}
    data_by_tag = get_table_data_per_tag(tags)
    enzyme_full_names = get_enzyme_types_with_full_names(tags)

    return render_template('edit_papers/high_importance_papers.html',
                           data_by_tag=data_by_tag,
                           tags=tags, tag_ids=tag_ids,
                           enzyme_full_names=enzyme_full_names)
