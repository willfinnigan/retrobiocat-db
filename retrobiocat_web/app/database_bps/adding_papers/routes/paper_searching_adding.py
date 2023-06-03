from flask import render_template, flash, redirect, url_for, session, request, jsonify, current_app
from flask_security import roles_required, current_user
from retrobiocat_web.analysis.enzyme_type_info_scores.paper_search import keyword_query_pubmed
from retrobiocat_web.analysis.uniprot_and_web.crossref_lookup import get_metadata_from_crossref
from retrobiocat_web.analysis.uniprot_and_web.pubmed import query_pubmed
from retrobiocat_web.app.database_bps.adding_papers import bp
from retrobiocat_web.app.database_bps.adding_papers.forms import PaperSearch
from retrobiocat_web.mongo.models.biocatdb_models import Paper, IrrelevantPaper, EnzymeType
from retrobiocat_web.mongo.models.user_models import User
import mongoengine as db

@bp.route('/query_papers_form', methods=['GET', 'POST'])
@roles_required('paper_adder')
def query_papers_form():
    form = PaperSearch()
    form.set_choices()

    if form.validate_on_submit() == True:
        keyword = form.keyword.data
        enzyme_type = form.enzyme_type.data
        return redirect(url_for("adding_papers.paper_search_results", enzyme_type=enzyme_type, keywords=keyword))

    return render_template('add_paper_workflow/query_papers_form.html', form=form)

@bp.route('/paper_search_results/<enzyme_type>/<keywords>/', methods=['GET', 'POST'])
@roles_required('paper_adder')
def paper_search_results(enzyme_type, keywords):
    papers_data = task_run_paper_query(keywords, enzyme_type)
    return render_template('add_paper_workflow/query_results_table.html',
                           papers_data=papers_data,
                           enzyme_type=enzyme_type)


@bp.route('/_add_paper_from_query', methods=['GET', 'POST'])
@roles_required('paper_adder')
def add_paper_from_query():
    doi = request.form['doi']
    tags = request.form['tags']
    title = request.form['title']
    journal = request.form['journal']

    doi = doi.replace(' ', '').lower()

    user = User.objects(id=current_user.id)[0]

    new_paper = Paper(doi=doi,
                      title=title,
                      html='https://doi.org/' + doi,
                      journal=journal,
                      tags=tags.split(', '),
                      added_by=user,
                      status='Data required')
    new_paper.save()

    current_app.db_queue.enqueue(task_update_paper_metadata_from_crossref, new_paper.id)

    result = {'new_id': str(new_paper.id)}
    return jsonify(result=result)

@bp.route('/_toggle_not_relevant', methods=['GET', 'POST'])
@roles_required('paper_adder')
def toggle_not_relevant():
    doi = request.form['doi']
    enzyme_type = request.form['enzyme_type']
    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    et_q = db.Q(enzyme_type=enzyme_type_obj)
    doi_q = db.Q(doi=doi)
    query = IrrelevantPaper.objects(et_q & doi_q)

    if len(query) == 0:
        irrelevant_paper = IrrelevantPaper(doi=doi,
                                        enzyme_type=enzyme_type_obj)
        irrelevant_paper.save()
        status = 'Not relevant'
    else:
        irrelevant_paper = query[0]
        irrelevant_paper.delete()
        status = 'Missing'

    result = {'status': status}
    return jsonify(result=result)

@bp.route('/_toggle_priority_paper', methods=['GET', 'POST'])
@roles_required('paper_adder')
def toggle_priority_paper():
    paper_id = request.form['paper_id']
    paper = Paper.objects(id=paper_id)[0]

    if paper.high_importance == True:
        paper.high_importance = False
        priority = ''
    else:
        paper.high_importance = True
        priority = 'Priority'
    paper.save()

    result = {'priority': priority}
    return jsonify(result=result)


def task_update_paper_metadata_from_crossref(paper_id):
    paper_query = Paper.objects(id=paper_id)
    if len(paper_query) == 0:
        return

    paper = paper_query[0]
    doi = paper.doi
    title, authors_list, journal, date, cite_mini = get_metadata_from_crossref(doi)
    if cite_mini == '':
        title, authors_list, journal, date, cite_mini = query_pubmed(doi)

    if cite_mini != '':
        paper.short_citation = cite_mini
        paper.title = title
        paper.journal = journal
        paper.date = date
        paper.authors = authors_list
        paper.save()

def get_papers_data(keyword):
    results = keyword_query_pubmed(keyword)
    for i, paper_dict in enumerate(results):
        results[i]['date'] = results[i]['date'].isoformat()
        results[i]['row_num'] = i

        for key, value in paper_dict.items():
            if value == None:
                results[i][key] = 'None'

    return results

def filter_by_fields(papers_data):
    fields = ['title', 'doi', 'authors', 'journal', 'date', 'abstract', 'row_num']
    for i, paper_dict in enumerate(papers_data):
        papers_data[i] = {x: paper_dict[x] for x in fields if x in paper_dict}
    return papers_data

def query_biocatdb(papers_data):
    list_dois = []
    for paper_dict in papers_data:
        doi = paper_dict['doi']
        list_dois.append(doi)

    biocatdb_query = Paper.objects(doi__in=list_dois).as_pymongo()
    doi_tags = {}
    for p in biocatdb_query:
        owner = "no"
        if 'owner' in p:
            owner = "yes"

        priority = ""
        if 'high_importance' in p:
            if p['high_importance'] == True:
                priority = "Priority paper"
        tags = p.get('tags', [])
        status = p.get('status', '')
        doi_tags[p['doi']] = [tags, status, owner, str(p['_id']), priority]

    for i, paper_dict in enumerate(papers_data):
        if paper_dict['doi'] in list(doi_tags.keys()):
            status = 'Present'
            tags = doi_tags[paper_dict['doi']][0]
            db_status = doi_tags[paper_dict['doi']][1]
            owner = doi_tags[paper_dict['doi']][2]
            id = str(doi_tags[paper_dict['doi']][3])
            priority = doi_tags[paper_dict['doi']][4]

        else:
            status = 'Missing'
            tags = []
            db_status = ''
            owner = "no"
            id = ""
            priority = ""
        papers_data[i]['biocatdb'] = status
        papers_data[i]['db_status'] = db_status
        papers_data[i]['tags'] = tags
        papers_data[i]['owner'] = owner
        papers_data[i]['id'] = id
        papers_data[i]['priority'] = priority

    return papers_data

def query_irrelevant_papers(papers_data, enzyme_type):
    list_dois = []
    for paper_dict in papers_data:
        doi = paper_dict['doi']
        list_dois.append(doi)

    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    q_enzyme = db.Q(enzyme_type=enzyme_type_obj)
    q_doi = db.Q(doi__in=list_dois)

    dois_not_relevant = list(IrrelevantPaper.objects(q_enzyme & q_doi).distinct('doi'))

    for i, paper_dict in enumerate(papers_data):
        if (paper_dict['doi'] in dois_not_relevant) and (paper_dict['biocatdb'] != 'Present'):
            papers_data[i]['biocatdb'] = 'Not relevant'

    return papers_data


def task_run_paper_query(keyword, enzyme_type):
    papers_data = get_papers_data(keyword)
    papers_data = filter_by_fields(papers_data)
    papers_data = query_biocatdb(papers_data)
    papers_data = query_irrelevant_papers(papers_data, enzyme_type)
    return papers_data





