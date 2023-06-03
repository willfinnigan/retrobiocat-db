from flask import request, jsonify
from flask_security import roles_required
from retrobiocat_web.analysis.uniprot_and_web import crossref_lookup, pubmed
from retrobiocat_web.app.database_bps.curation import bp
from retrobiocat_web.app.database_bps.curation.functions.list_parse import list_to_string
from retrobiocat_web.mongo.model_queries import paper_queries


def get_fields_from_pubmed(doi):

    title, authors_list, journal, date, cite_mini = pubmed.query_pubmed(doi)

    if cite_mini == '':
        result = {'status': 'danger',
                  'msg': 'DOI not found in pubmed',
                  'issues': []}
        return jsonify(result=result)

    paper_dict = {'short_cit': cite_mini,
                  'doi': doi,
                  'date': str(date),
                  'title': title,
                  'journal': journal,
                  'authors': list_to_string(authors_list)}

    result = {'status': 'success',
              'msg': 'Fields updated, click save to update the database',
              'issues': [],
              'paper': paper_dict}
    return jsonify(result=result)


def get_fields_from_crossref(doi):
    title, authors_list, journal, date, cite_mini = crossref_lookup.get_metadata_from_crossref(doi)

    if cite_mini == '':
        result = {'status': 'danger',
                  'msg': 'DOI not found in crossref',
                  'issues': []}
        return jsonify(result=result)

    paper_dict = {'short_cit': cite_mini,
                  'doi': doi,
                  'date': str(date),
                  'title': title,
                  'journal': journal,
                  'authors': list_to_string(authors_list)}

    result = {'status': 'success',
              'msg': 'Fields updated, click save to update the database',
              'issues': [],
              'paper': paper_dict}
    return jsonify(result=result)


@bp.route('/_query_pubmed_or_crossref', methods=['GET', 'POST'])
@roles_required('contributor')
def query_pubmed_or_crossref():
    service_name = request.form['service_name']
    paper = paper_queries.paper_from_id(request.form['paper_id'])
    doi = str(paper.doi).replace(' ', '')  # should be unnecessary, but just encase clean doi of spaces

    if service_name == 'pubmed':
        return get_fields_from_pubmed(doi)
    elif service_name == 'crossref':
        return get_fields_from_crossref(doi)
    else:
        result = {'status': 'danger',
                  'msg': 'Service must be pubmed or crossref',
                  'issues': []}
        return jsonify(result=result)
