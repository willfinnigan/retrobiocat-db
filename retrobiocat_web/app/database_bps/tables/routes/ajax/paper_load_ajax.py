from flask import request, jsonify

from retrobiocat_web.app.database_bps.tables import bp
from retrobiocat_web.mongo.model_queries import paper_queries

@bp.route('/_load_paper_data', methods=['GET', 'POST'])
def load_paper_data():

    paper = paper_queries.paper_from_id(request.form['paper_id'])

    result = {'title': paper.title,
              'authors': paper.authors_str(),
              'owner': paper.owner_text(),
              'journal': paper.journal,
              'date': paper.str_date(),
              'tags': paper.tags_str(),
              'status': paper.status,
              'doi': paper.doi,
              'paper_id': str(paper.id),
              'v_short_cit': f"<i>{paper.v_short_cit()}</i>"}

    return jsonify(result=result)