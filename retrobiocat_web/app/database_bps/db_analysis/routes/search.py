from flask import render_template

from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.mongo.model_queries import sequence_queries


@bp.route('/search/', methods=['GET'])
def search():
    return render_template('search.html')

@bp.route('/enzyme_name_search/', methods=['GET'])
def enzyme_name_search():
    enzyme_names = sequence_queries.seqs_of_type_with_other_names('All')
    return render_template('enzyme_name_search.html', enzyme_names=enzyme_names)
