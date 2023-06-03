from retrobiocat_web.app.database_bps.db_analysis import bp
from flask import render_template
from retrobiocat_web.mongo.models.biocatdb_models import Paper, Activity, Sequence


@bp.route('/db_stats', methods=['GET'])
def db_stats():

    num_activity = Activity.objects(reviewed=True).count()
    num_enzymes = Sequence.objects(reviewed=True).count()
    num_product_mols = len(Activity.objects(reviewed=True).distinct('product_1_smiles'))
    num_papers = Paper.objects(status='Complete').count()

    return render_template('db_stats.html',
                           num_activity=num_activity,
                           num_enzymes=num_enzymes,
                           num_papers=num_papers,
                           num_product_mols=num_product_mols)



