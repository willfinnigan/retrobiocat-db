from flask import render_template

from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
import time
from retrobiocat_web.mongo.model_queries.enzyme_type_queries import enzyme_types_ordered_by_database_score


@bp.route("/enzymes/", methods=["GET"])
def enzymes():
    t0 = time.time()
    enzyme_types = enzyme_types_ordered_by_database_score()

    enzymes = []
    enzyme_info = {}
    scores = []
    for et in enzyme_types:
        enzymes.append(et.enzyme_type)
        enzyme_info[et.enzyme_type] = {'full_name': et.full_name,
                                       'score': et.database_score,
                                       'score_dict': et.database_score_dict}
        scores.append(et.database_score)
        smi_rxn_svg = ''
        if et.rep_reaction != None:
            rxn_string = et.rep_reaction.example_rxn_string
            if rxn_string != None:
                smi_rxn_svg = smiles_rxn_to_svg(rxn_string, rxnSize=(350, 70))
        enzyme_info[et.enzyme_type]['svg'] = smi_rxn_svg

    normalised_scores = [round((score/max(scores))*100, 1) for score in scores]
    for i, et in enumerate(enzyme_types):
        enzyme_info[et.enzyme_type]['normalised_score'] = normalised_scores[i]

    t1 = time.time()
    print(f"Enzymes page loaded in {round(t1-t0, 1)} seconds")
    return render_template('enzyme_toolbox_launch_page.html', enzymes=enzymes, enzyme_info=enzyme_info)