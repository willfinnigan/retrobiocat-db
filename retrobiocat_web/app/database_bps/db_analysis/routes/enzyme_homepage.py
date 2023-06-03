from flask import render_template, request, jsonify

from retrobiocat_web.app.database_bps.curation.functions.check_permission import check_team_permission
from retrobiocat_web.app.database_bps.db_analysis import bp
from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, Activity, EnzymeTypeSearchHistory, SSN_record
from retrobiocat_web.mongo.models.reaction_models import Reaction
from retrobiocat_web.app.retrobiocat.functions.get_images import smiles_rxn_to_svg
import mongoengine as db
import time
from retrobiocat_web.mongo.models.user_models import User
from flask_security import current_user
from retrobiocat_web.analysis.enzyme_type_info_scores.enzyme_type_search_score import get_history_data
from datetime import datetime, timedelta
from retrobiocat_web.app.database_bps.db_analysis.functions.get_enzyme_team_information import get_team
from retrobiocat_web.analysis.enzyme_type_info_scores import enzyme_scores

# http://0.0.0.0:5000/heatmap/reaction=Aldehyde%20oxidation__enzyme_type=ADH__paperid=None__enzyme_names=None__moltype=substrate_1_smiles__onlyreviewed=True/
# "localhost:5000/heatmap?enzyme_type=CAR&paper_id=24934239"

def get_ssn_data(enzyme_type_obj):
    ssn_obj_q = SSN_record.objects(enzyme_type=enzyme_type_obj)
    if len(ssn_obj_q) == 0:
        return [], 'Not available', 0
    ssn_obj = ssn_obj_q[0]
    ssn_status = ssn_obj.status

    ssn_alignment_options = []
    score_idt_above_40 = None
    for score in ssn_obj.num_at_alignment_score:
        clusters = ssn_obj.num_at_alignment_score[score]
        idt = ssn_obj.identity_at_alignment_score[score]
        if idt[0] > 0.4 and score_idt_above_40 == None:
            score_idt_above_40 = score

        choice_text = f"{score}, {clusters} clusters, avg identity {idt[0]} ± {idt[1]}"
        ssn_alignment_options.append([score, choice_text])
    return ssn_alignment_options, ssn_status, score_idt_above_40



def check_enzyme_type_obj(enzyme_type):
    enzyme_type_obj_q = EnzymeType.objects(enzyme_type=enzyme_type)
    if len(enzyme_type_obj_q) == 0:
        return False
    return True


def get_colour(pc):
    if pc > 65:
        colour = "bg-success"
    elif pc > 30:
        colour = "bg-warning"
    else:
        colour = "bg-danger"
    return colour

def get_width(pc):
    return f"{pc}%"



@bp.route("/enzyme_homepage/<enzyme_type>", methods=["GET"])
def enzyme_homepage(enzyme_type):
    t0 = time.time()
    enzyme_type_obj = EnzymeType.objects(enzyme_type__iexact=enzyme_type)[0]
    fullname = enzyme_type_obj.full_name
    description = enzyme_type_obj.description

    reaction_objs = Reaction.objects(enzyme_types__iexact=enzyme_type)
    reaction_dict = {}
    for r in reaction_objs:
        paper_ids = Activity.objects(db.Q(reaction=r.name) & db.Q(enzyme_type=enzyme_type)).distinct('paper')
        num_papers = len(paper_ids)

        if r.example_rxn_string != None:
            example_rxn_string = r.example_rxn_string
        else:
            example_rxn_string = ">>"

        try:
            smi_rxn_svg = smiles_rxn_to_svg(example_rxn_string, rxnSize=(350, 100))
        except Exception as e:
            print(e)
            smi_rxn_svg = ''

        cofactors_used = ", ".join(r.cofactors[enzyme_type]['cofactors_minus'])
        cofactors_produced = ", ".join(r.cofactors[enzyme_type]['cofactors_plus'])

        reaction_dict[r.name] = {"cofactors_used": cofactors_used,
                                 'cofactors_produced': cofactors_produced,
                                 'svg': smi_rxn_svg,
                                 'num_papers': num_papers,
                                 'experimental': r.experimental,
                                 'multi_step': r.two_step,
                                 'requires_absence_of_water': r.requires_absence_of_water}

    reaction_dict = {k: v for k, v in sorted(reaction_dict.items(), key=lambda item: item[1]['num_papers'], reverse=True)}

    t1 = time.time()
    print(f"Loaded reactions for enzyme homepage in {round(t1 - t0, 2)} seconds")

    # class (colour), width, text
    data_entry_completion = int(enzyme_type_obj.database_score_dict['papers_complete']*100)
    paper_search_score = int(enzyme_type_obj.database_score_dict['search_score']*100)
    priority_completion = int(enzyme_type_obj.database_score_dict['priority_complete']*100)

    database_dict = {'data_entry_completion': [get_colour(data_entry_completion),
                                               get_width(data_entry_completion),
                                               enzyme_type_obj.database_score_dict['papers_complete_string']],

                     'paper_search': [get_colour(paper_search_score),
                                      get_width(paper_search_score),
                                      enzyme_type_obj.database_score_dict['search_score_string']],

                     'priority_completion': [get_colour(priority_completion),
                                             get_width(priority_completion),
                                             enzyme_type_obj.database_score_dict['priority_complete_string']]
                     }

    t2 = time.time()
    print(f"Loaded database info for enzyme homepage in {round(t2 - t1, 2)} seconds")

    team_info = get_team(enzyme_type_obj)

    t3 = time.time()
    print(f"Loaded team info for enzyme homepage in {round(t3 - t2, 2)} seconds")

    has_admin_rights = False
    user_in_team = False
    if current_user.is_authenticated:
        user = User.objects(id=current_user.id)[0]
        has_admin_rights = check_team_permission(user, enzyme_type_obj)
        if enzyme_type_obj in user.enzyme_teams:
            user_in_team = True

    paper_searches = get_history_data(enzyme_type)
    t4 = time.time()
    print(f"Loaded paper searchs for enzyme homepage in {round(t4 - t3, 2)} seconds")

    ssn_alignment_options, ssn_status, score_idt_above_40 = get_ssn_data(enzyme_type_obj)
    t5 = time.time()
    print(f"Loaded ssn info for enzyme homepage in {round(t5 - t4, 2)} seconds")

    print(f"Loaded enzyme homepage in {round(t5-t0,2)} seconds")

    return render_template('enzyme_homepage/enzyme_homepage.html',
                           enzyme_type=enzyme_type,
                           fullname=fullname,
                           reactions=reaction_dict,
                           description=description,
                           has_admin_rights=has_admin_rights,
                           database_dict=database_dict,
                           paper_searches=paper_searches,
                           team_info=team_info,
                           db_stats=enzyme_type_obj.database_score_dict,
                           ssn_alignment_options=ssn_alignment_options,
                           ssn_status=ssn_status,
                           score_idt_above_40=score_idt_above_40,
                           user_in_team=user_in_team
                           )

@bp.route("/_submit_papers_search_record", methods=["POST"])
def submit_papers_search_record():

    enzyme_type = request.form['enzyme_type']
    score = int(request.form['score'])
    user = User.objects(id=current_user.id)[0]

    if not check_enzyme_type_obj(enzyme_type):
        print(f'Enzyme type - {enzyme_type} - does not exist')
        return jsonify(result={})

    enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

    if not check_team_permission(user, enzyme_type_obj):
        print('No user access to update enzyme type')
        return jsonify(result={})

    date_from = datetime.utcnow() - timedelta(hours=24)
    query = EnzymeTypeSearchHistory.objects(db.Q(enzyme_type=enzyme_type_obj) &
                                            db.Q(user=user) &
                                            db.Q(date__gte=date_from))



    if len(query) == 0:
        new_history = EnzymeTypeSearchHistory(enzyme_type=enzyme_type_obj,
                                              user=user,
                                              score=score)
        new_history.save()
    else:
        history = query[0]
        history.score = score
        history.save()

    enzyme_scores.update_single_enzyme_score(enzyme_type_obj)

    return jsonify(result={})