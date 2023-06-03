from retrobiocat_web.mongo.models.biocatdb_models import EnzymeType, SSN_record, UniRef50, Sequence, EnzymeTypeSearchHistory, EnzymeType, Paper, Activity
import mongoengine as db
import time
from retrobiocat_web.analysis.enzyme_type_info_scores import enzyme_type_search_score

def get_db_stats(enzyme_type_obj):
    '''Function not in use'''

    q_rev = (db.Q(reviewed=True))
    q_et = db.Q(enzyme_type=enzyme_type_obj.enzyme_type)
    q_eto = db.Q(enzyme_type=enzyme_type_obj)
    q_has_seq = db.Q(sequence__ne=None)
    q_has_seq2 = db.Q(sequence__ne="")
    q_act = db.Q(binary=True)


    num_enzymes = Sequence.objects(q_rev & q_et).count()
    num_enzymes_w_seq = Sequence.objects(q_rev & q_et & q_has_seq & q_has_seq2).count()
    num_uniref50 = UniRef50.objects(q_eto).count()
    num_activity = Activity.objects(q_et & q_rev).count()
    num_act = Activity.objects(q_et & q_rev & q_act).count()
    num_inact = num_activity - num_act
    num_products = len(Activity.objects(q_et & q_rev).distinct('product_1_smiles'))

    db_stats = {'num_enzymes': num_enzymes,
                'num_enzymes_w_seqs': num_enzymes_w_seq,
                'num_uniref50': num_uniref50,
                'num_activity': num_activity,
                'pc_active': int((num_act/num_activity)*100),
                'pc_inactive': int((num_inact/num_activity)*100),
                'num_unique_products': num_products}

    return db_stats

def get_db_wide_stats():
    all_seqs = Sequence.objects(db.Q(reviewed=True) &
                                db.Q(sequence__ne=None) &
                                db.Q(sequence__ne="")).count()

    total_unique_products = len(Activity.objects(db.Q(reviewed=True) &
                                                 db.Q(auto_generated=False)).distinct('product_1_smiles'))

    all_activity = Activity.objects(db.Q(reviewed=True) &
                                    db.Q(auto_generated=False)).count()

    return all_seqs, total_unique_products, all_activity

def get_score_for_enzyme(enzyme_type, all_seqs, total_unique_products, all_activity,
                         weights=None , print_log=False):
    if print_log==True:
        print(f"--- {enzyme_type.enzyme_type} ---")

    paper_search_score, score_string = enzyme_type_search_score.calc_score(enzyme_type)

    search_score = round(paper_search_score / 100, 3)

    complete_papers = Paper.objects(db.Q(tags=enzyme_type.enzyme_type) &
                                    db.Q(status='Complete')).count()
    all_papers = Paper.objects(db.Q(tags=enzyme_type.enzyme_type)).count()

    if all_papers == 0:
        pc_complete = 0
    else:
        pc_complete = round((complete_papers / all_papers), 3)

    complete_priority = Paper.objects(db.Q(tags=enzyme_type.enzyme_type) &
                                      db.Q(status='Complete') &
                                      db.Q(high_importance=True)).count()
    num_priority = Paper.objects(db.Q(tags=enzyme_type.enzyme_type) &
                                 db.Q(high_importance=True)).count()

    if num_priority == 0:
        pc_priority = 0
    else:
        pc_priority = round((complete_priority / num_priority), 3)

    num_seqs = Sequence.objects(db.Q(reviewed=True) & db.Q(enzyme_type=enzyme_type.enzyme_type)).count()

    seqs_w_prot = Sequence.objects(db.Q(enzyme_type=enzyme_type.enzyme_type) &
                                   db.Q(reviewed=True) &
                                   db.Q(sequence__ne=None) &
                                   db.Q(sequence__ne="") &
                                   db.Q(sequence_unavailable__ne=True)).count()

    if all_seqs == 0:
        pc_protein = 0
    else:
        pc_protein = round((seqs_w_prot / all_seqs), 3)

    unique_products = len(Activity.objects(db.Q(reviewed=True) &
                                           db.Q(auto_generated=False) &
                                           db.Q(enzyme_type=enzyme_type.enzyme_type)).distinct('product_1_smiles'))

    pc_unique_products = round((unique_products / total_unique_products), 3)

    num_uniref50 = UniRef50.objects(db.Q(enzyme_type=enzyme_type)).count()

    num_activity = Activity.objects(db.Q(enzyme_type=enzyme_type.enzyme_type) &
                                    db.Q(reviewed=True) &
                                    db.Q(auto_generated=False)).count()
    num_act = Activity.objects(db.Q(enzyme_type=enzyme_type.enzyme_type) &
                               db.Q(reviewed=True) &
                               db.Q(auto_generated=False) &
                               db.Q(binary=True)).count()
    num_inact = num_activity - num_act

    if num_activity == 0:
        pc_active = 0
        pc_inactive = 0
        pc_of_all_activity = 0
    else:
        pc_active = round((num_act / num_activity), 3)
        pc_inactive = round((num_inact / num_activity), 3)
        pc_of_all_activity = round((num_activity / all_activity), 3)

    if weights == None:
        w = [1, 0.5, 1, 2, 2, 1]
    else:
        w = weights

    if print_log == True:
        print(
            f"search score = {search_score} %papers complete = {pc_complete}, %priority complete = {pc_priority}, %proteins /w seqs = {pc_protein}, %unique products = {pc_unique_products}, %all_activity = {pc_of_all_activity}")

    score = (w[0]*search_score) + (w[1]*pc_complete) + (w[2]*pc_priority) + (w[3]*pc_protein) + \
            (w[4]*pc_unique_products) + (w[5]*pc_of_all_activity)
    score = score / sum(w)
    score = round(score, 2)
    if print_log == True:
        print(f"score = {score}")
        print()

    return {'search_score': search_score,
           'search_score_string': score_string,
           'papers_complete_string': f"{complete_papers} out of {all_papers}",
           'papers_complete': pc_complete,
           'priority_complete': pc_priority,
           'num_priority': num_priority,
           'num_priority_complete': complete_priority,
           'priority_complete_string': f"{complete_priority} out of {num_priority}",
           '%_proteins': pc_protein,
           'num_unique_products': unique_products,
           '%_unique_products': pc_unique_products,
           '% all activity': pc_of_all_activity,
           'num_enzymes': num_seqs,
           'num_enzymes_w_protein': seqs_w_prot,
           'num_uniref50': num_uniref50,
           'num_activity': num_activity,
           'pc_active': pc_active,
           'pc_inactive': pc_inactive,

           'total_score': score}

def get_scores(queryset_enzyme_types, weights=None):
    '''
    paper search score
    %_complete_papers
    num seqs with protein / total num seqs with protein
    complete priority papers / total complete priority papers
    num unique products / total unique products (not auto)
    then normalise
    '''
    t0 = time.time()

    all_scores = {}
    all_seqs, total_unique_products, all_activity = get_db_wide_stats()
    for enzyme_type in queryset_enzyme_types:
        all_scores[enzyme_type.enzyme_type] = get_score_for_enzyme(enzyme_type, all_seqs, total_unique_products, all_activity, weights=weights)

    all_scores = {k: v for k, v in sorted(all_scores.items(), key=lambda item: item[1]['total_score'], reverse=True)}

    t1 = time.time()
    print(f"Scores calculated in {round(t1-t0, 1)} seconds")
    return all_scores

def update_scores():
    enzyme_types = EnzymeType.objects()

    all_scores_dict = get_scores(enzyme_types)

    for enzyme_type in enzyme_types:
        scores_dict = all_scores_dict[enzyme_type.enzyme_type]
        enzyme_type.database_score = scores_dict['total_score']
        enzyme_type.database_score_dict = scores_dict
        enzyme_type.save()

    return all_scores_dict

def update_single_enzyme_score(enzyme_type_object):
    all_seqs, total_unique_products, all_activity = get_db_wide_stats()
    scores_dict = get_score_for_enzyme(enzyme_type_object, all_seqs, total_unique_products, all_activity)
    enzyme_type_object.database_score = scores_dict['total_score']
    enzyme_type_object.database_score_dict = scores_dict
    enzyme_type_object.save()