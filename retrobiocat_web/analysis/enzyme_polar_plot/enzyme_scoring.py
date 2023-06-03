import pandas as pd

def get_act_df(activity_data):
    return pd.DataFrame(list(activity_data))

def score_enzyme(act_df, label_dict, label_reps, enzyme, smi_col, max_norm=5, show_log=False):
    enzyme_scores = {}
    for label, smis in label_dict.items():
        temp_df = act_df[(act_df['enzyme_name'] == enzyme) & (act_df['binary'] == 'True') & act_df[smi_col].isin(smis)]
        temp_df = temp_df.drop_duplicates(subset=[smi_col])
        num_active = len(temp_df.index)

        if show_log == True:
            print(f"Group {label} - {num_active} active")
            print(f"Representative smiles = {label_reps[label]}")

        norm_to = len(smis)
        if norm_to > max_norm:
            norm_to = max_norm
        score = round(num_active / norm_to, 2)
        if score > 1:
            score = 1
        enzyme_scores[label] = score

    total_score = round((sum(enzyme_scores.values()) / len(label_dict)) * 100, 1)

    return total_score, enzyme_scores

def score_all_enzymes(act_df, label_dict, label_reps, list_enzymes, smi_col, max_norm=5):
    all_total_scores = {}
    all_group_scores = {}
    for enzyme in list_enzymes:
        total_score, enzyme_scores = score_enzyme(act_df, label_dict, label_reps, enzyme, smi_col, max_norm=max_norm)
        all_total_scores[enzyme] = total_score
        all_group_scores[enzyme] = enzyme_scores

    return all_total_scores, all_group_scores

if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    from retrobiocat_web.analysis.data_query import get_data
    dq = get_data.DataQuery(enzyme_type='CAR', reaction="Carboxylic acid reduction", smi_col='substrate_1_smiles', log_level=1)
    unique_smis = dq.unique_smiles()
    unique_enzymes = dq.unique_enzymes()
    fps = get_data.get_substrates_fps(unique_smis)
    fp_dict = get_data.make_fp_dict(unique_smis, fps)

    from retrobiocat_web.analysis.enzyme_polar_plot.cluster_molecules import get_clusters
    labels, label_dict, label_reps, rep_labels = get_clusters(unique_smis, fps, 12, random_seed=0)

    act_df = get_act_df(dq.activity_data)
    all_total_scores, all_group_scores = score_all_enzymes(act_df, label_dict, label_reps, unique_enzymes, 'substrate_1_smiles', max_norm=5)
    print(all_total_scores['NbCAR'])
    print(all_group_scores['NbCAR'])

    nb_car_score, nb_car_dict = score_enzyme(act_df, label_dict, label_reps, 'NbCAR', 'substrate_1_smiles', max_norm=5, show_log=True)
    print(nb_car_dict)


