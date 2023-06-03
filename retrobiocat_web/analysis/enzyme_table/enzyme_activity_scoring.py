from retrobiocat_web.mongo.models.biocatdb_models import Sequence, Paper, Activity
from retrobiocat_web.analysis.data_query import get_data
import mongoengine as db
from sklearn.preprocessing import MinMaxScaler

# 1. num papers
# 2. num active
# 3. num high activity

def normalise(scores):
    scaler = MinMaxScaler()
    normalised = scaler.fit_transform(scores)
    return normalised

def score_enzymes(df, unique_enzymes):

    enzyme_scores = []
    enzyme_info_dict = {}
    for enz in unique_enzymes:
        enz_df = df[df['enzyme_name'] == enz]
        num_papers = len(enz_df['paper'].unique())

        true_df = enz_df[enz_df['binary'] == 'True']
        num_active = len(true_df['product_1_smiles'].unique())

        high_df = enz_df[enz_df['categorical'] == 'High']
        num_high = len(high_df['product_1_smiles'].unique())

        enzyme_scores.append([num_papers, num_active, num_high])
        enzyme_info_dict[enz] = [num_papers, num_active, num_high]

    # normalise each to between 0 and 1
    normalised_scores = normalise(enzyme_scores)

    # append total
    final_scores = []
    for norm_scores in normalised_scores:
        final_scores.append(list(norm_scores) + [round(sum(norm_scores)/3, 2)])

    # create dict of enzymes and sort by final score
    final_dict = {e: s for e, s in zip(unique_enzymes, final_scores)}
    final_dict = {k: v for k, v in sorted(final_dict.items(), key=lambda item: item[1][3], reverse=True)}

    return final_dict, enzyme_info_dict

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    enzyme_type = 'IRED'
    reaction = 'Reductive amination'

    dq = get_data.DataQuery(enzyme_type=enzyme_type, reaction=reaction)
    unique_enzymes = dq.unique_enzymes()
    activity_df = dq.activity_df()

    enzyme_scores, enzyme_info = score_enzymes(activity_df, unique_enzymes)
    print(enzyme_scores)

