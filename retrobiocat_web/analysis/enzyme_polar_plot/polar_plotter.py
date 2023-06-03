from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.analysis.enzyme_polar_plot.enzyme_scoring import score_all_enzymes
from retrobiocat_web.analysis.enzyme_polar_plot.cluster_molecules import get_clusters, order_smis
from retrobiocat_web.analysis.enzyme_polar_plot.draw_polar_plot import create_polar_plot
import pandas as pd
from bokeh.layouts import gridplot
from bokeh.io import output_file, show

class PolarPlotter(object):

    def __init__(self, enzyme_type, reaction, smi_col, n_clusters, max_norm=5, random_seed=0):
        self.random_seed = random_seed
        self.smi_col = smi_col
        dq = get_data.DataQuery(enzyme_type=enzyme_type, reaction=reaction, smi_col=smi_col, log_level=1)
        self.unique_smis = dq.unique_smiles()
        self.unique_enzymes = dq.unique_enzymes()
        self.fps = get_data.get_substrates_fps(self.unique_smis)
        self.fp_dict = get_data.make_fp_dict(self.unique_smis, self.fps)
        self.act_df = pd.DataFrame(list(dq.activity_data))
        self.labels, self.label_dict, self.label_reps, self.rep_labels = get_clusters(self.unique_smis, self.fps, n_clusters, random_seed=random_seed)
        self.substrates_ordered = order_smis(self.rep_labels)
        self.all_total_scores, self.all_group_scores = score_all_enzymes(self.act_df, self.label_dict, self.label_reps, self.unique_enzymes, smi_col, max_norm=max_norm)

        # sort by total scores
        self.all_total_scores = {k: v for k, v in sorted(self.all_total_scores.items(), key=lambda item: item[1], reverse=True)}
        self.all_group_scores = {k: self.all_group_scores[k] for k in self.all_total_scores}
        self.unique_enzymes = [k for k in self.all_total_scores]

    def make_polar_plot(self, enzyme_name, inner_radius=25, gap_size=0, num_axes=0,
                        mol_dims=(50, 50), mols_to_outer=False):

        score_dict = self.all_group_scores[enzyme_name]
        polar_dict = {}
        for smi in self.substrates_ordered:
            l = self.rep_labels[smi]
            polar_dict[l] = [score_dict[l] * 100, smi]

        polar_plot = create_polar_plot(polar_dict,
                                       inner_radius=inner_radius, gap_size=gap_size, num_axes=num_axes,
                                       mol_dims=mol_dims, mols_to_outer=mols_to_outer,
                                       enzyme_name=enzyme_name, min_y_value=1)

        return polar_plot

if __name__ == "__main__":
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection(host='138.68.135.53', database='test')

    pp = PolarPlotter('CAR', 'Carboxylic acid reduction', 'product_1_smiles', 12)
    plots = []
    for enzyme in pp.unique_enzymes:
        plots.append(pp.make_polar_plot(enzyme, num_axes=3, mol_dims=(45, 45), mols_to_outer=True))

    layout = gridplot(plots, ncols=2)

    output_file("polar_plots_test.html")
    show(layout)