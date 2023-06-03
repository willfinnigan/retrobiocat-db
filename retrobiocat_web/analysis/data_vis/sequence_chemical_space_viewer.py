from retrobiocat_web.analysis.data_query import get_data
from sklearn.decomposition import PCA
import pandas as pd

from bokeh.layouts import row
from bokeh.models import ColumnDataSource, CustomJS
from bokeh.transform import linear_cmap
from bokeh.models import HoverTool, LassoSelectTool
from bokeh import plotting
from bokeh.colors import RGB

class SpaceViewer():
    enzyme_col = 'enzyme_name'
    active_colour = RGB(44, 181, 0)
    inactive_colour = RGB(135, 12, 0)
    medium_colour = RGB(232, 236, 0)
    low_colour = RGB(203, 138, 38)
    empty_colour = RGB(191, 191, 191)

    def __init__(self, data_query, log_level=0, plot_mode='binary',
                 smi_col='product_1_smiles'):

        self.log_level = log_level
        self.plot_mode = plot_mode  # alternative is categorical

        self.smi_col = smi_col
        self.smi_map = {'product_1_smiles': 'Product',
                        'substrate_1_smiles': 'Substrate 1',
                        'substrate_2_smiles': 'Substrate 2'}

        self.log("HeatMapper initialised")

        self.activity_data = data_query.activity_data

    def pca(self, high_d):
        pca = PCA(n_components=2)
        crds = pca.fit_transform(high_d)
        pca_df = pd.DataFrame(crds, columns=["x", "y"])
        return pca_df

    def get_embeddings_pca(self, embeddings):
        emd_df = self.pca(embeddings)
        emb_xs, emb_ys = list(emd_df['x']), list(emd_df['y'])
        return emb_xs, emb_ys

    def get_fp_pca(self, substrates, fp_dict):
        fp_df = self.pca(list(fp_dict.values()))
        xs, ys = list(fp_df['x']), list(fp_df['y'])
        xy_dict = {}
        for i, smi in enumerate(substrates):
            xy_dict[smi] = (xs[i], ys[i])
        return xy_dict, xs, ys

    def activity_data_to_plot_format(self, activity_data):
        data = {'smiles': [],
                'enzyme_name': [],
                'active': [],
                'x': [],
                'y': [],
                'info': [],
                'id': []}

        for entry in activity_data:
            data['smiles'].append(entry[self.smi_col])
            data['enzyme_name'].append(entry[self.enzyme_col])

            if entry['binary'] == 'True':
                data['active'].append(1)
            elif entry['binary'] == 'False':
                data['active'].append(0)
            else:
                data['active'].append('nan')

            data['x'].append(entry['x'])
            data['y'].append(entry['y'])
            data['info'].append('info_goes_here')
            data['id'].append(entry['_id'])

        return data

    def create_plot(self, source_act, source_mol, source_seq):

        p_act = plotting.figure(height=500, width=500,
                                tools=['reset,box_zoom,wheel_zoom,zoom_in,zoom_out,pan,tap'],
                                title="Activity")
        p_act.sizing_mode = 'scale_height'
        p_seq = plotting.figure(height=500, width=500,
                                tools=['reset,box_zoom,wheel_zoom,zoom_in,zoom_out,pan,tap'],
                                title="Sequence")

        p_seq.sizing_mode = 'scale_height'

        mapper = linear_cmap(field_name='active', palette=['red', 'green'], low=0, high=1)

        p_act.circle('x', 'y', size=10, name='activity', source=source_mol, alpha=0.01,
                     nonselection_fill_alpha=0,
                     nonselection_line_alpha=0)

        p_act.circle('x', 'y', size=10, source=source_act, alpha=0.5,
                     color=mapper,
                     nonselection_fill_alpha=0,
                     nonselection_line_alpha=0)

        p_seq.circle('x', 'y', size=10, source=source_seq,
                     nonselection_fill_alpha=0,
                     nonselection_line_alpha=0,
                     name='enzymes')

        lasso_select1 = LassoSelectTool()
        lasso_select2 = LassoSelectTool()

        HT = HoverTool(names=['activity'])
        HT.tooltips = [("", "<div>@svg{safe}</div>")]
        HT.mode = 'mouse'
        p_act.add_tools(HT)
        p_act.add_tools(lasso_select1)


        HT2 = HoverTool(names=['enzymes'])
        HT2.tooltips = [("", "@enzyme_name")]
        HT2.mode = 'mouse'
        p_seq.add_tools(HT2)
        p_seq.add_tools(lasso_select2)

        source_seq.selected.js_on_change('indices',
                                         CustomJS(args=dict(s1=source_act, s2=source_seq, sm=source_mol), code="""

                    var seq_inds = cb_obj.indices;

                    var enzyme_names = []
                    var smi_names = []
                    s2.data.enzyme_name.forEach(function (name, index) {
                      if (seq_inds.includes(index)) {
                          enzyme_names.push(name)
                      }
                    });
                    var enzyme_names_set = new Set(enzyme_names);

                    var act_inds = []
                    s1.data.enzyme_name.forEach(function (name, index) {
                        if (enzyme_names_set.has(name)) {
                          act_inds.push(index)
                        }
                    })
                    
                    var smi_names = []
                    s1.data.smiles.forEach(function (smi, index) {
                      if (act_inds.includes(index)) {
                          smi_names.push(smi)
                      }
                    });
                    var smi_names_set = new Set(smi_names);

                    s1.selected.indices = act_inds;
                    s1.change.emit();
                    
                    var selected_seqs = new Set(enzyme_names)
                    window.selected_seqs = Array.from(selected_seqs);
                    window.selected_mols = Array.from(smi_names_set);
                    
                    if (window.selected_seqs.length == 0) {
                        window.selected_seqs = window.all_seqs
                    }
                    
                    if (window.selected_mols.length == 0) {
                        window.selected_mols = window.all_mols
                    }
                    
                     """))

        source_act.selected.js_on_change('indices',
                                         CustomJS(args=dict(s1=source_act, s2=source_seq, sm=source_mol), code="""

                    var act_inds = cb_obj.indices;

                    var enzyme_names = []
                    s1.data.enzyme_name.forEach(function (name, index) {
                      if (act_inds.includes(index)) {
                          enzyme_names.push(name)
                      }
                    });
                    var enzyme_names_set = new Set(enzyme_names);

                    var seq_inds = []
                    s2.data.enzyme_name.forEach(function (name, index) {
                        if (enzyme_names_set.has(name)) {
                          seq_inds.push(index)
                      }
                    })

                    s2.selected.indices = seq_inds;
                    s2.change.emit();

                    var smi_names = []
                    s1.data.smiles.forEach(function (smi, index) {
                      if (act_inds.includes(index)) {
                          smi_names.push(smi)
                      }
                    });
                    var smi_names_set = new Set(smi_names);

                    var smi_inds = []
                    sm.data.smiles.forEach(function (smi, index) {
                        if (smi_names_set.has(smi)) {
                          smi_inds.push(index)
                      }
                    })

                    sm.selected.indices = smi_inds;
                    sm.change.emit();
                    
                    var selected_seqs = new Set(enzyme_names)
                    window.selected_seqs = Array.from(selected_seqs);
                    window.selected_mols = Array.from(smi_names_set);
                    
                    if (window.selected_seqs.length == 0) {
                        window.selected_seqs = window.all_seqs
                    }
                    
                    if (window.selected_mols.length == 0) {
                        window.selected_mols = window.all_mols
                    }
                    
                    """))

        fig_layout = row(p_act, p_seq)
        fig_layout.sizing_mode = 'scale_height'

        p_act.toolbar.active_drag = lasso_select1
        p_seq.toolbar.active_drag = lasso_select2

        return fig_layout

    def create_figure(self):
        unique_smiles = get_data.get_unique_smiles(self.activity_data, self.smi_col)
        if len(unique_smiles) > 2:
            fps = get_data.get_substrates_fps(unique_smiles)
            fp_dict = get_data.make_fp_dict(unique_smiles, fps)
            xy_dict, fp_xs, fp_ys = self.get_fp_pca(unique_smiles, fp_dict)
        else:
            fps = []
            fp_dict = {smi: 'nan' for smi in unique_smiles}
            xy_dict = {smi: (0, 0) for smi in unique_smiles}
            fp_xs = [0 * len(unique_smiles)]
            fp_ys = [0 * len(unique_smiles)]

        enzymes = get_data.get_unique_enzymes(self.activity_data, self.enzyme_col)
        embeddings_map = get_data.get_enzyme_embeddings(enzymes)

        enz_with_embeddings = list(embeddings_map.keys())
        embeddings = list(embeddings_map.values())

        if len(enz_with_embeddings) > 2:
            emb_xs, emb_ys = self.get_embeddings_pca(embeddings)
        else:
            emb_xs, emb_ys = [0 * len(enz_with_embeddings)], [0 * len(enz_with_embeddings)]

        other_enzymes = [enz for enz in enzymes if enz not in enz_with_embeddings]

        for i, data in enumerate(self.activity_data):
            smi = data[self.smi_col]
            self.activity_data[i]['fp'] = fp_dict[smi]
            self.activity_data[i]['x'] = xy_dict[smi][0]
            self.activity_data[i]['y'] = xy_dict[smi][1]

        svg_map, svg_url_map = get_data.get_mol_svgs(unique_smiles, molSize=(250, 100))

        act_data = self.activity_data_to_plot_format(self.activity_data)

        source_act = ColumnDataSource(data=act_data)

        source_mol = ColumnDataSource(data={'smiles': unique_smiles,
                                            'svg': list(svg_map.values()),
                                            'x': fp_xs,
                                            'y': fp_ys})

        source_seq = ColumnDataSource(data={'enzyme_name': enz_with_embeddings,
                                            'x': emb_xs,
                                            'y': emb_ys})

        layout = self.create_plot(source_act, source_mol, source_seq)

        return layout

    def log(self, msg, level=1):
        if self.log_level >= level:
            print(f"SpaceViewer ({level}): {msg}")

if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction',
                            log_level=1)

    sv = SpaceViewer(dq, log_level=1)
    layout = sv.create_figure()

    from bokeh.plotting import figure, output_file, save
    output_file('space_viewer_test.html')
    save(layout)