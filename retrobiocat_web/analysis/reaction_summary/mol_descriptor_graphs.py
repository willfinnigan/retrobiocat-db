import numpy as np
from bokeh.io import output_file

from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.mongo.model_queries.mol_stats_queries import get_descriptors_for_list_smis

from bokeh.models import ColumnDataSource, HoverTool
from bokeh.plotting import figure
from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.logging import add_logger


class MolDescriptorHistrogramCreator():

    def __init__(self):
        self.bins = 'auto'
        self.colors = {'total_active': "#1613d1",
                       'inactive': "#a84232"}
        self.continuous_descriptors = ['mw', 'logp']
        self.discrete_descriptors = ['num_nh_or_oh', 'num_n_or_o', 'num_h_acceptors', 'num_h_donors', 'num_rings']

    def create_graphs(self, data_query):
        df = data_query.get_activity_df()
        mol_counts_dict = data_query.mol_counts()
        product_smis = list(df.product_1_smiles)
        descriptors_dict = get_descriptors_for_list_smis(product_smis)  # descriptors will be a dict {smi: descriptors, ..}

        graphs = []
        for descriptor_name in self.continuous_descriptors:
            graph = self._create_hist_graph(descriptors_dict, mol_counts_dict, descriptor_name)
            graphs.append(graph)

        for descriptor_name in self.discrete_descriptors:
            graph = self._create_bar_graph(descriptors_dict, mol_counts_dict, descriptor_name)
            graphs.append(graph)

        return graphs

    def _create_hist_graph(self, descriptors_dict, mol_counts_dict, name_descriptor):
        hist_data = self._get_hist_data(descriptors_dict, mol_counts_dict, name_descriptor)
        source = ColumnDataSource(hist_data)
        p = self._create_plot(name_descriptor)
        self._add_hist_bars(p, source, name_descriptor)

        return p

    def _create_bar_graph(self, descriptors_dict, mol_counts_dict, name_descriptor):
        bar_data = self._get_bar_data(descriptors_dict, mol_counts_dict, name_descriptor)
        source = ColumnDataSource(bar_data)
        p = self._create_plot(name_descriptor)
        self._add_bar_chart(p, source, name_descriptor)
        return p

    def _create_plot(self, descriptor_name):
        p = figure(
            title=f"{descriptor_name} product distribution",
            tools='', toolbar_location=None,
            frame_width=300, frame_height=300)

        p.y_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"
        p.xaxis.axis_label = f"{descriptor_name}"
        p.yaxis.axis_label = 'Count'
        return p

    def _add_bar_chart(self, plot, source, descriptor_name):
        plot.vbar(source=source, x='values', top='stacked_total_active', fill_color=self.colors['total_active'],
                  line_width=1, line_color='black', legend_label="Active")
        plot.vbar(source=source, x='values', top='inactive', fill_color=self.colors['inactive'], line_width=1,
                  line_color='black', legend_label="Inactive")
        plot.vbar(source=source, x='values', top='stacked_total_active', color=(0, 0, 0), alpha=0.0, hover_alpha=0.3,
                  name='hover_bar')

        HT = HoverTool(names=['hover_bar'])
        HT.tooltips = [(f'{descriptor_name}', '@values'),
                       ('Num. molecules', '@counts'),
                       ("High", '@high'),
                       ("Medium", '@medium'),
                       ("Low", '@low'),
                       ("Total active", '@total_active'),
                       ("Inactive", '@inactive')]
        HT.mode = 'mouse'
        plot.add_tools(HT)


    def _add_hist_bars(self, plot, source, descriptor_name):
        plot.quad(source=source, bottom='inactive', top='stacked_total_active', left='lefts', right='rights',
               fill_color=self.colors['total_active'], line_width=1, line_color='black', legend_label="Active")
        plot.quad(source=source, bottom=0, top='inactive', left='lefts', right='rights',
                  fill_color=self.colors['inactive'], line_width=1, line_color='black', legend_label="Inactive")
        plot.quad(source=source, bottom=0, top='stacked_total_active', left='lefts', right='rights',
               color=(0, 0, 0), alpha=0.0, hover_alpha=0.3, name='hover_bar')

        HT = HoverTool(names=['hover_bar'])
        HT.tooltips = [(f'{descriptor_name}', '@lefts to @rights'),
                       ('Num. molecules', '@counts'),
                       ("High", '@high'),
                       ("Medium", '@medium'),
                       ("Low", '@low'),
                       ("Total active", '@total_active'),
                       ("Inactive", '@inactive')]
        HT.mode = 'mouse'
        plot.add_tools(HT)

    def _get_bar_data(self, descriptors_dict, mol_counts_dict, name_descriptor):
        descriptor_values = {}
        for smi in descriptors_dict.keys():
            value = descriptors_dict[smi][name_descriptor]
            if value not in descriptor_values:
                descriptor_values[value] = {'count': 0,
                                            'high': 0,
                                            'medium': 0,
                                            'low': 0,
                                            'active': 0,
                                            'inactive': 0,
                                            'total_active': 0,
                                            'stacked_total_active': 0}

            descriptor_values[value]['count'] += 1
            descriptor_values[value]['high'] += mol_counts_dict[smi]['high']
            descriptor_values[value]['medium'] += mol_counts_dict[smi]['medium']
            descriptor_values[value]['low'] += mol_counts_dict[smi]['low']
            descriptor_values[value]['active'] += mol_counts_dict[smi]['active']
            descriptor_values[value]['inactive'] += mol_counts_dict[smi]['inactive']

            total_active = mol_counts_dict[smi]['high'] + mol_counts_dict[smi]['medium'] + mol_counts_dict[smi]['low'] + \
                           mol_counts_dict[smi]['active']
            descriptor_values[value]['total_active'] += total_active
            descriptor_values[value]['stacked_total_active'] += total_active + mol_counts_dict[smi]['inactive']


        values, high, medium, low, active, total_active, stacked_total_active, inactive, counts = [], [], [], [], [], [], [], [], []
        for k, v in descriptor_values.items():
            values.append(k)
            high.append(v['high'])
            medium.append(v['medium'])
            low.append(v['low'])
            active.append(v['active'])
            total_active.append(v['total_active'])
            stacked_total_active.append(v['stacked_total_active'])
            inactive.append(v['inactive'])
            counts.append(v['count'])

        data = {'values': values,
                'counts': counts,
                 'high': high,
                 'medium': medium,
                 'low': low,
                 'active': active,
                 'total_active': total_active,
                 'stacked_total_active': stacked_total_active,
                 'inactive': inactive
                 }

        return data

    def _get_hist_data(self, descriptors_dict, mol_counts_dict, name_descriptor):
        descriptor_values = []
        for smi in descriptors_dict.keys():
            value = descriptors_dict[smi][name_descriptor]
            descriptor_values.append(value)

        hist, edges = np.histogram(descriptor_values, bins=self.bins)

        hist_data = {'lefts': edges[:-1],
                     'rights': edges[1:],
                     'counts': hist,
                     'high': [0] * len(hist),
                     'medium': [0] * len(hist),
                     'low': [0] * len(hist),
                     'active': [0] * len(hist),
                     'total_active': [0] * len(hist),
                     'stacked_total_active': [0] * len(hist),
                     'inactive': [0] * len(hist)
                     }

        for smi in descriptors_dict:
            value = descriptors_dict[smi][name_descriptor]
            bin = self._get_bin(value, edges)
            hist_data['high'][bin] += mol_counts_dict[smi]['high']
            hist_data['medium'][bin] += mol_counts_dict[smi]['medium']
            hist_data['low'][bin] += mol_counts_dict[smi]['low']
            hist_data['active'][bin] += mol_counts_dict[smi]['active']
            hist_data['inactive'][bin] += mol_counts_dict[smi]['inactive']
            total_active = mol_counts_dict[smi]['high'] + mol_counts_dict[smi]['medium'] + mol_counts_dict[smi]['low'] + mol_counts_dict[smi]['active']
            hist_data['total_active'][bin] += total_active
            hist_data['stacked_total_active'][bin] += total_active + mol_counts_dict[smi]['inactive']

        return hist_data

    def _get_bin(self, value, edges):
        for i, edge in enumerate(edges[:-1]):
            if (value >= edges[i]) and (value <= edges[i+1]):
                return i


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from bokeh.plotting import output_file, save

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction', log_level=1)

    mol_des_hist_creator = MolDescriptorHistrogramCreator()
    descriptor_graphs = mol_des_hist_creator.create_graphs(dq)

    output_file('descriptor_test.html')
    save(descriptor_graphs[0])

