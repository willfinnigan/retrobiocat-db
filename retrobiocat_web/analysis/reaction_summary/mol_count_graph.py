import copy
import itertools
import math
from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool, ImageURL
from bokeh.palettes import Spectral6, viridis
from retrobiocat_web.analysis.data_query import get_data
from bokeh.plotting import figure

class MolCountGraphCreator():

    def __init__(self):
        pass

    def create_graph(self, data_query, only_top=30, min_categories=9):
        mol_counts = data_query.mol_counts(only_top=only_top)
        mol_names = list(mol_counts.keys())
        svg_map, svg_url_map = get_data.get_mol_svgs(mol_names, molSize=(250, 100))
        high = [c['high'] for c in mol_counts.values()]
        medium = [c['medium'] for c in mol_counts.values()]
        low = [c['low'] for c in mol_counts.values()]
        active = [c['active'] for c in mol_counts.values()]
        inactive = [c['inactive'] for c in mol_counts.values()]
        total = [c['total'] for c in mol_counts.values()]
        total_active = [c['high']+c['medium']+c['low']+c['active'] for c in mol_counts.values()]
        categories = ['active (no cat.)', 'high', 'medium', 'low', 'inactive']
        colors = ["#cdcfd1", "#32ab34", "#a6ba54",  "#edae61", "#a84232"]
        svgs = [svg_map.get(smi, smi) for smi in mol_names]
        svg_urls = [svg_url_map.get(smi, smi) for smi in mol_names]

        data = {'mols': mol_names,
                'svgs': svgs,
                'svg_urls': svg_urls,
                'high': high,
                'medium': medium,
                'low': low,
                'active (no cat.)': active,
                'total': total,
                'total_active': total_active,
                'inactive': inactive,
                'mol_x': [0 for mol in mol_names],
                'mol_y': [i+0.5 for i in range(len(mol_names))]}
        source = ColumnDataSource(data)

        # will pad categories to min size
        y_range = copy.copy(mol_names)
        if len(y_range) < min_categories:
            for i in range(min_categories - len(y_range)):
                j = i + 1
                y_range.append(' ' * j)

        mol_width = total[0]*0.11
        x_end = total[0] + total[0]*0.033

        p = figure(title=f"Number of entries for {data_query.smi_col} (top {len(mol_names)} of {len(data_query.unique_smiles())})",
                   y_range=y_range, x_range=(-mol_width, x_end),
                   toolbar_location=None, tools="",
                   frame_width=555, frame_height=int(len(y_range)*25))
        p.hbar_stack(categories, y='mols', height=0.9, source=source, color=colors,
                     legend_label=categories,
                     line_color='black', line_width=0.5)

        p.hbar(y='mols', right='total', height=0.9, source=source,
               color=(0, 0, 0), alpha=0.0, hover_alpha=0.3,
               name='mol_counts')


        image1 = ImageURL(url='svg_urls', x='mol_x', y='mol_y', anchor='center_right', h=1, w=mol_width)
        p.add_glyph(source, image1)

        p.line(x=[0, 0], y=[0, len(y_range)], line_color='black', line_width=1)
        p.line(x=[0, x_end], y=[len(y_range), len(y_range)], line_color='black', line_width=1)
        p.line(x=[x_end, x_end], y=[0, len(y_range)], line_color='black', line_width=1)

        p.xgrid.visible = False
        p.ygrid.visible = False
        p.yaxis.visible = False

        HT = HoverTool(names=['mol_counts'])
        HT.tooltips = [('Product mol', '@svgs{safe}'),
                       ("Total entries", '@total'),
                       ("High", '@high'),
                       ("Medium", '@medium'),
                       ("Low", '@low'),
                       ("Total active", '@total_active'),
                       ("Inactive", '@inactive')]
        HT.mode = 'mouse'
        p.add_tools(HT)

        callback = CustomJS(args=dict(source=source),
                            code=f"""var ind = cb_data.source.selected.indices[0];
                                     var smi = source.data.mols[ind]
                                     var svg = source.data.svgs[ind]
                                     molecule_bar_clicked(smi, svg)
                                    """)

        TT = TapTool(names=['mol_counts'])
        TT.callback = callback
        p.add_tools(TT)


        return p


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    from bokeh.plotting import output_file, save

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction', log_level=1)

    mol_count_graph_creator = MolCountGraphCreator()
    count_graph = mol_count_graph_creator.create_graph(dq)
    output_file('mol_count_graph_test.html')
    save(count_graph)
