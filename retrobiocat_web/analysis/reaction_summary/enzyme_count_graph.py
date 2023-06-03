import copy
import math

from bokeh.models import HoverTool, TapTool, CustomJS, Legend, FixedTicker, CategoricalTicker
from bokeh.models.sources import ColumnDataSource
from retrobiocat_web.analysis.data_query import get_data
from bokeh.plotting import figure

class EnzymeCountGraphCreator():

    def __init__(self):
        pass

    def create_graph(self, data_query, only_top=30, min_categories=9):
        enzyme_counts = data_query.enzyme_counts(only_top=only_top)
        enzyme_names = list(enzyme_counts.keys())
        high = [c['high'] for c in enzyme_counts.values()]
        medium = [c['medium'] for c in enzyme_counts.values()]
        low = [c['low'] for c in enzyme_counts.values()]
        active = [c['active'] for c in enzyme_counts.values()]
        inactive = [c['inactive'] for c in enzyme_counts.values()]
        total = [c['total'] for c in enzyme_counts.values()]
        total_active = [c['high'] + c['medium'] + c['low'] + c['active'] for c in enzyme_counts.values()]
        enzyme_types = [c['enzyme_type'] for c in enzyme_counts.values()]
        categories = ['active (no cat.)', 'high', 'medium', 'low', 'inactive']
        colors = ["#cdcfd1", "#32ab34", "#a6ba54",  "#edae61", "#a84232"]

        data = {'enzyme_names': enzyme_names,
                'enzyme_types': enzyme_types,
                'high': high,
                'medium': medium,
                'low': low,
                'active (no cat.)': active,
                'inactive': inactive,
                'total': total,
                'total_active': total_active}
        source = ColumnDataSource(data)

        # will pad categories to min size
        y_range = copy.copy(enzyme_names)
        if len(y_range) < min_categories:
            for i in range(min_categories-len(y_range)):
                j = i+1
                y_range.append(' '*j)

        p = figure(title=f"Number of entries per enzyme (top {len(enzyme_names)} of {len(data_query.unique_enzymes())})", y_range=y_range, toolbar_location=None, tools="",
                   frame_width=500, frame_height=int(len(y_range)*25))
        #p.yaxis.major_tick_line_color = None
        p.hbar_stack(categories, y='enzyme_names', height=0.9, source=source, color=colors,
                     legend_label=categories, line_color='black', line_width=0.5)

        p.hbar(y='enzyme_names', right='total', height=0.9, source=data,
               color=(0, 0, 0), alpha=0.0, hover_alpha=0.3,
               name='enzyme_counts')

        p.x_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"

        HT = HoverTool(names=['enzyme_counts'])
        HT.tooltips = [('Enzyme type', '@enzyme_types'),
                       ("Enzyme name", '@enzyme_names'),
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
                                     var enzyme_name = source.data.enzyme_names[ind]
                                    load_enzyme_modal(enzyme_name)
                                    """)

        TT = TapTool(names=['enzyme_counts'])
        TT.callback = callback
        p.add_tools(TT)

        return p


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    from bokeh.plotting import output_file, save

    dq = get_data.DataQuery(reaction='Primary alcohol oxidation', log_level=1)

    enzyme_count_graph_creator = EnzymeCountGraphCreator()
    count_graph = enzyme_count_graph_creator.create_graph(dq)
    output_file('enzyme_count_graph_test.html')
    save(count_graph)
