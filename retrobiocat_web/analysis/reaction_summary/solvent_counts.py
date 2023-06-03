import copy
import math

from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool
from retrobiocat_web.analysis.data_query import get_data
from bokeh.plotting import figure

class Solvent_Processor():

    def __init__(self):
        pass

    def parse_solvents(self, solvents):
        solvents = self._drop_non_strs(solvents)
        solvents = [s.lower() for s in solvents]
        solvents = self._drop_empty_strs(solvents)
        return solvents

    def _drop_non_strs(self, solvents):
        strs_only = []
        for sol in solvents:
            if isinstance(sol, str):
                strs_only.append(sol)
        return strs_only

    def _drop_empty_strs(self, solvents):
        not_empty = []
        for sol in solvents:
            if sol != '':
                not_empty.append(sol)
        return not_empty


class SolventCountGraphCreator():

    def __init__(self):
        pass

    def _get_solvents(self, df):
        if 'solvent' not in df.columns:
            return []
        return list(df.solvent)

    def create_graph(self, data_query):
        df = data_query.get_activity_df()
        solvents = self._get_solvents(df)
        parsed_solvents = Solvent_Processor().parse_solvents(solvents)

        solvent_count_dict = self._count_solvents(parsed_solvents)
        solvent_categories = list(solvent_count_dict.keys())
        solvent_counts = list(solvent_count_dict.values())
        solvent_categories, solvent_counts = self._sort_alphabetically(solvent_categories, solvent_counts)

        data = {'solvents': solvent_categories,
                'counts': solvent_counts}
        source = ColumnDataSource(data)

        max_width = 900
        min_width = 200
        frame_width = int(len(solvent_categories) * 25)
        if frame_width > max_width:
            frame_width = max_width
        if frame_width < min_width:
            frame_width = min_width


        p = figure(title=f"Solvent usage ({len(parsed_solvents)} out of {len(solvents)})",
                   x_range=solvent_categories, toolbar_location=None, tools="",
                   frame_width=frame_width, frame_height=400)

        p.vbar(x='solvents', top='counts', width=0.9, source=source,
               hover_alpha=0.65, name='solvent_counts',
               color="#8739db", line_width=1.5, line_color='black')

        p.xaxis.axis_label = "Solvent"
        p.yaxis.axis_label = "Count"
        p.y_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"

        p.xaxis.major_label_orientation = math.pi/2 - 0.7

        HT = HoverTool(names=['solvent_counts'])
        HT.tooltips = [('Solvent', '@solvents'),
                       ('Count', '@counts')]
        HT.mode = 'mouse'
        p.add_tools(HT)

        callback = CustomJS(args=dict(source=source),
                            code=f"""var ind = cb_data.source.selected.indices[0];
                                     var solvent = source.data.solvents[ind]
                                     launch_solvent_modal(solvent)
                                     """)

        TT = TapTool(names=['solvent_counts'])
        TT.callback = callback
        p.add_tools(TT)

        return p

    def _count_solvents(self, list_of_solvents):
        solvent_count_dict = {}
        for solvent in set(list_of_solvents):
            solvent_count_dict[solvent] = list_of_solvents.count(solvent)
        return solvent_count_dict

    def _sort_alphabetically(self, categories, counts):
        cat_counts = zip(categories, counts)
        sorted_cat_counts = sorted(cat_counts, key=lambda x: x[0], reverse=False)
        categories = [cat for cat, count in sorted_cat_counts]
        counts = [count for cat, count in sorted_cat_counts]
        return categories, counts


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()
    from bokeh.plotting import output_file, save

    dq = get_data.DataQuery(enzyme_type='CAR', reaction='Carboxylic acid reduction', log_level=1)

    solvent_graph_creator = SolventCountGraphCreator()
    count_graph = solvent_graph_creator.create_graph(dq)
    output_file('solvent_graph_test.html')
    save(count_graph)
