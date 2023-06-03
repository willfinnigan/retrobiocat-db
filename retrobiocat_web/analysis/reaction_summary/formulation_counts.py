import copy
import math

from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool
from retrobiocat_web.analysis.data_query import get_data
from bokeh.plotting import figure

class Formulation_Processor():

    def __init__(self):
        pass

    def parse_solvents(self, formulations):
        formulations = self._drop_non_strs(formulations)
        formulations = [s.lower() for s in formulations]
        formulations = self._drop_empty_strs(formulations)
        return formulations

    def _drop_non_strs(self, formulations):
        strs_only = []
        for sol in formulations:
            if isinstance(sol, str):
                strs_only.append(sol)
        return strs_only

    def _drop_empty_strs(self, formulations):
        not_empty = []
        for sol in formulations:
            if sol != '':
                not_empty.append(sol)
        return not_empty


class FormulationCountGraphCreator():

    def __init__(self):
        pass

    def _get_formulations(self, df):
        if 'formulation' not in df.columns:
            return []
        return list(df.formulation)

    def create_graph(self, data_query):
        df = data_query.get_activity_df()
        formulations = self._get_formulations(df)
        parsed_formulations = Formulation_Processor().parse_solvents(formulations)

        formulation_count_dict = self._count_formulations(parsed_formulations)
        formulation_categories = list(formulation_count_dict.keys())
        formulation_counts = list(formulation_count_dict.values())
        formulation_categories, formulation_counts = self._sort_alphabetically(formulation_categories, formulation_counts)

        data = {'formulations': formulation_categories,
                'counts': formulation_counts}
        source = ColumnDataSource(data)

        max_width = 900
        min_width = 200
        frame_width = int(len(formulation_categories) * 25)
        if frame_width > max_width:
            frame_width = max_width
        if frame_width < min_width:
            frame_width = min_width

        p = figure(title=f"Formulations ({len(parsed_formulations)} out of {len(formulations)})",
                   x_range=formulation_categories, toolbar_location=None, tools="",
                   frame_width=frame_width, frame_height=400)

        p.vbar(x='formulations', top='counts', width=0.9, source=source,
               hover_alpha=0.65, name='formulation_counts',
               color="#43c452", line_width=1.5, line_color='black')

        p.xaxis.axis_label = "Formulation"
        p.yaxis.axis_label = "Count"
        p.y_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"

        p.xaxis.major_label_orientation = math.pi / 2 - 0.7

        HT = HoverTool(names=['formulation_counts'])
        HT.tooltips = [('formulation', '@formulations'),
                       ('Count', '@counts')]
        HT.mode = 'mouse'
        p.add_tools(HT)

        callback = CustomJS(args=dict(source=source),
                            code=f"""var ind = cb_data.source.selected.indices[0];
                                               var formulation = source.data.formulations[ind]
                                               launch_formulation_modal(formulation)
                                               """)

        TT = TapTool(names=['formulation_counts'])
        TT.callback = callback
        p.add_tools(TT)

        return p

    def _count_formulations(self, list_of_formulations):
        formulation_count_dict = {}
        for formulation in set(list_of_formulations):
            formulation_count_dict[formulation] = list_of_formulations.count(formulation)
        return formulation_count_dict

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

    formulation_graph_creator = FormulationCountGraphCreator()
    count_graph = formulation_graph_creator.create_graph(dq)
    output_file('formulation_graph_test.html')
    save(count_graph)
