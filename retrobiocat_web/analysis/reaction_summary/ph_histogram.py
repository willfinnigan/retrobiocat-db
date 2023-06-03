import numpy as np
from bokeh.models import ColumnDataSource, HoverTool, CustomJS, TapTool
from bokeh.plotting import figure

from retrobiocat_web.analysis.data_query import get_data
from retrobiocat_web.logging import add_logger


class pH_Processor():

    def __init__(self, log_level='WARNING'):
        self.logger = add_logger('pH_Processor', level=log_level)

    def parse_phs(self, phs):
        """
        Parse the phs so that all numbers
        - remove empty entries ('')
        - all floats
        """

        parsed_phs = []
        for ph in phs:

            ph = self.parse_single_ph(ph)

            if ph is not None:
                parsed_phs.append(ph)

        return parsed_phs

    def parse_single_ph(self, ph):
        if isinstance(ph, str):
            ph = self._remove_spaces(ph)
            ph = self._empty_str_to_none(ph)

        ph = self._ph_to_float(ph)
        ph = self._is_not_a_nan(ph)
        return ph

    def _empty_str_to_none(self, ph):
        if ph == '':
            return None

        return ph

    def _ph_to_float(self, ph):
        """Converts ph to an float or returns a None"""
        if ph is None:
            return None

        try:
            return float(ph)
        except:
            self.logger.debug(f'{ph} could not be converted to a float')
            return None

    def _is_not_a_nan(self, ph):
        if ph is None:
            return None
        elif np.isnan(ph):
            return None

        return ph

    def _remove_spaces(self, ph):
        return ph.replace(' ', '')

class pHHistrogramCreator():

    def __init__(self):
        self.pH_processor = pH_Processor()

    def _get_phs(self, df):
        if 'ph' not in df.columns:
            return []

        return list(df.ph)

    def create_graph(self, data_query, bins=15):

        df = data_query.get_activity_df()
        phs = self._get_phs(df)

        try:
            parsed_phs = [ph for ph in phs if not np.isnan(ph)]
        except:
            print("ERROR PARSING PH DATA")
            parsed_phs = []

        #parsed_phs = self.pH_processor.parse_phs(phs)

        hist, edges = np.histogram(parsed_phs, bins=bins)

        data = {'counts': hist,
                'lefts': edges[:-1],
                'rights': edges[1:]}
        source = ColumnDataSource(data)

        p = figure(title=f'Reaction pH histogram ({len(parsed_phs)} out of {len(phs)} rows with valid data)', tools='', toolbar_location=None,
                   frame_width=500, frame_height=400)
        p.quad(source=source, top='counts', bottom=0, left='lefts', right='rights',
               fill_color="#572bd9", line_width=1.5, line_color='black', name='ph_hist',
               hover_alpha=0.65)

        p.xaxis.axis_label = r"$$pH$$"
        p.yaxis.axis_label = r'$$Count$$'
        p.y_range.start = 0
        p.xgrid.visible = False
        p.ygrid.visible = False
        p.outline_line_width = 1
        p.outline_line_color = "black"

        HT = HoverTool(names=['ph_hist'])
        HT.tooltips = [('Count', '@counts'),
                       ('pH', '@lefts to @rights')]
        HT.mode = 'mouse'
        p.add_tools(HT)

        callback = CustomJS(args=dict(source=source),
                            code=f"""var ind = cb_data.source.selected.indices[0];
                                     var ph_left = Math.round(source.data.lefts[ind] * 100) / 100
                                     var ph_right = Math.round(source.data.rights[ind] * 100) / 100
                                     launch_ph_modal(ph_left, ph_right)""")

        TT = TapTool(names=['ph_hist'])
        TT.callback = callback
        p.add_tools(TT)


        return p


if __name__ == '__main__':
    from retrobiocat_web.mongo.default_connection import make_default_connection

    make_default_connection()
    from bokeh.plotting import output_file, save, figure

    dq = get_data.DataQuery(enzyme_type='Lipase', reaction=None, log_level=1)
    ph_hist_maker = pHHistrogramCreator()
    ph_graph = ph_hist_maker.create_graph(dq)

    output_file('ph_graph_test.html')
    save(ph_graph)




